#include "io/sigmf.hpp"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace si {

Capture Capture::from_iq(cfvec iq, double fs_, double fc_, std::optional<UtcTime> t0_, std::string name_,
                         std::optional<double> gain, std::optional<double> bw) {
  Capture c;
  c.name = std::move(name_);
  c.fs = fs_;
  c.fc = fc_;
  c.t0 = t0_;
  c.gain_db = gain;
  c.bandwidth_hz = bw;
  c.meta = {{"global", {{"core:datatype", "cf32_le"}, {"core:sample_rate", fs_}}},
            {"captures", Json::array({{{"core:frequency", fc_}, {"core:sample_start", 0}}})}};
  if (gain) c.meta["global"]["signalvault:gain_db"] = *gain;
  if (bw) c.meta["global"]["signalvault:rf_bandwidth_hz"] = *bw;
  if (t0_) c.meta["captures"][0]["core:datetime"] = t0_->isoformat();
  c.samples = std::make_shared<const cfvec>(std::move(iq));
  return c;
}

cfvec convert_iq(const void* data, size_t n, const std::string& dt) {
  cfvec out(n);
  if (dt == "cf32_le" || dt == "cf32") {
    std::memcpy(out.data(), data, n * sizeof(cf));
  } else if (dt == "cf64_le" || dt == "cf64") {
    const double* p = static_cast<const double*>(data);
    for (size_t i = 0; i < n; ++i) out[i] = cf(static_cast<float>(p[2 * i]), static_cast<float>(p[2 * i + 1]));
  } else if (dt == "ci16_le" || dt == "ci16") {
    const int16_t* p = static_cast<const int16_t*>(data);
    for (size_t i = 0; i < n; ++i) out[i] = cf(static_cast<float>(p[2 * i]), static_cast<float>(p[2 * i + 1]));
  } else if (dt == "ci8") {
    const int8_t* p = static_cast<const int8_t*>(data);
    for (size_t i = 0; i < n; ++i) out[i] = cf(static_cast<float>(p[2 * i]), static_cast<float>(p[2 * i + 1]));
  } else if (dt == "cu8") {
    const uint8_t* p = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < n; ++i)
      out[i] = cf(static_cast<float>(p[2 * i]) - 127.5f, static_cast<float>(p[2 * i + 1]) - 127.5f);
  } else {
    throw std::invalid_argument("unsupported SigMF datatype '" + dt + "'");
  }
  return out;
}

namespace {

size_t sample_bytes(const std::string& dt) {
  if (dt == "cf32_le" || dt == "cf32") return 8;
  if (dt == "cf64_le" || dt == "cf64") return 16;
  if (dt == "ci16_le" || dt == "ci16") return 4;
  if (dt == "ci8" || dt == "cu8") return 2;
  throw std::invalid_argument("unsupported SigMF datatype '" + dt + "'");
}

std::string read_text(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + p.string());
  return std::string(std::istreambuf_iterator<char>(f), {});
}

fs::path meta_of(const std::string& path) {
  fs::path p(path);
  if (fs::is_directory(p)) {
    std::vector<fs::path> metas;
    for (const auto& e : fs::directory_iterator(p))
      if (e.path().extension() == ".sigmf-meta") metas.push_back(e.path());
    if (metas.size() != 1)
      throw std::runtime_error(p.string() + ": expected one .sigmf-meta, found " + std::to_string(metas.size()));
    return metas[0];
  }
  if (p.extension() == ".sigmf-data" || p.extension() == ".png") p.replace_extension(".sigmf-meta");
  return p;
}

}  // namespace

Capture load_sigmf(const std::string& path) {
  fs::path mp = meta_of(path);
  Capture c;
  c.meta = Json::parse(read_text(mp));
  c.meta_path = mp.string();
  std::string fn = mp.filename().string();
  c.name = fn.substr(0, fn.size() - std::string(".sigmf-meta").size());
  const Json& g = c.meta.at("global");
  const Json& cap = c.meta.at("captures").at(0);
  c.fc = cap.at("core:frequency").get<double>();
  c.fs = g.at("core:sample_rate").get<double>();
  c.datatype = g.at("core:datatype").get<std::string>();
  if (cap.contains("core:datetime") && cap["core:datetime"].is_string())
    c.t0 = UtcTime::parse(cap["core:datetime"].get<std::string>());
  if (g.contains("signalvault:gain_db") && g["signalvault:gain_db"].is_number())
    c.gain_db = g["signalvault:gain_db"].get<double>();
  if (g.contains("signalvault:rf_bandwidth_hz") && g["signalvault:rf_bandwidth_hz"].is_number() &&
      g["signalvault:rf_bandwidth_hz"].get<double>() > 0)
    c.bandwidth_hz = g["signalvault:rf_bandwidth_hz"].get<double>();

  fs::path dp = mp;
  dp.replace_extension(".sigmf-data");
  const size_t sb = sample_bytes(c.datatype);
  std::ifstream f(dp, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + dp.string());
  f.seekg(0, std::ios::end);
  const size_t bytes = static_cast<size_t>(f.tellg());
  f.seekg(0);
  const size_t n = bytes / sb;
  if (c.datatype == "cf32_le" || c.datatype == "cf32") {
    auto v = std::make_shared<cfvec>(n);
    f.read(reinterpret_cast<char*>(v->data()), static_cast<std::streamsize>(n * sb));
    c.samples = v;
  } else {
    std::vector<char> raw(n * sb);
    f.read(raw.data(), static_cast<std::streamsize>(raw.size()));
    c.samples = std::make_shared<const cfvec>(convert_iq(raw.data(), n, c.datatype));
  }
  return c;
}

std::vector<std::string> find_captures(const std::string& folder) {
  struct Item { double fc; std::string dt; std::string path; };
  std::vector<Item> items;
  for (const auto& e : fs::recursive_directory_iterator(folder)) {
    if (e.path().extension() != ".sigmf-meta") continue;
    Json m = Json::parse(read_text(e.path()));
    const Json& c = m.at("captures").at(0);
    items.push_back({c.at("core:frequency").get<double>(), c.value("core:datetime", std::string()),
                     e.path().string()});
  }
  std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
    if (a.fc != b.fc) return a.fc < b.fc;
    if (a.dt != b.dt) return a.dt < b.dt;
    return a.path < b.path;
  });
  std::vector<std::string> out;
  for (auto& i : items) out.push_back(i.path);
  return out;
}

}  // namespace si

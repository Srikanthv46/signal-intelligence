// si_analyze — run spectrum_intel on SigMF captures or raw IQ files.
//
//   si_analyze PATH [options]          PATH: .sigmf-meta/.sigmf-data, a folder with one, or a folder tree
//   si_analyze --raw FILE --format ci16 --fs 61.44e6 --fc 2.437e9 [--bw 56e6] [--start ISO] [options]
//
// options:
//   --out DIR|FILE   write <name>_response.json (or _report.json) there; default: stdout (one capture)
//   --report         write the per-block report (python report.writer format) instead of the response
//   --enable a,b     only these receivers ("none": energy detection and band rules only)
//   --disable a,b    all but these
//   --ppm X          the receiver's frequency error in ppm (if known)
//   --threads N      worker threads            --compact       response without provenance detail
//   --indent N       JSON indent (default 1)   --receivers     list receivers   --requirements  print them
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <spectrum_intel/spectrum_intel.hpp>

#include "analyze.hpp"
#include "core/parallel.hpp"
#include "decode/registry.hpp"
#include "io/sigmf.hpp"
#include "report/writer.hpp"

namespace fs = std::filesystem;

namespace {

std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string x;
  while (std::getline(ss, x, ','))
    if (!x.empty()) out.push_back(x);
  return out;
}

[[noreturn]] void usage(const char* msg) {
  std::fprintf(stderr, "si_analyze: %s\nsee the header of apps/si_analyze.cpp for usage\n", msg);
  std::exit(2);
}

si::IqFormat parse_format(const std::string& f) {
  if (f == "cf32" || f == "cf32_le") return si::IqFormat::CF32;
  if (f == "cf64" || f == "cf64_le") return si::IqFormat::CF64;
  if (f == "ci16" || f == "ci16_le") return si::IqFormat::CI16;
  if (f == "ci8") return si::IqFormat::CI8;
  if (f == "cu8") return si::IqFormat::CU8;
  usage(("unknown --format " + f).c_str());
}

const char* sigmf_datatype(si::IqFormat f) {
  switch (f) {
    case si::IqFormat::CF32: return "cf32_le";
    case si::IqFormat::CF64: return "cf64_le";
    case si::IqFormat::CI16: return "ci16_le";
    case si::IqFormat::CI8: return "ci8";
    case si::IqFormat::CU8: return "cu8";
  }
  return "cf32_le";
}

size_t sample_size(si::IqFormat f) {
  switch (f) {
    case si::IqFormat::CF32: return 8;
    case si::IqFormat::CF64: return 16;
    case si::IqFormat::CI16: return 4;
    default: return 2;
  }
}

void emit(const si::Json& j, const std::string& out, const std::string& name, const char* suffix, int indent,
          bool many) {
  const std::string text = j.dump(indent, ' ', false, si::Json::error_handler_t::replace);
  if (out.empty() || out == "-") {
    std::cout << text << "\n";
    return;
  }
  fs::path p(out);
  if (many || fs::is_directory(p) || out.back() == '/') {
    fs::create_directories(p);
    p /= name + suffix;
  }
  std::ofstream(p) << text << "\n";
  std::fprintf(stderr, "wrote %s\n", p.string().c_str());
}

}  // namespace

int main(int argc, char** argv) {
  std::string path, raw, format = "cf32", out, start;
  double fsr = 0, fc = 0, bw = 0;
  std::optional<double> ppm;
  bool report = false, compact = false;
  int threads = 0, indent = 1;
  std::vector<std::string> enable, disable;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) usage(("missing value after " + a).c_str());
      return argv[++i];
    };
    if (a == "--raw") raw = next();
    else if (a == "--format") format = next();
    else if (a == "--fs") fsr = std::atof(next().c_str());
    else if (a == "--fc") fc = std::atof(next().c_str());
    else if (a == "--bw") bw = std::atof(next().c_str());
    else if (a == "--start") start = next();
    else if (a == "--ppm") ppm = std::atof(next().c_str());
    else if (a == "--out") out = next();
    else if (a == "--report") report = true;
    else if (a == "--compact") compact = true;
    else if (a == "--enable") enable = split(next());
    else if (a == "--disable") disable = split(next());
    else if (a == "--threads") threads = std::atoi(next().c_str());
    else if (a == "--indent") indent = std::atoi(next().c_str());
    else if (a == "--receivers") {
      for (const auto& n : si::available_receivers()) std::cout << n << "\n";
      return 0;
    } else if (a == "--requirements") {
      std::cout << si::requirements().dump(1) << "\n";
      return 0;
    } else if (a == "-h" || a == "--help") usage("help");
    else if (!a.empty() && a[0] == '-') usage(("unknown option " + a).c_str());
    else path = a;
  }
  if (threads > 0) si::set_threads(threads);
  const bool no_receivers = enable.size() == 1 && enable[0] == "none";
  auto make_plugins = [&] { return no_receivers ? si::PluginList{} : si::plugins(enable, disable); };
  std::optional<si::Json> ppm_cal;
  if (ppm) ppm_cal = si::Json{{"ppm", *ppm}, {"source", "caller (IQBlock.receiver_ppm)"}};
  if (no_receivers && !report) usage("--enable none is for --report runs");

  si::EngineConfig cfg;
  cfg.enable = enable;
  cfg.disable = disable;
  cfg.threads = threads;
  cfg.compact = compact;
  try {
    if (!raw.empty()) {
      if (fsr <= 0) usage("--raw needs --fs and --fc");
      const si::IqFormat f = parse_format(format);
      std::ifstream in(raw, std::ios::binary);
      if (!in) usage(("cannot open " + raw).c_str());
      std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
      si::IQBlock b;
      b.data = bytes.data();
      b.n_samples = bytes.size() / sample_size(f);
      b.format = f;
      b.sample_rate_hz = fsr;
      b.center_frequency_hz = fc;
      if (bw > 0) b.bandwidth_hz = bw;
      if (!start.empty()) b.start_time_utc = start;
      b.receiver_ppm = ppm;
      b.source_id = fs::path(raw).stem().string();
      if (report) {
        si::cfvec x = si::convert_iq(b.data, b.n_samples, sigmf_datatype(f));
        auto cap = si::Capture::from_iq(std::move(x), fsr, fc, std::nullopt, *b.source_id, std::nullopt,
                                        bw > 0 ? std::optional<double>(bw) : std::nullopt);
        auto plugins = make_plugins();
        emit(si::build_report(si::analyze(cap, plugins, ppm_cal)), out, *b.source_id, "_report.json", indent, false);
      } else {
        si::Engine engine(cfg);
        emit(engine.analyze(b), out, *b.source_id, "_response.json", indent, false);
      }
      return 0;
    }
    if (path.empty()) usage("no input");
    std::vector<std::string> metas;
    if (fs::is_directory(path)) {
      metas = si::find_captures(path);
    } else {
      metas.push_back(path);
    }
    const bool many = metas.size() > 1;
    auto plugins = report ? make_plugins() : si::PluginList{};
    std::unique_ptr<si::Engine> engine_ptr;
    if (!report) engine_ptr = std::make_unique<si::Engine>(cfg);
    for (const auto& m : metas) {
      const auto t = std::chrono::steady_clock::now();
      si::Capture cap = si::load_sigmf(m);
      if (report) {
        emit(si::build_report(si::analyze(cap, plugins, ppm_cal)), out, cap.name, "_report.json", indent, many);
      } else {
        si::IQBlock b;
        b.data = cap.iq().data();
        b.n_samples = cap.n();
        b.format = si::IqFormat::CF32;
        b.sample_rate_hz = cap.fs;
        b.center_frequency_hz = cap.fc;
        b.bandwidth_hz = cap.bandwidth_hz;
        b.gain_db = cap.gain_db;
        b.receiver_ppm = ppm;
        if (cap.t0) b.start_time_utc = cap.t0->iso_z();
        b.source_id = cap.name;
        emit(engine_ptr->analyze(b), out, cap.name, "_response.json", indent, many);
      }
      std::fprintf(stderr, "%s: %.2f s\n", cap.name.c_str(),
                   std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count());
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "si_analyze: %s\n", e.what());
    return 1;
  }
  return 0;
}

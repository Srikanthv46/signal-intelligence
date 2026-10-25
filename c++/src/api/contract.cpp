#include "api/contract.hpp"

#include <algorithm>
#include <mutex>
#include <regex>
#include <stdexcept>

#include "core/embedded.hpp"

namespace si::contract {

const char* const kSchemaVersion = "1.4.0";

namespace {

const Json& load(const char* name) {
  // parsed once per file, on first use (thread-safe static initialisation)
  auto parse = [](const char* n) {
    const char* txt = embedded_file(n);
    if (!txt) throw std::runtime_error(std::string("schema file not compiled in: ") + n);
    return Json::parse(txt);
  };
  if (std::string(name) == "technologies.json") {
    static const Json j = parse(name);
    return j;
  }
  if (std::string(name) == "requirements.json") {
    static const Json j = parse(name);
    return j;
  }
  if (std::string(name) == "capabilities.json") {
    static const Json j = parse(name);
    return j;
  }
  if (std::string(name) == "descriptions.json") {
    static const Json j = parse(name);
    return j;
  }
  static const Json j = parse("response.schema.json");
  return j;
}

struct Rule {
  std::string tech, family;
  bool any_decoder;
  std::string decoder;
  bool any_label;
  std::regex rx;
};

const std::vector<Rule>& rules() {
  static const std::vector<Rule> r = [] {
    std::vector<Rule> out;
    for (const auto& j : technologies().at("rules")) {
      Rule x;
      x.tech = j.at("technology").get<std::string>();
      x.family = j.at("family").get<std::string>();
      x.any_decoder = j.at("decoder").is_null();
      if (!x.any_decoder) x.decoder = j.at("decoder").get<std::string>();
      const std::string lab = j.at("label").get<std::string>();
      x.any_label = lab.empty();
      if (!x.any_label) x.rx = std::regex(lab, std::regex::ECMAScript | std::regex::optimize);
      out.push_back(std::move(x));
    }
    return out;
  }();
  return r;
}

}  // namespace

const Json& technologies() { return load("technologies.json"); }
const Json& requirements() { return load("requirements.json"); }
const Json& response_schema() { return load("response.schema.json"); }
const Json& capabilities() { return load("capabilities.json"); }
const Json& descriptions() { return load("descriptions.json"); }

namespace {
bool empty_field(const Json& v) {
  return v.is_null() || ((v.is_array() || v.is_object()) && v.empty()) || (v.is_string() && v.get<std::string>().empty());
}
}  // namespace

Json validation(const std::string& tech, const std::string& decoder, const Json& fields) {
  const Json& reg = capabilities();
  Json out = Json::object();
  const Json& techs = reg.at("technologies");
  const bool known = techs.contains(tech);
  out["level"] = known ? techs.at(tech).at("level") : Json("research");
  if (known && techs.at(tech).contains("proof")) out["proof"] = techs.at(tech).at("proof");
  if (known && techs.at(tech).contains("sources")) out["sources"] = techs.at(tech).at("sources");
  Json rf = Json::array(), src = Json::array();
  for (const auto& c : reg.at("field_capabilities")) {
    if (c.at("level") != "research") continue;
    if (c.contains("decoder") && !c.at("decoder").is_null() && c.at("decoder").get<std::string>() != decoder) continue;
    bool any = false;
    for (const auto& k : c.at("fields")) {
      const std::string key = k.get<std::string>();
      if (fields.is_object() && fields.contains(key) && !empty_field(fields.at(key))) {
        any = true;
        if (std::find(rf.begin(), rf.end(), k) == rf.end()) rf.push_back(k);
      }
    }
    if (any && std::find(src.begin(), src.end(), c.at("sources")) == src.end()) src.push_back(c.at("sources"));
  }
  if (!rf.empty()) {
    out["research_fields"] = rf;
    out["research_sources"] = src;
  }
  return out;
}

std::pair<std::string, std::string> technology(const std::string& label, const std::string& decoder) {
  for (const auto& r : rules()) {
    if (!r.any_decoder && r.decoder != decoder) continue;
    if (!r.any_label && !std::regex_search(label, r.rx)) continue;
    return {r.tech, r.family};
  }
  return {"unknown", "unknown"};
}

const std::vector<std::string>& identifier_keys() {
  static const std::vector<std::string> k = technologies().at("identifier_keys").get<std::vector<std::string>>();
  return k;
}

namespace {
struct Desc {
  std::string tech;
  bool any_label;
  std::regex rx;
  std::string text;
};
const std::vector<Desc>& desc_rules() {
  static const std::vector<Desc> r = [] {
    std::vector<Desc> out;
    for (const auto& e : descriptions().at("entries")) {
      Desc d;
      d.tech = e.at("technology").get<std::string>();
      d.any_label = !e.contains("label");
      if (!d.any_label) d.rx = std::regex(e.at("label").get<std::string>(), std::regex::ECMAScript | std::regex::optimize);
      d.text = e.at("text").get<std::string>();
      out.push_back(std::move(d));
    }
    return out;
  }();
  return r;
}
}  // namespace

std::string description(const std::string& label, const std::string& decoder, const Json& flags) {
  auto [tech, fam] = technology(label, decoder);
  if (fam != "artifact" && flags.is_array()) {
    for (const char* x : kArtifactFlags)
      if (std::any_of(flags.begin(), flags.end(), [&](const Json& v) { return v == x; })) {
        tech = x;
        fam = "artifact";
        break;
      }
  }
  for (const auto& d : desc_rules())
    if (d.tech == tech && (d.any_label || std::regex_search(label, d.rx))) return d.text;
  const Json& fams = technologies().at("families");
  if (!fams.contains(fam)) return "";
  std::string f = fams.at(fam).get<std::string>();
  if (f.empty()) return "";
  if (f[0] >= 'a' && f[0] <= 'z') f[0] = static_cast<char>(f[0] - 'a' + 'A');
  return f + ".";
}

}  // namespace si::contract

#!/usr/bin/env python3
"""Loose coupling of the C++ library, checked (python/tests/test_modularity.py's twin):

  * a receiver's DSP files (src/receivers/<name>/*, except plugin.cpp) include only
    core/, io/ and their own folder;
  * plugin.cpp may also include decode/ (the plugin interface);
  * nothing outside src/receivers/ includes a receiver (the registry is generated);
  * core/ and io/ include nothing from the pipeline (detect, features, classify, decode,
    report, api, analyze) or from receivers;
  * every receiver folder has a plugin.cpp and a README.md.
Exit status 1 on a violation.
"""
import re
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
INC = re.compile(r'^\s*#\s*include\s+"([^"]+)"', re.M)
PIPELINE = ("detect/", "features/", "classify/", "decode/", "report/", "api/", "analyze")


def includes(p: Path) -> list[str]:
    return INC.findall(p.read_text(errors="replace"))


def main() -> int:
    bad = []
    rx_root = SRC / "receivers"
    for d in sorted(x for x in rx_root.iterdir() if x.is_dir()):
        if not (d / "plugin.cpp").is_file():
            bad.append(f"receivers/{d.name}: no plugin.cpp")
        if not (d / "README.md").is_file():
            bad.append(f"receivers/{d.name}: no README.md")
        for f in sorted(d.glob("*.[ch]pp")):
            for inc in includes(f):
                own = inc.startswith(f"receivers/{d.name}/") or "/" not in inc
                ok = own or inc.startswith(("core/", "io/")) or (f.name == "plugin.cpp" and inc.startswith("decode/"))
                if not ok:
                    bad.append(f"receivers/{d.name}/{f.name} includes {inc}")
    for f in sorted(SRC.rglob("*.[ch]pp")):
        rel = f.relative_to(SRC).as_posix()
        if rel.startswith("receivers/"):
            continue
        for inc in includes(f):
            if inc.startswith("receivers/"):
                bad.append(f"{rel} includes a receiver ({inc})")
            if rel.startswith(("core/", "io/")) and inc.startswith(PIPELINE):
                bad.append(f"{rel} (core) includes pipeline header {inc}")
    for b in bad:
        print("VIOLATION", b)
    print(f"{len(bad)} modularity violation(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

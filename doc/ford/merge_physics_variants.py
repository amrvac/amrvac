#!/usr/bin/env python3
"""
Auto-discovers which preprocessed source files differ across AGILE's five
physics modules (mhd/hd/ffhd/srhd/srmhd) and merges five per-physics staging
directories (produced by preprocessing every `.fpp` file once per physics
module) into one FORD source directory:

- Files that are byte-identical across all five passes -- the vast
  majority, anything with no `#:if PHYS == '...'`-guarded content -- are
  kept once, under their real name.
- Files that differ are kept once per physics module, suffixed `_<phys>`
  (both the filename and any `module <name>` / `end module <name>`
  statement inside it), so all five variants can coexist in one FORD site
  without name collisions.

This needs no hand-maintained list of "physics-sensitive" files: whatever
actually differs is discovered by comparison.
"""
import re
import sys
from pathlib import Path

PRIMARY = "mhd"
PHYS_LIST = ["mhd", "hd", "ffhd", "srhd", "srmhd"]

def rename_module(text, base, phys):
    escaped = re.escape(base)
    pattern = re.compile(
        rf"^(?P<indent1>[ \t]*module[ \t]+)(?P<name1>{escaped})(?P<tail1>[ \t]*)$"
        rf"|"
        rf"^(?P<indent2>[ \t]*end[ \t]+module[ \t]+)(?P<name2>{escaped})(?P<tail2>[ \t]*)$",
        re.IGNORECASE | re.MULTILINE,
    )

    def repl(m):
        if m.group("name1") is not None:
            return f"{m.group('indent1')}{base}_{phys}{m.group('tail1')}"
        return f"{m.group('indent2')}{base}_{phys}{m.group('tail2')}"

    return pattern.sub(repl, text)


def main():
    stage_root = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    stage_dirs = {p: stage_root / p for p in PHYS_LIST}
    names = sorted(f.name for f in stage_dirs[PRIMARY].glob("*.f90"))

    n_shared = 0
    n_variant = 0
    for name in names:
        base = name[: -len(".f90")]
        contents = {p: (stage_dirs[p] / name).read_text() for p in PHYS_LIST}
        if len(set(contents.values())) == 1:
            (out_dir / name).write_text(contents[PRIMARY])
            n_shared += 1
        else:
            for p in PHYS_LIST:
                (out_dir / f"{base}_{p}.f90").write_text(
                    rename_module(contents[p], base, p))
            n_variant += 1

    print(f"{n_shared} files identical across all {len(PHYS_LIST)} physics modules (kept once)")
    print(f"{n_variant} files differ by physics module "
          f"(kept as {n_variant * len(PHYS_LIST)} suffixed variants)")


if __name__ == "__main__":
    main()

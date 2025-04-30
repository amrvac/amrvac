"""
Reads a Fortran name list and extracts those parameters that are
relevant compile-time. Output is a Makefile.
"""

import sys
import f90nml
from typing import Any
from collections.abc  import Mapping
from enum import Enum


"""
The `Flag` class is an `Enum` containing all boolean
flags of components that can be enabled/disabled in
`amrvac.par`. Default values should be stored in the
`flag_default` dictionary.
"""
class Flag(Enum):
    GRAVITY = 1

flag_default: dict[Flag, bool] = {
    Flag.GRAVITY: False
}


"""
Tests that the `flag_default` dictionary is complete.
If it is not, the script will fail automatically.
"""
def test_flag_default_completeness():
    assert all(f in flag_default.keys() for f in Flag)


"""
Given a mapping of strings to bools (i.e. components  section
from `amrvac.par`, output all relevant configuration to Make.
"""
def handle_flag(cfg: Mapping[str, bool], flag: Flag):
    val = cfg.get(flag.name.lower(), flag_default[flag])
    if val:
        print(f"enabled += {flag.name}")
        print(f"fypp_flags += -D{flag.name}")


if __name__ == "__main__":
    test_flag_default_completeness()

    print("# Don't edit: this file was automatically generated")
    cfg = f90nml.read(sys.stdin)

    if "components" not in cfg.keys():
        sys.exit()

    components = cfg["components"]
    for f in Flag:
        handle_flag(components, f)


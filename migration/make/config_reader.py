"""
Reads a Fortran name list and extracts those parameters that are
relevant compile-time. Output is a Makefile.
"""

import sys
import f90nml
from typing import Any
from collections.abc  import Mapping
from enum import Enum
from dataclasses import dataclass
from pathlib import Path
import tomllib
import operator
from functools import reduce


@dataclass
class Flag:
    name: str
    flags: list[str]
    default: bool

    def handle(self, val):
        if val:
            if "hash" in self.flags:
                print(f"enabled += {self.name}")
            if "fypp" in self.flags:
                print(f"fypp_flags += -D{self.name}")


@dataclass
class Option:
    name: str
    flags: list[str]
    options: list[str]
    default: str

    def handle(self, val):
        assert(val in self.options)
        if "hash" in self.flags:
            print(f"enabled += \"{self.name}={val}\"")
        if "fypp" in self.flags:
            print(f"fypp_flags += -D{self.name}=\\'{val}\\'")


OPTION_TYPES = {
    "flag": Flag,
    "option": Option
}


def read_schema():
    filename = Path(__file__).parent / "config_schema.toml"
    raw_schema = tomllib.load(open(filename, "rb"))

    def convert(schema, path):
        if not isinstance(schema, dict):
            return {}

        if "type" in schema:
            t = schema["type"]
            if t in OPTION_TYPES:
                del schema["type"]
                return { ".".join(path): OPTION_TYPES[t](**schema) }
                
        return reduce(
            operator.or_,
            (convert(v, path + [k]) for k, v in schema.items()))

    return convert(raw_schema, [])


def apply_option(cfg, path, option):
    def get(cfg, path):
        if not path:
            return cfg
        if path[0] in cfg:
            return get(cfg[path[0]], path[1:])
        return option.default
    val = get(cfg, path.split("."))
    option.handle(val)


if __name__ == "__main__":
    print("# Don't edit: this file was automatically generated")
    schema = read_schema()
    cfg = f90nml.read(sys.stdin)

    for path, option in schema.items():
        apply_option(cfg, path, option)


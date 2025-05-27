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


class ConfigError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return f"Configuration error: {self.msg}"


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
        if val not in self.options:
            raise ConfigError(f"Option parameter {self.name} should be one of {self.options}, got: {val}")
        if "hash" in self.flags:
            print(f"enabled += \"{self.name}={val}\"")
        if "fypp" in self.flags:
            print(f"fypp_flags += -D{self.name}=\\'{val}\\'")


@dataclass
class Value:
    name: str
    flags: list[str]
    dtype: str
    default: Any

    def handle(self, val):
        if (self.dtype == "int" and not isinstance(val, int)) or \
           (self.dtype == "string" and not isinstance(val, string)):
            raise(ConfigError(f"Value parameter {self.name} should be {self.dtype}, got {type(val)}"))
        if self.dtype not in ["int", "string"]:
            raise(ConfigError(f"Unknown value type {self.dtype}"))

        if self.dtype == "string":
            val_str = f"\\'{val}\\'"
        else:
            val_str = f"{val}"

        if "hash" in self.flags:
            print(f"enabled += \"{self.name}={val}\"")
        if "fypp" in self.flags:
            print(f"fypp_flags += -D{self.name}={val_str}")


OPTION_TYPES = {
    "flag": Flag,
    "option": Option,
    "value": Value
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

    try:
        for path, option in schema.items():
            apply_option(cfg, path, option)
    except ConfigError as e:
        print(e, file=sys.stderr)
        sys.exit(1)


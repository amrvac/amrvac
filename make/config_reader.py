"""
Reads a Fortran name list and extracts those parameters that are
relevant compile-time. Output is a Makefile.
"""

import sys
import f90nml
from typing import Any
from collections.abc  import Mapping
from enum import Enum
from dataclasses import dataclass, field
from pathlib import Path
import tomllib
import operator
from functools import reduce


class ConfigError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return f"Configuration error: {self.msg}"


class PathDict(dict):
    """
    A little hackery around the Python `dict` to allow nested
    queries.

    Example:

        >>> d = PathDict({"a": {"b": 1}})
        >>> d["a"]
        {'b': 1}
        >>> d["a.b"]
        1
        >>> "a.b" in d
        True
        >>> "a.c" in d
        False
        >>> d.get("a.c", -1)
        -1
    """
    def __getitem__(self, key):
        path = key.split(".", 1)
        obj = dict.__getitem__(self, path[0])

        if len(path) == 1:
            return obj

        return PathDict(obj)[path[1]]

    def __contains__(self, key):
        path = key.split(".", 1)
        has = dict.__contains__(self, path[0])

        if not has:
            return False

        if len(path) == 1:
            return has

        obj = dict.__getitem__(self, path[0])
        return path[1] in PathDict(obj)

    def get(self, key, default=None):
        if key not in self:
            return default
        return self[key]


@dataclass
class Flag:
    name: str
    flags: list[str]
    default: bool
    implies: list[str] = field(default_factory=list)

    def handle(self, schema, val, visited):
        if not val:
            return

        if self.name in visited:
            return
        visited.add(self.name)

        for flag in self.implies:
            if flag not in schema:
                raise ConfigError(
                    f"flag `{flag}`, implied by `{self.name}` not "
                    f"present in schema")
            schema[flag].handle(schema, True, visited)

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

    def handle(self, schema, val, visited):
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

    def handle(self, schema, val, visited):
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


def apply_option(schema, cfg, path, option, visited):
    val = cfg.get(path, option.default)
    option.handle(schema, val, visited)


if __name__ == "__main__":
    print("# Don't edit: this file was automatically generated")
    schema = read_schema()
    visited = set()
    cfg = PathDict(f90nml.read(sys.stdin))

    try:
        for path, option in schema.items():
            apply_option(schema, cfg, path, option, visited)
    except ConfigError as e:
        print(e, file=sys.stderr)
        sys.exit(1)


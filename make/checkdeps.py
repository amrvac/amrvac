#!/usr/bin/env python3
from __future__ import annotations

# If you make changes to this file, please consider contributing
# to: https://github.com/jhidding/check-deps

import sys
import tomllib
from dataclasses import dataclass, field
from typing import Optional, List, Mapping, Tuple, Callable, TypeVar
from enum import Enum
import asyncio
import re
from contextlib import contextmanager, redirect_stdout
import textwrap
import io

assert sys.version_info[0] == 3, "This script only works with Python 3."


class ConfigError(Exception):
    pass


T = TypeVar("T")


class Relation(Enum):
    """Encodes ordinal relations among versions. Currently six operators are
    supported: `>=`, `<=`, `<`, `>`, `==`, `!=`."""

    GE = 1
    LE = 2
    LT = 3
    GT = 4
    EQ = 5
    NE = 6

    def __str__(self):
        return {
            Relation.GE: ">=",
            Relation.LE: "<=",
            Relation.LT: "<",
            Relation.GT: ">",
            Relation.EQ: "==",
            Relation.NE: "!=",
        }[self]


@dataclass
class Version:
    """Stores a version in the form of a tuple of ints and an optional string extension.
    This class supports (in)equality operators listed in `Relation`."""

    number: tuple[int, ...]
    extra: Optional[str]

    def __lt__(self, other):
        for n, m in zip(self.number, other.number):
            if n < m:
                return True
            elif n > m:
                return False
        return False

    def __gt__(self, other):
        return other < self

    def __le__(self, other):
        for n, m in zip(self.number, other.number):
            if n < m:
                return True
            elif n > m:
                return False
        return True

    def __ge__(self, other):
        return other <= self

    def __eq__(self, other):
        for n, m in zip(self.number, other.number):
            if n != m:
                return False
            return True

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return ".".join(map(str, self.number)) + (self.extra or "")


@dataclass
class VersionConstraint:
    """A VersionConstraint is a product of a `Version` and a `Relation`."""

    version: Version
    relation: Relation

    def __call__(self, other: Version) -> bool:
        method = f"__{self.relation.name}__".lower()
        return getattr(other, method)(self.version)

    def __str__(self):
        return f"{self.relation}{self.version}"


def split_at(split_chars: str, x: str) -> Tuple[str, str]:
    """Tries to split at character `x`. Returns a 2-tuple of the string
    before and after the given separator."""
    a = x.split(split_chars, maxsplit=1)
    if len(a) == 2:
        return a[0], a[1]
    else:
        return a[0], ""


def parse_split_f(split_chars: str, f: Callable[[str], T], x: str) -> Tuple[T, str]:
    """Given a string, splits at given character `x` and passes the left value
    through a function (probably a parser). The second half of the return tuple is the
    remainder of the string."""
    item, x = split_at(split_chars, x)
    val = f(item)
    return val, x


def parse_version(x: str) -> Tuple[Version, str]:
    """Parse a given string to a `Version`. A sequence of dot `.` separated integers
    is put into the numerical version component, while the remaining text ends up in
    the `extra` component."""
    _x = x
    number = []
    extra = None

    while True:
        try:
            n, _x = parse_split_f(".", int, _x)
            number.append(n)
        except ValueError:
            if len(x) > 0:
                m = re.match("([0-9]*)(.*)", _x)
                if lastn := m and m.group(1):
                    number.append(int(lastn))
                if suff := m and m.group(2):
                    extra = suff or None
                else:
                    extra = _x
            break

    if not number:
        raise ConfigError(f"A version needs a numeric component, got: {x}")

    return Version(tuple(number), extra), _x


def parse_relation(x: str) -> Tuple[Relation, str]:
    """Parses the operator of the version constraint."""
    op_map = {
        "<=": Relation.LE,
        ">=": Relation.GE,
        "<": Relation.LT,
        ">": Relation.GT,
        "==": Relation.EQ,
        "!=": Relation.NE,
    }
    for sym, op in op_map.items():
        if x.startswith(sym):
            return (op, x[len(sym) :])
    raise ConfigError(f"Not a comparison operator: {x}")


def parse_version_constraint(x: str) -> Tuple[VersionConstraint, str]:
    relation, x = parse_relation(x)
    version, x = parse_version(x)
    return VersionConstraint(version, relation), x


def async_cache(f):
    """Caches results from the `async` function `f`. This assumes `f` is a
    member of a class, where we have `_lock`, `_result` and `_done` members
    available."""

    async def g(self, *args, **kwargs):
        async with self._lock:
            if self._done:
                return self._result
            self._result = await f(self, *args, **kwargs)
            self._done = True
            return self._result

    return g


@dataclass
class Result:
    test: VersionTest
    success: bool
    failure_text: Optional[str] = None
    found_version: Optional[Version] = None

    def __bool__(self):
        return self.success


@dataclass
class VersionTest:
    name: str
    require: VersionConstraint
    get_version: str
    platform: Optional[str] = None
    pattern: Optional[str] = None
    suggestion_text: Optional[str] = None
    suggestion: Optional[str] = None
    depends: List[str] = field(default_factory=list)
    template: Optional[str] = None

    _lock: asyncio.Lock = field(default_factory=asyncio.Lock)
    _done: bool = False

    def print_formatted(self, msg):
        prefix = f"{self.name} {self.require}"
        print(f"{prefix:25}: {msg}")

    def print_not_found(self):
        self.print_formatted("not found")

    @async_cache
    async def run(self, recurse):
        for dep in self.depends:
            if not await recurse(dep):
                return Result(self, False, failure_text=f"Failed dependency: {dep}")

        proc = await asyncio.create_subprocess_shell(
            self.get_version,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        (stdout, stderr) = await proc.communicate()
        if proc.returncode != 0:
            self.print_not_found()
            return Result(
                self, success=False, failure_text=f"{stderr.decode().strip()}"
            )
        try:
            if self.pattern is not None:
                m = re.match(self.pattern, stdout.decode())
                if m is not None:
                    out, _ = parse_version(m.group(1).strip())
                else:
                    self.print_not_found()
                    msg = f"No regex match on pattern '{self.pattern}'"
                    return Result(self, False, failure_text=msg)
            else:
                out, _ = parse_version(stdout.decode().strip())
        except ConfigError as e:
            return Result(self, False, failure_text=str(e))

        if self.require(out):
            self.print_formatted(f"{str(out):10} Ok")
            return Result(self, True)
        else:
            self.print_formatted(f"{str(out):10} Fail")
            return Result(self, False, failure_text="Too old.", found_version=out)


def parse_config(name: str, config: Mapping[str, str], templates):
    if "template" in config.keys():
        _config = {}
        if config["template"] not in templates.keys():
            raise ConfigError(
                f"Template {config['template']} not found. Templates: {list(templates.keys())}"
            )

        for k, v in templates[config["template"]].items():
            if isinstance(v, str):
                _config[k] = v.format(name=name)
            else:
                _config[k] = v
        _config.update(config)
    else:
        _config = dict(config)

    _deps = map(str.strip, _config.get("depends", "").split(","))
    deps = list(filter(lambda x: x != "", _deps))

    assert "require" in _config, "Every item needs a `require` field"
    assert "get_version" in _config, "Every item needs a `get_version` field"

    require, _ = parse_version_constraint(_config["require"])

    return VersionTest(
        name=name,
        require=require,
        get_version=_config["get_version"],
        platform=_config.get("platform", None),
        pattern=_config.get("pattern", None),
        suggestion_text=_config.get("suggestion_text", None),
        suggestion=_config.get("suggestion", None),
        depends=deps,
        template=_config.get("template", None),
    )


@contextmanager
def indent(prefix: str):
    f = io.StringIO()
    with redirect_stdout(f):
        yield
    output = f.getvalue()
    print(textwrap.indent(output, prefix), end="")


async def main():
    config = tomllib.load(open("dependencies.toml", "rb"))

    if "template" in config.keys():
        templates = config["template"]
        del config["template"]
    else:
        templates = dict()

    try:
        tests = {
            name: parse_config(name, config[name], templates)
            for name in config
            if ":" not in name and name != "DEFAULT"
        }
    except (AssertionError, ConfigError) as e:
        print("Configuration error:", e)
        sys.exit(1)

    async def test_version(name: str):
        assert name in tests, f"unknown dependency {name}"
        x = await tests[name].run(test_version)
        return x

    result = await asyncio.gather(*(test_version(k) for k in tests))
    if all(r.success for r in result):
        print("Success")
        sys.exit(0)
    else:
        print("Failure")
        with indent("  |  "):
            for r in (r for r in result if not r.success):
                if r.failure_text:
                    print(f"{r.test.name}: {r.failure_text}")
                if r.found_version:
                    print(f"    found version {r.found_version}")
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())

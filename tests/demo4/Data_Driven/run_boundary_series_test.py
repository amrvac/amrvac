"""Generate three tiny boundary frames for test_boundary_series.f90."""

import struct
import sys
from pathlib import Path


target = Path(sys.argv[1])
target.mkdir(parents=True, exist_ok=True)
header = struct.Struct("=diidd")
for index, (time_value, field_value) in enumerate(
    ((0.0, 0.0), (10.0, 10.0), (20.0, 20.0)), 1
):
    path = target / "B_{:04d}.dat".format(index)
    with path.open("wb") as handle:
        handle.write(header.pack(time_value, 2, 2, 100.0, 100.0))
        handle.write(struct.pack("=12d", *([field_value] * 12)))

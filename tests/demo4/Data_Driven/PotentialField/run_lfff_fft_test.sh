#!/usr/bin/env bash
set -eu

if [ -z "${AMRVAC_DIR:-}" ]; then
  echo "AMRVAC_DIR is not set" >&2
  exit 2
fi

lfff_test_dir=$(mktemp -d)
trap 'rm -rf "${lfff_test_dir}"' EXIT
lfff_fc=${MPIFC:-mpif90}
lfff_lib_dir=${AMRVAC_DIR}/lib/3d3_default

if [ ! -f "${lfff_lib_dir}/mod_lfff.mod" ]; then
  echo "Build a 3D AMRVAC configuration before running this test" >&2
  exit 2
fi

"${lfff_fc}" -ffree-form -ffree-line-length-none \
  "${AMRVAC_DIR}/tests/demo4/Data_Driven/PotentialField/test_lfff_fft_modes.f90" \
  -I "${lfff_lib_dir}" -L "${lfff_lib_dir}" -lamrvac \
  -o "${lfff_test_dir}/test_lfff_fft_modes"
"${lfff_test_dir}/test_lfff_fft_modes"

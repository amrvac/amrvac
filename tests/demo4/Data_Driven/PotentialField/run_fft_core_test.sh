#!/usr/bin/env bash
set -eu

if [ -z "${AMRVAC_DIR:-}" ]; then
  echo "AMRVAC_DIR is not set" >&2
  exit 2
fi

fft_test_dir=$(mktemp -d)
trap 'rm -rf "${fft_test_dir}"' EXIT
fft_fc=${FC:-gfortran}

"${fft_fc}" -c -ffree-form -ffree-line-length-none -x f95 \
  "${AMRVAC_DIR}/src/physics/mod_fft.t" \
  -J "${fft_test_dir}" -o "${fft_test_dir}/mod_fft.o"
"${fft_fc}" -ffree-form -ffree-line-length-none \
  "${AMRVAC_DIR}/tests/demo4/Data_Driven/PotentialField/test_fft_core.f90" \
  "${fft_test_dir}/mod_fft.o" -I "${fft_test_dir}" \
  -o "${fft_test_dir}/test_fft_core"
"${fft_test_dir}/test_fft_core"

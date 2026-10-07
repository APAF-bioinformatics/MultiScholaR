#!/bin/sh
set -eu
source_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
output_dir=$1
mkdir -p "$output_dir"
cmake -S "$source_dir" -B "$output_dir/cmake" -DCMAKE_BUILD_TYPE=Release
cmake --build "$output_dir/cmake" --parallel 8
python3 "$source_dir/extension-ci-tools/scripts/append_extension_metadata.py" --library-file "$output_dir/cmake/libdia_rsum.so" --out-file "$output_dir/dia_rsum.duckdb_extension" --extension-name dia_rsum --duckdb-platform linux_amd64 --duckdb-version v1.2.0 --extension-version dia025-v1 --abi-type C_STRUCT

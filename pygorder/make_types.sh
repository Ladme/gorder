#!/bin/bash

cargo run --bin stub_gen

set -euo pipefail

sed -i '0,/import typing/{s/import typing/import typing\nimport gorder/}' python/gorder/__init__.pyi
sed -i '0,/from \. import results/{s/from \. import results/from . import results\nfrom . import exceptions/}' python/gorder/__init__.pyi
sed -i '0,/import typing/{s/import typing/import typing\nimport gorder/}' python/gorder/leaflets.pyi
sed -i '0,/import typing/{s/import typing/import typing\nimport gorder/}' python/gorder/results.pyi
sed -i '0,/import typing/{s/import typing/import typing\nfrom math import inf/}' python/gorder/geometry.pyi
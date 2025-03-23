# Python bindings for `gorder`

This directory contains python bindings for the `gorder` tool. See the [gorder manual](https://ladme.github.io/gorder-manual/) for more information.

## For developers

### Running tests

- Normally: `uv run pytest`
- If you change the core `gorder` code or anything inside `src`: `uv run --reinstall pytest`

### Generating documentation

- Activate a conda environment.
- Run maturin: `maturin develop --release`.
- Go to `docs` and run `make`: `cd docs && make clean html`.
- Open `build/index.html`.
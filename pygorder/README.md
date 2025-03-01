## Running tests

- normally: `uv run pytest`
- if you change the core `gorder` code or anything inside `src`: `uv run --reinstall pytest`

## Generating documentation

- Activate a conda environment.
- Run maturin: `maturin develop --release`.
- Go to `docs` and run `make`: `cd docs && make clean html`.
- Open `build/index.html`.
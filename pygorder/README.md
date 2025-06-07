# Python bindings for `gorder`

This project contains Python bindings for the `gorder` tool for calculating lipid order parameters from Gromacs simulations.

**See the [gorder manual](https://ladme.github.io/gorder-manual/) and [Python API documentation](https://ladme.github.io/pygorder-docs/) for more information.**

## Installation

To install `gorder` as a Python package, you can use `pip`:

```bash
pip install gorder
```

But you **should** use `uv` instead. Here is how you add `gorder` to your `uv` project:

```bash
uv add gorder
```

## Example

The following script calculates coarse-grained lipid order parameters using a TPR file `system.tpr` and a trajectory file `md.xtc`. The results are written into `order.yaml`.

```python
import gorder

analysis = gorder.Analysis(
    structure = "system.tpr",
    trajectory = "md.xtc",
    analysis_type = gorder.analysis_types.CGOrder("@membrane"),
    output_yaml = "order.yaml",
)

results = analysis.run()

results.write()
```
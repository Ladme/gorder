# gorder: Everything you will ever need for lipid order calculations

Command line tool for calculating atomistic or coarse-grained lipid order parameters from Gromacs simulations.

⚠️ `gorder` is currently in BETA. Therefore you will not be able to install it using the instructions below. If you want to try the program, use `cargo install gorder --version 0.2.0-beta`, but it might be better to wait for a bit more stable and validated version. ⚠️

## Installation

1. [Install Rust](https://www.rust-lang.org/tools/install)

2. Install gorder:
```bash
$ cargo install gorder
```

## Usage

1) Prepare a yaml file describing the parameters of your analysis. It's as simple as:

```yaml
structure: system.tpr
trajectory: md.xtc
analysis_type: !AAOrder
    heavy_atoms: "@membrane and element name carbon"
    hydrogens: "@membrane and element name hydrogen"
output: order.yaml
```

2) Run the `gorder` tool. The program will automatically find bonds between the selected atoms, classify lipid molecules, and calculate order parameters for them.

```
$ gorder YOUR_INPUT_YAML_FILE
```

3) See the results of the analysis in `order.yaml`.

**See the [gorder manual](https://ladme.github.io/gorder-manual/) for more information about using this program.**

## Features
- **Atomistic and coarse-grained systems.** `gorder` is able to calculate atomistic and coarse-grained order parameters for individual bonds of individual lipid types.
- **Powerful selection language.** `gorder` allows for simple yet powerful atom selection using a VMD-like selection language, supporting regular expressions and groups from NDX files.
- **Largely automatic.** `gorder` automatically recognizes bonds and classifies molecule types based on their topology.
- **Various output formats.** `gorder` can output results in YAML, XVG, CSV, and custom "table" format.
- **Leaflet-wise analysis.** `gorder` can perform scrambling-safe assignment of lipids to membrane leaflets using three different methods, and then calculate lipid order parameters for individual leaflets.
- **Order parameter maps.** `gorder` can construct 2D maps of order parameters, so you know what parts of the membrane are ordered and disordered.
- **Supports any force-field.** `gorder` is completely force-field agnostic. Martini? CHARMM? Slipids? Your own toy force-field? As long as your lipids have bonds, it will work.
- **Very fast with multithreading.** `gorder` is very fast (see below) through multithreading support.

## Planned
- [ ] Robust error estimation using cross-validation.
- [ ] Dynamic membrane normal calculation, supporting membrane vesicles.
- [ ] Improved multithreading (currently, multithreading is only implemented at the trajectory reading level).
- [ ] United-atom order parameters.
- [ ] Python API: using `gorder` as a Python library.

## Validation

### Atomistic order parameters
🏗️ Work in progress... 🏗️

### Coarse-grained order parameters
A Martini 3 simulation of a membrane consisting of 512 POPC lipids was used to validate the `gorder` program. In total, the system contained ~16,800 beads. The trajectory had a length of 1 μs and consisted of 10,000 frames. The following programs/libraries were used for validation:
- [Martini's `do-order` script](https://cgmartini.nl/docs/downloads/tools/other-tools.html#do-order)
- [`lipyphilic` library](https://lipyphilic.readthedocs.io/en/stable/index.html)
- [Ladme's `order` program](https://doi.org/10.5281/zenodo.8369479)

![Four completely overlapping curves, appearing as a single curve.](validation/cgorder_validation.png)

*Yes, all programs returned the same results. Minor variations due to the employed calculation approaches are too small to be visible in the chart.*

## Benchmarking
> **tl;dr `gorder` is fast. Like *really* fast.**

### Coarse-grained order parameters
Run time of the analyses performed in the Validation section by various programs:
![Bar chart showing the run time of various programs.](validation/cgorder_benchmark.png)

*Benchmarks were run on 8-core Intel Core i7-10700 CPU. Benchmarking of `order` and `gorder` was performed using [`hyperfine`](https://github.com/sharkdp/hyperfine). `do-order` and `lipyphilic` were really slow, so only an approximate value obtained using GNU's `time` is reported.*

<sup>a</sup> Note that the `do-order` script is not able to calculate order parameters for individual leaflets in the same run. In contrast, both `order` and `gorder` were run with this capability enabled.

<sup>b</sup> Note that the `lipyphilic` library is not able to calculate order parameters for the individual bonds simultaneously, requiring the analysis to be run multiple times. If you are only interested in average order parameters for one entire lipid tail, `lipyphilic` is quite fast (though still slower than `gorder`). Additionally, in this run, `lipyphilic` calculated order parameters only for the entire membrane, not for individual leaflets. In contrast, both `order` and `gorder` were run with this capability enabled.

## Rust API
`gorder` is also available as a Rust crate. See the [gorder manual](https://ladme.github.io/gorder-manual/) and the corresponding [docs.rs](https://docs.rs/gorder/latest/gorder) page for more information.

## License
The command line tool and the crate are both released under the MIT License.

## Limitations
Only works for simulations with orthogonal simulation boxes.
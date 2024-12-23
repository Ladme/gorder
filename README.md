# gorder: Everything you will ever need for lipid order calculations
[![DOI](https://zenodo.org/badge/848864524.svg)](https://doi.org/10.5281/zenodo.14391305) ![Crates.io Version](https://img.shields.io/crates/v/gorder)![Crates.io Downloads](https://img.shields.io/crates/d/gorder)![Crates.io License](https://img.shields.io/crates/l/gorder)

Command line tool for calculating atomistic and coarse-grained lipid order parameters from Gromacs simulations.

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
trajectory: md.xtc     # use your MD trajectory directly - no PBC handling or molecule fixing needed
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
- [ ] Dynamic selection of lipids for order parameter calculation based on geometric conditions (i.e., only calculating order parameters from a part of a membrane).
- [ ] Python API: using `gorder` as a Python library.
- [ ] United-atom order parameters.
- [ ] Improved multithreading (currently, multithreading is only implemented at the trajectory reading level).

## Validation

### Atomistic order parameters
A CHARMM36m simulation of a membrane consisting of 256 lipids was used to validate the calculation of atomistic order parameters by the `gorder` program. In total, the system contained ~64,500 atoms. The trajectory was 200 ns long and consisted of 10,000 frames. The following programs were used for validation:
- [`OrderParameter.py` by NMR Lipids](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/OrderParameter.py)
- [VMD's `calc_op.tcl` script](https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/att-14731/calc_op.tcl)

We also present a comparison of the calculated order parameters with the results from [`gmx order`](https://manual.gromacs.org/2021.4/onlinehelp/gmx-order.html) (version 2021.4). Note that `gmx order` actually calculates united atom order parameters, so it is not suitable for atomistic systems. (However, many users still utilize it.)

![Three completely overlapping curves, appearing as a single curve, and an additional dashed curve corresponding to `gmx order`.](validation/aaorder_validation.png)

*All programs, except for `gmx order`, produce the same results. Minor variations due to different calculation approaches are too small to be visible in the chart. `gmx order` returns slightly different values (and no value for carbon #16) because it calculates united atom, not atomistic, order parameters.*

### Coarse-grained order parameters
A Martini 3 simulation of a membrane consisting of 512 POPC lipids was used to validate the calculation of CG order parameters by the `gorder` program. In total, the system contained ~16,800 beads. The trajectory had a length of 1 μs and consisted of 10,000 frames. The following programs/libraries were used for validation:
- [Martini's `do-order` script](https://cgmartini.nl/docs/downloads/tools/other-tools.html#do-order)
- [`lipyphilic` library](https://lipyphilic.readthedocs.io/en/stable/index.html)
- [Ladme's `order` program](https://doi.org/10.5281/zenodo.8369479)

![Four completely overlapping curves, appearing as a single curve.](validation/cgorder_validation.png)

*All programs produce the same results. Minor variations due to the employed calculation approaches are too small to be visible in the chart.*

## Benchmarking
### Atomistic order parameters
Run time of the analyses performed in the Validation section by various programs:
![Bar chart showing the run time of various programs.](validation/aaorder_benchmark.png)

*Benchmarks were conducted on Debian 12 with an 8-core Intel Core i7-10700 CPU. Benchmarking of `gmx order` and `gorder` was performed using [`hyperfine`](https://github.com/sharkdp/hyperfine). The `NMR lipids` script and `calc_op.tcl` were dramatically slower, so only approximate values obtained using GNU's `time` are reported.*

<sup>a</sup> Note that, unlike `calc_op.tcl` and `gmx order`, the NMR lipids script provides information about the order of individual C-H bonds. `gorder` also provides this information.

<sup>b</sup> Note that `gmx order` calculates united atom order parameters, not atomistic order parameters. For saturated tails, this calculation is reasonably accurate; however, for unsaturated tails, it is **very** inaccurate [[1]](https://doi.org/10.1021/acs.jctc.7b00643). Additionally, it is slower and more tedious to use than `gorder`, so there is little justification for using it for atomistic systems.

### Coarse-grained order parameters
Run time of the analyses performed in the Validation section by various programs:
![Bar chart showing the run time of various programs.](validation/cgorder_benchmark.png)

*Benchmarks were conducted on Debian 12 with with an 8-core Intel Core i7-10700 CPU. Benchmarking of `order` and `gorder` was performed using [`hyperfine`](https://github.com/sharkdp/hyperfine). `do-order` and `lipyphilic` were dramatically slower, so only an approximate value obtained using GNU's `time` is reported.*

<sup>a</sup> Note that the `do-order` script is not able to calculate order parameters for individual leaflets in the same run. In contrast, both `order` and `gorder` were run with this capability enabled.

<sup>b</sup> Note that the `lipyphilic` library is not able to calculate order parameters for the individual bonds simultaneously, requiring the analysis to be run multiple times. If you are only interested in average order parameters for one entire lipid tail, `lipyphilic` is quite fast (though still slower than `gorder`). Additionally, in this run, `lipyphilic` calculated order parameters only for the entire membrane, not for individual leaflets. In contrast, both `order` and `gorder` were run with this capability enabled.

## Rust API
`gorder` is also available as a Rust crate. See the [gorder manual](https://ladme.github.io/gorder-manual/) and the corresponding [docs.rs](https://docs.rs/gorder/latest/gorder) page for more information.

## License
The command line tool and the crate are both released under the MIT License.

## Limitations
- `gorder` only works for simulations with orthogonal simulation boxes.
- `gorder` only supports TPR files generated with Gromacs 5.1 or newer.
- `gorder` is developed on Linux for Linux. While it should work on other operating systems, it is not guaranteed.
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
- Calculation of atomistic and coarse-grained order parameters for individual lipid types.
- Simple atom selection using a VMD-like selection language, supporting regular expressions and groups from NDX files.
- Automatic identification of bonds and topology-based classification of molecule types.
- Output in YAML, XVG, CSV, and custom "table" format.
- Scrambling-safe assignment of lipids to membrane leaflets using three different methods, and calculation of lipid order parameters for individual leaflets.
- Calculation of order parameter maps.
- Completely force-field agnostic.
- Robust: use any XTC trajectory even with molecules broken at periodic boundaries.
- Fast, multi-threaded analysis.

## Planned
- [ ] Robust error estimation using cross-validation.
- [ ] Dynamic membrane normal calculation, supporting membrane vesicles.
- [ ] Improved multithreading (currently, multithreading is only implemented at the trajectory reading level).
- [ ] United-atom order parameters.
- [ ] Python API: using `gorder` as a Python library.

## Validation

🏗️ Work in progress... 🏗️

## Rust API
`gorder` is also available as a Rust crate. See the [gorder manual](https://ladme.github.io/gorder-manual/) and the corresponding [docs.rs](https://docs.rs/gorder/latest/gorder) page for more information.

## License
The command line tool and the crate are both released under the MIT License.

## Limitations
Only works for simulations with orthogonal simulation boxes.
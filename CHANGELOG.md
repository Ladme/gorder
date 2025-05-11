## Version 0.7.2
- Changes to the Rust API: Made all user-provided parameters of all methods for leaflet assignment public.

## Version 0.7.1
- Fixed undefined behavior in the `groan_rs` library's C source code that caused TRR files to be read incorrectly when compiled with certain C compilers (e.g., clang).

## Version 0.7.0
- **Leaflet classification for curved membranes:** `gorder` can now use spectral clustering to classify lipids into leaflets for any membrane geometry, including buckled membranes and vesicles. See the [manual](https://ladme.github.io/gorder-manual/leaflets.html#clustering-method-for-leaflet-classification) for more information.
- **!!BREAKING CHANGE!!** **Removed support for some trajectory formats:** Removed support for PDB, Amber NetCDF, DCD, and LAMMPSTRJ trajectories. (PDB structure files are still supported.) The parsing of these trajectories relied on the `chemfiles` library, whose maintenance standards and capabilities just do not meet the quality requirements for `gorder`. If you need to use these trajectory formats, you can still use `gorder v0.6`.

***

## Version 0.6.0
- **United-atom order parameters:** `gorder` is now able to calculate order parameters in united-atom systems. See the [manual](https://ladme.github.io/gorder-manual/uaorder_basics.html) for more information.
- **Python API:** `gorder` is now available as a Python package. See the [manual](https://ladme.github.io/gorder-manual/python_api.html) for more information.
- **Trajectory concatenation:** You can now provide multiple trajectory files which will be all joined into one trajectory and analyzed. In case there are duplicate frames at trajectory boundaries, `gorder` will analyze only one of the duplicate frames. This feature is currently only supported for XTC and TRR files. See the [manual](https://ladme.github.io/gorder-manual/multiple_trajectories.html) for more information.
- **Manual membrane normals:** Membrane normals can be now assigned manually for each lipid molecule in each trajectory frame. See the [manual](https://ladme.github.io/gorder-manual/manual_normals.html) for more information.
- **Optimizations:** Local leaflet classification method and dynamic local membrane normal calculation are now much faster, especially for very large systems, through using cell lists. **Warning:** The changes may cause small differences in the calculated order parameters compared to version 0.5, especially for very short trajectories. These differences should be on the order of 0.0001 arb. u. or lower.

***

## Version 0.5.0
- **Dynamic membrane normal calculation:** Membrane normals can be now calculated dynamically from actual membrane shape which allows the calculation of order parameters for vesicles and similar systems (see the [manual](https://ladme.github.io/gorder-manual/membrane_normal.html)).
- **Ignoring PBC:** You can now choose to ignore periodic boundary conditions. This allows analyzing simulations with non-orthogonal simulation boxes with some small additional friction (making molecules whole). See the [manual](https://ladme.github.io/gorder-manual/no_pbc.html) for more information.
- **Reworked manual lipid assignment:** It is now possible to classify lipids into membrane leaflets using NDX files, enabling integration of `gorder` with `FATSLiM`. **BREAKING CHANGE** The keywords to request manual leaflet assignment using a leaflet assignment file have been changed, see the [manual](https://ladme.github.io/gorder-manual/manual_leaflets.html#assigning-lipids-using-a-leaflet-assignment-file) for more information.
- **Ordermaps visualization:** A python script is now generated inside any created `ordermaps` directory which can be used to easily plot the ordermaps. Changed the default range of the colorbar in ordermaps to more reasonable values.
- **More trajectory formats:** Added **experimental** support for more trajectory formats, namely TRR, GRO, PDB, Amber NetCDF, DCD, and LAMMPSTRJ. Always prefer using XTC trajectory as `gorder` is optimized to read it very fast. There are also some limitations connected with using different trajectory formats, see the [manual](https://ladme.github.io/gorder-manual/other_input.html#trajectory-file-formats).
- **Bug fixes and other changes:**
  - `gorder` now returns an error if the center of geometry calculation for leaflet classification is nonsensical (i.e., `nan`).
  - If molecule classification runs longer than expected, progress is logged. By default, progress output begins after 500 ms, but you can adjust this delay using the environment variable `GORDER_MOLECULE_CLASSIFICATION_TIME_LIMIT`.
  - Added more color to information written to standard output during analysis. Changed logging for output file writing. Changed logging for molecule types.
  - Reading an ndx file with invalid or duplicated group names no longer results in a hard error but instead raises a warning.

***

## Version 0.4.0
- **Geometry selection:** Added the ability to select a geometric region for analysis. Users can now specify cuboidal, spherical, or cylindrical regions, and order parameters will be calculated only for bonds located within the selected region.
- **Support for reading GRO, PDB, and PQR files:** These file formats are now supported as input structure files. In some cases, an additional "bonds" file specifying the system's connectivity may be required. Refer to the manual for more details.
- **Manual assignment of lipids to leaflets:** Lipids can now be manually assigned to leaflets using a provided leaflet assignment file. Refer to the manual for detailed instructions.
- **Calculating average results for the entire system:** YAML and TAB files now include information about the average order parameters calculated across all bonds and molecule types in the system. Additionally, ordermaps are generated for the entire system.

***

## Version 0.3.0
- **Error estimation and convergence analysis**: Implemented error estimation and convergence analysis. Refer to the corresponding section of the manual for more details.
- **Leaflet classification**: Leaflet classification can now also be performed either every N analyzed trajectory frames or only once at the start of the analysis.
- **Improved XTC file reading**: Switched to using the `molly` crate for reading XTC files. This allows reading only the parts of the XTC files that are needed, making `gorder` more than twice as fast compared to version 0.2.
- **Reworked Rust API**: The Rust API has been restructured, enabling access to results without the need to write output files.
- **Enhanced output files**:
  - YAML and TAB files now display the average order parameter calculated from all bonds of a single molecule type.
  - An average ordermap, summarizing order parameters collected from all bonds of a single molecule type, is now automatically generated during ordermap analysis.
- **YAML format updates**: The YAML format has been revised: molecule types are no longer stored as a list but are instead represented as a dictionary.
- **Export configuration**: Analysis parameters can now be exported into a YAML file using the `--export-config` argument.
- **Updated heavy atom order parameter calculation**: Adjusted the calculation of average order parameters for heavy atoms. This may result in order parameters being available for atoms even if all bond order parameters are reported as NaN. (And that's okay.)
- **Bug fixes**:
  - Fixed a minor issue causing occasional rounding discrepancies between YAML and other output formats, which could lead to slightly different results.
  - Corrected CG bond numbering in XVG files, ensuring they are properly numbered starting from 1.
- **Improved error messages**: Certain error messages have been clarified to enhance their readability and comprehension.

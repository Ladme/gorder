## Version 0.4.0
- **Geometry selection:** Implemented geometry selection. A cuboidal, spherical, or cylindrical region can be now specified. Order parameters will be then calculated only using bonds located in this region.
- **Support for reading GRO, PDB, and PQR files:** These file formats are now supported as input structure files. In some cases, an additional "bonds" file specifying the connectivity in the system must be supported. See the manual for more information.

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

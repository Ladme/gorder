## Version 0.3.0
- Implemented error estimation and convergence analysis (see the corresponding section of the manual for details).
- Leaflet classification can now also be performed every N analyzed trajectory frames or only once at that start of the analysis.
- Reworked the Rust API so you can now access the results without having to write output files.
- YAML and TAB files now display average order parameter calculated from all bonds of a single molecule type.
- Average ordermap displaying order parameters collected from all bonds of a single molecule type is now automatically written when performing ordermaps analysis.
- YAML format has been reworked, molecule types are no longer stored in a list but are instead stored in a dictionary.
- Changed how average order parameters for heavy atoms are calculated leading to situations where you may get order parameter for an atom even if all order parameters for its bonds are reported as NaN.
- Fixed a small bug which caused rounding errors to be different between YAML and the other formats occasionally leading to slightly different results.
- CG bonds in XVG files are now correctly numbered starting from 1.
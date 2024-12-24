## Version 0.3.0
- YAML and TAB files now display average order parameter calculated from all bonds of a single molecule type.
- Average ordermap displaying order parameters collected from all bonds of a single molecule type is now automatically written when performing ordermaps analysis.
- YAML format has been reworked, molecule types are no longer stored in a list but are instead stored in a dictionary.
- Changed how average order parameters for heavy atoms are calculated leading to situations where you may get order parameter for an atom even if all order parameters for its bonds are reported as NaN.
- Fixed a small bug which caused rounding errors to be different between YAML and the other formats occasionally leading to slightly different results.
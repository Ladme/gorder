// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::io::Write;
use std::{fs::File, io::BufWriter, path::Path};

use crate::analysis::molecule::MoleculeType;
use crate::analysis::topology::SystemTopology;
use crate::{
    analysis::{molecule::AtomType, ordermap::Map},
    errors::OrderMapWriteError,
    Leaflet,
};
use crate::{GORDER_VERSION, PANIC_MESSAGE};

impl Map {
    /// Write the map of order parameters for a single heavy atom into an output file.
    /// Leaflet `None` corresponds to ordermap for the full membrane.
    pub(crate) fn write_atom_map(
        &self,
        atom: &AtomType,
        leaflet: Option<Leaflet>,
    ) -> Result<(), OrderMapWriteError> {
        let filename = match leaflet {
            Some(Leaflet::Upper) => format!("ordermap_{}_upper.dat", atom),
            Some(Leaflet::Lower) => format!("ordermap_{}_lower.dat", atom),
            None => format!("ordermap_{}_total.dat", atom),
        };

        let comment = format!("# Map of average order parameters calculated for bonds involving atom type {}.\n# Calculated with 'gorder v{}'.", atom, GORDER_VERSION);

        self.write_ordermap(&filename, &comment)
    }

    /// Write the map of order parameters for a single bond into an output file.
    /// Leaflet `None` corresponds to ordermap for the full membrane.
    pub(crate) fn write_bond_map(
        &self,
        atom1: &AtomType,
        atom2: &AtomType,
        leaflet: Option<Leaflet>,
    ) -> Result<(), OrderMapWriteError> {
        let filename = match leaflet {
            Some(Leaflet::Upper) => format!("ordermap_{}--{}_upper.dat", atom1, atom2),
            Some(Leaflet::Lower) => format!("ordermap_{}--{}_lower.dat", atom1, atom2),
            None => format!("ordermap_{}--{}_total.dat", atom1, atom2),
        };

        let comment = format!("# Map of average order parameters calculated for bonds between atom types {} and {}.\n# Calculated with 'gorder v{}'.", atom1, atom2, GORDER_VERSION);

        self.write_ordermap(&filename, &comment)
    }

    fn write_ordermap(&self, filename: &str, comment: &str) -> Result<(), OrderMapWriteError> {
        // directory must already exist
        let directory = self.params().output_directory();

        let full_path = format!("{}/{}", directory, filename);
        let output_file = File::create(&full_path).map_err(|_| {
            OrderMapWriteError::CouldNotCreateFile(Box::from(Path::new(&full_path)))
        })?;
        let mut output = BufWriter::new(output_file);

        writeln!(output, "{}", comment)
            .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(Path::new(&full_path))))?;

        writeln!(
            output,
            "@ xlabel x-dimension [nm]\n@ ylabel y-dimension [nm]\n@ zrange -1 1 0.2"
        )
        .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(Path::new(&full_path))))?;

        writeln!(output, "$ type colorbar\n$ colormap seismic_r")
            .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(Path::new(&full_path))))?;

        // write the values
        for (value, samples) in self
            .values()
            .extract_convert()
            .zip(self.samples().extract_convert())
        {
            let average = if samples.2 < self.params().min_samples() {
                f32::NAN
            } else {
                value.2 / samples.2 as f32
            };

            writeln!(output, "{:.4} {:.4} {:.4}", value.0, value.1, average).map_err(|_| {
                OrderMapWriteError::CouldNotWriteLine(Box::from(Path::new(&full_path)))
            })?;
        }

        Ok(())
    }
}

impl MoleculeType {
    /// Write ordermaps constructed for all the bond types of this molecule type.
    fn write_ordermaps_bonds(&self) -> Result<(), OrderMapWriteError> {
        for bond in self.order_bonds().bond_types() {
            if let Some(map) = bond.total_map() {
                map.write_bond_map(bond.atom1(), bond.atom2(), None)?;
            }

            if let Some(map) = bond.upper_map() {
                map.write_bond_map(bond.atom1(), bond.atom2(), Some(Leaflet::Upper))?;
            }

            if let Some(map) = bond.lower_map() {
                map.write_bond_map(bond.atom1(), bond.atom2(), Some(Leaflet::Lower))?;
            }
        }

        Ok(())
    }

    fn write_ordermaps_atoms(&self) -> Result<(), OrderMapWriteError> {
        for heavy_atom in self.order_atoms().atoms() {
            let mut relevant_maps = Vec::new();

            for bond in self.order_bonds().bond_types() {
                if bond.contains(heavy_atom) {
                    relevant_maps.push((
                        bond.total_map().clone(),
                        bond.upper_map().clone(),
                        bond.lower_map().clone(),
                    ));
                }
            }

            let (total_total, total_upper, total_lower) = MoleculeType::merge_maps(relevant_maps);

            if let Some(map) = total_total {
                map.write_atom_map(heavy_atom, None)?;
            }

            if let Some(map) = total_upper {
                map.write_atom_map(heavy_atom, Some(Leaflet::Upper))?;
            }

            if let Some(map) = total_lower {
                map.write_atom_map(heavy_atom, Some(Leaflet::Lower))?;
            }
        }

        Ok(())
    }

    /// Merge order maps.
    fn merge_maps(
        maps: Vec<(Option<Map>, Option<Map>, Option<Map>)>,
    ) -> (Option<Map>, Option<Map>, Option<Map>) {
        let (mut total_total, mut total_upper, mut total_lower) = (None, None, None);

        for (total, upper, lower) in maps.into_iter() {
            if let Some(map) = total {
                total_total = Some(match total_total {
                    Some(acc) => acc + map,
                    None => map,
                });
            }
            if let Some(map) = upper {
                total_upper = Some(match total_upper {
                    Some(acc) => acc + map,
                    None => map,
                });
            }
            if let Some(map) = lower {
                total_lower = Some(match total_lower {
                    Some(acc) => acc + map,
                    None => map,
                });
            }
        }

        (total_total, total_upper, total_lower)
    }
}

impl SystemTopology {
    /// Write all ordermaps consturcted for all bonds of this topology.
    pub(crate) fn write_ordermaps_bonds(&self) -> Result<(), OrderMapWriteError> {
        for molecule in self.molecule_types() {
            molecule.write_ordermaps_bonds()?;
        }

        Ok(())
    }

    pub(crate) fn write_ordermaps_atoms(&self) -> Result<(), OrderMapWriteError> {
        for molecule in self.molecule_types() {
            molecule.write_ordermaps_atoms()?;
        }

        Ok(())
    }

    /// Create or back up the directory for the ordermaps.
    pub(crate) fn handle_ordermap_directory(
        &self,
        overwrite: bool,
    ) -> Result<(), OrderMapWriteError> {
        if let Some(map_unwrapped) = self
            .molecule_types()
            .iter()
            .filter_map(|m| m.order_bonds().bond_types().get(0))
            .filter_map(|bond| bond.total_map().as_ref())
            .next()
        {
            let directory = Path::new(map_unwrapped.params().output_directory());

            log::info!(
                "Writing ordermaps into a directory '{}'...",
                directory.to_str().expect(PANIC_MESSAGE)
            );

            log::logger().flush();

            // back up the directory if it already exists
            if directory.is_dir() {
                if !overwrite {
                    log::warn!(
                        "Output directory for ordermaps '{}' already exists. Backing it up.",
                        directory.to_str().expect(PANIC_MESSAGE)
                    );
                    backitup::backup(directory).map_err(|_| {
                        OrderMapWriteError::CouldNotBackupDirectory(directory.into())
                    })?;
                } else {
                    log::warn!("Output directory for ordermaps '{}' already exists. It will be overwritten as requested.", 
                        directory.to_str().expect(PANIC_MESSAGE)
                    );
                    std::fs::remove_dir_all(directory).map_err(|_| {
                        OrderMapWriteError::CouldNotRemoveDirectory(directory.into())
                    })?;
                }
            }

            // create a new directory
            std::fs::create_dir_all(directory)
                .map_err(|_| OrderMapWriteError::CouldNotCreateDirectory(directory.into()))?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use groan_rs::prelude::SimBox;
    use tempfile::TempDir;

    use crate::OrderMap;

    use super::*;

    #[test]
    fn test_write_map() {
        let directory = TempDir::new().unwrap();
        let path = directory.path().to_str().unwrap();

        let params = OrderMap::new()
            .output_directory(path)
            .bin_size_x(1.0)
            .bin_size_y(1.0)
            .min_samples(10)
            .build()
            .unwrap();

        let mut map = Map::new(params, &SimBox::from([3.0, 3.0, 10.0]));

        *(map.values_mut().get_mut_at(0.0, 0.0).unwrap()) = 2.618;
        *(map.samples_mut().get_mut_at(0.0, 0.0).unwrap()) = 17;

        *(map.values_mut().get_mut_at(0.0, 1.0).unwrap()) = -2.262;
        *(map.samples_mut().get_mut_at(0.0, 1.0).unwrap()) = 10;

        *(map.values_mut().get_mut_at(1.0, 0.0).unwrap()) = 0.364;
        *(map.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 2;

        *(map.values_mut().get_mut_at(2.0, 2.0).unwrap()) = 0.764;
        *(map.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 9;

        *(map.values_mut().get_mut_at(3.0, 1.0).unwrap()) = 8.826;
        *(map.samples_mut().get_mut_at(3.0, 1.0).unwrap()) = 15;

        let atom1 = AtomType::new_raw(4, "POPC", "C22");
        let atom2 = AtomType::new_raw(6, "POPC", "H22");

        map.write_bond_map(&atom1, &atom2, None).unwrap();
        map.write_bond_map(&atom1, &atom2, Some(Leaflet::Upper))
            .unwrap();
        map.write_bond_map(&atom1, &atom2, Some(Leaflet::Lower))
            .unwrap();

        map.write_atom_map(&atom1, None).unwrap();
        map.write_atom_map(&atom1, Some(Leaflet::Upper)).unwrap();
        map.write_atom_map(&atom1, Some(Leaflet::Lower)).unwrap();

        let results_bonds = [
            "ordermap_POPC-C22-4--POPC-H22-6_total.dat",
            "ordermap_POPC-C22-4--POPC-H22-6_upper.dat",
            "ordermap_POPC-C22-4--POPC-H22-6_lower.dat",
        ];

        let results_atom = [
            "ordermap_POPC-C22-4_total.dat",
            "ordermap_POPC-C22-4_upper.dat",
            "ordermap_POPC-C22-4_lower.dat",
        ];

        for result_file in results_bonds.iter() {
            let full_path = format!("{}/{}", path, result_file);
            let mut result = File::open(&full_path).unwrap();
            let mut expected = File::open("tests/files/ordermap_bonds_expected.dat").unwrap();

            assert!(file_diff::diff_files(&mut result, &mut expected));
        }

        for result_file in results_atom.iter() {
            let full_path = format!("{}/{}", path, result_file);
            let mut result = File::open(&full_path).unwrap();
            let mut expected = File::open("tests/files/ordermap_atom_expected.dat").unwrap();

            assert!(file_diff::diff_files(&mut result, &mut expected));
        }
    }

    #[test]
    fn test_merge_maps() {
        let params = OrderMap::new()
            .output_directory(".")
            .bin_size_x(1.0)
            .bin_size_y(1.0)
            .min_samples(10)
            .build()
            .unwrap();

        let mut map1 = Map::new(params.clone(), &SimBox::from([3.0, 3.0, 10.0]));
        *(map1.values_mut().get_mut_at(0.0, 0.0).unwrap()) = 2.618;
        *(map1.samples_mut().get_mut_at(0.0, 0.0).unwrap()) = 17;

        *(map1.values_mut().get_mut_at(0.0, 1.0).unwrap()) = -2.262;
        *(map1.samples_mut().get_mut_at(0.0, 1.0).unwrap()) = 10;

        *(map1.values_mut().get_mut_at(1.0, 0.0).unwrap()) = 0.364;
        *(map1.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 2;

        *(map1.values_mut().get_mut_at(2.0, 2.0).unwrap()) = 0.764;
        *(map1.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 9;

        *(map1.values_mut().get_mut_at(3.0, 1.0).unwrap()) = 8.826;
        *(map1.samples_mut().get_mut_at(3.0, 1.0).unwrap()) = 15;

        let mut map2 = Map::new(params.clone(), &SimBox::from([3.0, 3.0, 10.0]));
        *(map2.values_mut().get_mut_at(0.0, 0.0).unwrap()) = 4.174;
        *(map2.samples_mut().get_mut_at(0.0, 0.0).unwrap()) = 35;

        *(map2.values_mut().get_mut_at(1.0, 1.0).unwrap()) = 1.654;
        *(map2.samples_mut().get_mut_at(1.0, 1.0).unwrap()) = 4;

        *(map2.values_mut().get_mut_at(1.0, 0.0).unwrap()) = 0.549;
        *(map2.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 5;

        *(map2.values_mut().get_mut_at(2.0, 2.0).unwrap()) = 1.677;
        *(map2.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 12;

        *(map2.values_mut().get_mut_at(3.0, 2.0).unwrap()) = -3.453;
        *(map2.samples_mut().get_mut_at(3.0, 2.0).unwrap()) = 15;

        let mut map3 = Map::new(params.clone(), &SimBox::from([3.0, 3.0, 10.0]));
        *(map3.values_mut().get_mut_at(1.0, 0.0).unwrap()) = 2.764;
        *(map3.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 14;

        *(map3.values_mut().get_mut_at(1.0, 1.0).unwrap()) = 0.874;
        *(map3.samples_mut().get_mut_at(1.0, 1.0).unwrap()) = 8;

        *(map3.values_mut().get_mut_at(2.0, 2.0).unwrap()) = -0.123;
        *(map3.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 41;

        *(map3.values_mut().get_mut_at(3.0, 0.0).unwrap()) = 0.434;
        *(map3.samples_mut().get_mut_at(3.0, 0.0).unwrap()) = 3;

        // total_a = map1 + map2 + map3
        // total_b = map3
        // total_c = map2 + map3

        let maps = vec![
            (Some(map1.clone()), None, Some(map2.clone())),
            (Some(map2.clone()), Some(map3.clone()), None),
            (Some(map3.clone()), None, Some(map3.clone())),
        ];

        let (total_a, total_b, total_c) = MoleculeType::merge_maps(maps);

        let values_a = total_a.as_ref().unwrap().values();
        let expected_values = [
            ((0.0, 0.0), 6.792),
            ((0.0, 1.0), -2.262),
            ((0.0, 2.0), 0.0),
            ((1.0, 0.0), 3.677),
            ((1.0, 1.0), 2.528),
            ((1.0, 2.0), 0.0),
            ((2.0, 0.0), 0.0),
            ((2.0, 1.0), 0.0),
            ((2.0, 2.0), 2.318),
            ((3.0, 0.0), 0.434),
            ((3.0, 1.0), 8.826),
            ((3.0, 2.0), -3.453),
        ];

        for &((x, y), expected) in &expected_values {
            assert_relative_eq!(*values_a.get_at(x, y).unwrap(), expected);
        }

        let samples_a = total_a.as_ref().unwrap().samples();
        let expected_samples = [
            ((0.0, 0.0), 52),
            ((0.0, 1.0), 10),
            ((0.0, 2.0), 0),
            ((1.0, 0.0), 21),
            ((1.0, 1.0), 12),
            ((1.0, 2.0), 0),
            ((2.0, 0.0), 0),
            ((2.0, 1.0), 0),
            ((2.0, 2.0), 62),
            ((3.0, 0.0), 3),
            ((3.0, 1.0), 15),
            ((3.0, 2.0), 15),
        ];

        for &((x, y), expected) in &expected_samples {
            assert_eq!(*samples_a.get_at(x, y).unwrap(), expected);
        }

        let values_b = total_b.as_ref().unwrap().values();
        let expected_values = [
            ((0.0, 0.0), 0.0),
            ((0.0, 1.0), 0.0),
            ((0.0, 2.0), 0.0),
            ((1.0, 0.0), 2.764),
            ((1.0, 1.0), 0.874),
            ((1.0, 2.0), 0.0),
            ((2.0, 0.0), 0.0),
            ((2.0, 1.0), 0.0),
            ((2.0, 2.0), -0.123),
            ((3.0, 0.0), 0.434),
            ((3.0, 1.0), 0.0),
            ((3.0, 2.0), 0.0),
        ];

        for &((x, y), expected) in &expected_values {
            assert_relative_eq!(*values_b.get_at(x, y).unwrap(), expected);
        }

        let samples_b = total_b.as_ref().unwrap().samples();
        let expected_samples = [
            ((0.0, 0.0), 0),
            ((0.0, 1.0), 0),
            ((0.0, 2.0), 0),
            ((1.0, 0.0), 14),
            ((1.0, 1.0), 8),
            ((1.0, 2.0), 0),
            ((2.0, 0.0), 0),
            ((2.0, 1.0), 0),
            ((2.0, 2.0), 41),
            ((3.0, 0.0), 3),
            ((3.0, 1.0), 0),
            ((3.0, 2.0), 0),
        ];

        for &((x, y), expected) in &expected_samples {
            assert_eq!(*samples_b.get_at(x, y).unwrap(), expected);
        }

        let values_c = total_c.as_ref().unwrap().values();
        let expected_values = [
            ((0.0, 0.0), 4.174),
            ((0.0, 1.0), 0.0),
            ((0.0, 2.0), 0.0),
            ((1.0, 0.0), 3.313),
            ((1.0, 1.0), 2.528),
            ((1.0, 2.0), 0.0),
            ((2.0, 0.0), 0.0),
            ((2.0, 1.0), 0.0),
            ((2.0, 2.0), 1.554),
            ((3.0, 0.0), 0.434),
            ((3.0, 1.0), 0.0),
            ((3.0, 2.0), -3.453),
        ];

        for &((x, y), expected) in &expected_values {
            assert_relative_eq!(*values_c.get_at(x, y).unwrap(), expected);
        }

        let samples_c = total_c.as_ref().unwrap().samples();
        let expected_samples = [
            ((0.0, 0.0), 35),
            ((0.0, 1.0), 0),
            ((0.0, 2.0), 0),
            ((1.0, 0.0), 19),
            ((1.0, 1.0), 12),
            ((1.0, 2.0), 0),
            ((2.0, 0.0), 0),
            ((2.0, 1.0), 0),
            ((2.0, 2.0), 53),
            ((3.0, 0.0), 3),
            ((3.0, 1.0), 0),
            ((3.0, 2.0), 15),
        ];

        for &((x, y), expected) in &expected_samples {
            assert_eq!(*samples_c.get_at(x, y).unwrap(), expected);
        }
    }

    #[test]
    fn test_merge_maps_none() {
        let (a, b, c) = MoleculeType::merge_maps(vec![
            (None, None, None),
            (None, None, None),
            (None, None, None),
        ]);

        assert!(a.is_none());
        assert!(b.is_none());
        assert!(c.is_none());
    }
}

// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Methods for writing maps of order parameters.

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

use super::OrderType;

impl Map {
    /// Write the map of order parameters for a single heavy atom type into an output file.
    /// Leaflet `None` corresponds to ordermap for the full membrane.
    #[inline]
    pub(crate) fn write_atom_map<O: OrderType>(
        &self,
        atom: &AtomType,
        leaflet: Option<Leaflet>,
        molname: &str,
    ) -> Result<(), OrderMapWriteError> {
        let filename = match leaflet {
            Some(Leaflet::Upper) => format!("ordermap_{}_upper.dat", atom),
            Some(Leaflet::Lower) => format!("ordermap_{}_lower.dat", atom),
            None => format!("ordermap_{}_full.dat", atom),
        };

        let comment = format!("# Map of average order parameters calculated for the atom type {}.\n# Calculated with 'gorder v{}'.", atom, GORDER_VERSION);

        self.write_ordermap::<O>(&filename, molname, &comment)
    }

    /// Write the map of order parameters for a single bond type into an output file.
    /// Leaflet `None` corresponds to ordermap for the full membrane.
    #[inline]
    pub(crate) fn write_bond_map<O: OrderType>(
        &self,
        atom1: &AtomType,
        atom2: &AtomType,
        leaflet: Option<Leaflet>,
        molname: &str,
    ) -> Result<(), OrderMapWriteError> {
        let filename = match leaflet {
            Some(Leaflet::Upper) => format!("ordermap_{}--{}_upper.dat", atom1, atom2),
            Some(Leaflet::Lower) => format!("ordermap_{}--{}_lower.dat", atom1, atom2),
            None => format!("ordermap_{}--{}_full.dat", atom1, atom2),
        };

        let comment = format!("# Map of average order parameters calculated for bonds between atom types {} and {}.\n# Calculated with 'gorder v{}'.", atom1, atom2, GORDER_VERSION);

        self.write_ordermap::<O>(&filename, molname, &comment)
    }

    /// Write the map of order parameters collected across all bonds.
    /// Leaflet `None` corresponds to ordermap for the full membarne.
    #[inline]
    fn write_average_map<O: OrderType>(
        &self,
        leaflet: Option<Leaflet>,
        molname: &str,
    ) -> Result<(), OrderMapWriteError> {
        let filename = match leaflet {
            Some(Leaflet::Upper) => "ordermap_average_upper.dat",
            Some(Leaflet::Lower) => "ordermap_average_lower.dat",
            None => "ordermap_average_full.dat",
        };

        let comment = format!(
            "# Map of average order parameters calculated for molecule type '{}'.\n# Calculated with 'gorder v{}'.",
            molname, GORDER_VERSION
        );

        self.write_ordermap::<O>(filename, molname, &comment)
    }

    /// Write the map of order parameters for a bond type or an atom type into an output file.
    fn write_ordermap<O: OrderType>(
        &self,
        filename: &str,
        molname: &str,
        comment: &str,
    ) -> Result<(), OrderMapWriteError> {
        // directory must already exist
        let directory = self.params().output_directory();

        let full_path = format!("{}/{}/{}", directory, molname, filename);
        let output_file = File::create(&full_path).map_err(|_| {
            OrderMapWriteError::CouldNotCreateFile(Box::from(Path::new(&full_path)))
        })?;
        let mut output = BufWriter::new(output_file);

        writeln!(output, "{}", comment)
            .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(Path::new(&full_path))))?;

        let (label_x, label_y) = self.params().plane().expect(PANIC_MESSAGE).get_labels();

        let label_z = O::zlabel();

        writeln!(
            output,
            "@ xlabel {label_x}-dimension [nm]\n@ ylabel {label_y}-dimension [nm]\n@ zlabel {label_z}\n@ zrange -1 1 0.2"
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
                O::convert(value.2 / samples.2 as f32, None).value
            };

            writeln!(output, "{:.4} {:.4} {:.4}", value.0, value.1, average).map_err(|_| {
                OrderMapWriteError::CouldNotWriteLine(Box::from(Path::new(&full_path)))
            })?;
        }

        Ok(())
    }
}

impl MoleculeType {
    /// Create output directory for ordermaps of a single molecule.
    fn create_molecule_dir(&self) -> Result<(), OrderMapWriteError> {
        if let Some(map) = self
            .order_bonds()
            .bond_types()
            .first()
            .and_then(|bond| bond.total_map().as_ref())
        {
            let dirname = format!("{}/{}", map.params().output_directory(), self.name());
            let path = Path::new(&dirname);

            std::fs::create_dir(path).map_err(|_| {
                OrderMapWriteError::CouldNotCreateDirectory(Box::from(Path::new(path)))
            })?;
        }

        Ok(())
    }

    /// Write ordermaps constructed for all the bond types of this molecule type.
    fn write_ordermaps_bonds<O: OrderType>(&self, molname: &str) -> Result<(), OrderMapWriteError> {
        let mut all_maps = Vec::new();
        for bond in self.order_bonds().bond_types() {
            all_maps.push((
                bond.total_map().clone(),
                bond.upper_map().clone(),
                bond.lower_map().clone(),
            ));

            if let Some(map) = bond.total_map() {
                map.write_bond_map::<O>(bond.atom1(), bond.atom2(), None, molname)?;
            }

            if let Some(map) = bond.upper_map() {
                map.write_bond_map::<O>(bond.atom1(), bond.atom2(), Some(Leaflet::Upper), molname)?;
            }

            if let Some(map) = bond.lower_map() {
                map.write_bond_map::<O>(bond.atom1(), bond.atom2(), Some(Leaflet::Lower), molname)?;
            }
        }

        let (average_total, average_upper, average_lower) = MoleculeType::merge_maps(all_maps);

        if let Some(map) = average_total {
            map.write_average_map::<O>(None, molname)?;
        }

        if let Some(map) = average_upper {
            map.write_average_map::<O>(Some(Leaflet::Upper), molname)?;
        }

        if let Some(map) = average_lower {
            map.write_average_map::<O>(Some(Leaflet::Lower), molname)?;
        }

        Ok(())
    }

    /// Write ordermaps for all heavy atom types of this molecule type.
    /// Ordermap for an atom type is constructed by merging maps for bond types the atom type is involved in.
    fn write_ordermaps_atoms<O: OrderType>(&self, molname: &str) -> Result<(), OrderMapWriteError> {
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
                map.write_atom_map::<O>(heavy_atom, None, molname)?;
            }

            if let Some(map) = total_upper {
                map.write_atom_map::<O>(heavy_atom, Some(Leaflet::Upper), molname)?;
            }

            if let Some(map) = total_lower {
                map.write_atom_map::<O>(heavy_atom, Some(Leaflet::Lower), molname)?;
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
    /// Prepare output directories for the individual molecules.
    #[inline(always)]
    pub(crate) fn prepare_directories(&self) -> Result<(), OrderMapWriteError> {
        self.molecule_types()
            .iter()
            .try_for_each(|mol| mol.create_molecule_dir())
    }

    /// Write all ordermaps constructed for all bonds of this topology.
    #[inline(always)]
    pub(crate) fn write_ordermaps_bonds<O: OrderType>(&self) -> Result<(), OrderMapWriteError> {
        self.molecule_types()
            .iter()
            .try_for_each(|mol| mol.write_ordermaps_bonds::<O>(mol.name()))
    }

    /// Write all ordermaps constructed for all atoms of this topology.
    #[inline(always)]
    pub(crate) fn write_ordermaps_atoms<O: OrderType>(&self) -> Result<(), OrderMapWriteError> {
        self.molecule_types()
            .iter()
            .try_for_each(|mol| mol.write_ordermaps_atoms::<O>(mol.name()))
    }

    /// Create or back up the directory for the ordermaps.
    pub(crate) fn handle_ordermap_directory(
        &self,
        overwrite: bool,
    ) -> Result<(), OrderMapWriteError> {
        if let Some(map_unwrapped) = self
            .molecule_types()
            .iter()
            .filter_map(|m| m.order_bonds().bond_types().first())
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
    use std::io::{BufRead, BufReader};

    use approx::assert_relative_eq;
    use groan_rs::prelude::SimBox;
    use tempfile::TempDir;

    use crate::{
        analysis::order::OrderValue, input::ordermap::Plane, input::OrderMap, presentation::CGOrder,
    };

    use super::*;

    /// Test utility. Diff the contents of two files without the first `skip` lines.
    pub(crate) fn diff_files_ignore_first(file1: &str, file2: &str, skip: usize) -> bool {
        let content1 = read_file_without_first_lines(file1, skip);
        let content2 = read_file_without_first_lines(file2, skip);
        content1 == content2
    }

    fn read_file_without_first_lines(file: &str, skip: usize) -> Vec<String> {
        let reader = BufReader::new(File::open(file).unwrap());
        reader
            .lines()
            .skip(skip) // Skip the first line
            .map(|line| line.unwrap())
            .collect()
    }

    #[test]
    fn test_write_map() {
        let directory = TempDir::new().unwrap();
        let path = directory.path().to_str().unwrap();
        let path_to_inner = format!("{}/{}", path, "molecule");
        std::fs::create_dir(&path_to_inner).unwrap();

        let params = OrderMap::new()
            .output_directory(path)
            .bin_size([1.0, 1.0])
            .min_samples(10)
            .plane(Plane::XY)
            .build()
            .unwrap();

        let mut map = Map::new(params, &SimBox::from([3.0, 3.0, 10.0])).unwrap();

        *(map.values_mut().get_mut_at(0.0, 0.0).unwrap()) = OrderValue::from(2.618);
        *(map.samples_mut().get_mut_at(0.0, 0.0).unwrap()) = 17;

        *(map.values_mut().get_mut_at(0.0, 1.0).unwrap()) = OrderValue::from(-2.262);
        *(map.samples_mut().get_mut_at(0.0, 1.0).unwrap()) = 10;

        *(map.values_mut().get_mut_at(1.0, 0.0).unwrap()) = OrderValue::from(0.364);
        *(map.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 2;

        *(map.values_mut().get_mut_at(2.0, 2.0).unwrap()) = OrderValue::from(0.764);
        *(map.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 9;

        *(map.values_mut().get_mut_at(3.0, 1.0).unwrap()) = OrderValue::from(8.826);
        *(map.samples_mut().get_mut_at(3.0, 1.0).unwrap()) = 15;

        let atom1 = AtomType::new_raw(4, "POPC", "C22");
        let atom2 = AtomType::new_raw(6, "POPC", "H22");

        map.write_bond_map::<CGOrder>(&atom1, &atom2, None, "molecule")
            .unwrap();
        map.write_bond_map::<CGOrder>(&atom1, &atom2, Some(Leaflet::Upper), "molecule")
            .unwrap();
        map.write_bond_map::<CGOrder>(&atom1, &atom2, Some(Leaflet::Lower), "molecule")
            .unwrap();

        map.write_atom_map::<CGOrder>(&atom1, None, "molecule")
            .unwrap();
        map.write_atom_map::<CGOrder>(&atom1, Some(Leaflet::Upper), "molecule")
            .unwrap();
        map.write_atom_map::<CGOrder>(&atom1, Some(Leaflet::Lower), "molecule")
            .unwrap();

        let results_bonds = [
            "ordermap_POPC-C22-4--POPC-H22-6_full.dat",
            "ordermap_POPC-C22-4--POPC-H22-6_upper.dat",
            "ordermap_POPC-C22-4--POPC-H22-6_lower.dat",
        ];

        let results_atom = [
            "ordermap_POPC-C22-4_full.dat",
            "ordermap_POPC-C22-4_upper.dat",
            "ordermap_POPC-C22-4_lower.dat",
        ];

        for result_file in results_bonds.iter() {
            let full_path = format!("{}/{}", path_to_inner, result_file);

            assert!(diff_files_ignore_first(
                &full_path,
                "tests/files/ordermap_bonds_expected.dat",
                2
            ));
        }

        for result_file in results_atom.iter() {
            let full_path = format!("{}/{}", path_to_inner, result_file);

            assert!(diff_files_ignore_first(
                &full_path,
                "tests/files/ordermap_atom_expected.dat",
                2
            ));
        }
    }

    #[test]
    fn test_merge_maps() {
        let params = OrderMap::new()
            .output_directory(".")
            .bin_size([1.0, 1.0])
            .min_samples(10)
            .plane(Plane::XY)
            .build()
            .unwrap();

        let mut map1 = Map::new(params.clone(), &SimBox::from([3.0, 3.0, 10.0])).unwrap();
        *(map1.values_mut().get_mut_at(0.0, 0.0).unwrap()) = OrderValue::from(2.618);
        *(map1.samples_mut().get_mut_at(0.0, 0.0).unwrap()) = 17;

        *(map1.values_mut().get_mut_at(0.0, 1.0).unwrap()) = OrderValue::from(-2.262);
        *(map1.samples_mut().get_mut_at(0.0, 1.0).unwrap()) = 10;

        *(map1.values_mut().get_mut_at(1.0, 0.0).unwrap()) = OrderValue::from(0.364);
        *(map1.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 2;

        *(map1.values_mut().get_mut_at(2.0, 2.0).unwrap()) = OrderValue::from(0.764);
        *(map1.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 9;

        *(map1.values_mut().get_mut_at(3.0, 1.0).unwrap()) = OrderValue::from(8.826);
        *(map1.samples_mut().get_mut_at(3.0, 1.0).unwrap()) = 15;

        let mut map2 = Map::new(params.clone(), &SimBox::from([3.0, 3.0, 10.0])).unwrap();
        *(map2.values_mut().get_mut_at(0.0, 0.0).unwrap()) = OrderValue::from(4.174);
        *(map2.samples_mut().get_mut_at(0.0, 0.0).unwrap()) = 35;

        *(map2.values_mut().get_mut_at(1.0, 1.0).unwrap()) = OrderValue::from(1.654);
        *(map2.samples_mut().get_mut_at(1.0, 1.0).unwrap()) = 4;

        *(map2.values_mut().get_mut_at(1.0, 0.0).unwrap()) = OrderValue::from(0.549);
        *(map2.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 5;

        *(map2.values_mut().get_mut_at(2.0, 2.0).unwrap()) = OrderValue::from(1.677);
        *(map2.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 12;

        *(map2.values_mut().get_mut_at(3.0, 2.0).unwrap()) = OrderValue::from(-3.453);
        *(map2.samples_mut().get_mut_at(3.0, 2.0).unwrap()) = 15;

        let mut map3 = Map::new(params.clone(), &SimBox::from([3.0, 3.0, 10.0])).unwrap();
        *(map3.values_mut().get_mut_at(1.0, 0.0).unwrap()) = OrderValue::from(2.764);
        *(map3.samples_mut().get_mut_at(1.0, 0.0).unwrap()) = 14;

        *(map3.values_mut().get_mut_at(1.0, 1.0).unwrap()) = OrderValue::from(0.874);
        *(map3.samples_mut().get_mut_at(1.0, 1.0).unwrap()) = 8;

        *(map3.values_mut().get_mut_at(2.0, 2.0).unwrap()) = OrderValue::from(-0.123);
        *(map3.samples_mut().get_mut_at(2.0, 2.0).unwrap()) = 41;

        *(map3.values_mut().get_mut_at(3.0, 0.0).unwrap()) = OrderValue::from(0.434);
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
            assert_relative_eq!(values_a.get_at_convert(x, y).unwrap(), expected);
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
            assert_relative_eq!(values_b.get_at_convert(x, y).unwrap(), expected);
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
            assert_relative_eq!(values_c.get_at_convert(x, y).unwrap(), expected);
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

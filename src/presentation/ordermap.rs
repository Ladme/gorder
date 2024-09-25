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
    /// Write ordermaps constructed for all the bonds of this molecule type.
    pub(crate) fn write_ordermaps_bonds(&self) -> Result<(), OrderMapWriteError> {
        for bond in self.order_bonds().bonds() {
            if let Some(map) = bond.total_map() {
                map.write_bond_map(bond.bond_type().atom1(), bond.bond_type().atom2(), None)?;
            }

            if let Some(map) = bond.upper_map() {
                map.write_bond_map(
                    bond.bond_type().atom1(),
                    bond.bond_type().atom2(),
                    Some(Leaflet::Upper),
                )?;
            }

            if let Some(map) = bond.lower_map() {
                map.write_bond_map(
                    bond.bond_type().atom1(),
                    bond.bond_type().atom2(),
                    Some(Leaflet::Lower),
                )?;
            }
        }

        Ok(())
    }
}

impl SystemTopology {
    /// Write all ordermaps consturcted for all bonds of this topology.
    pub(crate) fn write_ordermaps_bonds(&self) -> Result<(), OrderMapWriteError> {
        for molecule in self.molecules() {
            molecule.write_ordermaps_bonds()?;
        }

        Ok(())
    }

    /// Create or back up the directory for the ordermaps.
    pub(crate) fn handle_ordermap_directory(
        &self,
        overwrite: bool,
    ) -> Result<(), OrderMapWriteError> {
        if let Some(map_unwrapped) = self
            .molecules()
            .iter()
            .filter_map(|m| m.order_bonds().bonds().get(0))
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
}

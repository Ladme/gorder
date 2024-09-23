// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

use std::io::Write;
use std::{fs::File, io::BufWriter, path::Path};

use crate::analysis::molecule::MoleculeType;
use crate::analysis::topology::SystemTopology;
use crate::GORDER_VERSION;
use crate::{
    analysis::{molecule::AtomType, ordermap::Map},
    errors::OrderMapWriteError,
    Leaflet,
};

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

        let comment = format!("# Map of average order parameters calculated for bonds involving atom type {}.\n# Calculated with 'gorder v{}.", atom, GORDER_VERSION);

        self.write_ordermap(&filename, &comment)
    }

    /// Write the map of order parameters for a single bond into an output file.
    /// /// Leaflet `None` corresponds to ordermap for the full membrane.
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

        let comment = format!("# Map of average order parameters calculated for bonds between atom types {} and {}.\n# Calculated with 'gorder v{}.'", atom1, atom2, GORDER_VERSION);

        self.write_ordermap(&filename, &comment)
    }

    fn write_ordermap(&self, filename: &str, comment: &str) -> Result<(), OrderMapWriteError> {
        let directory = self.params().output_directory();

        // create the directory if it does not exist
        if !Path::new(directory).is_dir() {
            std::fs::create_dir_all(directory).map_err(|_| {
                OrderMapWriteError::CouldNotCreateDirectory(Box::from(Path::new(directory)))
            })?;
        }

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
}

// Released under MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! Structures and methods for writing ordermaps.

use crate::errors::{OrderMapWriteError, WriteError};
use crate::input::Plane;
use crate::presentation::aaresults::{AAAtomResults, AAMoleculeResults};
use crate::presentation::cgresults::CGMoleculeResults;
use crate::presentation::{
    BondResults, GridMapF32, MoleculeResults, OrderResults, OrderType, OutputFormat, Presenter,
    PresenterProperties,
};
use crate::{GORDER_VERSION, PANIC_MESSAGE};
use getset::Getters;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Collection of (up to) 3 order maps: for the full membrane, the upper leaflet,
/// and the lower leaflet.
#[derive(Debug, Clone, Default, Getters)]
pub struct OrderMapsCollection {
    #[getset(get = "pub")]
    total: Option<GridMapF32>,
    #[getset(get = "pub")]
    upper: Option<GridMapF32>,
    #[getset(get = "pub")]
    lower: Option<GridMapF32>,
}

impl OrderMapsCollection {
    pub(super) fn new(
        total: Option<GridMapF32>,
        upper: Option<GridMapF32>,
        lower: Option<GridMapF32>,
    ) -> Self {
        Self {
            total,
            upper,
            lower,
        }
    }
}

/// Structure handling the writing of ordermaps.
#[derive(Debug, Clone)]
pub(super) struct OrderMapPresenter<'a, R: OrderResults> {
    /// Results of the analysis.
    results: &'a R,
    /// Parameters necessary for the writing of the ordermaps.
    properties: OrderMapProperties,
}

/// Structure containing parameters necessary for the writing of the xvg file.
#[derive(Debug, Clone)]
pub(crate) struct OrderMapProperties {
    /// Plane in which the ordermaps are constructed.
    plane: Plane,
}

/// Trait implemented by all structures from which ordermaps can be written.
pub(crate) trait MapWrite {
    /// Write the ordermaps that are stored in the structure.
    fn write_map<O: OrderType>(
        &self,
        path: impl AsRef<Path>,
        properties: &OrderMapProperties,
    ) -> Result<(), OrderMapWriteError>;
}

impl PresenterProperties for OrderMapProperties {
    #[allow(unused)]
    fn leaflets(&self) -> bool {
        panic!("FATAL GORDER ERROR | OrderMapProperties::leaflets | This method should never be called. {}", PANIC_MESSAGE);
    }
}

impl OrderMapProperties {
    /// Create new structure capturing the properties of the ordermaps.
    pub(super) fn new(plane: Plane) -> Self {
        Self { plane }
    }
}

impl<'a, R: OrderResults> Presenter<'a, R> for OrderMapPresenter<'a, R> {
    type Properties = OrderMapProperties;

    fn new(results: &'a R, properties: OrderMapProperties) -> Self {
        Self {
            results,
            properties,
        }
    }

    fn file_format(&self) -> OutputFormat {
        OutputFormat::MAP
    }

    #[allow(unused)]
    fn write_results(&self, writer: &mut impl Write) -> Result<(), WriteError> {
        panic!("FATAL GORDER ERROR | OrderMapPresenter::write_results | This method should never be called. {}", PANIC_MESSAGE);
    }

    /// Does nothing, empty map is just not written.
    fn write_empty_order(
        _writer: &mut impl Write,
        _properties: &OrderMapProperties,
    ) -> Result<(), WriteError> {
        Ok(())
    }

    /// Write all ordermaps. Handles backing of the directory.
    fn write(&self, directory: impl AsRef<Path>, overwrite: bool) -> Result<(), WriteError> {
        prepare_ordermaps_directory(&directory, overwrite)
            .map_err(|e| WriteError::CouldNotWriteOrderMap(e))?;
        for molecule in self.results.molecules() {
            molecule
                .write_map::<R::OrderType>(&directory, &self.properties)
                .map_err(|e| WriteError::CouldNotWriteOrderMap(e))?;
        }

        Ok(())
    }
}

impl OrderMapsCollection {
    /// Write maps from a single OrderMapsCollection.
    fn write_maps<O: OrderType>(
        &self,
        directory: &impl AsRef<Path>,
        name_pattern: &str,
        plane: Plane,
        comment: &str,
    ) -> Result<(), OrderMapWriteError> {
        for (map_type, leaflet) in [self.total(), self.upper(), self.lower()]
            .into_iter()
            .zip(["full", "upper", "lower"])
        {
            if let Some(map) = map_type {
                let filename = format!("{}_{}.dat", name_pattern, leaflet);
                let path = directory.as_ref().join(Path::new(&filename));
                write_ordermap::<O>(map, path, plane, comment)?;
            }
        }

        Ok(())
    }
}

impl MapWrite for BondResults {
    fn write_map<O: OrderType>(
        &self,
        directory: impl AsRef<Path>,
        properties: &OrderMapProperties,
    ) -> Result<(), OrderMapWriteError> {
        let atom1 = self.bond().atom1();
        let atom2 = self.bond().atom2();
        let name = format!("ordermap_{}--{}", atom1, atom2);
        let comment = format!("# Map of average order parameters calculated for bonds between atom types {} and {} of molecule type {}.\n# Calculated with 'gorder v{}'.",
                              atom1, atom2, self.molecule(), GORDER_VERSION);

        self.ordermaps()
            .write_maps::<O>(&directory, &name, properties.plane, &comment)
    }
}

impl MapWrite for AAMoleculeResults {
    fn write_map<O: OrderType>(
        &self,
        path: impl AsRef<Path>,
        properties: &OrderMapProperties,
    ) -> Result<(), OrderMapWriteError> {
        let name = self.molecule();
        let directory = path.as_ref().join(Path::new(name));

        std::fs::create_dir(&directory)
            .map_err(|_| OrderMapWriteError::CouldNotCreateDirectory(Box::from(path.as_ref())))?;

        // write average ordermaps
        let comment = format!("# Map of average order parameters calculated for molecule type {}.\n# Calculated with 'gorder v{}'.",
                              self.molecule(), GORDER_VERSION);
        let name = "ordermap_average";
        self.average_ordermaps()
            .write_maps::<O>(&directory, &name, properties.plane, &comment)?;

        // write ordermaps for atoms
        self.order()
            .values()
            .try_for_each(|atom| atom.write_map::<O>(&directory, properties))
    }
}

impl MapWrite for AAAtomResults {
    fn write_map<O: OrderType>(
        &self,
        directory: impl AsRef<Path>,
        properties: &OrderMapProperties,
    ) -> Result<(), OrderMapWriteError> {
        let name = format!("ordermap_{}", self.atom());
        let comment = format!("# Map of average order parameters calculated for atom type {} of molecule type {}.\n# Calculated with 'gorder v{}'.",
                              self.atom(), self.molecule(), GORDER_VERSION);

        // write the ordermaps for the atom itself
        self.ordermaps()
            .write_maps::<O>(&directory, &name, properties.plane, &comment)?;

        // write ordermaps for the individual bonds
        self.bonds()
            .values()
            .try_for_each(|bond| bond.write_map::<O>(&directory, properties))
    }
}

impl MapWrite for CGMoleculeResults {
    fn write_map<O: OrderType>(
        &self,
        path: impl AsRef<Path>,
        properties: &OrderMapProperties,
    ) -> Result<(), OrderMapWriteError> {
        let name = self.molecule();
        let directory = path.as_ref().join(Path::new(name));

        std::fs::create_dir(&directory)
            .map_err(|_| OrderMapWriteError::CouldNotCreateDirectory(Box::from(path.as_ref())))?;

        // write average ordermaps
        let comment = format!("# Map of average order parameters calculated for molecule type {}.\n# Calculated with 'gorder v{}'.",
                              self.molecule(), GORDER_VERSION);
        let name = "ordermap_average";
        self.average_ordermaps()
            .write_maps::<O>(&directory, &name, properties.plane, &comment)?;

        // write ordermaps for individual bonds
        self.order()
            .values()
            .try_for_each(|bond| bond.write_map::<O>(&directory, properties))
    }
}

fn prepare_ordermaps_directory(
    directory: &impl AsRef<Path>,
    overwrite: bool,
) -> Result<(), OrderMapWriteError> {
    log::info!(
        "Writing ordermaps into a directory '{}'...",
        directory.as_ref().to_str().expect(PANIC_MESSAGE)
    );

    log::logger().flush();

    if directory.as_ref().is_dir() {
        if !overwrite {
            log::warn!(
                "Output directory for ordermaps '{}' already exists. Backing it up.",
                directory.as_ref().to_str().expect(PANIC_MESSAGE)
            );
            backitup::backup(directory).map_err(|_| {
                OrderMapWriteError::CouldNotBackupDirectory(Box::from(directory.as_ref()))
            })?;
        } else {
            log::warn!("Output directory for ordermaps '{}' already exists. It will be overwritten as requested.",
                        directory.as_ref().to_str().expect(PANIC_MESSAGE)
                    );
            std::fs::remove_dir_all(directory).map_err(|_| {
                OrderMapWriteError::CouldNotRemoveDirectory(Box::from(directory.as_ref()))
            })?;
        }
    }

    // create a new directory
    std::fs::create_dir_all(directory)
        .map_err(|_| OrderMapWriteError::CouldNotCreateDirectory(Box::from(directory.as_ref())))?;

    Ok(())
}

/// Write the map of order parameters for a bond type or an atom type into an output file.
fn write_ordermap<O: OrderType>(
    map: &GridMapF32,
    full_path: impl AsRef<Path>,
    plane: Plane,
    comment: &str,
) -> Result<(), OrderMapWriteError> {
    // directory must already exist
    let output_file = File::create(&full_path)
        .map_err(|_| OrderMapWriteError::CouldNotCreateFile(Box::from(full_path.as_ref())))?;
    let mut output = BufWriter::new(output_file);

    writeln!(output, "{}", comment)
        .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(full_path.as_ref())))?;

    let (label_x, label_y) = plane.get_labels();
    let label_z = O::zlabel();

    writeln!(
        output,
        "@ xlabel {label_x}-dimension [nm]\n@ ylabel {label_y}-dimension [nm]\n@ zlabel {label_z}\n@ zrange -1 1 0.2"
    )
        .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(full_path.as_ref())))?;

    writeln!(output, "$ type colorbar\n$ colormap seismic_r")
        .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(full_path.as_ref())))?;

    for (x, y, z) in map.extract_convert() {
        writeln!(output, "{:.4} {:.4} {:.4}", x, y, z)
            .map_err(|_| OrderMapWriteError::CouldNotWriteLine(Box::from(full_path.as_ref())))?;
    }

    Ok(())
}

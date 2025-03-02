// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Structures and methods for calculations of order parameters from united-atom simulations.

use std::{f32::consts::PI, ops::Deref};

use groan_rs::{prelude::Vector3D, system::System};
use nalgebra::{Rotation3, Unit};

use crate::{errors::AnalysisError, prelude::AtomType, Leaflet, PANIC_MESSAGE};

use super::{
    calc_sch,
    geometry::GeometrySelection,
    leaflets::MoleculeLeafletClassification,
    molecule::{BondLike, BondType},
    normal::MoleculeMembraneNormal,
    order::{AnalysisOrder, OrderValue},
    ordermap::Map,
    pbc::PBCHandler,
};

/// Tetrahedral angle.
const TETRAHEDRAL_ANGLE: f32 = 1.910633;
/// Half of the tetrahedral angle.
const TETRAHEDRAL_ANGLE_HALF: f32 = 0.9553165;
/// Length of the C-H bond in nm.
const BOND_LENGTH: f32 = 0.109;
/// Rotation angle used for predicting hydrogens in the CH3 group.
const CH3_ANGLE: f32 = 2.0943952;

/// Enum representing a type of united atom.
/// Contains indices of all atoms of this type and its helper atoms.
#[derive(Debug, Clone)]
pub(super) enum UAOrderAtomType {
    CH3(CH3Atom),
    CH2(CH2Atom),
    CH1(CH1Atom),
}

impl UAOrderAtomType {
    fn analyze_frame<'a, Geom: GeometrySelection>(
        &mut self,
        frame: &'a System,
        leaflet_classification: &mut Option<MoleculeLeafletClassification>,
        pbc_handler: &'a impl PBCHandler<'a>,
        membrane_normal: &mut MoleculeMembraneNormal,
        frame_index: usize,
        geometry: &Geom,
    ) -> Result<(), AnalysisError> {
        macro_rules! impl_analyze_frame {
            ($enum:ident => $($variant:ident),+) => {
                match self {
                    $(
                        Self::$variant(x) => x.analyze_frame(
                            frame,
                            leaflet_classification,
                            pbc_handler,
                            membrane_normal,
                            frame_index,
                            geometry,
                        ),
                    )+
                }
            };
        }

        impl_analyze_frame!(Self => CH3, CH2, CH1)
    }
}

trait UAAtom<const N: usize> {
    /// Predict positions of N hydrogens for molecule with the given index.
    fn predict_hydrogens<'a>(
        triplet: &AtomTriplet,
        frame: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; N], AnalysisError>;

    /// Iterate over the atom triplets of the UA atom.
    fn triplets_iter(&self) -> impl Iterator<Item = &AtomTriplet>;

    /// Get mutable reference to UA bond type with the given index.
    /// Should panic if the index is out of range.
    fn get_bond(&mut self, index: usize) -> &mut VirtualBondType;

    /// Use the predicted hydrogen positions to calculate order parameters and bond positions.
    /// If a bond is not inside the specified geometry, return `None` for it.
    fn calculate_sch<'a, Geom: GeometrySelection>(
        target: &Vector3D,
        hydrogens: [Vector3D; N],
        pbc: &'a impl PBCHandler<'a>,
        normal: &Vector3D,
        geometry: &Geom,
    ) -> ([Option<f32>; N], [Vector3D; N]) {
        let bond_positions = hydrogens.clone().map(|h| {
            let vec = pbc.vector_to(target, &h);
            h + (&vec / 2.0)
        });

        let mut sch = [None; N];
        for (i, (h, b)) in hydrogens.iter().zip(bond_positions.iter()).enumerate() {
            sch[i] = if pbc.inside(b, geometry) {
                Some(calc_sch(h, normal))
            } else {
                None
            };
        }

        (sch, bond_positions)
    }

    /// Calculate order parameters for a single atom type in a single trajectory frame.
    fn analyze_frame<'a, Geom: GeometrySelection>(
        &mut self,
        frame: &'a System,
        leaflet_classification: &mut Option<MoleculeLeafletClassification>,
        pbc_handler: &'a impl PBCHandler<'a>,
        membrane_normal: &mut MoleculeMembraneNormal,
        frame_index: usize,
        geometry: &Geom,
    ) -> Result<(), AnalysisError> {
        let self_ptr = self as *mut Self;

        for (molecule_index, triplet) in self.triplets_iter().enumerate() {
            let normal =
                membrane_normal.get_normal(frame_index, molecule_index, frame, pbc_handler)?;

            let hydrogens = Self::predict_hydrogens(triplet, frame, pbc_handler)?;
            let target_pos = triplet.get_target(frame)?;
            let (sch, bond_pos) =
                Self::calculate_sch(target_pos, hydrogens, pbc_handler, normal, geometry);

            for (i, (order, pos)) in sch.into_iter().zip(bond_pos.into_iter()).enumerate() {
                match order {
                    // safety: self lives the entire method
                    // we are modifying different part of the structure than is iterated by `triplets_iter`
                    Some(x) => unsafe { &mut *self_ptr }.get_bond(i).add_order(
                        molecule_index,
                        x,
                        &pos,
                        leaflet_classification,
                        frame_index,
                    ),
                    None => continue,
                }
            }
        }

        Ok(())
    }
}

/// Methyl atom type.
#[derive(Debug, Clone)]
struct CH3Atom {
    atom_type: AtomType,
    atoms: Vec<AtomTriplet>,
    bonds: [VirtualBondType; 3],
}

impl UAAtom<3> for CH3Atom {
    #[inline(always)]
    fn triplets_iter(&self) -> impl Iterator<Item = &AtomTriplet> {
        self.atoms.iter()
    }

    #[inline(always)]
    fn get_bond(&mut self, index: usize) -> &mut VirtualBondType {
        self.bonds.get_mut(index).unwrap_or_else(|| {
            panic!(
                "FATAL GORDER ERROR | CH3Atom::get_bond | Bond index `{}` out of range. {}",
                index, PANIC_MESSAGE
            )
        })
    }

    #[inline(always)]
    fn predict_hydrogens<'a>(
        triplet: &AtomTriplet,
        frame: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; 3], AnalysisError> {
        triplet.predict_hydrogens_ch3(frame, pbc)
    }
}
/// Methylene atom type.
#[derive(Debug, Clone)]
struct CH2Atom {
    atom_type: AtomType,
    atoms: Vec<AtomTriplet>,
    bonds: [VirtualBondType; 2],
}

impl UAAtom<2> for CH2Atom {
    #[inline(always)]
    fn triplets_iter(&self) -> impl Iterator<Item = &AtomTriplet> {
        self.atoms.iter()
    }

    #[inline(always)]
    fn get_bond(&mut self, index: usize) -> &mut VirtualBondType {
        self.bonds.get_mut(index).unwrap_or_else(|| {
            panic!(
                "FATAL GORDER ERROR | CH2Atom::get_bond | Bond index `{}` out of range. {}",
                index, PANIC_MESSAGE
            )
        })
    }

    #[inline(always)]
    fn predict_hydrogens<'a>(
        triplet: &AtomTriplet,
        frame: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; 2], AnalysisError> {
        triplet.predict_hydrogens_ch2(frame, pbc)
    }
}

/// Methine atom type.
#[derive(Debug, Clone)]
struct CH1Atom {
    atom_type: AtomType,
    atoms: Vec<AtomTriplet>,
    bond: VirtualBondType,
}

impl UAAtom<1> for CH1Atom {
    #[inline(always)]
    fn triplets_iter(&self) -> impl Iterator<Item = &AtomTriplet> {
        self.atoms.iter()
    }

    #[inline(always)]
    fn get_bond(&mut self, index: usize) -> &mut VirtualBondType {
        if index != 0 {
            panic!(
                "FATAL GORDER ERROR | CH1Atom::get_bond | Bond index `{}` out of range. {}",
                index, PANIC_MESSAGE
            )
        }

        &mut self.bond
    }

    #[inline(always)]
    fn predict_hydrogens<'a>(
        triplet: &AtomTriplet,
        frame: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; 1], AnalysisError> {
        triplet.predict_hydrogen_ch1(frame, pbc)
    }
}

/// Structure storing indices of two helper atoms and a target atom.
#[derive(Debug, Clone)]
struct AtomTriplet {
    helper1: usize,
    target: usize,
    helper2: usize,
}

impl AtomTriplet {
    /// Get positions of the atoms which indices are stored in the structure.
    fn unpack2pos<'a>(
        &self,
        system: &'a System,
    ) -> Result<(&'a Vector3D, &'a Vector3D, &'a Vector3D), AnalysisError> {
        let positions: Result<Vec<_>, _> = [self.helper1, self.target, self.helper2]
            .into_iter()
            .map(|index| {
                // SAFETY: indices must be valid
                let atom = unsafe { system.get_atom_unchecked(index) };
                atom.get_position()
                    .ok_or_else(|| AnalysisError::UndefinedPosition(index))
            })
            .collect();

        match positions {
            Ok(positions) => Ok((positions[0], positions[1], positions[2])),
            Err(e) => Err(e),
        }
    }

    /// Get the position of the target atom.
    #[inline(always)]
    fn get_target<'a>(&self, system: &'a System) -> Result<&'a Vector3D, AnalysisError> {
        // SAFETY: index must be valid
        unsafe { system.get_atom_unchecked(self.target) }
            .get_position()
            .ok_or_else(|| AnalysisError::UndefinedPosition(self.target))
    }

    /// Predict positions of hydrogens of a methyl group.
    /// Adapted from https://buildh.readthedocs.io/en/latest/algorithms_Hbuilding.html#building-ch3.
    fn predict_hydrogens_ch3<'a>(
        &self,
        system: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; 3], AnalysisError> {
        let (helper1, target, helper2) = self.unpack2pos(system)?;

        let th1 = pbc.vector_to(&target, &helper1).to_unit();
        let th2 = pbc.vector_to(&target, &helper2).to_unit();

        let rot_axis = th2.cross(&th1).normalize();
        let rotation1 =
            Rotation3::from_axis_angle(&Unit::new_normalize(rot_axis.clone()), TETRAHEDRAL_ANGLE);

        let hydrogen1_vec = th1.clone().rotate(&rotation1.into());
        let mut hydrogen1 = target.clone();
        hydrogen1.shift(hydrogen1_vec.clone(), BOND_LENGTH);
        pbc.wrap(&mut hydrogen1);

        let normalized_th1 = Unit::new_normalize(th1.deref().clone());

        let rotation2 = Rotation3::from_axis_angle(&normalized_th1, CH3_ANGLE);
        let hydrogen2_vec = hydrogen1_vec.clone().rotate(&rotation2.into());
        let mut hydrogen2 = target.clone();
        hydrogen2.shift(hydrogen2_vec, BOND_LENGTH);
        pbc.wrap(&mut hydrogen2);

        let rotation3 = Rotation3::from_axis_angle(&normalized_th1, -CH3_ANGLE);
        let hydrogen3_vec = hydrogen1_vec.clone().rotate(&rotation3.into());
        let mut hydrogen3 = target.clone();
        hydrogen3.shift(hydrogen3_vec, BOND_LENGTH);
        pbc.wrap(&mut hydrogen3);

        Ok([hydrogen1, hydrogen2, hydrogen3])
    }

    /// Predict positions of hydrogens of a methylene group.
    /// Adapted from https://buildh.readthedocs.io/en/latest/algorithms_Hbuilding.html#building-ch2.
    fn predict_hydrogens_ch2<'a>(
        &self,
        system: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; 2], AnalysisError> {
        let (helper1, target, helper2) = self.unpack2pos(system)?;

        let th1 = pbc.vector_to(&target, &helper1).to_unit();
        let th2 = pbc.vector_to(&target, &helper2).to_unit();
        let plane_normal = th2.cross(&th1).normalize();
        let rot_axis = (th1 - th2).to_unit();
        let rot_vec = Vector3D::from(plane_normal.cross(&rot_axis).normalize());

        let rotation_positive = Rotation3::from_axis_angle(
            &Unit::new_normalize(rot_axis.deref().clone()),
            TETRAHEDRAL_ANGLE_HALF,
        );

        let rotation_negative = Rotation3::from_axis_angle(
            &Unit::new_normalize(rot_axis.deref().clone()),
            -TETRAHEDRAL_ANGLE_HALF,
        );

        let hydrogen1_vec = rot_vec.clone().rotate(&rotation_positive.into());
        let hydrogen2_vec = rot_vec.rotate(&rotation_negative.into());

        let mut hydrogen1 = target.clone();
        hydrogen1.shift(hydrogen1_vec, BOND_LENGTH);
        pbc.wrap(&mut hydrogen1);

        let mut hydrogen2 = target.clone();
        hydrogen2.shift(hydrogen2_vec, BOND_LENGTH);
        pbc.wrap(&mut hydrogen2);

        Ok([hydrogen1, hydrogen2])
    }

    /// Predict positions of hydrogens of a methine group.
    /// Adapted from https://buildh.readthedocs.io/en/latest/algorithms_Hbuilding.html#building-ch-on-a-double-bond.
    fn predict_hydrogen_ch1<'a>(
        &self,
        system: &'a System,
        pbc: &'a impl PBCHandler<'a>,
    ) -> Result<[Vector3D; 1], AnalysisError> {
        let (helper1, target, helper2) = self.unpack2pos(system)?;

        let th1 = pbc.vector_to(&target, &helper1).to_unit();
        let th2 = pbc.vector_to(&target, &helper2).to_unit();
        let gamma = th1.angle(&th2);
        let rot_axis = th1.cross(&th2).normalize();

        let rotation =
            Rotation3::from_axis_angle(&Unit::new_normalize(rot_axis.clone()), PI - (gamma / 2.0));

        let hydrogen_vec = th2.clone().rotate(&rotation.into());
        let mut hydrogen = target.clone();
        hydrogen.shift(hydrogen_vec, BOND_LENGTH);
        pbc.wrap(&mut hydrogen);

        Ok([hydrogen])
    }
}

/// Order parameters and ordermaps calculated for a single bond between a real atom and a virtual atom.
#[derive(Debug, Clone)]
struct VirtualBondType {
    total: AnalysisOrder<super::timewise::AddExtend>,
    upper: Option<AnalysisOrder<super::timewise::AddExtend>>,
    lower: Option<AnalysisOrder<super::timewise::AddExtend>>,
    total_map: Option<Map>,
    upper_map: Option<Map>,
    lower_map: Option<Map>,
}

impl BondLike for VirtualBondType {
    #[inline(always)]
    fn get_total(&mut self) -> &mut AnalysisOrder<super::timewise::AddExtend> {
        &mut self.total
    }

    #[inline(always)]
    fn get_upper(&mut self) -> Option<&mut AnalysisOrder<super::timewise::AddExtend>> {
        self.upper.as_mut()
    }

    #[inline(always)]
    fn get_lower(&mut self) -> Option<&mut AnalysisOrder<super::timewise::AddExtend>> {
        self.lower.as_mut()
    }

    #[inline(always)]
    fn get_total_map(&mut self) -> Option<&mut Map> {
        self.total_map.as_mut()
    }

    #[inline(always)]
    fn get_upper_map(&mut self) -> Option<&mut Map> {
        self.upper_map.as_mut()
    }

    #[inline(always)]
    fn get_lower_map(&mut self) -> Option<&mut Map> {
        self.lower_map.as_mut()
    }
}

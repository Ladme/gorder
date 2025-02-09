// Released under MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! Contains methods for working with and without periodic boundary conditions.

use groan_rs::{
    errors::{AtomError, GroupError},
    prelude::{
        AtomIterable, AtomIteratorWithBox, Cylinder, Dimension, ImmutableAtomIterable, NaiveShape, Shape, SimBox, Vector3D
    },
    system::System,
};

use crate::{errors::AnalysisError, PANIC_MESSAGE};

use super::geometry::GeometrySelection;

/// Trait implemented by structures that should handle PBC.
pub(crate) trait PBCHandler {
    /// Get the geometric center of a group.
    fn group_get_center(&self, system: &System, group: &str) -> Result<Vector3D, GroupError>;

    /// Take atoms of a group in a specified geometry and collect their positions.
    fn group_filter_geometry_get_pos<S: Shape + NaiveShape>(&self, system: &System, group: &str, shape: S) -> Result<Vec<Vector3D>, AnalysisError>;

    /// Take atoms of a group in a specified geometry and calculate their center of geometry.
    fn group_filter_geometry_get_center<S: Shape + NaiveShape>(
        &self,
        system: &System,
        group: &str,
        shape: S,
    ) -> Result<Vector3D, AtomError>;

    /// Calculate distance between two points in the specified dimensions.
    fn distance(&self, point1: &Vector3D, point2: &Vector3D, dim: Dimension) -> f32;

    /// Calculate distance between two atoms of target indices.
    fn atoms_distance(
        &self,
        system: &System,
        index1: usize,
        index2: usize,
        dim: Dimension,
    ) -> Result<f32, AtomError>;

    /// Check if a point is inside a geometric shape.
    fn inside<Geom: GeometrySelection>(&self, point: &Vector3D, shape: &Geom) -> bool;

    /// Calculate shortest vector connecting point1 with point2.
    fn vector_to(&self, point1: &Vector3D, point2: &Vector3D) -> Vector3D;

    /// Wrap the point into the simulation box. Or not, if PBC are ignored.
    fn wrap(&self, point: &mut Vector3D);

    /// Get the reference position and length in the corresponding dimension for infinite shape.
    fn get_infinite_span(&self, pos: &mut f32) -> f32;

    /// Get the simulation box center or PANIC if PBC are ignored.
    fn get_box_center(&self) -> Vector3D;

    /// Get refrence to the simulation box. Returns None, if PBC are ignored.
    fn get_simbox(&self) -> Option<&SimBox>;

    /// Construct a cylinder for local center of geometry calculation.
    fn cylinder_for_local_center(&self, x: f32, y: f32, radius: f32, orientation: Dimension) -> Cylinder;
}

/// PBCHandler that ignores all periodic boundary conditions.
#[derive(Debug, Clone)]
pub(crate) struct NoPBC;

impl PBCHandler for NoPBC {
    #[inline(always)]
    fn group_get_center(&self, system: &System, group: &str) -> Result<Vector3D, GroupError> {
        system.group_get_center_naive(group)
    }

    #[inline(always)]
    fn group_filter_geometry_get_center<S: Shape + NaiveShape>(
        &self,
        system: &System,
        group: &str,
        shape: S,
    ) -> Result<Vector3D, AtomError> {
        system
            .group_iter(group)
            .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | NoPBC::group_filter_geometry_get_center | Unknown group `{}`. {}", group, PANIC_MESSAGE))
            .filter_geometry_naive(shape)
            .get_center_naive()
    }

    fn group_filter_geometry_get_pos<S: Shape + NaiveShape>(
        &self,
        system: &System,
        group: &str,
        shape: S,
    ) -> Result<Vec<Vector3D>, AnalysisError> {
        system
            .group_iter(group)
            .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | NoPBC::group_filter_geometry_get_pos | Unknown group `{}`. {}", group, PANIC_MESSAGE))
            .filter_geometry_naive(shape)
            .map(|atom| atom
                .get_position()
                .ok_or_else(|| AnalysisError::UndefinedPosition(atom.get_index()))
                .cloned()
            )
            .collect::<Result<Vec<_>, _>>()
    }

    #[inline(always)]
    fn distance(&self, point1: &Vector3D, point2: &Vector3D, dim: Dimension) -> f32 {
        point1.distance_naive(point2, dim)
    }

    #[inline(always)]
    fn atoms_distance(
        &self,
        system: &System,
        index1: usize,
        index2: usize,
        dim: Dimension,
    ) -> Result<f32, AtomError> {
        let atom1 = system.get_atom(index1)?;
        let atom2 = system.get_atom(index2)?;

        atom1.distance_naive(atom2, dim)
    }

    #[inline(always)]
    fn inside<Geom: GeometrySelection>(&self, point: &Vector3D, shape: &Geom) -> bool {
        shape.inside_naive(point)
    }

    #[inline(always)]
    fn vector_to(&self, point1: &Vector3D, point2: &Vector3D) -> Vector3D {
        point2 - point1
    }

    #[inline(always)]
    fn wrap(&self, _point: &mut Vector3D) {} // do nothing; no wrapping is performed if PBC are ignored

    #[inline(always)]
    fn get_infinite_span(&self, pos: &mut f32) -> f32 {
        *pos = f32::MIN;
        f32::INFINITY
    }

    #[inline(always)]
    fn get_box_center(&self) -> Vector3D {
        panic!("FATAL GORDER ERROR | NoPBC::get_box_center | PBC are ignored. Can't get box center. {}", PANIC_MESSAGE);
    }

    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        None
    }

    #[inline(always)]
    fn cylinder_for_local_center(&self, x: f32, y: f32, radius: f32, orientation: Dimension) -> Cylinder {
        Cylinder::new(
            Vector3D::new(x, y, f32::MIN),
            radius,
            f32::INFINITY,
            orientation,
        )
    }
}

/// PBCHandler that assumes periodic boundary conditions in all three dimensions.
#[derive(Debug, Clone)]
pub(crate) struct PBC3D<'a>(&'a SimBox);

impl<'a> PBCHandler for PBC3D<'a> {
    #[inline(always)]
    fn group_get_center(&self, system: &System, group: &str) -> Result<Vector3D, GroupError> {
        system.group_get_center(group)
    }

    #[inline(always)]
    fn group_filter_geometry_get_center<S: Shape + NaiveShape>(
        &self,
        system: &System,
        group: &str,
        shape: S,
    ) -> Result<Vector3D, AtomError> {
        system
            .group_iter(group)
            .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | PBC3D::group_filter_geometry_get_center | Unknown group `{}`. {}", group, PANIC_MESSAGE))
            .filter_geometry(shape)
            .get_center()
    }

    fn group_filter_geometry_get_pos<S: Shape + NaiveShape>(
        &self,
        system: &System,
        group: &str,
        shape: S,
    ) -> Result<Vec<Vector3D>, AnalysisError> {
        system
            .group_iter(group)
            .unwrap_or_else(|_| panic!("FATAL GORDER ERROR | PBC3D::group_filter_geometry_get_pos | Unknown group `{}`. {}", group, PANIC_MESSAGE))
            .filter_geometry(shape)
            .map(|atom| atom
                .get_position()
                .ok_or_else(|| AnalysisError::UndefinedPosition(atom.get_index()))
                .cloned()
            )
            .collect::<Result<Vec<_>, _>>()
    }

    #[inline(always)]
    fn distance(&self, point1: &Vector3D, point2: &Vector3D, dim: Dimension) -> f32 {
        point1.distance(point2, dim, self.0)
    }

    #[inline(always)]
    fn atoms_distance(
        &self,
        system: &System,
        index1: usize,
        index2: usize,
        dim: Dimension,
    ) -> Result<f32, AtomError> {
        let atom1 = system.get_atom(index1)?;
        let atom2 = system.get_atom(index2)?;

        atom1.distance(atom2, dim, self.0)
    }

    #[inline(always)]
    fn inside<Geom: GeometrySelection>(&self, point: &Vector3D, shape: &Geom) -> bool {
        shape.inside(point, self.0)
    }

    #[inline(always)]
    fn vector_to(&self, point1: &Vector3D, point2: &Vector3D) -> Vector3D {
        // this calculation introduces minor numerical errors compared to the NoPBC version
        // as a result, the results when using PBC3D and NoPBC handlers will differ,
        // even when calculating order parameters for lipids near the box center
        // which are not affected by numerical errors introduced by making the molecules whole
        // these discrepancies become more noticeable with shorter trajectories
        point1.vector_to(point2, self.0)
    }

    #[inline(always)]
    fn wrap(&self, point: &mut Vector3D) {
        point.wrap(self.0)
    }

    #[inline(always)]
    fn get_infinite_span(&self, pos: &mut f32) -> f32 {
        *pos = 0.0;
        f32::INFINITY
    }

    #[inline(always)]
    fn get_box_center(&self) -> Vector3D {
        Vector3D::new(self.0.x / 2.0f32, self.0.y / 2.0f32, self.0.z / 2.0f32)
    }

    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        Some(self.0)
    }

    #[inline(always)]
    fn cylinder_for_local_center(&self, x: f32, y: f32, radius: f32, orientation: Dimension) -> Cylinder {
        Cylinder::new(
            Vector3D::new(x, y, 0.0),
            radius,
            f32::INFINITY,
            orientation,
        )
    }
}

impl<'a> PBC3D<'a> {
    pub(crate) fn new(simbox: &'a SimBox) -> Self {
        Self(simbox)
    }

    /// Create a periodic boundary handler for a system.
    /// Assumes that the simulation box is defined and orthogonal.
    pub(super) fn from_system(system: &'a System) -> Self {
        Self(
            system.get_box()
            .unwrap_or_else(|| 
                panic!("FATAL GORDER ERROR | PBC3D::from_system | Simulation box is undefined but this should have been handled before.")
            )
        )
    }
}

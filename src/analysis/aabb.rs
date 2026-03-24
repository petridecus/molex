//! Axis-aligned bounding box.

use glam::Vec3;

/// Axis-aligned bounding box (AABB).
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner of the bounding box.
    pub min: Vec3,
    /// Maximum corner of the bounding box.
    pub max: Vec3,
}

impl Aabb {
    /// Geometric center of the box.
    #[must_use]
    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }

    /// Size along each axis (max - min).
    #[must_use]
    pub fn extents(&self) -> Vec3 {
        self.max - self.min
    }

    /// Half-diagonal length (bounding sphere radius from center).
    #[must_use]
    pub fn radius(&self) -> f32 {
        self.extents().length() * 0.5
    }

    /// Merge two AABBs into one that contains both.
    #[must_use]
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }

    /// Build AABB from positions. Returns `None` if the slice is empty.
    #[must_use]
    pub fn from_positions(positions: &[Vec3]) -> Option<Aabb> {
        let first = *positions.first()?;
        let mut min = first;
        let mut max = first;
        for &p in &positions[1..] {
            min = min.min(p);
            max = max.max(p);
        }
        Some(Aabb { min, max })
    }

    /// Build unified AABB from multiple AABBs.
    #[must_use]
    pub fn from_aabbs(aabbs: &[Aabb]) -> Option<Aabb> {
        aabbs.iter().copied().reduce(|a, b| a.union(&b))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_from_positions() {
        let positions = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, 5.0, 0.0),
            Vec3::new(3.0, -2.0, 7.0),
        ];
        let aabb = Aabb::from_positions(&positions).unwrap();
        assert_eq!(aabb.min, Vec3::new(-1.0, -2.0, 0.0));
        assert_eq!(aabb.max, Vec3::new(3.0, 5.0, 7.0));
    }

    #[test]
    fn test_aabb_empty() {
        assert!(Aabb::from_positions(&[]).is_none());
    }

    #[test]
    fn test_aabb_union() {
        let a = Aabb {
            min: Vec3::new(0.0, 0.0, 0.0),
            max: Vec3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Vec3::new(-1.0, 2.0, -3.0),
            max: Vec3::new(0.5, 4.0, 0.5),
        };
        let merged = a.union(&b);
        assert_eq!(merged.min, Vec3::new(-1.0, 0.0, -3.0));
        assert_eq!(merged.max, Vec3::new(1.0, 4.0, 1.0));
    }

    #[test]
    fn test_aabb_from_aabbs() {
        let aabbs = vec![
            Aabb {
                min: Vec3::ZERO,
                max: Vec3::ONE,
            },
            Aabb {
                min: Vec3::splat(2.0),
                max: Vec3::splat(3.0),
            },
        ];
        let merged = Aabb::from_aabbs(&aabbs).unwrap();
        assert_eq!(merged.min, Vec3::ZERO);
        assert_eq!(merged.max, Vec3::splat(3.0));
        assert!(Aabb::from_aabbs(&[]).is_none());
    }

    #[test]
    fn test_aabb_center_extents_radius() {
        let aabb = Aabb {
            min: Vec3::ZERO,
            max: Vec3::new(4.0, 6.0, 8.0),
        };
        assert_eq!(aabb.center(), Vec3::new(2.0, 3.0, 4.0));
        assert_eq!(aabb.extents(), Vec3::new(4.0, 6.0, 8.0));
        let expected_radius = Vec3::new(4.0, 6.0, 8.0).length() * 0.5;
        assert!((aabb.radius() - expected_radius).abs() < 1e-6);
    }
}

use std::ptr::eq;

use super::*;

#[derive(Debug, Clone, Copy)]
pub enum InteractionKind<'a> {
    Clash {
        k2: f32,
        atom3: AtomRef<'a>,
        sq_dist: f32,
    },
    Stretch {
        tolerance: f32,
    },
    Bend {
        k2: f32,
        atom3: AtomRef<'a>,
    },
    Inversion {
        atom3: AtomRef<'a>,
        atom4: AtomRef<'a>,
        is_z: bool,
        force_move: bool,
    },
}

#[derive(Debug)]
pub struct Interaction<'a> {
    pub kind: InteractionKind<'a>,
    pub atom1: AtomRef<'a>,
    pub atom2: AtomRef<'a>,
    pub rest_v: f32,
    pub k: f32,
}
impl Interaction<'_> {
    #[inline]
    fn energy(&mut self, score: &mut f32) {
        match &mut self.kind {
            InteractionKind::Clash { k2, atom3, sq_dist } => {
                *sq_dist = math::sq_dist_point_segment(
                    self.atom2.borrow().coordinates,
                    self.atom1.borrow().coordinates,
                    atom3.borrow().coordinates,
                )
                .0;
                let dr = self.rest_v - *sq_dist;
                if dr > 0.0 {
                    *score += 0.5 * self.k * *k2 * dr;
                }
            }
            InteractionKind::Bend { k2, atom3, .. } => {
                let da = math::unsigned_angle(
                    self.atom1.borrow().coordinates,
                    self.atom2.borrow().coordinates,
                    atom3.borrow().coordinates,
                );
                *score += self.k * *k2 * da * da * 5.0;
            }
            InteractionKind::Inversion {
                atom3, atom4, is_z, ..
            } => {
                if math::same_side(
                    self.atom1.borrow().coordinates,
                    atom4.borrow().coordinates,
                    self.atom2.borrow().coordinates,
                    atom3.borrow().coordinates,
                ) != *is_z
                {
                    *score += 5000.0;
                }
            }
            _ => {
                let l = self.atom1.borrow().coordinates - self.atom2.borrow().coordinates;
                *score += 0.5 * self.k * l.sq_length();
            }
        }
    }
    pub fn score(&mut self, score: &mut f32) {
        self.energy(score);
        if eq(self.atom1, self.atom2) {
            return;
        }
        match self.kind {
            InteractionKind::Clash { k2, atom3, sq_dist } => {
                if sq_dist > self.rest_v {
                    return;
                }
                if eq(self.atom1, atom3) || eq(self.atom2, atom3) {
                    return;
                }
                let mut atom1 = self.atom1.borrow_mut();
                let mut atom2 = self.atom2.borrow_mut();
                let mut atom3 = atom3.borrow_mut();
                let proj =
                    math::project_on_line(atom2.coordinates, atom1.coordinates, atom3.coordinates);
                let mut f = atom2.coordinates - proj;
                f.normalize();
                f *= (self.rest_v - sq_dist) * self.k * k2;
                atom2.force += f;
                atom1.force -= f * 0.5;
                atom3.force -= f * 0.5;
            }
            InteractionKind::Stretch { tolerance } => {
                let mut atom1 = self.atom1.borrow_mut();
                let mut atom2 = self.atom2.borrow_mut();
                let mut l = atom2.coordinates - atom1.coordinates;
                let m = l.length();
                let dr = if m < self.rest_v - tolerance {
                    self.rest_v - tolerance - m
                } else if m > self.rest_v + tolerance {
                    self.rest_v + tolerance - m
                } else {
                    return;
                };
                let short = self.rest_v * 0.4;
                let short_penalty = (short - m).min(0.0);
                if m > EPSILON {
                    l /= m;
                }
                l *= self.k * dr + short_penalty * 10.0;
                atom1.force += l;
                atom2.force -= l;
            }
            InteractionKind::Bend { k2, atom3, .. } => {
                let mut atom1 = self.atom1.borrow_mut();
                let mut atom2 = self.atom2.borrow_mut();
                let mut atom3 = atom3.borrow_mut();
                let p1 = atom1.coordinates;
                let p2 = atom2.coordinates;
                let p3 = atom3.coordinates;
                let a = math::unsigned_angle(p1, p2, p3);
                let target = self.rest_v.min(360.0 - self.rest_v);
                let da = target - a;
                let v1 = p1 - p2;
                let v2 = p3 - p2;
                let v3 = p3 - p1;
                let mut n1 = PointF(v1.1, -v1.0);
                let mut n2 = PointF(v2.1, -v2.0);
                if n1.dot(v3) > 0.0 {
                    n1 *= -1.0;
                }
                if n2.dot(v3) > 0.0 {
                    n2 *= -1.0;
                }
                n1.normalize();
                n2.normalize();
                n1 *= self.k * k2 * da;
                n2 *= self.k * k2 * da;
                atom1.force += n1;
                atom3.force += n2;
                atom2.force -= n1 + n2;
            }
            InteractionKind::Inversion {
                atom3,
                atom4,
                is_z,
                force_move,
            } => {
                if eq(self.atom1, atom3)
                    || eq(self.atom1, atom4)
                    || eq(self.atom2, atom3)
                    || eq(self.atom2, atom4)
                    || eq(atom3, atom4)
                {
                    return;
                }
                let mut atom1 = self.atom1.borrow_mut();
                let mut atom2 = self.atom2.borrow_mut();
                let mut atom3 = atom3.borrow_mut();
                let mut atom4 = atom4.borrow_mut();
                let p1 = atom1.coordinates;
                let p2 = atom2.coordinates;
                let p3 = atom3.coordinates;
                let p4 = atom4.coordinates;
                if math::same_side(p1, p4, p2, p3) == is_z {
                    return;
                }
                let proj1 = math::project_on_line(p1, p2, p3);
                let proj2 = math::project_on_line(p4, p1, p2);
                let mut side = &mut atom1;
                let mut bonded = &mut atom2;
                let mut proj = proj1;
                if (p1 - proj1).sq_length() > (p4 - proj2).sq_length() {
                    side = &mut atom4;
                    bonded = &mut atom3;
                    proj = proj2;
                }
                let mut force = proj - side.coordinates;
                if force_move {
                    side.coordinates += force;
                    bonded.coordinates -= force;
                    side.force = PointF::default();
                    bonded.force = PointF::default();
                } else {
                    force.normalize();
                    force *= 10.0;
                    side.force += force;
                    bonded.force -= force;
                }
            }
        }
    }
}

pub type InteractionRef<'a> = &'a RefCell<Interaction<'a>>;

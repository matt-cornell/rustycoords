use super::*;
use std::cell::RefCell;
use std::ptr::eq;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Stereo {
    Cis,
    Trans,
    #[default]
    Unspecified,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct StereoInfo<'a> {
    pub atom1: Option<AtomRef<'a>>,
    pub atom2: Option<AtomRef<'a>>,
    pub stereo: Stereo,
}

#[derive(Debug, Clone)]
pub struct Bond<'a> {
    pub start: AtomRef<'a>,
    pub end: AtomRef<'a>,
    pub bond_order: u8,
    pub skip: bool,
    pub stereo: StereoInfo<'a>,
    pub(crate) is_reversed: bool,
    pub(crate) is_ze_active: bool,
    pub(crate) is_z: bool,
    pub(crate) ignore_ze: bool,
    pub(crate) sssr_visited: bool,
    pub(crate) rings: Vec<RingRef<'a>>,
    pub(crate) sssr_parent: Option<BondRef<'a>>,
}
impl<'a> Bond<'a> {
    pub const fn new(start: AtomRef<'a>, end: AtomRef<'a>) -> Self {
        Self {
            start,
            end,
            bond_order: 1,
            is_reversed: false,
            skip: false,
            is_ze_active: false,
            is_z: false,
            ignore_ze: false,
            sssr_visited: false,
            stereo: StereoInfo {
                atom1: None,
                atom2: None,
                stereo: Stereo::Unspecified,
            },
            rings: Vec::new(),
            sssr_parent: None,
        }
    }
    pub(crate) fn set_absolute_stereo(&mut self) {
        if self.is_stereo() {
            let start = self.start_first_cip_neighbor().unwrap();
            let end = self.end_first_cip_neighbor().unwrap();
            let mut inv = false;
            if self
                .stereo
                .atom1
                .map_or(true, |a| !eq(a, start) && !eq(a, end))
            {
                inv = true;
            }
            if self
                .stereo
                .atom2
                .map_or(true, |a| !eq(a, start) && !eq(a, end))
            {
                inv = !inv;
            }
            self.is_z = (self.stereo.stereo == Stereo::Cis) ^ inv;
        }
        if self.stereo.stereo == Stereo::Unspecified {
            self.ignore_ze = true;
        }
    }
    pub(crate) fn check_stereo(&self) -> bool {
        if self.is_stereo() {
            return true;
        }
        if self.is_in_small_ring() {
            return true;
        }
        let Some(start) = self.start_first_cip_neighbor() else {
            return true;
        };
        let Some(end) = self.start_first_cip_neighbor() else {
            return true;
        };
        let res = math::same_side(
            start.borrow().coordinates,
            end.borrow().coordinates,
            self.start.borrow().coordinates,
            self.end.borrow().coordinates,
        );
        res == self.is_z
    }
    pub fn is_in_small_ring(&self) -> bool {
        self.rings.iter().any(|r| !r.borrow().is_macrocycle())
    }
    pub fn is_in_macrocycle(&self) -> bool {
        self.rings.iter().any(|r| r.borrow().is_macrocycle())
    }
    pub fn is_terminal(&self) -> bool {
        self.start.borrow().bonds.len() == 1 || self.end.borrow().bonds.len() == 1
    }
    pub(crate) fn is_inter_fragment(&self) -> bool {
        if Atom::share_a_ring(self.start, self.end).is_some() {
            return false;
        };
        !self.is_stereo()
    }
    pub fn is_stereo(&self) -> bool {
        if self.bond_order != 2 {
            return false;
        }
        if self.ignore_ze {
            return false;
        }
        let ring = Atom::share_a_ring(self.start, self.end);
        ring.map_or(true, |r| r.borrow().is_macrocycle())
    }
    pub(crate) fn marked_as_cis(&self, a1: &Atom<'a>, a2: &Atom<'a>) -> bool {
        let Some(n1) = self.start_first_cip_neighbor() else {
            return false;
        };
        let Some(n2) = self.end_first_cip_neighbor() else {
            return false;
        };
        let n1 = n1.borrow();
        let n2 = n2.borrow();
        let mut cis = self.is_z;
        if eq(&*n1, a1) && eq(&*n2, a2) {
            cis = !cis;
        }
        if eq(&*n1, a2) && eq(&*n2, a1) {
            cis = !cis;
        }
        cis
    }
    pub(crate) fn flip(&mut self) {
        todo!()
    }
    pub(crate) fn start_first_cip_neighbor(&self) -> Option<AtomRef<'a>> {
        if self.bond_order < 2 {
            return None;
        }
        let a = self.start.borrow();
        match a.neighbors.len() {
            2 => {
                if eq(a.neighbors[0], self.end) {
                    Some(a.neighbors[1])
                } else {
                    Some(a.neighbors[0])
                }
            }
            3 => {
                let mut ats = a.neighbors.iter().filter(|&&i| !eq(i, self.end));
                let first = ats.next()?;
                let second = ats.next()?;
                if ats.next().is_some() {
                    None
                } else {
                    Some(Atom::cip_priority(first, second, self.start))
                }
            }
            _ => None,
        }
    }
    pub(crate) fn end_first_cip_neighbor(&self) -> Option<AtomRef<'a>> {
        if self.bond_order < 2 {
            return None;
        }
        let a = self.end.borrow();
        match a.neighbors.len() {
            2 => {
                if eq(a.neighbors[0], self.start) {
                    Some(a.neighbors[1])
                } else {
                    Some(a.neighbors[0])
                }
            }
            3 => {
                let mut ats = a.neighbors.iter().filter(|&&i| !eq(i, self.start));
                let first = ats.next()?;
                let second = ats.next()?;
                if ats.next().is_some() {
                    None
                } else {
                    Some(Atom::cip_priority(first, second, self.end))
                }
            }
            _ => None,
        }
    }
}

pub type BondRef<'a> = &'a RefCell<Bond<'a>>;

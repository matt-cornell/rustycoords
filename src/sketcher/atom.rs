use super::*;
use std::cell::RefCell;

#[derive(Debug)]
pub struct Atom<'a> {
    pub cross_layout: bool,
    pub fixed: bool,
    pub constrained: bool,
    pub rigid: bool,
    pub is_shared_and_inner: bool,
    pub atom_number: u8,
    pub charge: i8,
    pub(crate) valence: u8,
    pub(crate) general_use_n: usize,
    pub(crate) general_use_n2: usize,
    pub(crate) chm_n: usize,
    pub(crate) general_use_visited: bool,
    pub(crate) general_use_visited2: bool,
    pub(crate) clockwise_invert: bool,
    pub(crate) ignore_ring_chirality: bool,
    pub(crate) rs_priorities: Vec<u8>,
    pub(crate) implicit_h: u8,
    pub(crate) molecule: Option<MoleculeRef<'a>>,
    pub(crate) fragment: Option<FragmentRef<'a>>,
    pub(crate) neighbors: Vec<AtomRef<'a>>,
    pub(crate) bonds: Vec<BondRef<'a>>,
    pub(crate) rings: Vec<RingRef<'a>>,
    pub coordinates: PointF,
    pub template: PointF,
    pub force: PointF,
    pub coordinates_set: bool,
}
impl Default for Atom<'_> {
    fn default() -> Self {
        Self {
            cross_layout: false,
            fixed: false,
            constrained: false,
            rigid: false,
            is_shared_and_inner: false,
            atom_number: 0,
            charge: 0,
            valence: u8::MAX,
            general_use_n: 0,
            general_use_n2: 0,
            chm_n: 0,
            general_use_visited: false,
            general_use_visited2: false,
            clockwise_invert: false,
            ignore_ring_chirality: false,
            rs_priorities: Vec::new(),
            implicit_h: u8::MAX,
            molecule: None,
            fragment: None,
            neighbors: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
            coordinates: PointF::default(),
            template: PointF::default(),
            force: PointF::default(),
            coordinates_set: false,
        }
    }
}
impl<'a> Atom<'a> {
    pub fn new(atom_number: u8) -> Self {
        Self {
            atom_number,
            ..Self::default()
        }
    }
    pub(crate) fn cip_priority(
        a1: AtomRef<'a>,
        a2: AtomRef<'a>,
        start: AtomRef<'a>,
    ) -> AtomRef<'a> {
        todo!()
    }
    pub(crate) fn shares_a_ring_with(&self, a2: &Atom<'a>) -> Option<RingRef<'a>> {
        if self.rings.is_empty() {
            return None;
        }
        if a2.rings.is_empty() {
            return None;
        }
        let mut r: Option<RingRef<'a>> = None;
        for &r1 in &self.rings {
            let ring = r1.borrow();
            if ring.is_macrocycle() {
                continue;
            }
            for &r2 in &a2.rings {
                if !std::ptr::eq(r1, r2) {
                    continue;
                }
                if let Some(old) = r {
                    if ring.atoms.len() > old.borrow().atoms.len() {
                        r = Some(r1);
                    }
                } else {
                    r = Some(r1);
                }
            }
        }
        for &r1 in &self.rings {
            let ring = r1.borrow();
            for &r2 in &a2.rings {
                if !std::ptr::eq(r1, r2) {
                    continue;
                }
                if let Some(old) = r {
                    if ring.atoms.len() > old.borrow().atoms.len() {
                        r = Some(r1);
                    }
                } else {
                    r = Some(r1);
                }
            }
        }
        r
    }
    pub(crate) fn shares_a_ring_with_both(
        &self,
        a2: &Atom<'a>,
        a3: &Atom<'a>,
    ) -> Option<RingRef<'a>> {
        if self.rings.is_empty() {
            return None;
        }
        if a2.rings.is_empty() {
            return None;
        }
        if a3.rings.is_empty() {
            return None;
        }
        let mut r: Option<RingRef<'a>> = None;
        for &r1 in &self.rings {
            let ring = r1.borrow();
            if ring.is_macrocycle() {
                continue;
            }
            for &r2 in &a2.rings {
                if !std::ptr::eq(r1, r2) {
                    continue;
                }
                for &r3 in &a3.rings {
                    if !std::ptr::eq(r1, r3) {
                        continue;
                    }
                    if let Some(old) = r {
                        if ring.atoms.len() > old.borrow().atoms.len() {
                            r = Some(r1);
                        }
                    } else {
                        r = Some(r1);
                    }
                }
            }
        }
        for &r1 in &self.rings {
            let ring = r1.borrow();
            for &r2 in &a2.rings {
                if !std::ptr::eq(r1, r2) {
                    continue;
                }
                for &r3 in &a3.rings {
                    if !std::ptr::eq(r1, r3) {
                        continue;
                    }
                    if let Some(old) = r {
                        if ring.atoms.len() > old.borrow().atoms.len() {
                            r = Some(r1);
                        }
                    } else {
                        r = Some(r1);
                    }
                }
            }
        }
        r
    }
    pub(crate) fn clockwise_ordered_neighbors(&self, out: &mut Vec<AtomRef<'a>>) {
        out.clone_from(&self.neighbors); // TODO: actually order this
    }
    pub fn is_metal(atom_number: u8) -> bool {
        [3, 4, 11, 12, 13, 31, 49, 32, 50, 51].contains(&atom_number)
            || [19..=30, 37..=48, 55..=84, 87..=112]
                .iter()
                .any(|r| r.contains(&atom_number))
    }
    pub fn ex_valence(atom_number: u8) -> u8 {
        match atom_number {
            1 | 9 | 17 | 35 => 1,
            5 => 3,
            6 | 14 => 4,
            7 | 15 => 3,
            8 | 16 | 34 | 53 => 2,
            _ => 0,
        }
    }
    pub fn find_h_number(&self) -> u8 {
        let mut valence = self.valence as i16;
        if valence == 255 {
            valence = Self::ex_valence(self.atom_number) as i16;
        }
        let bond_orders: i16 = self
            .bonds
            .iter()
            .map(|b| b.borrow().bond_order as i16)
            .sum();
        if self.atom_number == 15 || self.atom_number == 16 {
            let mut it = self
                .neighbors
                .iter()
                .zip(&self.bonds)
                .filter(|(n, b)| n.borrow().atom_number == 8 && b.borrow().bond_order == 2);
            if it.next().is_some() {
                if it.next().is_some() {
                    if it.next().is_none() {
                        valence += 2;
                    }
                } else {
                    valence += 2;
                }
            }
        }
        (valence - bond_orders + self.charge as i16).clamp(0, 4) as u8
    }
    pub fn set_coords(&mut self, coords: PointF) {
        self.coordinates = coords;
        self.coordinates.round(2);
        self.coordinates_set = true;
    }
    #[inline(always)]
    pub fn set_coords_to_template(&mut self) {
        self.set_coords(self.template);
    }
}

pub type AtomRef<'a> = &'a RefCell<Atom<'a>>;

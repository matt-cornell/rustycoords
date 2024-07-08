use super::*;
use std::cell::RefCell;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Molecule<'a> {
    pub atoms: Vec<AtomRef<'a>>,
    pub bonds: Vec<BondRef<'a>>,
    pub(crate) main_fragment: Option<FragmentRef<'a>>,
    pub(crate) rings: Vec<RingRef<'a>>,
    pub(crate) fragments: Vec<FragmentRef<'a>>,
    pub(crate) proximity_relations: Vec<BondRef<'a>>,
    pub(crate) has_fixed_frags: bool,
    pub(crate) has_constrained_frags: bool,
}
impl<'a> Molecule<'a> {
    pub fn force_update_struct(&mut self) {
        for a in &self.atoms {
            let mut a = a.borrow_mut();
            a.bonds.clear();
            a.neighbors.clear();
            a.rings.clear();
        }
        for bond in &self.bonds {
            let b = bond.borrow();
            let mut start = b.start.borrow_mut();
            let mut end = b.end.borrow_mut();
            start.bonds.push(bond);
            end.bonds.push(bond);
            start.neighbors.push(b.end);
            end.neighbors.push(b.start);
        }
        for a in &self.atoms {
            let mut a = a.borrow_mut();
            if a.implicit_h == u8::MAX {
                a.implicit_h = a.find_h_number();
            }
        }
    }
}
pub type MoleculeRef<'a> = &'a RefCell<Molecule<'a>>;

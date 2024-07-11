use super::*;

#[derive(Debug, Default)]
pub struct Molecule<'a> {
    pub(crate) atoms: Vec<AtomRef<'a>>,
    pub(crate) bonds: Vec<BondRef<'a>>,
    pub(crate) main_fragment: Option<FragmentRef<'a>>,
    pub(crate) rings: Vec<RingRef<'a>>,
    pub(crate) fragments: Vec<FragmentRef<'a>>,
    pub(crate) proximity_relations: Vec<BondRef<'a>>,
    pub(crate) has_fixed_frags: bool,
    pub(crate) has_constrained_frags: bool,
}
impl<'a> Molecule<'a> {
    pub fn add_atom(&mut self, atom: Atom<'a>, intern: &'a dyn Interner) -> AtomRef<'a> {
        let a = intern.intern_atom(atom);
        self.atoms.push(a);
        a
    }
    pub fn add_new_atom(&mut self, atom_number: u8, intern: &'a dyn Interner) -> AtomRef<'a> {
        let a = intern.intern_atom(Atom::new(atom_number));
        self.atoms.push(a);
        a
    }
    pub fn add_bond(&mut self, bond: Bond<'a>, intern: &'a dyn Interner) -> BondRef<'a> {
        let b = intern.intern_bond(bond);
        self.bonds.push(b);
        b
    }
    pub fn add_new_bond(&mut self, start: AtomRef<'a>, end: AtomRef<'a>, order: u8, intern: &'a dyn Interner) -> BondRef<'a> {
        let mut bond = Bond::new(start, end);
        bond.bond_order = order;
        let b = intern.intern_bond(bond);
        self.bonds.push(b);
        b
    }

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

use super::*;

pub struct Builder<'a> {
    pub intern: &'a dyn Interner,
    mol: Molecule<'a>,
}
impl<'a> Builder<'a> {
    pub fn new(intern: &'a dyn Interner) -> Self {
        Self {
            intern,
            mol: Molecule::default(),
        }
    }
    pub fn add_atom(&mut self, atom_number: u8) -> AtomRef<'a> {
        self.mol.add_new_atom(atom_number, self.intern)
    }
    pub fn add_bond(&mut self, start: AtomRef<'a>, end: AtomRef<'a>, order: u8) -> BondRef<'a> {
        self.mol.add_new_bond(start, end, order, self.intern)
    }
    pub fn finish(self) -> MoleculeRef<'a> {
        self.intern.intern_molecule(self.mol)
    }
}

use std::cell::RefCell;

use super::*;

pub trait Interner {
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a>;
}
pub struct Leak;
impl Interner for Leak {
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a> {
        Box::leak(Box::new(RefCell::new(mol)))
    }
}

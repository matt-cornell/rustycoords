use std::cell::RefCell;

use super::*;

pub trait Interner {
    fn intern_fragment<'a>(&'a self, frg: Fragment<'a>) -> FragmentRef<'a>;
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a>;
    fn intern_frag_dof<'a>(&'a self, dof: FragmentDof<'a>) -> FragmentDofRef<'a>;
}

/// Interner that just leaks everything, not recommended for obvious reasons
pub struct Leak;
impl Interner for Leak {
    fn intern_fragment<'a>(&'a self, frg: Fragment<'a>) -> FragmentRef<'a> {
        Box::leak(Box::new(RefCell::new(frg)))
    }
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a> {
        Box::leak(Box::new(RefCell::new(mol)))
    }
    fn intern_frag_dof<'a>(&'a self, dof: FragmentDof<'a>) -> FragmentDofRef<'a> {
        Box::leak(Box::new(RefCell::new(dof)))
    }
}

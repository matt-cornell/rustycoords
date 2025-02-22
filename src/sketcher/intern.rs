use super::*;
use std::cell::RefCell;
use std::mem::transmute;

pub(crate) trait InternerImpl {
    fn intern_atom<'a>(&'a self, atom: Atom<'a>) -> AtomRef<'a>;
    fn intern_bond<'a>(&'a self, bond: Bond<'a>) -> BondRef<'a>;
    fn intern_ring<'a>(&'a self, ring: Ring<'a>) -> RingRef<'a>;
    fn intern_fragment<'a>(&'a self, frg: Fragment<'a>) -> FragmentRef<'a>;
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a>;
    fn intern_frag_dof<'a>(&'a self, dof: FragmentDof<'a>) -> FragmentDofRef<'a>;
    fn intern_interaction<'a>(&'a self, int: Interaction<'a>) -> InteractionRef<'a>;
}

#[allow(private_bounds)]
pub trait Interner: InternerImpl {}
impl<T: InternerImpl> Interner for T {}

/// Interner that just leaks everything, not recommended for obvious reasons
pub struct Leak;
impl InternerImpl for Leak {
    fn intern_atom<'a>(&'a self, atom: Atom<'a>) -> AtomRef<'a> {
        Box::leak(Box::new(RefCell::new(atom)))
    }
    fn intern_bond<'a>(&'a self, bond: Bond<'a>) -> BondRef<'a> {
        Box::leak(Box::new(RefCell::new(bond)))
    }
    fn intern_ring<'a>(&'a self, ring: Ring<'a>) -> RingRef<'a> {
        Box::leak(Box::new(RefCell::new(ring)))
    }
    fn intern_fragment<'a>(&'a self, frg: Fragment<'a>) -> FragmentRef<'a> {
        Box::leak(Box::new(RefCell::new(frg)))
    }
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a> {
        Box::leak(Box::new(RefCell::new(mol)))
    }
    fn intern_frag_dof<'a>(&'a self, dof: FragmentDof<'a>) -> FragmentDofRef<'a> {
        Box::leak(Box::new(RefCell::new(dof)))
    }
    fn intern_interaction<'a>(&'a self, int: Interaction<'a>) -> InteractionRef<'a> {
        Box::leak(Box::new(RefCell::new(int)))
    }
}

#[derive(Debug)]
enum InternedKind<'a> {
    Atom(RefCell<Atom<'a>>),
    Bond(RefCell<Bond<'a>>),
    Ring(RefCell<Ring<'a>>),
    Fragment(RefCell<Fragment<'a>>),
    Molecule(RefCell<Molecule<'a>>),
    FragmentDof(RefCell<FragmentDof<'a>>),
    Interaction(RefCell<Interaction<'a>>),
}

#[derive(Debug, Default)]
pub struct Intern {
    inner: orx_imp_vec::ImpVec<InternedKind<'static>>,
}
impl Intern {
    pub fn new() -> Self {
        Self::default()
    }
}
impl InternerImpl for Intern {
    fn intern_atom<'a>(&'a self, atom: Atom<'a>) -> AtomRef<'a> {
        unsafe {
            let atom = transmute::<Atom<'a>, Atom<'static>>(atom);
            self.inner.imp_push(InternedKind::Atom(RefCell::new(atom)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::Atom(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<Atom<'static>>, AtomRef<'a>>(out)
        }
    }
    fn intern_bond<'a>(&'a self, bond: Bond<'a>) -> BondRef<'a> {
        unsafe {
            let bond = transmute::<Bond<'a>, Bond<'static>>(bond);
            self.inner.imp_push(InternedKind::Bond(RefCell::new(bond)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::Bond(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<Bond<'static>>, BondRef<'a>>(out)
        }
    }
    fn intern_ring<'a>(&'a self, ring: Ring<'a>) -> RingRef<'a> {
        unsafe {
            let ring = transmute::<Ring<'a>, Ring<'static>>(ring);
            self.inner.imp_push(InternedKind::Ring(RefCell::new(ring)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::Ring(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<Ring<'static>>, RingRef<'a>>(out)
        }
    }
    fn intern_fragment<'a>(&'a self, frg: Fragment<'a>) -> FragmentRef<'a> {
        unsafe {
            let frg = transmute::<Fragment<'a>, Fragment<'static>>(frg);
            self.inner
                .imp_push(InternedKind::Fragment(RefCell::new(frg)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::Fragment(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<Fragment<'static>>, FragmentRef<'a>>(out)
        }
    }
    fn intern_molecule<'a>(&'a self, mol: Molecule<'a>) -> MoleculeRef<'a> {
        unsafe {
            let mol = transmute::<Molecule<'a>, Molecule<'static>>(mol);
            self.inner
                .imp_push(InternedKind::Molecule(RefCell::new(mol)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::Molecule(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<Molecule<'static>>, MoleculeRef<'a>>(out)
        }
    }
    fn intern_frag_dof<'a>(&'a self, dof: FragmentDof<'a>) -> FragmentDofRef<'a> {
        unsafe {
            let dof = transmute::<FragmentDof<'a>, FragmentDof<'static>>(dof);
            self.inner
                .imp_push(InternedKind::FragmentDof(RefCell::new(dof)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::FragmentDof(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<FragmentDof<'static>>, FragmentDofRef<'a>>(out)
        }
    }
    fn intern_interaction<'a>(&'a self, int: Interaction<'a>) -> InteractionRef<'a> {
        unsafe {
            let int = transmute::<Interaction<'a>, Interaction<'static>>(int);
            self.inner
                .imp_push(InternedKind::Interaction(RefCell::new(int)));
            let last = self
                .inner
                .fragments()
                .iter()
                .rev()
                .find_map(|f| f.last())
                .unwrap();
            let InternedKind::Interaction(out) = last else {
                unreachable!()
            };
            transmute::<&'a RefCell<Interaction<'static>>, InteractionRef<'a>>(out)
        }
    }
}

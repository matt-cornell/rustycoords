use super::*;
use ahash::AHashMap;
use std::cell::RefCell;

pub trait FragmentDof {
    fn num_states(&self) -> u8;
    fn tier(&self) -> u8;
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Fragment<'a> {
    pub atoms: Vec<AtomRef<'a>>,
    pub bonds: Vec<BondRef<'a>>,
    pub rings: Vec<RingRef<'a>>,
    pub inter_frags: Vec<BondRef<'a>>,
    pub children: Vec<FragmentRef<'a>>,
    pub parent: Option<(FragmentRef<'a>, BondRef<'a>)>,
    pub num_children: usize,
    pub child_rank: f32,
    pub longest_chain: f32,
    pub is_chain: bool,
    pub fixed: bool,
    pub constrained: bool,
    pub constrained_flip: bool,
    pub coords: AHashMap<*const RefCell<Atom<'a>>, PointF>,
}
impl<'a> Fragment<'a> {
    pub fn count_fixed(&self) -> usize {
        if !self.fixed {
            return 0;
        }
        self.atoms.iter().filter(|a| a.borrow().fixed).count()
    }
    pub fn count_constrained(&self) -> usize {
        if !self.constrained {
            return 0;
        }
        self.atoms.iter().filter(|a| a.borrow().constrained).count()
    }
    pub fn count_heavy(&self) -> usize {
        self.atoms
            .iter()
            .filter(|a| a.borrow().atom_number != 6)
            .count()
    }
    pub fn total_weight(&self) -> usize {
        self.atoms
            .iter()
            .map(|a| {
                let a = a.borrow();
                (a.atom_number + a.implicit_h) as usize
            })
            .sum()
    }
    pub fn count_double_bonds(&self) -> usize {
        self.bonds
            .iter()
            .filter(|b| b.borrow().bond_order == 2)
            .count()
    }
}
pub type FragmentRef<'a> = &'a RefCell<Fragment<'a>>;

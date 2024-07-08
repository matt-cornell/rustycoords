use super::*;
use ahash::AHashMap;
use std::cell::RefCell;

pub trait FragmentDof {
    fn num_states(&self) -> u8;
    fn tier(&self) -> u8;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Fragment<'a> {
    pub atoms: Vec<AtomRef<'a>>,
    pub rings: Vec<RingRef<'a>>,
    pub children: Vec<FragmentRef<'a>>,
    pub parent: Option<(FragmentRef<'a>, BondRef<'a>)>,
    pub num_children: usize,
    pub child_rank: f32,
    pub longest_chain: f32,
    pub coords: AHashMap<*const RefCell<Atom<'a>>, PointF>,
}
pub type FragmentRef<'a> = &'a RefCell<Fragment<'a>>;

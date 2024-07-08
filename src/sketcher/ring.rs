use super::*;
use std::cell::RefCell;

#[derive(Debug, Clone, PartialEq)]
pub struct Ring<'a> {
    pub atoms: Vec<AtomRef<'a>>,
    pub fused_with: Vec<(RingRef<'a>, Vec<AtomRef<'a>>)>,
}
impl<'a> Ring<'a> {
    pub fn is_macrocycle(&self) -> bool {
        todo!()
    }
    pub fn find_center(&self) -> PointF {
        todo!()
    }
}

pub type RingRef<'a> = &'a RefCell<Ring<'a>>;

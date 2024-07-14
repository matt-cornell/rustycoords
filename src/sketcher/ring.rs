use super::*;

#[derive(Debug, Clone)]
pub struct Ring<'a> {
    pub atoms: Vec<AtomRef<'a>>,
    pub fused_with: Vec<(RingRef<'a>, Vec<AtomRef<'a>>)>,
    pub fusion_bonds: Vec<BondRef<'a>>,
}
impl<'a> Ring<'a> {
    pub fn is_macrocycle(&self) -> bool {
        todo!()
    }
    pub fn find_center(&self) -> PointF {
        self.atoms
            .iter()
            .map(|a| a.borrow().coordinates)
            .sum::<PointF>()
            / self.atoms.len() as f32
    }
}

pub type RingRef<'a> = &'a RefCell<Ring<'a>>;

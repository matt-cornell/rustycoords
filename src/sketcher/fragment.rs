use super::*;
use ahash::AHashMap;

#[derive(Debug, Default)]
pub struct Fragment<'a> {
    pub atoms: Vec<AtomRef<'a>>,
    pub bonds: Vec<BondRef<'a>>,
    pub rings: Vec<RingRef<'a>>,
    pub inter_frags: Vec<BondRef<'a>>,
    pub children: Vec<FragmentRef<'a>>,
    pub dofs: Vec<FragmentDofRef<'a>>,
    pub parent: Option<(FragmentRef<'a>, BondRef<'a>)>,
    pub num_children: usize,
    pub child_rank: f32,
    pub longest_chain: f32,
    pub is_chain: bool,
    pub fixed: bool,
    pub constrained: bool,
    pub constrained_flip: bool,
    pub templated: bool,
    pub coords: AHashMap<*const RefCell<Atom<'a>>, PointF>,
    pub atom_dofs: AHashMap<*const RefCell<Atom<'a>>, Vec<FragmentDofRef<'a>>>,
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
    pub(crate) fn store_coordinate_information(&mut self) {
        self.coords.clear();
        let origin;
        let angle;
        if let Some((_, bond)) = self.parent {
            let bond = bond.borrow();
            let start = bond.start.borrow().coordinates;
            origin = bond.end.borrow().coordinates;
            let p = start - origin;
            angle = p.1.atan2(-p.0);
        } else {
            angle = 0.0;
            if self.constrained || self.fixed {
                origin = PointF::default();
            } else {
                origin = self.atoms[0].borrow().coordinates;
            }
        }
        let (sin, cos) = angle.sin_cos();
        for a in self.atoms.iter().copied().chain(
            self.children
                .iter()
                .map(|c| c.borrow().parent.unwrap().1.borrow().end),
        ) {
            let mut c = a.borrow().coordinates - origin;
            c.rotate(-sin, cos);
            self.coords.insert(a as *const _, c);
        }
    }
}
pub type FragmentRef<'a> = &'a RefCell<Fragment<'a>>;

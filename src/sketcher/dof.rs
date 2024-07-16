use super::*;

#[derive(Debug, Clone, Copy)]
pub enum FragmentDofKind<'a> {
    RotateFrag,
    FlipFrag,
    ScaleAtom(AtomRef<'a>),
    #[allow(dead_code)]
    ScaleFrag,
    ChangeParentBond,
    InvertBond {
        pivot: AtomRef<'a>,
        bound: AtomRef<'a>,
    },
    #[allow(dead_code)]
    FlipRing {
        pivot1: AtomRef<'a>,
        pivot2: AtomRef<'a>,
        penalty: f32,
    },
}

#[derive(Debug)]
pub struct FragmentDof<'a> {
    pub kind: FragmentDofKind<'a>,
    pub atoms: Vec<AtomRef<'a>>,
    pub frag: FragmentRef<'a>,
    pub current_state: usize,
}
impl<'a> FragmentDof<'a> {
    pub fn num_states(&self) -> usize {
        match self.kind {
            FragmentDofKind::RotateFrag => {
                if self.frag.borrow().parent.is_some() {
                    5
                } else {
                    1
                }
            }
            FragmentDofKind::FlipFrag => {
                if self.frag.borrow().parent.is_some() {
                    2
                } else {
                    1
                }
            }
            FragmentDofKind::ScaleAtom(_) => 2,
            FragmentDofKind::ScaleFrag => {
                if self.frag.borrow().rings.is_empty() {
                    1
                } else {
                    5
                }
            }
            FragmentDofKind::ChangeParentBond { .. } => 7,
            FragmentDofKind::InvertBond { .. } => 2,
            FragmentDofKind::FlipRing { .. } => 2,
        }
    }
    pub fn change_state(&mut self) {
        self.current_state += 1;
        self.current_state %= self.num_states();
    }
    pub fn current_penalty(&self, frag: &Fragment<'a>) -> f32 {
        match self.kind {
            FragmentDofKind::RotateFrag => {
                if self.current_state == 0 {
                    0.0
                } else {
                    (self.current_state as f32 + 1.0) * 200.0
                }
            }
            FragmentDofKind::FlipFrag => {
                let mut penalty = 0.0;
                if self.current_state != 0 && frag.constrained_flip {
                    penalty += 1000.0;
                }
                if frag.is_chain && frag.parent.map_or(false, |p| p.0.borrow().is_chain) {
                    penalty += 10.0;
                }
                penalty
            }
            FragmentDofKind::ScaleAtom(_) => {
                if self.current_state == 0 {
                    0.0
                } else {
                    50.0 * self.atoms.len() as f32
                }
            }
            FragmentDofKind::ScaleFrag => {
                if self.current_state == 0 {
                    0.0
                } else {
                    (self.current_state as f32 + 1.0) * 250.0
                }
            }
            FragmentDofKind::ChangeParentBond => {
                if self.current_state == 0 {
                    0.0
                } else {
                    (self.current_state as f32 + 1.0) * 100.0
                }
            }
            FragmentDofKind::InvertBond { .. } => {
                if self.current_state == 0 {
                    0.0
                } else {
                    100.0
                }
            }
            FragmentDofKind::FlipRing { penalty, .. } => {
                if self.current_state == 0 {
                    0.0
                } else {
                    penalty * 200.0
                }
            }
        }
    }
    pub fn tier(&self) -> u8 {
        match self.kind {
            FragmentDofKind::RotateFrag => 3,
            FragmentDofKind::FlipFrag => 0,
            FragmentDofKind::ScaleFrag => 5,
            FragmentDofKind::ScaleAtom(_) => 4,
            FragmentDofKind::ChangeParentBond => 2,
            FragmentDofKind::InvertBond { .. } => 1,
            FragmentDofKind::FlipRing { .. } => 1,
        }
    }
    pub fn apply(&self, frag: &Fragment<'a>) {
        if self.current_state != 0 {
            match self.kind {
                FragmentDofKind::RotateFrag => {
                    let mut angle = (((self.current_state + 1) / 2) as f32 * 15.0).to_radians();
                    if self.current_state & 1 == 0 {
                        angle = -angle;
                    }
                    let (sin, cos) = angle.sin_cos();
                    let origin = PointF(-BOND_LENGTH, 0.0);
                    for atom in frag.coords.keys() {
                        let mut atom = unsafe { (**atom).borrow_mut() };
                        let mut coords = atom.coordinates - origin;
                        coords.rotate(sin, cos);
                        atom.set_coords(coords + origin);
                    }
                }
                FragmentDofKind::FlipFrag => {
                    for atom in frag.coords.keys() {
                        unsafe { (**atom).borrow_mut() }.coordinates *= -1.0;
                    }
                }
                FragmentDofKind::ScaleFrag => {
                    let mut scale = 1.4f32.powi((self.current_state as i32 + 1) / 2);
                    if self.current_state & 1 == 0 {
                        scale = scale.recip();
                    }
                    for atom in frag.coords.keys() {
                        unsafe { (**atom).borrow_mut() }.coordinates *= scale;
                    }
                }
                FragmentDofKind::ScaleAtom(pivot) => {
                    let pivot = pivot.borrow().coordinates;
                    for atom in &self.atoms {
                        let mut atom = atom.borrow_mut();
                        let dist = atom.coordinates - pivot;
                        atom.set_coords(dist * 0.4 + pivot);
                    }
                }
                FragmentDofKind::ChangeParentBond => {
                    let mut scale = 1.6f32.powi((self.current_state as i32 + 1) / 2);
                    if self.current_state & 1 == 0 {
                        scale = scale.recip();
                    }
                    let move_by = BOND_LENGTH * (scale - 1.0);
                    for atom in frag.coords.keys() {
                        unsafe { (**atom).borrow_mut() }.coordinates.0 += move_by;
                    }
                }
                FragmentDofKind::InvertBond { pivot, bound } => {
                    let pivot = pivot.borrow().coordinates;
                    let bond_dir = bound.borrow().coordinates - pivot;
                    let normal = PointF(bond_dir.1, -bond_dir.0);
                    let point1 = pivot + normal;
                    let point2 = pivot - normal;
                    for atom in &self.atoms {
                        let mut atom = atom.borrow_mut();
                        atom.coordinates = math::mirror_point(atom.coordinates, point1, point2);
                    }
                }
                FragmentDofKind::FlipRing { pivot1, pivot2, .. } => {
                    let point1 = pivot1.borrow().coordinates;
                    let point2 = pivot2.borrow().coordinates;
                    for atom in &self.atoms {
                        let mut atom = atom.borrow_mut();
                        atom.coordinates = math::mirror_point(atom.coordinates, point1, point2);
                    }
                }
            }
        }
    }
}

pub type FragmentDofRef<'a> = &'a RefCell<FragmentDof<'a>>;

use super::*;

#[derive(Debug, Clone, Copy)]
pub enum FragmentDofKind<'a> {
    RotateFrag,
    FlipFrag,
    ScaleAtom(AtomRef<'a>),
    ScaleFrag,
    ChangeParentBond,
    InvertBond {
        pivot: AtomRef<'a>,
        bound: AtomRef<'a>,
    },
    FlipRing {
        pivot1: AtomRef<'a>,
        pivot2: AtomRef<'a>,
    },
}

#[derive(Debug)]
pub struct FragmentDof<'a> {
    pub kind: FragmentDofKind<'a>,
    pub atoms: Vec<AtomRef<'a>>,
    pub frag: FragmentRef<'a>,
}

pub type FragmentDofRef<'a> = &'a RefCell<FragmentDof<'a>>;

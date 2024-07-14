use super::*;

#[derive(Debug)]
pub enum InteractionKind<'a> {
    Clash {
        k2: f32,
        atom3: AtomRef<'a>,
        sq_dist: f32,
    },
    Stretch {
        tolerance: f32,
    },
    Bend {
        k2: f32,
        atom3: AtomRef<'a>,
        is_ring: bool,
    },
    Inversion {
        atom3: AtomRef<'a>,
        atom4: AtomRef<'a>,
        is_z: bool,
        force: bool,
    },
}

#[derive(Debug)]
pub struct Interaction<'a> {
    pub kind: InteractionKind<'a>,
    pub atom1: AtomRef<'a>,
    pub atom2: AtomRef<'a>,
    pub rest_v: f32,
    pub k: f32,
}
impl Interaction<'_> {
    pub fn score(&self, score: &mut f32) {
        todo!()
    }
}

pub type InteractionRef<'a> = &'a RefCell<Interaction<'a>>;

pub mod coordgen;
pub mod sketcher;

const BOND_LENGTH: f32 = 50.0;
const EPSILON: f32 = 0.001;

pub mod prelude {
    pub use super::sketcher::{AtomRef, BondRef, Builder, MoleculeRef, Sketcher, Unsync};
}

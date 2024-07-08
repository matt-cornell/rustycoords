use super::*;

#[derive(Debug)]
pub struct FragmentBuilder {
    pub precision: f32,
}
impl FragmentBuilder {
    pub const fn new() -> Self {
        Self { precision: 0.01 }
    }
    pub fn initialize_coords(&self, f: sketcher::FragmentRef<'_>) {
        todo!()
    }
}

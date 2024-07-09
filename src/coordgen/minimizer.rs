use super::*;

#[derive(Debug, Default)]
pub struct Minimizer<'a> {
    pub atoms: Vec<AtomRef<'a>>,
}
impl<'a> Minimizer<'a> {
    pub const fn new() -> Self {
        Self { atoms: Vec::new() }
    }
    pub fn build_from_fragments(&mut self, b: bool) {
        todo!()
    }
    pub fn avoid_clashes(&mut self) -> bool {
        todo!()
    }
}

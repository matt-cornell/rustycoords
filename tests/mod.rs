use coordgen_rs::prelude::*;
use coordgen_rs::sketcher::intern::Leak;

#[test]
fn test_simple() {
    let intern = Leak;
    let mut sketch = Sketcher::new(&intern);
    let mut mol = Molecule::default();
    let center_c = mol.add_new_atom(6, &intern);
    let chain_c = mol.add_new_atom(6, &intern);
    let keto_o = mol.add_new_atom(8, &intern);
    let hydroxy_o = mol.add_new_atom(8, &intern);
    mol.add_new_bond(center_c, chain_c, 1, &intern);
    mol.add_new_bond(center_c, keto_o, 2, &intern);
    mol.add_new_bond(center_c, hydroxy_o, 1, &intern);
    sketch.initialize(mol);
    sketch.run_generate_coordinates();
}
use coordgen_rs::prelude::*;
use coordgen_rs::sketcher::intern::Leak;

#[test]
fn test_simple() {
    let intern = Leak;
    let mut mol = Builder::new(&intern);
    let center_c = mol.add_atom(6);
    let chain_c = mol.add_atom(6);
    let keto_o = mol.add_atom(8);
    let hydroxy_o = mol.add_atom(8);
    mol.add_bond(center_c, chain_c, 1);
    mol.add_bond(center_c, keto_o, 2);
    mol.add_bond(center_c, hydroxy_o, 1);
    Sketcher::new(&intern).generate(mol.finish());
}

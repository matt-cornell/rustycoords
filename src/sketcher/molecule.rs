use super::*;
use std::collections::VecDeque;
use std::ptr::eq;

#[derive(Debug, Default)]
pub struct Molecule<'a> {
    pub(crate) atoms: Vec<AtomRef<'a>>,
    pub(crate) bonds: Vec<BondRef<'a>>,
    pub(crate) main_fragment: Option<FragmentRef<'a>>,
    pub(crate) rings: Vec<RingRef<'a>>,
    pub(crate) fragments: Vec<FragmentRef<'a>>,
    pub(crate) proximity_relations: Vec<BondRef<'a>>,
    pub(crate) has_fixed_frags: bool,
    pub(crate) has_constrained_frags: bool,
    pub(crate) requires_minimization: bool,
}
impl<'a> Molecule<'a> {
    pub fn add_atom(&mut self, atom: Atom<'a>, intern: &'a dyn Interner) -> AtomRef<'a> {
        let a = intern.intern_atom(atom);
        self.atoms.push(a);
        a
    }
    pub fn add_new_atom(&mut self, atom_number: u8, intern: &'a dyn Interner) -> AtomRef<'a> {
        let a = intern.intern_atom(Atom::new(atom_number));
        self.atoms.push(a);
        a
    }
    pub fn add_bond(&mut self, bond: Bond<'a>, intern: &'a dyn Interner) -> BondRef<'a> {
        let b = intern.intern_bond(bond);
        self.bonds.push(b);
        b
    }
    pub fn add_new_bond(
        &mut self,
        start: AtomRef<'a>,
        end: AtomRef<'a>,
        order: u8,
        intern: &'a dyn Interner,
    ) -> BondRef<'a> {
        let mut bond = Bond::new(start, end);
        bond.bond_order = order;
        let b = intern.intern_bond(bond);
        self.bonds.push(b);
        b
    }

    pub fn force_update_struct(&mut self, intern: &'a dyn Interner) {
        // assign bonds and neighbors
        for a in &self.atoms {
            let mut a = a.borrow_mut();
            a.bonds.clear();
            a.neighbors.clear();
            a.rings.clear();
        }
        for bond in &self.bonds {
            let b = bond.borrow();
            {
                let mut start = b.start.borrow_mut();
                start.bonds.push(bond);
                start.neighbors.push(b.end);
            }
            {
                let mut end = b.end.borrow_mut();
                end.bonds.push(bond);
                end.neighbors.push(b.start);
            }
        }
        for a in &self.atoms {
            let mut a = a.borrow_mut();
            if a.implicit_h == u8::MAX {
                a.implicit_h = a.find_h_number();
            }
        }
        // find rings
        self.rings.clear();
        let mut queue = VecDeque::new();
        for &b in &self.bonds {
            for b in &self.bonds {
                let mut b = b.borrow_mut();
                b.sssr_visited = false;
                b.sssr_parent = None;
                b.sssr_parent_at_start = true;
            }
            b.borrow_mut().sssr_visited = true;
            queue.clear();
            queue.push_back(b);
            'bfs: while let Some(l) = queue.pop_front() {
                let last = l.borrow();
                let pivot = if !last.sssr_parent_at_start {
                    last.start
                } else {
                    last.end
                };
                for &n in &pivot.borrow().bonds {
                    if eq(l, n) {
                        continue;
                    }
                    let mut next = n.borrow_mut();
                    if next.sssr_visited {
                        if eq(n, b) {
                            drop(next);
                            let mut ring = Ring::default();
                            let mut last = Some(n);
                            while let Some(bond) = last {
                                ring.bonds.push(bond);
                                last = bond.borrow().sssr_parent;
                            }
                            let mut found = false;
                            'rings: for r in &self.rings {
                                let r = r.borrow();
                                if r.bonds.len() != ring.bonds.len() {
                                    continue;
                                }
                                for &b in &ring.bonds {
                                    if !r.bonds.iter().any(|&b2| eq(b, b2)) {
                                        continue 'rings;
                                    }
                                }
                                found = true;
                            }
                            if !found {
                                let ring = intern.intern_ring(ring);
                                self.rings.push(ring);
                            }
                            break 'bfs;
                        }
                    } else {
                        next.sssr_parent_at_start = !eq(next.end, pivot);
                        next.sssr_parent = Some(b);
                        queue.push_back(n);
                    }
                }
            }
        }
        for ring in &self.rings {
            for b in &ring.borrow().bonds {
                b.borrow_mut().rings.push(ring);
            }
        }
        for b in &self.bonds {
            let bond = b.borrow();
            if eq(bond.start, bond.end) {
                continue;
            }
            let mut start = bond.start.borrow_mut();
            let mut end = bond.end.borrow_mut();
            for &r in &bond.rings {
                if !start.rings.iter().any(|&r2| eq(r, r2)) {
                    start.rings.push(r);
                }
                if !end.rings.iter().any(|&r2| eq(r, r2)) {
                    end.rings.push(r);
                }
            }
        }
        for a in &self.atoms {
            for r in &a.borrow().rings {
                r.borrow_mut().atoms.push(a);
            }
        }
    }
}
pub type MoleculeRef<'a> = &'a RefCell<Molecule<'a>>;

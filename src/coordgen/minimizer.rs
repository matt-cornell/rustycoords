use super::*;
use std::ptr::eq;

struct Solutions<'a, 'b> {
    min: &'b mut Minimizer<'a>,
    sketcher: &'b Sketcher<'a>,
    mol: &'b Molecule<'a>,
    dofs: &'b [FragmentDofRef<'a>],
}
impl<'a> Solutions<'a, '_> {
    pub fn run_search(&mut self) -> bool {
        todo!()
    }
}

struct BendStub<'a> {
    atom1: AtomRef<'a>,
    atom2: AtomRef<'a>,
    atom3: AtomRef<'a>,
    rest_v: f32,
    k: f32,
    is_ring: bool,
}

pub struct Minimizer<'a> {
    interactions: Vec<InteractionRef<'a>>,
    intra_inters: Vec<InteractionRef<'a>>,
}
impl<'a> Minimizer<'a> {
    pub fn new() -> Self {
        Self {
            interactions: Vec::new(),
            intra_inters: Vec::new(),
        }
    }
    pub fn build_from_fragments(&mut self, first: bool) {
        todo!()
    }
    fn build_mol_from_fragments(&mut self, mol: &Molecule<'a>, first: bool) {
        todo!()
    }
    pub fn avoid_clashes(&mut self, sketcher: &Sketcher<'a>) -> bool {
        let mut all_clean = true;
        let mut dofs = Vec::new();
        let mut ring_ints = Vec::new();
        let mut non_ring_ints = Vec::new();
        let mut atom_buf = Vec::new();
        let mut lowest_energy_coords = Vec::new();
        let mut local_energy_list = Vec::new();
        let mut old_pos = Vec::new();
        for m in &sketcher.molecules {
            dofs.clear();
            // avoid clashes of molecule
            let mol = &mut m.borrow_mut();
            self.add_mol_clashes(mol, false, sketcher);
            let clash_e = self.score_clashes(mol);
            if clash_e < 10.0 {
                continue;
            }
            // flip fragments
            let mut it = mol.fragments.iter().rev().peekable();
            while let Some(frag) = it.next() {
                let f = frag.borrow();
                if f.fixed {
                    continue;
                }
                for &dof in &f.dofs {
                    if dof.borrow().num_states() > 1 {
                        dofs.push(dof);
                    }
                }
            }
            let clean_pose = Solutions {
                min: self,
                sketcher,
                mol: &mol,
                dofs: &dofs,
            }
            .run_search();
            self.build_from_fragments(true);
            if !clean_pose {
                all_clean = false;
                if clash_e >= 0.1 {
                    for &bond in &mol.bonds {
                        let b = bond.borrow();
                        if !b.is_terminal() {
                            continue;
                        }
                        let (mut term, coords) = if eq(b.start, b.end) {
                            let t = b.start.borrow_mut();
                            let c = t.coordinates;
                            (t, c)
                        } else {
                            let mut root = b.start.borrow_mut();
                            let mut term = b.end.borrow_mut();
                            if term.neighbors.len() != 1 {
                                std::mem::swap(&mut root, &mut term);
                            }
                            (term, root.coordinates)
                        };
                        if term.fixed {
                            continue;
                        }
                        for &b2 in &mol.bonds {
                            if !eq(bond, b2) && Self::bonds_clash(&bond.borrow(), &b2.borrow()) {
                                let c = term.coordinates;
                                term.set_coords(coords + (c - coords) * 0.1);
                            }
                        }
                    }
                    self.score_clashes(&mol);
                    mol.requires_minimization = true;
                }
            }
            if mol.requires_minimization {
                // minimize molecule
                old_pos.clear();
                old_pos.extend(mol.atoms.iter().map(|a| a.borrow().coordinates));
                self.interactions.clear();
                self.intra_inters.clear();
                self.add_mol_clashes(mol, false, sketcher);

                // add stretch interactions
                for b in &mol.bonds {
                    let b = b.borrow();
                    let at1 = b.start.borrow();
                    let at2 = b.end.borrow();
                    let rest_v = if at1.rigid && at2.rigid {
                        (at2.coordinates - at1.coordinates).length()
                    } else {
                        BOND_LENGTH
                    };
                    let k = at1.shares_a_ring_with(&at2).map_or(0.1, |r| {
                        if r.borrow().is_macrocycle() {
                            5.0
                        } else {
                            0.1
                        }
                    });
                    let int = Interaction {
                        kind: InteractionKind::Stretch { tolerance: 0.0 },
                        atom1: b.start,
                        atom2: b.end,
                        rest_v,
                        k,
                    };
                    let int = sketcher.interner.intern_interaction(int);
                    self.interactions.push(int);
                }
                // add bend interactions
                for a in &mol.atoms {
                    let at = a.borrow();
                    let nbonds = at.neighbors.len();
                    let mut inverted = false;
                    if nbonds > 1 {
                        at.clockwise_ordered_neighbors(&mut atom_buf);
                        let rest_v = if let [b0, b1] = at.bonds.as_slice() {
                            if b0.borrow().bond_order + b1.borrow().bond_order > 3 {
                                180.0
                            } else {
                                120.0
                            }
                        } else {
                            360.0 / nbonds as f32
                        };
                        for i in 0..nbonds {
                            let atom3 = atom_buf[i];
                            let atom1 = atom_buf[(i + 1) % nbonds];
                            let r = at.shares_a_ring_with_both(&atom1.borrow(), &atom3.borrow());
                            if let Some(r) = r {
                                let r = r.borrow();
                                if !r.is_macrocycle() {
                                    let mut extras = 0;
                                    for (fr, fa) in &r.fused_with {
                                        let fr = fr.borrow();
                                        if !fr.is_macrocycle() {
                                            extras += fr.atoms.len() - fa.len();
                                        }
                                    }
                                    let rest_v = 180.0 - (360.0 / (r.atoms.len() + extras) as f32);
                                    ring_ints.push(BendStub {
                                        atom1,
                                        atom3,
                                        rest_v,
                                        k: 100.0,
                                        atom2: a,
                                        is_ring: true,
                                    });
                                } else {
                                    if nbonds == 3 && !inverted {
                                        let at1 = atom1.borrow();
                                        let at3 = atom3.borrow();
                                        let other = atom_buf
                                            .iter()
                                            .find(|&&i| !(eq(atom1, i) || eq(atom3, i)));
                                        if let Some(other) = other {
                                            if math::same_side(
                                                at3.coordinates,
                                                other.borrow().coordinates,
                                                at1.coordinates,
                                                at.coordinates,
                                            ) {
                                                inverted = true;
                                            }
                                        }
                                    }
                                    let fused = inverted
                                        || if atom_buf.len() > 2 {
                                            true
                                        } else {
                                            atom_buf.iter().all(|n| {
                                                at.shares_a_ring_with(&n.borrow()).is_some()
                                            })
                                        };
                                    let stub = BendStub {
                                        atom1,
                                        atom3,
                                        rest_v,
                                        k: 1.0,
                                        atom2: a,
                                        is_ring: false,
                                    };
                                    if fused {
                                        ring_ints.push(stub);
                                    } else {
                                        non_ring_ints.push(stub);
                                    }
                                }
                            } else {
                                non_ring_ints.push(BendStub {
                                    atom1,
                                    atom3,
                                    rest_v,
                                    k: 1.0,
                                    atom2: a,
                                    is_ring: false,
                                });
                            }
                        }
                    }
                    if ring_ints.len() != 1 || non_ring_ints.len() != 2 {
                        inverted = false;
                    }
                    if !ring_ints.is_empty() {
                        let mut total: f32 = ring_ints.iter().map(|i| i.rest_v).sum();
                        if inverted {
                            total = 360.0 - total;
                        }
                        let rest_v = (360.0 - total) / non_ring_ints.len() as f32;
                        for i in &mut non_ring_ints {
                            i.rest_v = rest_v;
                        }
                    }
                    if let [nr1, nr2, nr3, nr4] = non_ring_ints.as_mut_slice() {
                        if at.cross_layout || sketcher.even_angles {
                            nr1.rest_v = 90.0;
                            nr2.rest_v = 90.0;
                            nr3.rest_v = 90.0;
                            nr4.rest_v = 90.0;
                        } else {
                            let idx = non_ring_ints
                                .iter()
                                .enumerate()
                                .max_by_key(|(_, i)| {
                                    (math::unsigned_angle(
                                        i.atom1.borrow().coordinates,
                                        i.atom2.borrow().coordinates,
                                        i.atom3.borrow().coordinates,
                                    ) * 10000000.0) as u32
                                })
                                .unwrap()
                                .0;
                            non_ring_ints[idx].rest_v = 120.0;
                            non_ring_ints[(idx + 1) % 4].rest_v = 90.0;
                            non_ring_ints[(idx + 2) % 4].rest_v = 60.0;
                            non_ring_ints[(idx + 3) % 4].rest_v = 90.0;
                        }
                    } else if non_ring_ints.len() > 4 {
                        let rest_v = 360.0 / non_ring_ints.len() as f32;
                        for i in &mut non_ring_ints {
                            i.rest_v = rest_v;
                        }
                    }
                    for &BendStub {
                        atom1,
                        atom2,
                        atom3,
                        rest_v,
                        k,
                        is_ring,
                    } in ring_ints.iter().chain(&non_ring_ints)
                    {
                        if at.fixed && atom1.borrow().fixed && atom3.borrow().fixed {
                            continue;
                        }
                        let int = Interaction {
                            kind: InteractionKind::Bend {
                                k2: 0.05,
                                atom3,
                                is_ring,
                            },
                            atom1,
                            atom2,
                            rest_v,
                            k,
                        };
                        let int = sketcher.interner.intern_interaction(int);
                        self.interactions.push(int);
                    }
                }
                // add chiral inversion constraints
                for ring in &mol.rings {
                    let r = ring.borrow();
                    if !r.is_macrocycle() {
                        continue;
                    }
                    atom_buf.clone_from(&r.atoms);
                    let start = atom_buf[0];
                    FragmentBuilder::order_chain_of_atoms(&mut atom_buf, start);
                    let len = atom_buf.len();
                    for i in 0..len {
                        let a1 = atom_buf[(i + len - 1) % len];
                        let a11 = atom_buf[(i + len - 2) % len];
                        let a2 = atom_buf[(i + 1) % len];
                        let ai = atom_buf[i];
                        let at1 = a1.borrow();
                        let bond = at1
                            .neighbors
                            .iter()
                            .zip(&at1.bonds)
                            .find(|(n, _)| eq(**n, ai))
                            .unwrap()
                            .1
                            .borrow();
                        if bond.is_stereo() {
                            let is_z = bond.marked_as_cis(&a11.borrow(), &a2.borrow());
                            let int = Interaction {
                                kind: InteractionKind::Inversion {
                                    atom3: a11,
                                    atom4: a1,
                                    is_z,
                                    force: false,
                                },
                                atom1: ai,
                                atom2: a2,
                                k: 1.0,
                                rest_v: 50.0,
                            };
                            let int = sketcher.interner.intern_interaction(int);
                            self.interactions.push(int);
                        }
                    }
                }
                // run minimization
                lowest_energy_coords.resize(mol.atoms.len(), PointF::default());
                local_energy_list.clear();
                let mut min = f32::INFINITY;
                for iter in 0..sketcher.max_iterations {
                    let mut total = 0.0;
                    for i in &self.interactions {
                        i.borrow().score(&mut total);
                    }
                    local_energy_list.push(total);
                    if total < min {
                        for (coord, atom) in lowest_energy_coords.iter_mut().zip(&mol.atoms) {
                            *coord = atom.borrow().coordinates;
                        }
                        min = total;
                    }
                    // apply forces
                    let mut distance = 0.0;
                    for atom in &mol.atoms {
                        let mut atom = atom.borrow_mut();
                        if atom.fixed {
                            continue;
                        }
                        let mut disp = atom.force * 0.3;
                        if disp.0.is_nan() || disp.1.is_nan() {
                            disp = PointF::default();
                        }
                        let dist = disp.sq_length();
                        let dsq = dist.max(EPSILON);
                        if dsq > 0.01 {
                            disp *= 0.1 * dsq.sqrt();
                        }
                        atom.coordinates += disp;
                        atom.force = PointF::default();
                        distance += dist;
                    }
                    if distance < 0.001 {
                        break;
                    }
                    if iter < 200 {
                        continue;
                    }
                    if local_energy_list[iter - 100] - total < 20.0 {
                        break;
                    }
                }
                if min < f32::INFINITY {
                    for (atom, coord) in mol.atoms.iter().zip(&lowest_energy_coords) {
                        atom.borrow_mut().coordinates = *coord;
                    }
                }
                // final stereo check
                if !mol.bonds.iter().all(|b| b.borrow().check_stereo()) {
                    for (atom, coords) in mol.atoms.iter().zip(&old_pos) {
                        atom.borrow_mut().coordinates = *coords;
                    }
                }
            }
        }
        all_clean
    }
    fn add_mol_clashes(&mut self, mol: &Molecule<'a>, intra: bool, sketcher: &Sketcher<'a>) {
        if mol.atoms.len() <= 1 {
            return;
        }
        for &atom in &mol.atoms {
            let at2 = atom.borrow();
            if at2.fixed || at2.rigid {
                continue;
            }
            let a2f = at2.fragment.unwrap();
            let dofs2 = !intra && a2f.borrow().atom_dofs.contains_key(&(atom as *const _));
            for bond in &mol.bonds {
                let b = bond.borrow();
                if eq(atom, b.start) || eq(atom, b.end) || eq(b.start, b.end) {
                    continue;
                }
                let at1 = b.start.borrow();
                let at3 = b.end.borrow();
                if at1.fixed || at3.fixed || at1.rigid || at3.rigid {
                    continue;
                }
                let a1f = at1.fragment.unwrap();
                let a3f = at3.fragment.unwrap();
                if dofs2
                    && a1f.borrow().atom_dofs.contains_key(&(b.start as *const _))
                    && a3f.borrow().atom_dofs.contains_key(&(b.end as *const _))
                {
                    if eq(a1f, a2f) || eq(a3f, a2f) {
                        continue;
                    }
                }
                if at2.neighbors.iter().any(|&n| {
                    eq(n, b.start)
                        || eq(n, b.end)
                        || at1.neighbors.iter().any(|&n2| eq(n, n2))
                        || at3.neighbors.iter().any(|&n2| eq(n, n2))
                }) {
                    continue;
                }
                let mut rest_vk = 0.8;
                if at2.atom_number == 6 && at2.charge == 0 {
                    rest_vk -= 0.1;
                }
                if at1.atom_number == 6
                    && at1.charge == 0
                    && at2.atom_number == 6
                    && at2.charge == 0
                {
                    rest_vk -= 0.1;
                }
                let rest_v = rest_vk * rest_vk * BOND_LENGTH * BOND_LENGTH;
                let interaction = Interaction {
                    kind: InteractionKind::Clash {
                        k2: 0.1,
                        atom3: b.end,
                        sq_dist: 0.0,
                    },
                    atom1: b.start,
                    atom2: atom,
                    rest_v,
                    k: 900.0,
                };
                let int = sketcher.interner.intern_interaction(interaction);
                self.intra_inters.push(int);
                self.interactions.push(int);
            }
        }
    }
    fn score_clashes(&mut self, mol: &Molecule<'a>) -> f32 {
        let mut e = 0.0;
        for i in &self.intra_inters {
            i.borrow().score(&mut e);
        }
        for f in &mol.fragments {
            for dof in &f.borrow().dofs {
                e += dof.borrow().current_penalty();
            }
        }
        if mol.bonds.len() > 2 {
            for (n, b1) in mol.bonds[..(mol.bonds.len() - 1)].iter().enumerate() {
                let b1 = b1.borrow();
                for b2 in &mol.bonds[(n + 1)..] {
                    let b2 = b2.borrow();
                    if !Self::bonds_clash(&b1, &b2) {
                        continue;
                    }
                    let mut penalty = 2500.0;
                    if b1.is_terminal() || b2.is_terminal() {
                        penalty *= 0.5;
                    }
                    if b1.is_in_macrocycle() || b2.is_in_macrocycle() {
                        penalty *= 8.0;
                    }
                    if b1.is_in_small_ring() || b2.is_in_small_ring() {
                        penalty *= 2.0;
                    }
                    e += penalty;
                }
            }
        }
        for r in &mol.rings {
            let r = r.borrow();
            if r.atoms.len() > 9 {
                continue;
            }
            if r.atoms.len() < 3 {
                continue;
            }
            let c = r.find_center();
            let frag = r.atoms[0].borrow().fragment.unwrap();
            for a in &mol.atoms {
                let a = a.borrow();
                if eq(a.fragment.unwrap(), frag) {
                    continue;
                }
                let d = c - a.coordinates;
                if d.0.abs() > BOND_LENGTH || d.1.abs() > BOND_LENGTH {
                    continue;
                }
                let sq = d.sq_length();
                if sq > BOND_LENGTH * BOND_LENGTH {
                    continue;
                }
                let dist = sq.sqrt();
                if dist < BOND_LENGTH {
                    e += 50.0 + 100.0 * (1.0 - dist - BOND_LENGTH);
                }
            }
        }
        e
    }
    fn bonds_clash(b1: &Bond<'a>, b2: &Bond<'a>) -> bool {
        if eq(b1.start, b2.start)
            || eq(b1.start, b2.end)
            || eq(b1.end, b2.start)
            || eq(b1.end, b2.end)
        {
            return false;
        }
        let start1 = b1.start.borrow().coordinates;
        let start2 = b2.start.borrow().coordinates;
        let end1 = b1.end.borrow().coordinates;
        let end2 = b2.end.borrow().coordinates;
        if start1.0.max(end1.0) < start2.0.min(end1.0)
            || start1.1.max(end1.1) < start2.1.min(end1.1)
            || start1.0.min(end1.0) < start2.0.max(end1.0)
            || start1.1.min(end1.1) < start2.1.max(end1.1)
        {
            return false;
        }
        if math::points_coincide(start1, start2)
            || math::points_coincide(start1, end2)
            || math::points_coincide(end1, start2)
            || math::points_coincide(end1, end2)
        {
            return true;
        }
        math::intersection_of_segments(start1, end1, start2, end2).is_some()
    }
}

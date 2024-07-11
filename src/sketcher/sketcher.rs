use super::*;
use ahash::AHashMap;
use std::collections::VecDeque;
use std::f32::consts::{FRAC_PI_2, FRAC_PI_6, PI};
use std::ptr::eq;

pub struct Sketcher<'a> {
    pub interner: &'a dyn Interner,
    fragments: Vec<FragmentRef<'a>>,
    minimizer: coordgen::Minimizer<'a>,
    frag_builder: coordgen::FragmentBuilder,
    atoms: Vec<AtomRef<'a>>,
    bonds: Vec<BondRef<'a>>,
    ref_atoms: Vec<AtomRef<'a>>,
    ref_bonds: Vec<BondRef<'a>>,
    proximity_relations: Vec<BondRef<'a>>,
    molecules: Vec<MoleculeRef<'a>>,
    sin_flip: f32,
    cos_flip: f32,
    center: PointF,
}
impl<'a> Sketcher<'a> {
    pub fn new(interner: &'a dyn Interner) -> Self {
        Self {
            interner,
            fragments: Vec::new(),
            minimizer: coordgen::Minimizer::new(),
            frag_builder: coordgen::FragmentBuilder::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            ref_atoms: Vec::new(),
            ref_bonds: Vec::new(),
            proximity_relations: Vec::new(),
            molecules: vec![interner.intern_molecule(Molecule::default())],
            sin_flip: 0.0,
            cos_flip: 0.0,
            center: PointF(0.0, 0.0),
        }
    }
    pub fn atoms(&self) -> &[AtomRef<'a>] {
        &self.atoms
    }
    pub fn bonds(&self) -> &[BondRef<'a>] {
        &self.bonds
    }
    pub fn generate(&mut self, mol: MoleculeRef<'a>) -> bool {
        self.initialize(mol);
        self.run_generate_coordinates()
    }
    fn initialize(&mut self, mol: MoleculeRef<'a>) {
        let mut mol = mol.borrow_mut();
        self.ref_atoms.clone_from(&mol.atoms);
        self.ref_bonds.clone_from(&mol.bonds);
        let mut bonds_to_atom = AHashMap::with_capacity(self.ref_atoms.len());
        for &b in &self.ref_bonds {
            let b = b.borrow();
            bonds_to_atom
                .entry(b.start as *const _)
                .and_modify(|e| *e += 1)
                .or_insert(1usize);
            bonds_to_atom
                .entry(b.end as *const _)
                .and_modify(|e| *e += 1)
                .or_insert(1usize);
        }
        for &b in &self.ref_bonds {
            let mut b = b.borrow_mut();
            if b.skip {
                continue;
            }
            if b.bond_order == 1 || b.bond_order == 2 {
                let term = bonds_to_atom[&(b.start as *const _)] == 1
                    || bonds_to_atom[&(b.end as *const _)] == 1;
                if !term
                    && (Atom::is_metal(b.start.borrow().atom_number)
                        || Atom::is_metal(b.end.borrow().atom_number))
                {
                    b.bond_order = 0;
                }
            }
        }
        mol.bonds.retain(|i| {
            let b = i.borrow();
            !b.skip
                && (b.bond_order > 0 || {
                    self.proximity_relations.push(*i);
                    false
                })
        });
        self.canonical_orderng(&mut mol);
        self.atoms.clone_from(&mol.atoms);
        self.bonds.clone_from(&mol.bonds);
        mol.force_update_struct();
        self.split_molecules(&mut mol);
        for &bref in &self.proximity_relations {
            let b = bref.borrow();
            let start = b.start.borrow();
            if let Some(mol) = start.molecule {
                mol.borrow_mut().proximity_relations.push(bref);
            }
            if !eq(b.start, b.end) {
                let end = b.end.borrow();
                if let Some(mol) = end.molecule {
                    mol.borrow_mut().proximity_relations.push(bref);
                }
            }
        }
        self.flag_cross_atoms();
        // TODO: assign coordgen stuff
    }
    fn run_generate_coordinates(&mut self) -> bool {
        self.find_fragments();
        // self.minimizer.build_from_fragments(true);
        // let clean_pose = self.minimizer.avoid_clashes();
        let clean_pose = true;
        self.best_rotation();
        // self.maybe_flip();
        // self.arrange_multiple_molecules();
        clean_pose
    }
    /// take self so we can use ref buffers
    fn canonical_orderng(&mut self, mol: &mut Molecule<'a>) {
        self.ref_atoms.clear();
        self.ref_bonds.clear();
        let mut scores = vec![0; mol.atoms.len()];
        for (n, a) in mol.atoms.iter().enumerate() {
            let mut a = a.borrow_mut();
            a.general_use_n = n;
            a.general_use_visited = false;
        }
        Self::morgan_scores(&mol.atoms, &mol.bonds, &mut scores);
        for (s, a) in scores.iter_mut().zip(&mol.atoms) {
            *s *= 100;
            *s += a.borrow().atom_number as usize;
        }
        for a in &mol.atoms {
            let mut a = a.borrow_mut();
            a.neighbors.clear();
            a.bonds.clear();
        }
        for bref in &mol.bonds {
            let mut b = bref.borrow_mut();
            {
                let mut start = b.start.borrow_mut();
                start.bonds.push(bref);
                start.neighbors.push(b.end);
            }
            {
                let mut end = b.end.borrow_mut();
                end.neighbors.push(b.start);
                end.bonds.push(bref);
            }
            b.sssr_visited = false;
        }
        let mut queue = VecDeque::new();
        self.ref_atoms.clear();
        self.ref_bonds.clear();
        // skip assignment here, nothing should be able to modify it since the last time?
        loop {
            let mut score_max_i = usize::MAX;
            for i in 0..scores.len() {
                if mol.atoms[i].borrow().general_use_visited {
                    continue;
                }
                if score_max_i == usize::MAX || scores[i] > scores[score_max_i] {
                    score_max_i = i;
                }
            }
            if score_max_i == usize::MAX {
                break;
            }
            queue.push_back(mol.atoms[score_max_i]);
            mol.atoms[score_max_i].borrow_mut().general_use_visited = true;
            while let Some(at) = queue.pop_front() {
                self.ref_atoms.push(at);
                loop {
                    let a = at.borrow();
                    let res = a.bonds.iter().zip(&a.neighbors).enumerate().fold(
                        None,
                        |r, (i, (b, n))| {
                            if b.borrow().sssr_visited {
                                r
                            } else if let Some((last, _)) = r {
                                let gn = n.borrow().general_use_n;
                                if last < gn {
                                    Some((gn, i))
                                } else {
                                    r
                                }
                            } else {
                                Some((n.borrow().general_use_n, i))
                            }
                        },
                    );
                    let Some((_, neighbor)) = res else { break };
                    a.bonds[neighbor].borrow_mut().sssr_visited = true;
                    let n = a.neighbors[neighbor];
                    if eq(n, at) {
                        continue;
                    }
                    let mut n = n.borrow_mut();
                    if !n.general_use_visited {
                        n.general_use_visited = true;
                        queue.push_back(a.neighbors[neighbor]);
                    }
                    let b = a.bonds[neighbor];
                    b.borrow_mut().sssr_visited = true;
                    self.ref_bonds.push(b);
                }
            }
        }
        mol.atoms.clone_from(&self.ref_atoms);
        mol.bonds.clone_from(&self.ref_bonds);
    }
    fn split_molecules(&mut self, mol: &mut Molecule<'a>) {
        if mol.atoms.is_empty() {
            return;
        }
        for a in &mol.atoms {
            a.borrow_mut().general_use_visited = false;
        }
        let mut q = VecDeque::new();
        while let Some(&first) = mol.atoms.first() {
            q.push_back(first);
            while let Some(a) = q.pop_front() {
                let mut a = a.borrow_mut();
                a.general_use_visited = true;
                q.extend(
                    a.neighbors
                        .iter()
                        .filter_map(|n| (!n.borrow().general_use_visited).then_some(*n)),
                );
            }
            let mut out = Molecule::default();
            let mut i = 0;
            while i < mol.rings.len() {
                if mol.rings[i].borrow().atoms[0].borrow().general_use_visited {
                    out.rings.push(mol.rings.remove(i));
                } else {
                    i += 1;
                }
            }
            i = 0;
            while i < mol.bonds.len() {
                if mol.bonds[i].borrow().start.borrow().general_use_visited {
                    out.bonds.push(mol.bonds.remove(i));
                } else {
                    i += 1;
                }
            }
            i = 0;
            while i < mol.atoms.len() {
                if mol.atoms[i].borrow().general_use_visited {
                    out.atoms.push(mol.atoms.remove(i));
                } else {
                    i += 1;
                }
            }
            self.molecules.push(self.interner.intern_molecule(out));
        }
    }
    fn flag_cross_atoms(&mut self) {
        for a in &self.atoms {
            let mut a = a.borrow_mut();
            if a.atom_number == 15 || a.atom_number == 16 {
                a.cross_layout = true;
            }
        }
        for a in &self.atoms {
            let mut a = a.borrow_mut();
            if a.cross_layout {
                continue;
            }
            a.cross_layout = a
                .neighbors
                .iter()
                .filter(|n| n.borrow().neighbors.len() > 3)
                .skip(3)
                .next()
                .is_some();
        }
    }
    fn find_fragments(&mut self) {
        let mut independent = Vec::new();
        for &mol in &self.molecules {
            coordgen::fragmenter::split_into_fragments(mol, self.interner);
            let mol = mol.borrow();
            if mol.fragments.is_empty() {
                continue;
            }
            self.fragments.extend(mol.fragments.iter().copied());
            independent.extend(mol.main_fragment);
        }
        if self.fragments.is_empty() {
            return;
        }
        for f in &independent {
            self.assign_num_children(&mut f.borrow_mut());
        }
        for &f in &self.fragments {
            self.frag_builder.initialize_coords(f, self.interner);
        }
        for f in &independent {
            self.assign_longest_chain(&mut f.borrow_mut(), None);
        }
    }
    fn assign_num_children(&self, frag: &mut Fragment<'a>) {
        let mut cumatoms = 0;
        let mut cumranks = 0.0;
        let mut childatoms = 0;
        for &child in &frag.children {
            let mut child = child.borrow_mut();
            self.assign_num_children(&mut child);
            cumatoms += child.num_children;
            cumranks += child.child_rank;
            childatoms += child.atoms.len();
        }
        frag.num_children = cumatoms + childatoms;
        frag.child_rank = cumranks * 0.01 + childatoms as f32;
    }
    fn assign_longest_chain(
        &self,
        frag: &mut Fragment<'a>,
        parent: Option<&AHashMap<*const RefCell<Atom<'a>>, PointF>>,
    ) {
        frag.longest_chain = 0.0;
        for &child in &frag.children {
            let mut child = child.borrow_mut();
            self.assign_longest_chain(&mut child, Some(&frag.coords));
            if child.longest_chain > frag.longest_chain {
                frag.longest_chain = child.longest_chain;
            }
        }
        if let Some((f, bond)) = frag.parent {
            let key = &(bond.borrow().end as *const _);
            frag.longest_chain += parent
                .map_or_else(|| f.borrow().coords[key], |p| p[key])
                .length();
        }
    }
    fn add_to_vector(weight: f32, angle: f32, angles: &mut Vec<(f32, f32)>) {
        let angle = ((angle * 100.0).floor() * 0.01).rem_euclid(PI);
        let len = angles.len();
        for (n, a) in angles.iter_mut().enumerate() {
            if a.1 < angle - EPSILON {
                if n == len - 1 {
                    angles.push((weight, angle));
                    break;
                }
            } else if (a.1 - angle).abs() < EPSILON {
                a.0 += weight;
                break;
            } else {
                angles.insert(n, (weight, angle));
                break;
            }
        }
        if angles.is_empty() {
            angles.push((weight, angle));
        }
    }
    fn best_rotation(&mut self) {
        let mut angles = Vec::new();
        for &mol in &self.molecules {
            angles.clear();
            let mol = mol.borrow_mut();
            for a in &mol.atoms {
                let a = a.borrow();
                if a.rings.is_empty() {
                    continue;
                }
                if a.neighbors.len() > 1 {
                    for (n, i) in a.neighbors[..(a.neighbors.len() - 1)].iter().enumerate() {
                        let i = i.borrow();
                        let mut w_base = 6.0;
                        if i.neighbors.len() != 1 {
                            w_base += 2.0;
                        }
                        // this should be j according the source, but it looks like a typo
                        if i.atom_number == 6 {
                            w_base += 1.0;
                        }
                        if i.charge == 0 {
                            w_base += 1.0;
                        }
                        for j in &a.neighbors[(n + 1)..] {
                            let j = j.borrow();
                            let mut weight = w_base;
                            if j.neighbors.len() != 1 {
                                weight += 2.0;
                            }
                            if j.atom_number == 6 {
                                weight += 1.0;
                            }
                            if j.charge == 0 {
                                weight += 1.0;
                            }
                            let p = i.coordinates - j.coordinates;
                            let angle = (-p.1).atan2(p.1);
                            Self::add_to_vector(weight, angle, &mut angles);
                        }
                    }
                }
            }
            for b in &mol.bonds {
                let b = b.borrow();
                let start = b.start.borrow();
                let end = b.end.borrow();
                let p = end.coordinates - start.coordinates;
                let angle = (((-p.0).atan2(p.1) * 100.0).floor() * 0.01).rem_euclid(PI);
                let mut last_angle = angle;
                for i in 0..6 {
                    let mut weight = match i {
                        1 | 5 => 5.0,
                        0 | 3 => 1.5,
                        _ => 1.0,
                    };
                    if b.bond_order == 2
                        && i == 3
                        && (start.neighbors.len() == 1 || end.neighbors.len() == 1)
                    {
                        weight += 1.5;
                    }
                    if start.neighbors.len() == 1 && end.neighbors.len() == 1 && i == 0 {
                        weight += 10.0;
                    }
                    Self::add_to_vector(weight, angle, &mut angles);
                    last_angle += FRAC_PI_6;
                    if last_angle >= PI {
                        last_angle -= PI;
                    }
                }
            }
            for f in &mol.fragments {
                let f = f.borrow();
                match f.rings.as_slice() {
                    [r1, r2] => {
                        let p = r2.borrow().find_center() - r1.borrow().find_center();
                        let angle = (-p.1).atan2(p.0);
                        Self::add_to_vector(25.0, angle, &mut angles);
                    }
                    [r1, r2, r3] => {
                        let mut r1p = r1.borrow().find_center();
                        let mut r2p = r2.borrow().find_center();
                        for r in [r1, r2, r3] {
                            let r = r.borrow();
                            if let [(_, f1a), (_, f2a)] = r.fused_with.as_slice() {
                                if let ([f1a1, f1a2], [f2a1, f2a2]) =
                                    (f1a.as_slice(), f2a.as_slice())
                                {
                                    r1p = (f1a1.borrow().coordinates + f1a2.borrow().coordinates)
                                        * 0.5;
                                    r2p = (f2a1.borrow().coordinates + f2a2.borrow().coordinates)
                                        * 0.5;
                                    break;
                                }
                            }
                        }
                        let p = r2p - r1p;
                        let angle = (-p.1).atan2(p.0);
                        Self::add_to_vector(50.0, angle, &mut angles);
                    }
                    _ => {
                        for ring in &f.rings {
                            let ring = ring.borrow();
                            if ring.atoms.len() != 6 {
                                continue;
                            }
                            for (_, fas) in &ring.fused_with {
                                if let [fa1, fa2] = fas.as_slice() {
                                    let fa1 = fa1.borrow();
                                    let fa2 = fa2.borrow();
                                    let p = fa1.coordinates - fa2.coordinates;
                                    let angle = (-p.1).atan2(p.0) - FRAC_PI_2;
                                    Self::add_to_vector(25.0, angle, &mut angles)
                                }
                            }
                        }
                    }
                }
            }
            if angles.len() > 1 {
                let last = angles.last().unwrap();
                if last.1 - angles[0].1 >= FRAC_PI_2 - 2.0 * EPSILON {
                    angles[0].0 += last.0;
                    angles.pop();
                }
            }
            if !angles.is_empty() {
                let best = angles.iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();
                let (s, c) = best.1.sin_cos();
                let mut center: PointF = mol.atoms.iter().map(|a| a.borrow().coordinates).sum();
                if !mol.atoms.is_empty() {
                    center /= mol.atoms.len() as f32;
                }
                for a in &mol.atoms {
                    let mut a = a.borrow_mut();
                    let mut v = a.coordinates - center;
                    v.rotate(-s, c);
                    a.coordinates = v + center;
                }
                self.sin_flip = -s;
                self.cos_flip = c;
                self.center = center;
            }
        }
    }
    #[allow(dead_code)]
    fn maybe_flip(&mut self) {
        // this can't be that important, right?
        todo!()
    }
    #[allow(dead_code)]
    fn arrange_multiple_molecules(&mut self) {
        // also probably fine
        todo!()
    }
    fn morgan_scores(atoms: &[AtomRef<'a>], bonds: &[BondRef<'a>], scores: &mut [usize]) {
        if atoms.len() < 2 {
            return;
        }
        scores.fill(1);
        let mut new_scores = vec![0; scores.len()];
        let mut ordered = Vec::new();
        let mut idx1;
        let mut idx2;
        let mut old_ties = atoms.len();
        loop {
            for b in bonds {
                let b = b.borrow();
                idx1 = b.start.borrow().general_use_n;
                idx2 = b.end.borrow().general_use_n;
                new_scores[idx1] += scores[idx2];
                new_scores[idx2] += scores[idx1];
            }
            ordered.clone_from(&new_scores);
            ordered.sort();
            let new_ties = ordered.windows(2).filter(|i| i[0] == i[1]).count();
            if new_ties >= old_ties {
                break;
            }
            old_ties = new_ties;
            scores.copy_from_slice(&ordered);
        }
    }
    pub(crate) fn alignment_matrix(refs: &[PointF], points: &[PointF]) -> [f32; 4] {
        debug_assert_eq!(refs.len(), points.len());
        let mut a = [0.0; 4];
        for (r, p) in refs.iter().zip(points) {
            a[0] += r.0 + p.0;
            a[1] += r.1 + p.0;
            a[2] += r.0 + p.1;
            a[3] += r.1 + p.1;
        }
        let su = [
            a[0] * a[2] + a[1] * a[1],
            a[0] * a[2] + a[1] * a[3],
            a[2] * a[0] + a[3] * a[1],
            a[2] * a[2] + a[3] * a[3],
        ];
        let phi = (su[1] + su[2]).atan2(su[0] - su[3]) * 0.5;
        let (sphi, cphi) = phi.sin_cos();
        let u = [-cphi, -sphi, -sphi, cphi];
        let sw = [
            a[0] * a[0] + a[2] * a[2],
            a[0] * a[1] + a[2] * a[3],
            a[1] * a[0] + a[3] * a[2],
            a[1] * a[1] + a[3] * a[3],
        ];
        let theta = (sw[1] + sw[2]).atan2(sw[0] - sw[3]) * 0.5;
        let (stheta, ctheta) = theta.sin_cos();
        let w = [ctheta, -stheta, stheta, ctheta];
        // let susum = su[0] + su[3];
        // let sudif = ((su[0] - su[3]) * (su[0] - su[3]) + 4.0 * su[1] * su[2]).sqrt();
        // let sig = [
        //     ((susum + sudif) * 0.5).sqrt(),
        //     0.0,
        //     0.0,
        //     ((susum - sudif) * 0.5).sqrt(),
        // ];
        let u1a = [
            u[0] * a[0] + u[2] * a[2],
            u[0] * a[1] + u[2] * a[3],
            u[1] * a[0] + u[3] * a[2],
            u[1] * u[1] + u[3] * a[3],
        ];
        let s = [
            u1a[0] * w[0] + u1a[1] * w[2],
            u1a[0] * w[1] + u1a[1] * w[3],
            u1a[2] * w[0] + u1a[3] * w[2],
            u1a[2] * w[1] + u1a[3] * w[3],
        ];
        let c = [s[0].signum(), 0.0, 0.0, s[3].signum()];
        let v = [
            w[0] * c[0] + w[1] * c[2],
            w[0] * c[1] + w[1] * c[3],
            w[2] * c[0] + w[3] * c[2],
            w[2] * c[1] + w[3] * c[3],
        ];
        [
            v[0] * u[0] + v[1] * c[2],
            v[0] * c[1] + v[1] * c[3],
            v[2] * c[0] + v[3] * c[2],
            v[2] * c[1] + v[3] * c[3],
        ]
    }
}

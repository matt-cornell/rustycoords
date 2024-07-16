use super::*;
use ahash::AHashSet;
use std::collections::VecDeque;
use std::f32::consts::{FRAC_PI_2, FRAC_PI_3, TAU};
use std::ptr::eq;

pub(crate) fn initialize_coords<'a>(f: FragmentRef<'a>, sketcher: &Sketcher<'a>) {
    let mut frag = f.borrow_mut();
    // build rings
    if !frag.rings.is_empty() {
        // initialize fused ring information
        if frag.rings.len() > 1 {
            for atom in &frag.atoms {
                let a = atom.borrow();
                if a.rings.len() < 2 {
                    continue;
                }
                for (n, &iref) in a.rings[..(a.rings.len() - 1)].iter().enumerate() {
                    let mut i = iref.borrow_mut();
                    for &jref in &a.rings[(n + 1)..] {
                        if eq(iref, jref) {
                            continue;
                        }
                        let mut j = jref.borrow_mut();
                        let mut found = false;
                        for (r, ats) in &mut i.fused_with {
                            if !eq(*r, jref) {
                                continue;
                            }
                            ats.push(atom);
                            for (rr, rats) in &mut j.fused_with {
                                if !eq(*rr, iref) {
                                    continue;
                                }
                                rats.push(atom);
                                break;
                            }
                            found = true;
                            break;
                        }
                        if !found {
                            i.fused_with.push((jref, vec![atom]));
                            j.fused_with.push((iref, vec![atom]));
                        }
                    }
                }
            }
            for bond in &frag.bonds {
                let b = bond.borrow();
                if b.bond_order != 1 {
                    continue;
                }
                if eq(b.start, b.end) {
                    continue;
                }
                if b.start
                    .borrow()
                    .shares_a_ring_with(&b.end.borrow())
                    .is_some()
                {
                    continue;
                }
                let start = b.start.borrow();
                let end = b.end.borrow();
                if start.rings.is_empty() || end.rings.is_empty() {
                    continue;
                }
                for &ring1 in &start.rings {
                    let mut r = ring1.borrow_mut();
                    for &ring2 in &end.rings {
                        if eq(ring1, ring2) {
                            continue; // should be unreachable
                        }
                        let mut r2 = ring2.borrow_mut();
                        r.fused_with.push((ring2, Vec::new()));
                        r2.fused_with.push((ring1, Vec::new()));
                        r.fusion_bonds.push(bond);
                        r2.fusion_bonds.push(bond);
                    }
                }
            }
            for ring in &frag.rings {
                let mut r = ring.borrow_mut();
                for (_, fa) in &mut r.fused_with {
                    if fa.len() <= 2 {
                        continue;
                    }
                    let &start = fa
                        .iter()
                        .find(|a| {
                            let a = a.borrow();
                            let mut it = a
                                .neighbors
                                .iter()
                                .filter(|&&n| fa.iter().any(|&f| eq(f, n)));
                            it.next().is_some() && it.next().is_none()
                        })
                        .unwrap();
                    order_chain_of_atoms(fa, start);
                }
            }
            for &ring in &frag.rings {
                for (fr, fa) in &ring.borrow().fused_with {
                    if fa.len() <= 2 {
                        continue;
                    }
                    for &atom in fa {
                        let mut found = false;
                        let mut a = atom.borrow_mut();
                        for &n in &a.neighbors {
                            if eq(atom, n) {
                                continue;
                            }
                            let n = n.borrow();
                            let mut it = n.rings.iter().filter(|&&nr| eq(nr, ring) || eq(nr, *fr));
                            if it.next().is_some() && it.next().is_none() {
                                found = true;
                                break;
                            }
                        }
                        a.is_shared_and_inner = !found;
                    }
                }
            }
        }
        if let [ring] = frag.rings.as_slice() {
            build_ring(ring);
        } else {
            todo!() // CoordgenFragmentBuilder.cpp:1085
        }
    }
    // build non ring atoms
    {
        let mut visited = AHashSet::new();
        let mut queue = VecDeque::new();
        for &atom in &frag.atoms {
            if !atom.borrow().rings.is_empty() {
                visited.insert(atom as *const _);
                queue.push_back(atom);
            }
        }
        if queue.is_empty() {
            let start = if let Some((_, bond)) = frag.parent {
                let bond = bond.borrow();
                bond.start
                    .borrow_mut()
                    .set_coords(PointF(-BOND_LENGTH, 0.0));
                bond.end
            } else {
                frag.atoms[0]
            };
            start.borrow_mut().set_coords(PointF::default());
            queue.push_back(start);
            visited.insert(start as *const _);
        }
        let mut neighbors = Vec::new();
        let mut angles = Vec::new();
        while let Some(at) = queue.pop_front() {
            let atom = at.borrow();
            let mut init = PointF(-BOND_LENGTH, 0.0);
            neighbors.clear();
            angles.clear();
            // init vars for coords
            {
                if atom.rings.is_empty() {
                    if atom.neighbors.len() != 4 {
                        neighbors.clone_from(&atom.neighbors);
                    } else {
                        todo!() // CoordgenFragmentBuilder.cpp:640
                    }
                } else {
                    todo!() // CoordgenFragmentBuilder.cpp:667
                }
                let n = neighbors
                    .iter()
                    .position(|a| {
                        visited.contains(&(*a as *const _)) && {
                            init = a.borrow().coordinates - atom.coordinates;
                            true
                        }
                    })
                    .unwrap_or(0);
                neighbors.rotate_left(n);
                // neighbors angles at center
                let mut div = atom.neighbors.len();
                if !sketcher.even_angles {
                    if let [n1, n2] = atom.neighbors.as_slice() {
                        if atom.atom_number == 6
                            && !(n1.borrow().cross_layout || n2.borrow().cross_layout)
                        {
                            div = 3;
                        }
                        if atom.bonds[0].borrow().bond_order + atom.bonds[1].borrow().bond_order
                            >= 4
                        {
                            div = 2;
                        }
                    } else if atom.neighbors.len() == 4 && !atom.cross_layout {
                        angles.extend_from_slice(&[
                            FRAC_PI_3,
                            FRAC_PI_2,
                            FRAC_PI_3 * 2.0,
                            FRAC_PI_2,
                        ]);
                    }
                }
                if angles.is_empty() {
                    angles.extend(std::iter::repeat(TAU / div as f32).take(div));
                }
            }
            for (&n, a) in neighbors.iter().zip(&angles) {
                if eq(n, at) {
                    continue;
                }
                if visited.contains(&(n as *const _)) {
                    continue;
                }
                let (s, c) = a.sin_cos();
                init.rotate(s, c);
                let mut neigh = n.borrow_mut();
                neigh.set_coords(atom.coordinates + init);
                if eq(neigh.fragment.unwrap(), f) {
                    queue.push_back(n);
                    visited.insert(n as *const _);
                } else {
                    neigh.coordinates_set = false;
                }
                // TODO: CoordgenFragmentBuilder.cpp:772
                if let Some(dofs) = frag.atom_dofs.get(&(at as *const _)) {
                    for dof in dofs {
                        let mut dof = dof.borrow_mut();
                        if eq(dof.frag, f) {
                            dof.atoms.push(at);
                        }
                    }
                }
            }
            // avoid Z/E inversions
            // maybe add macrocycle
            'mac: {
                if let [ring] = atom.rings.as_slice() {
                    let ring = ring.borrow();
                    if !ring.is_macrocycle() && atom.neighbors.len() == 3 {
                        break 'mac;
                    }
                    for bond in &atom.bonds {
                        let bond = bond.borrow();
                        if bond.is_stereo() && !bond.is_terminal() {
                            break 'mac;
                        }
                    }
                    for &n in &atom.neighbors {
                        if eq(n, at) {
                            continue;
                        }
                        let neigh = n.borrow_mut();
                        if atom.shares_a_ring_with(&neigh).is_none() {
                            let dof = FragmentDof {
                                kind: FragmentDofKind::InvertBond {
                                    pivot: at,
                                    bound: n,
                                },
                                atoms: vec![n],
                                frag: f,
                                current_state: 0,
                            };
                            let dof = sketcher.interner.intern_frag_dof(dof);
                            frag.atom_dofs.entry(n as *const _).or_default().push(dof);
                            frag.dofs.push(dof);
                        }
                    }
                }
            }
            for &n in &atom.neighbors {
                let nr = n.borrow();
                if atom.shares_a_ring_with(&nr).is_none() && eq(nr.fragment.unwrap(), f) {
                    let dof = FragmentDof {
                        kind: FragmentDofKind::ScaleAtom(at),
                        atoms: vec![],
                        frag: f,
                        current_state: 0,
                    };
                    let dof = sketcher.interner.intern_frag_dof(dof);
                    frag.dofs.push(dof);
                }
            }
        }
    }
    // avoid interal clashes
    for &a1 in &frag.atoms {
        let mut a = a1.borrow_mut();
        if a.neighbors.len() != 1 {
            continue;
        }
        if a.fixed {
            continue;
        }
        if frag.atom_dofs.contains_key(&(a1 as *const _)) {
            continue;
        }
        for &a2 in &frag.atoms {
            if eq(a1, a2) {
                continue;
            }
            if frag.atom_dofs.contains_key(&(a2 as *const _)) {
                continue;
            }
            if a.neighbors.iter().any(|n| eq(*n, a2)) {
                continue;
            }
            let mut a2 = a2.borrow_mut();
            let d = a.coordinates - a2.coordinates;
            if d.0 > BOND_LENGTH * 0.5 {
                continue;
            }
            if d.1 > BOND_LENGTH * 0.5 {
                continue;
            }
            if d.sq_length() > BOND_LENGTH * BOND_LENGTH * 0.25 {
                continue;
            }
            let vec = a.coordinates
                - a.neighbors[0]
                    .try_borrow()
                    .map_or(a2.coordinates, |a| a.coordinates);
            a.coordinates -= vec * 0.3;
            if a2.neighbors.len() == 1 {
                a2.coordinates += vec;
                a2.coordinates.round(2);
            }
        }
    }
    if frag.parent.is_none() && frag.constrained {
        if !(frag.fixed || frag.constrained || frag.templated) {
            let mut coc = PointF::default();
            let mut cnc = PointF::default();
            let constrained = frag
                .atoms
                .iter()
                .copied()
                .chain(
                    frag.children
                        .iter()
                        .map(|c| c.borrow().parent.unwrap().1.borrow().end),
                )
                .filter(|a| a.borrow().constrained)
                .collect::<Vec<_>>();
            for a in &constrained {
                let a = a.borrow();
                coc += a.template;
                cnc += a.coordinates;
            }
            if !constrained.is_empty() {
                let len = constrained.len() as f32;
                coc /= len;
                cnc /= len;
            }
            let (v1, v2): (Vec<_>, Vec<_>) = constrained
                .iter()
                .map(|a| {
                    let a = a.borrow();
                    (a.template - coc, a.coordinates - cnc)
                })
                .unzip();
            let rot_mat = Sketcher::alignment_matrix(&v1, &v2);
            for (a, set) in frag.atoms.iter().map(|a| (*a, true)).chain(
                frag.children
                    .iter()
                    .map(|c| (c.borrow().parent.unwrap().1.borrow().end, false)),
            ) {
                let mut a = a.borrow_mut();
                let mut v = a.coordinates - cnc;
                v = PointF(
                    v.0 * rot_mat[0] + v.1 * rot_mat[1],
                    v.0 * rot_mat[2] + v.1 * rot_mat[3],
                );
                v.round(2);
                a.coordinates = v;
                a.coordinates_set = set;
            }
        }
    }
    if frag.fixed {
        for atom in &frag.atoms {
            atom.borrow_mut().set_coords_to_template();
        }
        if let Some((_, b)) = frag.parent {
            let b = b.borrow();
            b.start.borrow_mut().set_coords_to_template();
            b.end.borrow_mut().set_coords_to_template();
        }
        for c in &frag.children {
            let b = c.borrow().parent.unwrap().1.borrow();
            b.start.borrow_mut().set_coords_to_template();
            b.end.borrow_mut().set_coords_to_template();
        }
    }
    frag.store_coordinate_information();
}
pub(crate) fn order_chain_of_atoms<'a>(chain: &mut [AtomRef<'a>], start: AtomRef<'a>) {
    todo!()
}
fn build_ring(ring: RingRef<'_>) {
    todo!()
}

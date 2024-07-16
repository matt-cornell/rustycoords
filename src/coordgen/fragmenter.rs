use super::*;
use std::collections::VecDeque;
use std::ptr::eq;

pub fn split_into_fragments<'a>(mol: MoleculeRef<'a>, intern: &'a dyn Interner) {
    let mut m = mol.borrow_mut();
    let m = &mut *m;
    m.fragments.clear();
    for a in &m.atoms {
        a.borrow_mut().fragment = None;
    }
    if let [only] = m.atoms.as_slice() {
        let mut frag = Fragment::default();
        frag.atoms.push(*only);
        let interned = intern.intern_fragment(frag);
        Fragment::finalize(interned, intern);
        only.borrow_mut().fragment = Some(interned);
        m.fragments.push(interned);
    }
    for &bond in &m.bonds {
        let b = bond.borrow();
        if b.is_inter_fragment() {
            {
                let mut start = b.start.borrow_mut();
                if start.fragment.is_none() {
                    let mut frag = Fragment::default();
                    frag.atoms.push(b.start);
                    let interned = intern.intern_fragment(frag);
                    Fragment::finalize(interned, intern);
                    start.fragment = Some(interned);
                    m.fragments.push(interned);
                }
            }
            {
                let mut end = b.end.borrow_mut();
                if end.fragment.is_none() {
                    let mut frag = Fragment::default();
                    frag.atoms.push(b.end);
                    let interned = intern.intern_fragment(frag);
                    Fragment::finalize(interned, intern);
                    end.fragment = Some(interned);
                    m.fragments.push(interned);
                }
            }
        } else {
            if eq(b.start, b.end) {
                let mut a = b.start.borrow_mut();
                if a.fragment.is_none() {
                    let mut frag = Fragment::default();
                    frag.atoms.push(b.start);
                    let interned = intern.intern_fragment(frag);
                    Fragment::finalize(interned, intern);
                    a.fragment = Some(interned);
                    m.fragments.push(interned);
                }
                continue;
            }
            let mut start = b.start.borrow_mut();
            let mut end = b.end.borrow_mut();
            match (start.fragment, end.fragment) {
                (None, None) => {
                    let mut frag = Fragment::default();
                    frag.atoms.extend_from_slice(&[b.start, b.end]);
                    let interned = intern.intern_fragment(frag);
                    Fragment::finalize(interned, intern);
                    start.fragment = Some(interned);
                    end.fragment = Some(interned);
                    m.fragments.push(interned);
                }
                (Some(f), None) => {
                    f.borrow_mut().atoms.push(b.end);
                    end.fragment = Some(f);
                }
                (None, Some(f)) => {
                    f.borrow_mut().atoms.push(b.start);
                    start.fragment = Some(f);
                }
                (Some(far), Some(fbr)) => {
                    if eq(far, fbr) {
                        continue;
                    }
                    let mut fa = far.borrow_mut();
                    let mut fb = fbr.borrow_mut();
                    let old = fa.atoms.len();
                    fa.atoms.append(&mut fb.atoms);
                    for a in &fa.atoms[old..] {
                        a.borrow_mut().fragment = Some(far);
                    }
                    let i = m.fragments.iter().position(|f| eq(*f, fbr)).unwrap();
                    m.fragments.swap_remove(i);
                }
            }
        }
    }
    // initialize information
    if !m.fragments.is_empty() {
        for bond in &m.bonds {
            let b = bond.borrow();
            let start = b.start.borrow();
            let end = b.end.borrow();
            let sf = start.fragment.unwrap();
            let ef = end.fragment.unwrap();
            if eq(sf, ef) {
                sf.borrow_mut().bonds.push(bond);
            } else {
                sf.borrow_mut().inter_frags.push(bond);
                ef.borrow_mut().inter_frags.push(bond);
            }
        }
        for ring in &m.rings {
            ring.borrow().atoms[0]
                .borrow()
                .fragment
                .unwrap()
                .borrow_mut()
                .rings
                .push(ring);
        }
        for frag in &m.fragments {
            let mut frag = frag.borrow_mut();
            let is_chain = frag.atoms.len() <= 3
                && frag.atoms.iter().all(|a| {
                    let a = a.borrow();
                    a.bonds.len() <= 3 && a.rings.is_empty()
                })
                && frag.bonds.iter().all(|b| b.borrow().bond_order <= 2);
            frag.is_chain = is_chain;
        }
        let mut fixed = false;
        let mut constrained = false;
        for frag in &m.fragments {
            let mut frag = frag.borrow_mut();
            let mut ffixed = false;
            let mut fconst = false;
            for a in &frag.atoms {
                let a = a.borrow();
                if a.fixed {
                    ffixed = true;
                    if fconst {
                        break;
                    }
                }
                if a.constrained {
                    fconst = true;
                    if ffixed {
                        break;
                    }
                }
            }
            frag.fixed = ffixed;
            frag.constrained = fconst;
            fixed |= ffixed;
            constrained |= fconst;
        }
        m.has_fixed_frags = fixed;
        m.has_constrained_frags = constrained;
        // find main fragment
        let mut main_fragment = m
            .fragments
            .iter()
            .copied()
            .min_by(|l, r| {
                let l = l.borrow();
                let r = r.borrow();
                l.count_fixed()
                    .cmp(&r.count_fixed())
                    .then_with(|| l.count_constrained().cmp(&r.count_constrained()))
                    .then(l.rings.len().cmp(&r.rings.len()))
                    .then(l.atoms.len().cmp(&r.atoms.len()))
                    .then(l.inter_frags.len().cmp(&r.inter_frags.len()))
                    .then_with(|| l.count_heavy().cmp(&r.count_heavy()))
                    .then_with(|| l.count_double_bonds().cmp(&r.count_double_bonds()))
            })
            .unwrap();
        let mut queue = VecDeque::new();
        // consider chains
        'chains: {
            for f in &m.fragments {
                let f = f.borrow();
                if f.fixed || f.constrained {
                    break 'chains;
                }
            }
            let mut longest_start = None;
            let mut longest_len = 0;
            for &f in &m.fragments {
                let frag = f.borrow();
                if !frag.is_chain {
                    continue;
                }
                let mut chain_it = frag.inter_frags.iter().filter(|b| {
                    let b = b.borrow();
                    let start = b.start.borrow();
                    let end = b.end.borrow();
                    if eq(start.fragment.unwrap(), f) {
                        end.fragment.unwrap().borrow().is_chain
                    } else {
                        start.fragment.unwrap().borrow().is_chain
                    }
                });
                if chain_it.next().is_some() && chain_it.next().is_some() {
                    continue;
                }
                queue.push_back(f);
                let mut rec = f;
                while let Some(last) = queue.pop_front() {
                    rec = last;
                    let l = last.borrow();
                    for bond in &l.inter_frags {
                        let b = bond.borrow();
                        let child = {
                            let start = b.start.borrow();
                            let end = b.end.borrow();
                            if eq(start.fragment.unwrap(), last) {
                                end.fragment.unwrap()
                            } else {
                                start.fragment.unwrap()
                            }
                        };
                        if eq(child, f) {
                            continue;
                        }
                        if l.parent.map_or(false, |p| eq(p.0, child)) {
                            continue;
                        }
                        let mut c = child.borrow_mut();
                        if !c.is_chain {
                            continue;
                        }
                        if c.parent.is_none() {
                            c.parent = Some((last, bond));
                        }
                        queue.push_back(child);
                    }
                }
                let mut it = std::iter::successors(Some(rec), |r| r.borrow().parent.map(|r| r.0));
                let Some(first) = it.next() else { continue };
                let len = it.count() + 1;
                if len > longest_len {
                    longest_len = len;
                    longest_start = Some(first);
                }
            }
            let acceptable_len = match main_fragment.borrow().rings.len() {
                0 => 1,
                1 => 5,
                2 => 8,
                3 => 10,
                _ => 12,
            };
            if longest_len > acceptable_len {
                main_fragment = longest_start.unwrap();
            }
        }
        m.main_fragment = Some(main_fragment);
        // add parent info
        for f in &m.fragments {
            f.borrow_mut().parent = None;
        }
        queue.push_back(main_fragment);
        while let Some(last) = queue.pop_front() {
            let mut l = last.borrow_mut();
            let l = &mut *l;
            for bond in &l.inter_frags {
                let mut b = bond.borrow_mut();
                let start = b.start.borrow();
                let end = b.end.borrow();
                let child = if eq(start.fragment.unwrap(), last) {
                    end.fragment.unwrap()
                } else {
                    let child = start.fragment.unwrap();
                    let b = &mut *b;
                    std::mem::swap(&mut b.start, &mut b.end);
                    b.is_reversed = !b.is_reversed;
                    child
                };
                if eq(child, main_fragment) {
                    continue;
                }
                let mut c = child.borrow_mut();
                if c.parent.is_none() {
                    c.parent = Some((last, bond));
                    l.children.push(child);
                    queue.push_back(child);
                }
            }
        }
        // order fragments
        let old_len = m.fragments.len();
        m.fragments.clear();
        queue.push_back(main_fragment);
        while let Some(frag) = queue.pop_front() {
            m.fragments.push(frag);
            queue.extend(frag.borrow().children.iter().copied());
        }
        assert_eq!(m.fragments.len(), old_len);
        for frag in &m.fragments {
            let mut frag = frag.borrow_mut();
            frag.constrained_flip = false;
            let mut it = frag.atoms.iter().filter(|a| a.borrow().constrained);
            if it.next().is_none() {
                continue;
            }
            if frag.atoms.len() == 1 {
                frag.constrained_flip = frag.children.iter().any(|f| f.borrow().constrained);
            } else {
                frag.constrained_flip = it.next().is_some();
            }
        }
    }
}

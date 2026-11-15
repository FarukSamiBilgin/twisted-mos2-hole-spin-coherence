# -*- coding: utf-8 -*-
"""kpts.out + proj.out ayristirici: bant sirasi (Gamma vs K VBM; K vs Q CBM)
ve Gamma VBM orbital bilesimi (Mo d_z2 / S p_z agirliklari)."""
import re, json, sys, os

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'results', 'dft')  # QE outputs live in results/dft
NELEC = 52          # 2x Mo(14) + 4x S(6)
VB = NELEC // 2     # 26. bant = VBM, 27. = CBM
KLABELS = ['Gamma', 'Q(0.167)', 'Q(0.200)', 'Q(0.250)', 'K']

def parse_kpts(path):
    txt = open(path, encoding='utf-8', errors='ignore').read()
    blocks = re.findall(r'k =\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*?bands \(ev\):\s*\n\n((?:\s+[-\d.]+)+)',
                        txt, re.S)
    out = []
    for kx, ky, kz, bandtxt in blocks:
        ev = [float(x) for x in bandtxt.split()]
        out.append(dict(k=(float(kx), float(ky)), ev=ev))
    return out

def parse_proj_states(path):
    """state #   N: atom   A (El ), wfc  W (l=L m= M) -> harita"""
    states = {}
    for m in re.finditer(r'state #\s*(\d+):\s*atom\s+(\d+)\s*\((\w+)\s*\),\s*wfc\s+(\d+)\s*\(l=(\d)\s*m=\s*(\d+)\)',
                         open(path, encoding='utf-8', errors='ignore').read()):
        n, atom, el, wfc, l, mm = m.groups()
        states[int(n)] = dict(atom=int(atom), el=el, l=int(l), m=int(mm), wfc=int(wfc))
    return states

def parse_proj_band(path, kindex, band):
    """kindex. k-noktasinda (1-tabanli) 'band'. bandin psi = ... satirini topla."""
    txt = open(path, encoding='utf-8', errors='ignore').read()
    kblocks = txt.split(' k =')
    if len(kblocks) <= kindex:
        raise RuntimeError('k blogu yok')
    blk = kblocks[kindex]
    m = re.search(r'==== e\(\s*%d\)\s*=\s*([-\d.]+)\s*eV ====\s*\n(.*?)\|psi\|\^2' % band, blk, re.S)
    if not m:
        raise RuntimeError('bant %d bulunamadi' % band)
    ene = float(m.group(1))
    coefs = re.findall(r'([\d.]+)\*\[#\s*(\d+)\]', m.group(2))
    return ene, [(float(c), int(s)) for c, s in coefs]

def orbital_summary(coefs, states):
    agg = {}
    for c, s in coefs:
        st = states[s]
        if st['el'].startswith('Mo') and st['l'] == 2:
            key = 'Mo_d_z2' if st['m'] == 1 else 'Mo_d_diger'
        elif st['el'].startswith('S') and st['l'] == 1:
            key = 'S_p_z' if st['m'] == 1 else 'S_p_xy'
        else:
            key = st['el'] + '_l%d' % st['l']
        agg[key] = agg.get(key, 0.0) + c
    return agg

if __name__ == '__main__':
    kp = parse_kpts(os.path.join(D, 'kpts.out'))
    print('=== Bant sirasi (eV) ===')
    res = {}
    for lab, d in zip(KLABELS, kp):
        vbm, cbm = d['ev'][VB-1], d['ev'][VB]
        res[lab] = dict(VBM=vbm, CBM=cbm)
        print('%-9s VBM(b%d) = %9.4f   CBM(b%d) = %9.4f' % (lab, VB, vbm, VB+1, cbm))
    dVBM = res['Gamma']['VBM'] - res['K']['VBM']
    cbm_min = min(res, key=lambda L: res[L]['CBM'])
    print()
    print('VBM(Gamma) - VBM(K) = %+.3f eV  (pozitifse Gamma kazanir)' % dVBM)
    print('CBM minimumu: %s' % cbm_min)
    gap_ind = res[cbm_min]['CBM'] - res['Gamma']['VBM']
    gap_dir = res['K']['CBM'] - res['K']['VBM']
    print('Indirekt aralik Gamma->%s = %.3f eV; direkt K->K = %.3f eV' % (cbm_min, gap_ind, gap_dir))

    states = parse_proj_states(os.path.join(D, 'proj.out'))
    ene, coefs = parse_proj_band(os.path.join(D, 'proj.out'), 1, VB)
    agg = orbital_summary(coefs, states)
    tot = sum(agg.values())
    print()
    print('=== Gamma VBM (E=%.3f eV) orbital bilesimi ===' % ene)
    for k in sorted(agg, key=agg.get, reverse=True):
        print('  %-10s %6.3f  (%5.1f%%)' % (k, agg[k], 100*agg[k]/tot))
    wd = agg.get('Mo_d_z2', 0.0) / tot
    wp = agg.get('S_p_z', 0.0) / tot
    print()
    print('Normalize: w_d(Mo d_z2) = %.3f, w_p(S p_z) = %.3f, diger = %.3f'
          % (wd, wp, 1-wd-wp))
    json.dump(dict(bands=res, dVBM_eV=dVBM, cbm_at=cbm_min,
                   gap_indirect_eV=gap_ind, gap_direct_eV=gap_dir,
                   gamma_vbm_composition={k: v/tot for k, v in agg.items()},
                   w_d=wd, w_p=wp),
              open(os.path.join(D, 'dft_sonuc.json'), 'w'), indent=2)
    print('dft_sonuc.json yazildi.')

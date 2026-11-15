# -*- coding: utf-8 -*-
"""
KOL B - Izotop-saf bukulmus MoS2: hyperfine-sinirli T2* hesabi
================================================================
Boru hatti:
  1. Gamma-vadisi continuum modeli (Angeli-MacDonald tarzi):
     H = -hbar^2 k^2 / 2m* + Delta(r),  Delta(r) = 2V sum_{i=1,3,5} cos(g_i.r + psi)
     Duzlem-dalga koseginlestirme -> moire bantlari, kuyu derinligi, seviye araligi.
  2. Lokalize taban durumu zarfi -> etkin cekirdek sayisi N.
  3. Overhauser dagilimi sigma_B = A_eff * sqrt(f*N)/N  (gorev tanimindaki form),
     T2* = hbar / (g mu_B sigma_B_alan)  [enerji cinsinden: T2* = hbar/sigma_E].
  4. Kiyas tablosu: dogal MoS2 / 98Mo32S2 / GaAs / dogal Si / 28Si.

Parametre blogundaki degerler literatur ankrajli; kaynaklar yanlarinda.
"""
import numpy as np
from scipy.linalg import eigh
import json, csv, io, sys

# ----------------------------------------------------------------------
# Sabitler (SI + pratik birimler)
# ----------------------------------------------------------------------
hbar   = 1.054571817e-34   # J s
me     = 9.1093837015e-31  # kg
muB    = 9.2740100783e-24  # J/T
eV     = 1.602176634e-19   # J
meV    = 1e-3 * eV
ueV    = 1e-6 * eV
nm     = 1e-9

# ----------------------------------------------------------------------
# PARAMETRE BLOGU  (literatur guncellemesi icin tek nokta)
# ----------------------------------------------------------------------
P = dict(
    # --- MoS2 orgu / moire geometrisi ---
    a0        = 0.3182,    # nm (Angeli-MacDonald Tablo I, dogrulandi)
    theta_deg = 2.0,       # TASARIM ACISI: dot rejimi (s-bant genisligi ~7 meV,
                           # s-p araligi ~32 meV, l~1.2 nm); 2-4 tarama main() icinde

    # --- Gamma-vadisi continuum modeli (Angeli-MacDonald, PNAS 118, e2021826118;
    #     arXiv:2008.01735v2 Tablo I — birincil kaynaktan dogrulandi) ---
    # H = hbar^2 k^2/2m* + U_delik(r);  U_delik = -Delta(r),
    # Delta(r) = sum_s sum_{j=1..6} V_s exp(i g_j^s.r + phi_s), phi_s = pi (tum kabuklar)
    # -> delik potansiyeli Fourier katsayisi her kabuk vektorunde +V_s.
    # Delikler AB/BA (Bernal) bolgelerinde lokalize -> bal petegi orgusu (AM ana sonucu;
    # Naik-Jain PRL 121, 266401 DFT ile uyumlu).
    m_eff     = 0.9,       # m0 (AM Tablo I, dogrulandi; monolayer degeri ~2.6 DEGIL —
                           #     bilayer Gamma antibonding dali daha hafif)
    V_shells  = (39.45, 6.5, 10.0),   # meV, V1..V3 (AM Tablo I)
    psi_deg   = 180.0,     # phi = pi, tum kabuklar (AM Tablo I)

    # --- hyperfine (Avdeev & Smirnov, Nanoscale Adv. 1, 2624 (2019); arXiv:1810.06449) ---
    # Element olcek sabitleri (dipolar, Fermi-kontakt ~0; s-karisimi ihmal):
    #   A_M(Mo) = 0.57 ueV (95Mo ve 97Mo ortak), A_X(S) = 0.75 ueV
    # Bant-cozumlu A_zz (Tablo II-III carpimlari):
    #   K iletim bandi (Mo d_z2):        A_zz = (4/7)A_M = 0.326;  S: (12/5)A_X = 1.8 (isaret -)
    #   K deger bandi  (Mo d_x2y2/dxy):  A_zz = (24/7)A_M = 1.954 (SAF ISING); S: (8/5)A_X = 1.2 (ISING)
    # GAMMA deger bandi: yayimlanmis deger YOK. Kestirim (bu calisma):
    #   Gamma VBM ~ w_d Mo d_z2 + w_p S p_z, L_z=0 -> yalnizca spin-dipolar.
    #   d_z2 icin (4/7)A_M = 0.326 ueV; p_z icin (4/5)A_X = 0.60 ueV.
    A_scale_Mo = 0.57,     # ueV (dogrulandi)
    A_scale_S  = 0.75,     # ueV (dogrulandi)
    # w_d/w_p: bu calismanin QE-PBE hesabi (kol_B/dft/, bilayer 2H MoS2, PAW,
    # projwfc Gamma VBM): Mo d_z2 %56.7, S p_z %39.6, diger %3.7 -> DOGRULANDI
    w_d        = 0.567,    # Gamma VBM Mo d_z2 agirligi (DFT, bu calisma)
    w_p        = 0.396,    # Gamma VBM S p_z agirligi (DFT, bu calisma)
    g_factor   = 2.0,      # Gamma deligi: L_z=0, vadi yok -> spin-only g~2

    # --- izotop bilesimleri (Avdeev-Smirnov Tablo IV ile dogrulandi) ---
    f_Mo95    = 0.1592,    # I=5/2
    f_Mo97    = 0.0955,    # I=5/2
    f_S33     = 0.0076,    # I=3/2
    # saflastirilmis hedef (zenginlestirme varsayimlari; duyarlilik taranir)
    f_Mo_pure = 2e-4,      # 98Mo zenginlestirme sonrasi kalinti spinli Mo (%0.02 varsayim)
    f_S_pure  = 1e-5,      # 32S zenginlestirme sonrasi kalinti 33S
)

def gamma_A_effective(P, senaryo='merkez'):
    """Gamma-vadisi deligi icin etkin A_zz (ueV) — kestirim + duyarlilik bandi.
    merkez   : orbital-agirlikli d_z2/p_z spin-dipolar kestirimi
    kotumser : K-deger bandi Ising degerleri (ust sinir; Gamma'da d_x2y2 yok ama
               kestirim belirsizligini kapsamak icin)
    iyimser  : merkez / 2 (radyal integral ve agirlik belirsizligi alt banti)"""
    A_Mo_c = P['w_d'] * (4.0/7.0) * P['A_scale_Mo']
    A_S_c  = P['w_p'] * (4.0/5.0) * P['A_scale_S']
    if senaryo == 'merkez':
        return A_Mo_c, A_S_c
    if senaryo == 'kotumser':
        return (24.0/7.0)*P['A_scale_Mo'], (8.0/5.0)*P['A_scale_S']
    if senaryo == 'iyimser':
        return A_Mo_c/2.0, A_S_c/2.0
    raise ValueError(senaryo)

# ----------------------------------------------------------------------
# 1) Gamma-vadisi continuum modeli — duzlem dalga koseginlestirme
# ----------------------------------------------------------------------
def moire_geometry(a0_nm, theta_deg):
    theta = np.deg2rad(theta_deg)
    aM = a0_nm / (2*np.sin(theta/2))          # moire periyodu (nm)
    gM = 4*np.pi / (np.sqrt(3)*aM)            # moire ters orgu buyuklugu (1/nm)
    # ilk kabuk g vektorleri (6 adet)
    gs = np.array([[gM*np.cos(a), gM*np.sin(a)]
                   for a in (np.pi/6 + np.arange(6)*np.pi/3)])
    return aM, gM, gs

def shell_vectors(gs_first):
    """Uc kabugun 6'sar vektoru: |g| = gM, sqrt(3) gM, 2 gM."""
    b1, b2 = gs_first[0], gs_first[1]
    gM = np.linalg.norm(b1)
    shells = {1: [], 2: [], 3: []}
    for i in range(-3, 4):
        for j in range(-3, 4):
            if i == 0 and j == 0:
                continue
            v = i*b1 + j*b2
            L = np.linalg.norm(v)/gM
            if abs(L-1) < 1e-6:            shells[1].append(v)
            elif abs(L-np.sqrt(3)) < 1e-6: shells[2].append(v)
            elif abs(L-2) < 1e-6:          shells[3].append(v)
    return {s: np.array(vs) for s, vs in shells.items()}

def hole_fourier_coeffs(P, gs_first):
    """Delik potansiyeli U(r) = -Delta(r) Fourier katsayilari.
    Delta katsayisi kabuk s'in her vektorunde V_s e^{i phi_s} (AM Es. 2);
    delik icin isaret ters: c_g = -V_s e^{i phi_s}. phi = pi -> c_g = +V_s."""
    phi = np.deg2rad(P['psi_deg'])
    shells = shell_vectors(gs_first)
    Vc = {}
    for s, V in zip((1, 2, 3), P['V_shells']):
        for g in shells[s]:
            Vc[tuple(np.round(g, 6))] = -V * np.exp(1j*phi)
    return Vc

def hole_potential_real(P, gs_first, X, Y):
    """U(r) gercek uzayda (meV)."""
    phi = np.deg2rad(P['psi_deg'])
    shells = shell_vectors(gs_first)
    U = np.zeros_like(X)
    for s, V in zip((1, 2, 3), P['V_shells']):
        for g in shells[s]:
            U += -V * np.cos(g[0]*X + g[1]*Y + phi)
    return U

def build_hamiltonian(kvec, gs_first, m_eff, Vc, nshell=6):
    """Duzlem-dalga tabani: moire ters orgu {G}; delik: +k^2/2m + U(r)."""
    b1, b2 = gs_first[0], gs_first[1]
    G = []
    for i in range(-nshell, nshell+1):
        for j in range(-nshell, nshell+1):
            G.append(i*b1 + j*b2)
    G = np.array(G)
    NG = len(G)
    kin = np.array([np.sum((kvec + g)**2) for g in G])  # 1/nm^2
    kin_meV = (hbar**2 * (kin * 1e18) / (2*m_eff*me)) / meV
    H = np.diag(kin_meV).astype(complex)
    for a in range(NG):
        for b in range(NG):
            key = tuple(np.round(G[a] - G[b], 6))
            if key in Vc:
                H[a, b] += Vc[key]
    return H, G

def solve_moire(P, klist, nshell=6):
    aM, gM, gs = moire_geometry(P['a0'], P['theta_deg'])
    Vc = hole_fourier_coeffs(P, gs)
    bands = []
    for k in klist:
        H, G = build_hamiltonian(k, gs, P['m_eff'], Vc, nshell)
        w = eigh(H, eigvals_only=True)
        bands.append(np.sort(w)[:8])
    return aM, gM, np.array(bands)

def k_path(gs, npts=30):
    """gamma_M -> K_M -> M_M -> gamma_M moire BZ yolu."""
    gM = np.linalg.norm(gs[0])
    K = np.array([gM/np.sqrt(3), 0.0])
    M = np.array([gM/2*np.cos(np.pi/6), gM/2*np.sin(np.pi/6)])
    Gm = np.array([0.0, 0.0])
    seg = lambda p, q: [p + (q-p)*t for t in np.linspace(0, 1, npts, endpoint=False)]
    pts = seg(Gm, K) + seg(K, M) + seg(M, Gm) + [Gm]
    return np.array(pts)

# ----------------------------------------------------------------------
# 2) Kuyu analizi: derinlik, harmonik seviye araligi, zarf genisligi
# ----------------------------------------------------------------------
def well_analysis(P):
    """Delik potansiyeli U(r) uzerinde: Bernal minimumu (bal petegi alt-orgusu),
    AA tepesine gore toplam derinlik, komsu Bernal kuyusuna eyer bariyeri,
    minimum etrafinda harmonik seviye araligi ve zarf uzunlugu."""
    aM, gM, gs = moire_geometry(P['a0'], P['theta_deg'])
    L = aM * 1.6
    n = 320
    xs = np.linspace(-L, L, n); ys = np.linspace(-L, L, n)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    U = hole_potential_real(P, gs, X, Y)
    depth_AA = U.max() - U.min()                # AA tepesi - Bernal dibi
    imin = np.unravel_index(np.argmin(U), U.shape)
    r_min = np.array([X[imin], Y[imin]])
    # komsu Bernal minimumu: bal petegi bag uzunlugu aM/sqrt(3), en yakin minimumu ara
    dist = np.sqrt((X-r_min[0])**2 + (Y-r_min[1])**2)
    mask = (dist > 0.4*aM/np.sqrt(3)) & (dist < 1.4*aM/np.sqrt(3))
    Umask = np.where(mask, U, np.inf)
    inb = np.unravel_index(np.argmin(Umask), U.shape)
    r_nb = np.array([X[inb], Y[inb]])
    # eyer bariyeri: iki minimum arasi dogru boyunca maksimum (yaklasik)
    ts = np.linspace(0, 1, 200)
    line = np.array([r_min + t*(r_nb - r_min) for t in ts])
    gsh = shell_vectors(gs)
    phi = np.deg2rad(P['psi_deg'])
    Uline = np.zeros(len(ts))
    for s, V in zip((1, 2, 3), P['V_shells']):
        for g in gsh[s]:
            Uline += -V*np.cos(line[:, 0]*g[0] + line[:, 1]*g[1] + phi)
    barrier = Uline.max() - U.min()             # komsu kuyuya tunel bariyeri
    # harmonik acilim (iki yonde sayisal egrilik, ortalama)
    dx = xs[1]-xs[0]
    d2x = (U[imin[0]+1, imin[1]] - 2*U[imin] + U[imin[0]-1, imin[1]]) / dx**2
    d2y = (U[imin[0], imin[1]+1] - 2*U[imin] + U[imin[0], imin[1]-1]) / dx**2
    kappa = 0.5*(d2x + d2y)                     # meV/nm^2
    m = P['m_eff']*me
    omega = np.sqrt(max(kappa, 1e-12)*meV*1e18/m)
    hw = hbar*omega/meV
    l_env = np.sqrt(hbar/(m*omega))/nm
    return dict(aM=aM, depth_meV=depth_AA, barrier_meV=barrier,
                hbar_omega_meV=hw, l_env_nm=l_env, kappa=kappa,
                Umin=U.min(), Umax=U.max(), U_bernal=float(U[imin]))

def ground_state_participation(P, nshell=6):
    """k=Gamma_M taban durumunun gercek-uzay zarfi ve katilim-orani N_eff.
    N_eff = 1 / (A_cell * int |psi|^4 d2r), psi moire hucresine normalize.
    Ayrica zarf RMS yaricapi dondurur."""
    aM, gM, gs = moire_geometry(P['a0'], P['theta_deg'])
    Vc = hole_fourier_coeffs(P, gs)
    H, G = build_hamiltonian(np.zeros(2), gs, P['m_eff'], Vc, nshell)
    w, v = eigh(H)
    c = v[:, 0]                                   # taban durum duzlem-dalga katsayilari
    # gercek uzay orneklemesi (moire birim hucresi): a1, a2 orgu vektorleri
    # ters orgu b1,b2 = gs[0],gs[1] -> a_i . b_j = 2pi delta_ij
    b1, b2 = gs[0], gs[1]
    Bm = np.array([b1, b2]).T                     # sutunlar b1, b2
    Am = 2*np.pi*np.linalg.inv(Bm)                # satirlar a1, a2 (a_i.b_j = 2pi delta)
    a1, a2 = Am[0], Am[1]
    n = 90
    u = np.linspace(0, 1, n, endpoint=False)
    density = np.zeros((n, n))
    XY = np.zeros((n, n, 2))
    for i, ui in enumerate(u):
        for j, uj in enumerate(u):
            r = ui*a1 + uj*a2
            XY[i, j] = r
            phase = G @ r
            amp = np.sum(c*np.exp(1j*phase))
            density[i, j] = np.abs(amp)**2
    A_uc = abs(a1[0]*a2[1] - a1[1]*a2[0])         # moire hucre alani nm^2
    dA = A_uc / n**2
    density /= density.sum()*dA                   # int |psi|^2 = 1
    ipr = np.sum(density**2)*dA                   # int |psi|^4
    A_cell = np.sqrt(3)/2 * P['a0']**2            # atomik hucre alani
    N_eff = 1.0/(A_cell*ipr)                      # etkin hucre sayisi (Bloch durumu)
    # bal petegi: Bloch durumu hucre basina 2 Bernal kuyusuna yayilir;
    # tepe sayisini yogunluk esiklemesiyle say (yerel maksimumlar, 0.5*max esik)
    pk = density > 0.5*density.max()
    from scipy import ndimage
    n_peaks = int(ndimage.label(pk)[1])
    return dict(N_eff_cells=N_eff, E0_meV=w[0],
                gap_sp_meV=w[1]-w[0], A_uc_nm2=A_uc, n_peaks=n_peaks)

# ----------------------------------------------------------------------
# 3) Hyperfine: N, sigma, T2*
# ----------------------------------------------------------------------
def hyperfine_T2(A_ueV_list, f_list, I_list, N_nuc_per_species, g_factor):
    """
    sigma_E^2 = sum_species  f * N * (A/N)^2 * I(I+1)/3   [zarf-ortalama, uniform yaklasim]
    Gorevdeki form sigma = A sqrt(fN)/N tek tur icin; burada I(I+1)/3 agirlikli genel hali.
    T2* = hbar / sigma_E  (gorevdeki tanim; sqrt(2) konvansiyon farklari not edilir)
    """
    var = 0.0
    for A, f, I, N in zip(A_ueV_list, f_list, I_list, N_nuc_per_species):
        if N <= 0: continue
        var += f * N * (A/N)**2 * (I*(I+1)/3.0)
    sigma_ueV = np.sqrt(var)
    if sigma_ueV == 0:
        return dict(sigma_ueV=0.0, sigma_mT=0.0, T2_ns=np.inf)
    sigma_J = sigma_ueV * ueV
    T2 = hbar / sigma_J
    sigma_mT = sigma_J / (g_factor*muB) * 1e3
    return dict(sigma_ueV=sigma_ueV, sigma_mT=sigma_mT, T2_ns=T2*1e9)

def N_nuclei(P, l_env_nm):
    """2D Gauss zarf |psi|^2 = exp(-r^2/l^2)/(pi l^2); etkin N = 2*pi*l^2/A_cell
    (katilim orani tanimiyla N_eff = 1/(A_cell * int|psi|^4) = 2 pi l^2 / A_cell)."""
    A_cell = np.sqrt(3)/2 * P['a0']**2   # nm^2, birim hucre (1 Mo + 2 S)
    N_cell = 2*np.pi*l_env_nm**2 / A_cell
    return dict(N_Mo=N_cell, N_S=2*N_cell, N_cell=N_cell)

# ----------------------------------------------------------------------
# hBN kapsul cekirdeklerinin uzaydan (dipolar) hiperfini — sizintisiz taban
# ----------------------------------------------------------------------
# Delik spin yogunlugu MoS2 bilayer duzlemlerinde Gauss zarf (l), hBN cekirdekleri
# surekli tabaka (n_s = 1/hucre, her tur icin). Secular dipol cekirdegi
#   B||z : K(q) = -2 pi q e^{-q h}       ;  B||x : K(q) = 2 pi (q_x^2/q) e^{-q h}
# (in-plane: |K|^2 acisal ortalama 3/8 (2 pi q)^2 e^{-2qh}).
# sigma^2 = mu^2 * sum_tur w_tur * n_s * sum_katman int d^2q/(2pi)^2 |rho(q) K(q)|^2
# rho(q) = exp(-q^2 l^2/4) * sum_duzlem w_p e^{-q (z_n - z_p)}.
# Bagimsiz uygulama; MIMARI/sorun_cozum_v1/S6_spin_kontrol/s13_ek.py ile capraz kontrol.
HBN_NUC = {'11B': (0.801, 1.5, 13.6629e6), '10B': (0.199, 3.0, 4.5751e6),
           '14N': (0.99636, 1.0, 3.0777e6), '15N': (0.00364, 0.5, -4.3173e6)}
A_HBN = 0.2504                      # nm, hBN orgu sabiti
# bilayer 2H-MoS2 spin-yogunlugu duzlemleri (z [A], agirlik): Mo d_z2 w_d/2 iki katman,
# S p_z w_p/2 ic + dis S; dis S'ye 0.099 (yari-yariya ic/dis) — duyarlilik dusuk
def _mos2_planes(P, w_outerS=0.099):
    wd, wp = P['w_d'], P['w_p']
    wi = wp / 2 - w_outerS
    return [(3.07, wd / 2), (-3.07, wd / 2), (1.49, wi), (-1.49, wi),
            (4.65, w_outerS), (-4.65, w_outerS)]

def hbn_dipolar_T2(P, l_env_nm, field='z', spacer_nm=0.0, nlayers=6, sides=2):
    """hBN kapsul dipolar Overhauser tabani. Donus: dict(sigma_ueV, T2_ns)."""
    from scipy.integrate import quad
    h_pl = 6.62607015e-34
    mu = 1e-7 * P['g_factor'] * muB                          # (mu0/4pi) g muB
    wsum = sum(f * I * (I + 1) / 3 * gm ** 2 for f, I, gm in HBN_NUC.values())
    n_s = 1.0 / (np.sqrt(3) / 2 * A_HBN ** 2)                 # site/nm^2 (her tur)
    ang = 1.0 if field == 'z' else 3.0 / 8.0
    planes = _mos2_planes(P)
    l = l_env_nm
    var = 0.0
    for k in range(nlayers):
        zn = (4.65 + 3.3 + 3.33 * k) * 0.1 + spacer_nm        # nm (MoS2 merkezinden)
        def f(q):
            s = sum(w * np.exp(-q * (zn - zp * 0.1)) for zp, w in planes)
            return q / (2 * np.pi) * ang * (2 * np.pi * q) ** 2 * \
                np.exp(-q * q * l * l / 2) * s * s
        var += quad(f, 0, 60 / l, limit=300)[0] * n_s * 1e54   # nm^-6 -> m^-6
    sigma_Hz = np.sqrt(sides * var * mu ** 2 * wsum)
    sigma_ueV = sigma_Hz * h_pl / eV * 1e6
    return dict(sigma_ueV=sigma_ueV, T2_ns=hbar / (sigma_ueV * ueV) * 1e9)

def combine_T2(*T2s_ns):
    """Bagimsiz quasi-static Gauss banyolar: 1/T2*^2 toplanir."""
    return 1.0 / np.sqrt(sum(1.0 / t ** 2 for t in T2s_ns))

# ----------------------------------------------------------------------
# Gorev 5: kalan dekoherans tabani (hyperfine olduktan sonra)
# ----------------------------------------------------------------------
def kalan_dekoherans(B_list=(0.1, 0.5, 1.0, 2.0, 5.0)):
    """Fonon T1: Pearce & Burkard, 2D Mater. 4, 025114 (2017), MoS2 QD R=50nm,
    admixture mekanizmasi, Gamma_ph ~ B^6; ankraj T1(1T) ~ 1 s (Fig. 3 okumasi).
    Gamma-vadisi deligi icin bu ust sinirdir: Gamma'da SOC ~0 ve vadi yok ->
    admixture kanali daha da zayif, gercek T1 daha uzun beklenir.
    Yuk gurultusu: TMD QD icin olcum yok (literatur taramasi). Parametrizasyon:
      iyimser  : Si-esdeger stack, 0.3 ueV/sqrt(Hz) @1Hz  -> Yoneda-2018 tipi
                 EDSR mimarisinde T2*_yuk ~ 20 us olcekle
      kotumser : 10x gurultulu stack -> T2*_yuk ~ 2 us
    (dogrudan olcum yerine olcekleme; belirsizlik acikca isaretlenir)"""
    fonon = [dict(B_T=B, T1_s=1.0*(1.0/B)**6) for B in B_list]
    yuk = dict(iyimser_us=20.0, kotumser_us=2.0,
               kaynak_iyimser='Yoneda et al. Nat.Nano 13,102 (2018): 28Si EDSR, '
                              'T2*=20 us, yuk-gurultusu-sinirli',
               kaynak_kotumser='ayni mimari, 10x gurultu varsayimi (TMD stack olgunlugu)')
    return dict(fonon=fonon, yuk=yuk)

def zenginlestirme_taramasi(P, N_Mo, N_S, senaryo='merkez'):
    """T2* vs kalinti spinli izotop orani f_res (Mo ve S ayni carpanla seyreltilir)."""
    A_Mo, A_S = gamma_A_effective(P, senaryo)
    f_Mo_nat = P['f_Mo95'] + P['f_Mo97']
    out = []
    for f_res in np.logspace(-7, np.log10(f_Mo_nat), 40):
        olcek = f_res / f_Mo_nat
        r = hyperfine_T2([A_Mo, A_S], [f_res, P['f_S33']*olcek],
                         [2.5, 1.5], [N_Mo, N_S], P['g_factor'])
        out.append(dict(f_res=f_res, T2_ns=r['T2_ns'], sigma_ueV=r['sigma_ueV']))
    return out

# ----------------------------------------------------------------------
# Ana akis
# ----------------------------------------------------------------------
def main(P, theta_scan=(2.0, 2.5, 3.0, 3.5, 4.0)):
    out = {}
    # --- bant yapisi calisma acisinda ---
    aM, gM, gs = moire_geometry(P['a0'], P['theta_deg'])
    kpts = k_path(gs, npts=24)
    _, _, bands = solve_moire(P, kpts)
    flat_bw = bands[:, 0].max() - bands[:, 0].min()
    gap_01  = bands[:, 1].min() - bands[:, 0].max()
    out['bands'] = dict(flat_bw_meV=flat_bw, gap01_meV=gap_01,
                        e0_min=bands[:,0].min(), aM_nm=aM)
    # --- kuyu ---
    w = well_analysis(P)
    out['well'] = w
    # --- taban durum katilim orani (capraz kontrol) ---
    gp = ground_state_participation(P)
    out['ground_state'] = gp
    # Birincil N: tek-kuyu harmonik zarf (qubit tek Bernal kuyusunda tanimli).
    # Capraz kontrol: Bloch katilim N_eff / tepe sayisi ~ tek-kuyu N olmali.
    A_cell = np.sqrt(3)/2*P['a0']**2
    N_harm = 2*np.pi*w['l_env_nm']**2/A_cell
    N_bloch_tek = gp['N_eff_cells']/max(gp['n_peaks'], 1)
    out['N_kontrol'] = dict(N_harmonik=N_harm, N_bloch_tek_kuyu=N_bloch_tek,
                            oran=N_bloch_tek/N_harm)
    # --- aci taramasi: kuyu, zarf, bant genisligi, yorunge araligi ---
    scan = []
    for th in theta_scan:
        Q = dict(P); Q['theta_deg'] = th
        ws = well_analysis(Q)
        Ns = N_nuclei(Q, ws['l_env_nm'])
        _, _, bd = solve_moire(Q, k_path(moire_geometry(Q['a0'], th)[2], npts=16))
        s_pair_bw = bd[:, :2].max() - bd[:, :2].min()   # s-cifti tam genisligi
        orb_gap = bd[:, 2].min() - bd[:, :2].max()      # s-cifti -> p-kabugu araligi
        scan.append(dict(theta=th, aM=ws['aM'], depth=ws['depth_meV'],
                         barrier=ws['barrier_meV'], hw=ws['hbar_omega_meV'],
                         l=ws['l_env_nm'], N_Mo=Ns['N_Mo'],
                         s_bant_genisligi=s_pair_bw, s_p_araligi=orb_gap))
    out['scan'] = scan
    # --- hyperfine kiyas (uc senaryo: iyimser/merkez/kotumser A kestirimi) ---
    Ns = N_nuclei(P, w['l_env_nm'])
    f_Mo_nat = P['f_Mo95'] + P['f_Mo97']
    rows = []
    for senaryo in ('iyimser', 'merkez', 'kotumser'):
        A_Mo, A_S = gamma_A_effective(P, senaryo)
        r_nat = hyperfine_T2([A_Mo, A_S], [f_Mo_nat, P['f_S33']],
                             [2.5, 1.5], [Ns['N_Mo'], Ns['N_S']], P['g_factor'])
        r_pur = hyperfine_T2([A_Mo, A_S], [P['f_Mo_pure'], P['f_S_pure']],
                             [2.5, 1.5], [Ns['N_Mo'], Ns['N_S']], P['g_factor'])
        rows.append((f'dogal MoS2 [{senaryo}]', f_Mo_nat, Ns['N_Mo'], r_nat))
        rows.append((f'98Mo32S2 [{senaryo}]', P['f_Mo_pure'], Ns['N_Mo'], r_pur))
    out['hyperfine_rows'] = rows
    out['N'] = Ns
    # --- GaAs dogrulama satiri: ayni formul deneysel ~10 ns'yi tutturuyor mu? ---
    # A ~ 90 ueV (69/71Ga + 75As ortalama; Coish-Baugh 2009), N ~ 1e6 (lateral QD)
    out['gaas_check'] = hyperfine_T2([90.0], [1.0], [1.5], [1.0e6], 2.0)
    # --- gorev 5: kalan dekoherans + zenginlestirme duyarliligi ---
    out['kalan'] = kalan_dekoherans()
    out['zengin_tarama'] = zenginlestirme_taramasi(P, Ns['N_Mo'], Ns['N_S'])
    return out

if __name__ == '__main__':
    res = main(P)
    print(json.dumps({k: v for k, v in res.items() if k != 'hyperfine_rows'},
                     default=lambda o: float(o) if isinstance(o, np.floating) else str(o),
                     indent=2))
    for name, f, N, r in res['hyperfine_rows']:
        print(f"{name:24s} f={f:.4g} N_Mo={N:.3g} sigma={r['sigma_ueV']:.4g} ueV "
              f"({r['sigma_mT']:.4g} mT)  T2*={r['T2_ns']:.4g} ns")
    g = res['gaas_check']
    print(f"GaAs dogrulama: sigma={g['sigma_ueV']:.3g} ueV T2*={g['T2_ns']:.3g} ns "
          f"(deney ~10 ns, Petta 2005)")

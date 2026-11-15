# -*- coding: utf-8 -*-
"""
Kol B preprint — şekil üretici
==============================
Girdi : t2_hesap.py (same folder), results/dft_sonuc.json
Çıktı : figures/fig1..fig4 (PDF + PNG) + results/preprint_sayilar.json
        (main.tex'teki her sayının kaynağı; statik çapraz kontrol için)

Çalıştırma:  python fig_uret.py   (numpy, scipy, matplotlib yeter; QE/WSL gerekmez)
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

import t2_hesap as T
from t2_hesap import (P, moire_geometry, hole_potential_real, well_analysis,
                      hyperfine_T2, N_nuclei, gamma_A_effective,
                      zenginlestirme_taramasi, solve_moire, k_path,
                      ground_state_participation, hole_fourier_coeffs,
                      build_hamiltonian)
from scipy.linalg import eigh

plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "legend.fontsize": 7.2, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "figure.dpi": 200, "savefig.bbox": "tight",
    "axes.linewidth": 0.7,
})

C_V = "#b2182b"   # değerlik / kırmızı
C_C = "#2166ac"   # iletim / mavi
C_G = "#4d4d4d"   # gri
C_OK = "#1b7837"  # yeşil
C_BAND = "#fddbc7"

OUT = {}  # preprint_sayilar.json içeriği


def save(fig, name):
    fig.savefig(f"figures/{name}.pdf")
    fig.savefig(f"figures/{name}.png", dpi=300)
    plt.close(fig)
    print(f"  figures/{name}.pdf/png yazildi")


# ======================================================================
# Şekil 1 — Karanlık hüküm: momentum aritmetiği + kendi DFT seviyeleri
# ======================================================================
def fig1():
    dft = json.load(open("results/dft/dft_sonuc.json"))
    a0 = P["a0"]
    GK = 4 * np.pi / (3 * a0)                      # Γ–K mesafesi (1/nm)
    aM2, gM2, _ = moire_geometry(a0, 2.0)
    aM4, gM4, _ = moire_geometry(a0, 4.0)
    n_umk2 = GK / gM2
    OUT["momentum"] = dict(GK_invnm=GK, gM_2deg=gM2, gM_4deg=gM4,
                           aM_2deg_nm=aM2, n_umklapp_2deg=n_umk2)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.9),
                                 gridspec_kw=dict(width_ratios=[1.15, 1.0]))

    # --- (a) BZ ve moiré umklapp zinciri ---
    th = np.deg2rad(np.arange(7) * 60 + 30)
    R = GK
    hexx, hexy = R * np.cos(th), R * np.sin(th)
    ax.plot(hexx, hexy, color=C_G, lw=1.0)
    ax.plot(0, 0, "o", color=C_V, ms=6, zorder=5)
    ax.annotate(r"$\Gamma$ (VBM)", (0, 0), (-0.18 * R, 0.10 * R), color=C_V,
                fontsize=9, fontweight="bold")
    Kpt = np.array([R * np.cos(np.pi / 6), R * np.sin(np.pi / 6)])
    for a in th[:-1]:
        ax.plot(R * np.cos(a), R * np.sin(a), "s", color=C_C, ms=5, zorder=5)
    ax.annotate(r"K (CBM)", Kpt, Kpt + np.array([0.06 * R, 0.14 * R]),
                color=C_C, fontsize=9, fontweight="bold", ha="right")
    Qpt = 0.55 * Kpt
    ax.plot(*Qpt, "D", color=C_C, ms=4, mfc="none", zorder=5)
    ax.annotate(r"Q", Qpt, Qpt + np.array([0.02 * R, -0.13 * R]), color=C_C,
                fontsize=8)
    # umklapp zinciri: Γ→K yolu boyunca |g_M| adımları
    uhat = Kpt / np.linalg.norm(Kpt)
    nstep = int(np.floor(n_umk2))
    for i in range(nstep):
        p0, p1 = i * gM2 * uhat, (i + 1) * gM2 * uhat
        ax.annotate("", xy=p1, xytext=p0,
                    arrowprops=dict(arrowstyle="-|>", color=C_OK, lw=0.9,
                                    mutation_scale=7))
    ax.annotate(rf"$|g_M|={gM2:.2f}$ nm$^{{-1}}$ ($\theta=2^\circ$)"
                + "\n" + rf"$\Gamma\!\to\!$K $= {GK:.1f}$ nm$^{{-1}}$"
                + rf" $\approx {int(n_umk2)}\,|g_M|$",
                (0.46 * R, -0.33 * R), fontsize=8, color=C_OK, ha="center")
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(a) moiré umklapp arithmetic", loc="left")

    # --- (b) kendi DFT seviye hizası ---
    lv = dft["bands"]
    xpos = {"Gamma": 0.0, "Q(0.167)": 1.0, "K": 2.0}
    lbl = {"Gamma": r"$\Gamma$", "Q(0.167)": "Q", "K": "K"}
    for k, x in xpos.items():
        for band, col in (("VBM", C_V), ("CBM", C_C)):
            E = lv[k][band]
            bx.hlines(E, x - 0.28, x + 0.28, color=col, lw=2.4)
    # indirekt / direkt oklar
    Ev_G, Ec_K, Ev_K = lv["Gamma"]["VBM"], lv["K"]["CBM"], lv["K"]["VBM"]
    bx.annotate("", xy=(xpos["K"] - 0.1, Ec_K), xytext=(xpos["Gamma"] + 0.1, Ev_G),
                arrowprops=dict(arrowstyle="-|>", color=C_OK, lw=1.4))
    bx.text(0.95, 0.5 * (Ev_G + Ec_K) + 0.1, f"indirect\n{dft['gap_indirect_eV']:.2f} eV",
            fontsize=8, color=C_OK, ha="center")
    bx.annotate("", xy=(xpos["K"] + 0.22, Ec_K), xytext=(xpos["K"] + 0.22, Ev_K),
                arrowprops=dict(arrowstyle="<|-|>", color=C_G, lw=1.0))
    bx.text(xpos["K"] + 0.30, 0.5 * (Ev_K + Ec_K), f"direct\n{dft['gap_direct_eV']:.2f} eV",
            fontsize=7.5, color=C_G, va="center")
    bx.annotate(rf"$\Delta$VBM$(\Gamma-$K$)=+{dft['dVBM_eV']:.3f}$ eV",
                (0.03, 0.97), xycoords="axes fraction", fontsize=8, color=C_V,
                va="top")
    bx.annotate(r"$\Gamma$ VBM: 56.7% Mo $d_{z^2}$ + 39.6% S $p_z$",
                (0.03, 0.88), xycoords="axes fraction", fontsize=7.5, color=C_G,
                va="top")
    bx.set_xticks(list(xpos.values()), [lbl[k] for k in xpos])
    bx.set_xlim(-0.55, 2.75)
    bx.set_ylim(0.9, 5.6)
    bx.set_ylabel("energy (eV, cell reference)")
    bx.set_title("(b) bilayer 2H MoS$_2$, QE-PBE (this work)", loc="left")
    save(fig, "fig1_karanlik_hukum")


# ======================================================================
# Şekil 2 — Γ-vadisi moiré kuyusu + açı taraması (dot rejimi)
# ======================================================================
def fig2():
    aM, gM, gs = moire_geometry(P["a0"], P["theta_deg"])
    w = well_analysis(P)
    gp = ground_state_participation(P)
    OUT["well_2deg"] = {k: float(v) for k, v in w.items()}
    OUT["ground_state_2deg"] = {k: float(v) for k, v in gp.items()}

    L = 1.05 * aM
    n = 240
    xs = np.linspace(-L, L, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    U = hole_potential_real(P, gs, X, Y)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.9),
                                 gridspec_kw=dict(width_ratios=[1.0, 1.15]))
    im = ax.pcolormesh(X, Y, U, cmap="RdBu_r", shading="auto",
                       rasterized=True)
    ax.contour(X, Y, U, levels=8, colors="k", linewidths=0.25, alpha=0.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
    cb.set_label(r"$U_h(\mathbf{r})$ (meV)", labelpad=1)
    fig.subplots_adjust(wspace=0.42)
    # bal peteği minimumları işaretle
    thr = U.min() + 0.08 * (U.max() - U.min())
    from scipy import ndimage
    lab, nl = ndimage.label(U < thr)
    for i in range(1, nl + 1):
        cy = np.mean(X[lab == i]); cx = np.mean(Y[lab == i])
        ax.plot(cy, cx, "o", mfc="none", mec=C_OK, ms=9, mew=1.2)
    ax.annotate("", xy=(aM / 2, -0.88 * L), xytext=(-aM / 2, -0.88 * L),
                arrowprops=dict(arrowstyle="<|-|>", color="k", lw=1.0))
    ax.text(0, -0.80 * L, rf"$a_M={aM:.1f}$ nm", ha="center", fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
    ax.set_title(r"(a) hole potential, $\theta=2^\circ$", loc="left")

    # --- (b) açı taraması ---
    thetas = np.array([2.0, 2.5, 3.0, 3.5, 4.0])
    bars, hws, ts, Ns, orb = [], [], [], [], []
    scan_rows = []
    for th in thetas:
        Q = dict(P); Q["theta_deg"] = th
        ws = well_analysis(Q)
        Nn = N_nuclei(Q, ws["l_env_nm"])
        _, _, bd = solve_moire(Q, k_path(moire_geometry(Q["a0"], th)[2], npts=16))
        s_bw = bd[:, :2].max() - bd[:, :2].min()
        orb_gap = bd[:, 2].min() - bd[:, :2].max()
        bars.append(ws["barrier_meV"]); hws.append(ws["hbar_omega_meV"])
        ts.append(s_bw / 6.0); orb.append(orb_gap)
        Ns.append(Nn["N_Mo"])
        scan_rows.append(dict(theta=th, aM=ws["aM"], barrier=ws["barrier_meV"],
                              hw=ws["hbar_omega_meV"], l=ws["l_env_nm"],
                              t=s_bw / 6.0, sp_gap=orb_gap, N_Mo=Nn["N_Mo"]))
    OUT["aci_taramasi"] = scan_rows
    bx.semilogy(thetas, bars, "o-", color=C_G, label="barrier to neighbor well")
    bx.semilogy(thetas, orb, "s-", color=C_C, label=r"$s\!\to\!p$ orbital gap")
    bx.semilogy(thetas, ts, "d-", color=C_V, label=r"tunneling $t$")
    bx.axvspan(1.9, 2.1, color=C_OK, alpha=0.15)
    bx.set_ylim(0.9, 300)
    bx.text(2.07, 1.25, "dot\nregime", ha="left", fontsize=7.5, color=C_OK)
    bx.set_xlabel(r"twist angle $\theta$ (deg)")
    bx.set_ylabel("energy (meV)")
    bx.set_title("(b) confinement vs. tunneling", loc="left")
    bx.legend(loc="lower right", framealpha=0.9)
    bx.grid(alpha=0.25, which="both", lw=0.3)
    save(fig, "fig2_moire_kuyusu")


# ======================================================================
# Şekil 3 — Koherans bütçesi: T2* vs kalıntı spinli izotop oranı
# ======================================================================
def fig3():
    w = well_analysis(P)
    Ns = N_nuclei(P, w["l_env_nm"])
    f_nat = P["f_Mo95"] + P["f_Mo97"]

    fig, ax = plt.subplots(figsize=(3.5, 3.1))
    taramalar = {}
    for sen, stil in (("iyimser", ":"), ("merkez", "-"), ("kotumser", "--")):
        tr = zenginlestirme_taramasi(P, Ns["N_Mo"], Ns["N_S"], senaryo=sen)
        f = np.array([r["f_res"] for r in tr])
        t2 = np.array([r["T2_ns"] for r in tr]) * 1e-3  # µs
        taramalar[sen] = (f, t2)
    f_c, t_c = taramalar["merkez"]
    ax.fill_between(taramalar["iyimser"][0], taramalar["kotumser"][1],
                    taramalar["iyimser"][1], color=C_BAND, alpha=0.8,
                    label="$A$-estimate band")
    ax.loglog(f_c, t_c, "-", color=C_V, lw=1.6, label="central scenario")

    # Mo-only zenginleştirme (S doğal %0.76 ³³S): ~0.41 µs'de doyar
    A_Mo0, A_S0 = gamma_A_effective(P, "merkez")
    t2_moonly = np.array([
        hyperfine_T2([A_Mo0, A_S0], [fv, P["f_S33"]], [2.5, 1.5],
                     [Ns["N_Mo"], Ns["N_S"]], P["g_factor"])["T2_ns"] * 1e-3
        for fv in f_c])
    ax.loglog(f_c, t2_moonly, "-.", color="k", lw=1.0,
              label="Mo-only enrichment (natural S)")

    # hBN kapsul dipolar tabani (sizintisiz, surekli tabaka) ve toplam
    hz = T.hbn_dipolar_T2(P, w["l_env_nm"], "z")["T2_ns"] * 1e-3
    hx = T.hbn_dipolar_T2(P, w["l_env_nm"], "x")["T2_ns"] * 1e-3
    ax.axhspan(hz, hx, color="0.55", alpha=0.25, lw=0)
    ax.text(2.2e-7, 1.04 * hx, "hBN-encapsulation floor", fontsize=6.4,
            color="0.25")
    tot_z = np.array([T.combine_T2(t * 1e3, hz * 1e3) for t in t_c]) * 1e-3
    ax.loglog(f_c, tot_z, "-", color="0.15", lw=1.3,
              label=r"total with hBN ($B\perp$)")
    OUT["hbn_taban"] = dict(T2_hBN_perp_us=hz, T2_hBN_inplane_us=hx,
                            toplam_2em4_perp_us=T.combine_T2(1514.1, hz * 1e3) * 1e-3)

    # yük gürültüsü tabanı
    ax.axhspan(2, 20, color=C_C, alpha=0.12)
    ax.text(1.5e-3, 3.4, "charge-noise floor (Si-scaled)", fontsize=7,
            color=C_C)

    # bölgeler
    ax.axvspan(1e-4, 1e-3, color="0.85", alpha=0.6)
    ax.text(1.5e-4, 1.3e2, "high-grade\nenrichment", ha="center", fontsize=6.5,
            color="0.35")
    ax.axvspan(1e-6, 1e-5, color=C_OK, alpha=0.12)
    ax.text(3.2e-6, 1.4e2, r"$^{28}$Si-class", ha="center", fontsize=6.5,
            color=C_OK)

    # işaret noktaları (merkez senaryo)
    def at(f0):
        return np.interp(np.log10(f0), np.log10(f_c), t_c)
    pts = [(f_nat, "natural\n42 ns"), (2e-4, r"1.5 µs"), (1e-6, r"22 µs")]
    for f0, lab in pts:
        ax.plot(f0, at(f0), "o", color=C_V, ms=4)
    ax.annotate("natural: 42 ns", (f_nat, at(f_nat)), (2.0e-2, 5e-3),
                fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.annotate("1.5 µs", (2e-4, at(2e-4)), (1.2e-3, 0.6), fontsize=7,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.annotate("22 µs", (1e-6, at(1e-6)), (1.0e-5, 60), fontsize=7,
                arrowprops=dict(arrowstyle="-", lw=0.6))

    # deneysel referanslar (farklı malzemeler; bağlam için)
    refs = [(1.0, 0.010, "GaAs (expt)"), (0.047, 1.0, "natural Si"),
            (8e-4, 120.0, r"$^{28}$Si MOS"), (8e-4, 20.0, r"$^{28}$Si EDSR")]
    for f0, t0, lab in refs:
        ax.plot(f0, t0, "^", color=C_G, ms=4)
        ax.annotate(lab, (f0, t0), (f0 * 1.6, t0 * 1.25), fontsize=6.2,
                    color=C_G)

    OUT["butce"] = dict(f_dogal=f_nat,
                        T2_dogal_ns=at(f_nat) * 1e3,
                        T2_2em4_us=at(2e-4), T2_1em5_us=at(1e-5),
                        T2_1em6_us=at(1e-6),
                        kazanc_2em4=at(2e-4) * 1e3 / (at(f_nat) * 1e3),
                        mo_only_doygunluk_us=float(t2_moonly[0]))
    ax.set_xlabel(r"residual spinful-isotope fraction $f_{\rm res}$")
    ax.set_ylabel(r"$T_2^{*}$ ($\mu$s)")
    ax.set_xlim(1.2e-7, 1.6)
    ax.set_ylim(1e-4, 4e2)
    ax.legend(loc="lower left", framealpha=0.9)
    ax.grid(alpha=0.25, which="both", lw=0.3)
    save(fig, "fig3_koherans_butcesi")


# ======================================================================
# Şekil 4 — A-kestirimi duyarlılık analizi (tornado)
# ======================================================================
def fig4():
    w = well_analysis(P)
    Ns = N_nuclei(P, w["l_env_nm"])
    f_pur = (P["f_Mo_pure"], P["f_S_pure"])
    I = [2.5, 1.5]

    def T2_of(A_Mo, A_S, NMo=None, NS=None):
        NMo = Ns["N_Mo"] if NMo is None else NMo
        NS = Ns["N_S"] if NS is None else NS
        r = hyperfine_T2([A_Mo, A_S], list(f_pur), I, [NMo, NS], P["g_factor"])
        return r["T2_ns"] * 1e-3  # µs

    A_Mo0, A_S0 = gamma_A_effective(P, "merkez")
    base = T2_of(A_Mo0, A_S0)

    satirlar = []
    # 1) radyal integral aktarımı (Avdeev-Smirnov ölçek sabitlerinin Γ'ya taşınması)
    satirlar.append((r"radial transfer $\xi\in[0.5,2]$",
                     T2_of(2 * A_Mo0, 2 * A_S0), T2_of(0.5 * A_Mo0, 0.5 * A_S0)))
    # 2) üst sınır: K-değerlik Ising katsayıları (kötümser senaryo tanımı)
    A_Mo_k = (24.0 / 7.0) * P["A_scale_Mo"]
    A_S_k = (8.0 / 5.0) * P["A_scale_S"]
    satirlar.append((r"upper bound: K-band $A_{zz}$ coeffs.",
                     T2_of(A_Mo_k, A_S_k), base))
    # 3) orbital ağırlıklar: atanmamış %3.7'nin tamamı d'ye veya p'ye
    wd, wp = P["w_d"], P["w_p"]
    rest = 1.0 - wd - wp
    hi = T2_of((wd + rest) * (4 / 7) * P["A_scale_Mo"], wp * (4 / 5) * P["A_scale_S"])
    lo = T2_of(wd * (4 / 7) * P["A_scale_Mo"], (wp + rest) * (4 / 5) * P["A_scale_S"])
    satirlar.append((r"DFT weights: residual 3.7% share", min(hi, lo), max(hi, lo)))
    # 4) zarf çekirdek sayısı: harmonik 106 vs Bloch-katılım 152 (tek kuyu)
    satirlar.append((r"envelope $N$: 106 (harm.) vs 152 (Bloch)",
                     T2_of(A_Mo0, A_S0, 152.0, 2 * 152.0), base))
    # 4b) zarf: pc-AFM ölçülen ~2.0 nm eksiton konfinmanı (Westenberg 2026)
    N2 = 2 * np.pi * 2.0 ** 2 / (np.sqrt(3) / 2 * P["a0"] ** 2)
    satirlar.append((r"envelope $\ell$: 1.22 vs 2.0 nm (pc-AFM)",
                     T2_of(A_Mo0, A_S0, N2, 2 * N2), base))
    # 5) T2* konvansiyonu ħ/σ vs √2ħ/σ
    satirlar.append((r"convention: $\hbar/\sigma$ vs $\sqrt{2}\hbar/\sigma$",
                     base * np.sqrt(2), base))

    OUT["duyarlilik"] = dict(
        taban_us=base,
        satirlar=[dict(ad=s[0], yuksek=max(s[1], s[2]), dusuk=min(s[1], s[2]))
                  for s in satirlar])

    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    ypos = np.arange(len(satirlar))[::-1]
    for y, (ad, a, b) in zip(ypos, satirlar):
        lo, hi = min(a, b), max(a, b)
        ax.barh(y, np.log10(hi) - np.log10(lo), left=np.log10(lo),
                color=C_BAND, edgecolor=C_V, height=0.55, lw=0.8)
        ax.text(np.log10(hi) + 0.06, y, f"{lo:.2g}–{hi:.2g}", va="center",
                fontsize=6.4, color=C_G)
    ax.axvline(np.log10(base), color=C_V, lw=1.3)
    ax.text(np.log10(base), len(satirlar) - 0.25,
            rf"central {base:.1f} $\mu$s", fontsize=7, color=C_V, ha="center")
    ax.set_yticks(ypos, [s[0] for s in satirlar], fontsize=7)
    xt = [-1, 0, 1]
    ax.set_xticks(xt, [f"$10^{{{v}}}$" for v in xt])
    ax.set_xlim(-1.2, 1.3)
    ax.set_xlabel(r"$T_2^{*}$ ($\mu$s), $^{98}$Mo$^{32}$S$_2$ at $f_{\rm res}=2\times10^{-4}$")
    ax.grid(axis="x", alpha=0.25, lw=0.3)
    save(fig, "fig4_duyarlilik")


# ======================================================================
if __name__ == "__main__":
    import os
    os.makedirs("figures", exist_ok=True)

    # çapraz kontrol: merkez senaryo ana satırları yeniden üret
    w = well_analysis(P)
    Ns = N_nuclei(P, w["l_env_nm"])
    f_nat = P["f_Mo95"] + P["f_Mo97"]
    A_Mo, A_S = gamma_A_effective(P, "merkez")
    nat = hyperfine_T2([A_Mo, A_S], [f_nat, P["f_S33"]], [2.5, 1.5],
                       [Ns["N_Mo"], Ns["N_S"]], P["g_factor"])
    pur = hyperfine_T2([A_Mo, A_S], [P["f_Mo_pure"], P["f_S_pure"]],
                       [2.5, 1.5], [Ns["N_Mo"], Ns["N_S"]], P["g_factor"])
    gaas = hyperfine_T2([90.0], [1.0], [1.5], [1.0e6], 2.0)
    OUT["ana_satirlar"] = dict(
        N_Mo=Ns["N_Mo"], N_S=Ns["N_S"],
        sigma_dogal_ueV=nat["sigma_ueV"], T2_dogal_ns=nat["T2_ns"],
        sigma_pur_ueV=pur["sigma_ueV"], T2_pur_ns=pur["T2_ns"],
        A_Mo_merkez_ueV=A_Mo, A_S_merkez_ueV=A_S,
        gaas_T2_ns=gaas["T2_ns"])
    print("Ana satırlar:", json.dumps(OUT["ana_satirlar"], indent=2))

    # kapalı formlar (tekrar-üretilebilirlik; makale Sec. III ile birebir)
    V1, V2, V3 = P["V_shells"]
    aM, gM, gs = moire_geometry(P["a0"], P["theta_deg"])
    OUT["kapali_formlar"] = dict(
        U_AA=6 * (V1 + V2 + V3), U_B=-3 * V1 + 6 * V2 - 3 * V3,
        U_SP=-2 * V1 - 2 * V2 + 6 * V3,
        bariyer_V1_8V2_9V3=V1 - 8 * V2 + 9 * V3,
        kappa_meV_nm2=(1.5 * V1 - 9 * V2 + 6 * V3) * gM ** 2,
        kappa_kabuk_paylari=[1.5 * V1, -9 * V2, 6 * V3])
    # Γ_M tam spektrum (s-çifti + p-kabuğu)
    Vc = hole_fourier_coeffs(P, gs)
    H, _ = build_hamiltonian(np.zeros(2), gs, P["m_eff"], Vc, nshell=6)
    ww = np.sort(eigh(H, eigvals_only=True))[:6]
    OUT["gammaM_spektrum_meV"] = [float(x) for x in (ww - ww[0])]
    print("Gamma_M spektrum (meV):", np.round(ww - ww[0], 2))
    print("kapali formlar:", OUT["kapali_formlar"])

    fig1(); fig2(); fig3(); fig4()

    with open("results/preprint_sayilar.json", "w") as f:
        json.dump(OUT, f, indent=2, default=float)
    print("results/preprint_sayilar.json yazildi")

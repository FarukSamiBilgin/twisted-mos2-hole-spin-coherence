# Bant Karakteri Kararı — Bükülmüş Bilayer MoS₂ (2–4°)

## HÜKÜM: KARANLIK

Küçük açılı bükülmüş bilayer MoS₂'nin en düşük enerjili arabant geçişi **momentum-indirekt,
dipol-yasaklı ve fonon-destekli**dir; malzeme fotolüminesansta karanlıktır. Bu sınıflandırma
bir "belirsizlik" değil, birbirinden bağımsız dört kanıt hattının kesişimidir. Görev 4–5,
görev tanımının öngördüğü gibi **elektriksel okuma varsayımıyla** sürdürülmüştür.

---

## Kanıt hattı 1 — Değerlik maksimumu Γ'dadır (moiré düz bantları Γ-vadisinden türer)

- Naik & Jain, PRL 121, 266401 (2018), DFT (doğrudan alıntı): "The VBM in all BLM stackings
  is at the Γ point." Düz değerlik bantları Γ-vadisi kökenli; K'daki 150 meV SOC yarılması
  Γ-türevli düz bantları etkilemiyor (Γ–K enerji ayrımı büyük).
- Angeli & MacDonald, PNAS 118, e2021826118 (2021): Γ-vadisi moiré bantları için continuum
  modeli; delikler AB/BA (Bernal) bölgelerinde lokalize olup **bal peteği örgüsü** kurar.
  Model parametreleri (Tablo I, birincil kaynaktan doğrulandı): V₁=39.45, V₂=6.5,
  V₃=10.0 meV; tüm fazlar ψ=π; m*=0.9 m₀; a₀=3.182 Å.
- Kendi continuum hesabımız (bu depo, `t2_hesap.py`): aynı parametrelerle 3°'de en alt iki
  moiré bandı dejenere s-çifti (bal peteği), bant genişliği ~11 meV — Naik-Jain'in 3.5°
  DFT değeri (23 meV) ile mertebe uyumlu. Model, iki bağımsız DFT çalışmasıyla tutarlı.

## Kanıt hattı 2 — İletim minimumu Γ'da DEĞİL: geçiş momentum-indirekt

- Bilayer 2H MoS₂ indirekt aralıklıdır: Mak et al., PRL 105, 136805 (2010): bulk indirekt
  1.29 eV; bilayer indirekt ~1.6 eV, direkt K–K ~1.88 eV. Cheiwchanchamnangij & Lambrecht,
  PRB 85, 205302 (2012), GW: bilayerde VBM Γ'da, CBM K'da (Γ→K indirekt; deney 1.6 eV).
- CBM'nin K mı Q/Λ mı olduğu literatürde bilayer için milimetrik yarış içindedir
  (Mak 2010 "Γ–K doğrultusu boyunca minimum" der; Lambrecht K der). **Karar için fark
  etmez:** iki aday da Γ'dan ~13 nm⁻¹ uzaktadır; geçiş her durumda indirekt.
- Momentum aritmetiği: Γ–K mesafesi 4π/(3a₀) ≈ 13.2 nm⁻¹. Moiré ters örgü vektörü
  2°'de |g_M| ≈ 0.81 nm⁻¹ (Westenberg 2025'in ölçtüğü λ_M ≈ 9 nm ile bağımsız doğrulandı),
  4°'de ≈1.6 nm⁻¹. K(c)'yi Γ(v) üzerine katlamak ~16 adet moiré umklapp ister; arabant
  tünelleme matris elemanları yalnızca tek monolayer ters örgü vektörü bağlar
  (Wang et al., PRB 95, 115429). **Moiré katlama geçişi direkte çeviremez.**

## Kanıt hattı 3 — Deney: her bükülme açısında karanlık

- Mak et al. (2010): PL kuantum verimi monolayer ~4×10⁻³, bilayer 10⁻⁵–10⁻⁶
  (bilayer, monolayerden ~400–4000× sönük).
- Liu et al., Nat. Commun. 5, 4966 (2014), 44 örnek, θ=0–60° taraması: **her açıda**
  indirekt pik (1.4–1.7 eV) mevcut; direkt K–K piki her açıda monolayerin ~1/50'sine
  bastırılmış. Direkt taban durumunu geri getiren hiçbir açı yok.
- Westenberg et al., arXiv:2511.20398 (2025), 2° bükülmüş MoS₂ homobilayer (tam bizim
  rejim): indirekt eksitonlar (Γ–K arakatman ~1.46 eV; karışık Γ–K/Γ–Λ 1.6–1.7 eV)
  **fotolüminesansta karanlık**, yalnızca fotoakım-AFM ile tespit edilebiliyor
  (doğrudan alıntı: "Unlike photoluminescence, where indirect excitons are typically
  dark..."). Moiré hibridizasyonu parlaklığı en alt duruma taşımıyor; en alt durumlar
  indirekt/karanlık kalıyor.

## Kanıt hattı 3.5 — Kendi DFT çaprazımız (bu çalışma, `kol_B/dft/`)

QE 7.5 (PBE, PAW, 12×12×1) ile bilayer 2H MoS₂ (deneysel geometri, a=3.182 Å,
katmanlar arası Mo-Mo 6.15 Å):

| Nicelik | Bu çalışma (PBE) | Literatür |
|---|---|---|
| VBM(Γ) − VBM(K) | **+0.477 eV** (Γ kazanır) | Naik-Jain, AM ile uyumlu |
| CBM konumu | **K** (Q'dan 84 meV aşağıda) | Lambrecht GW ile uyumlu |
| İndirekt Γ→K aralık | 1.169 eV | deney 1.6 eV (PBE alçaltması tipik) |
| Direkt K→K aralık | 1.647 eV | deney 1.88 eV |
| Γ VBM bileşimi | **Mo d_z² %56.7 + S p_z %39.6** | "d_z²+p_z" nitel literatürle uyumlu; yüzdeler bu çalışmanın katkısı |

Not: SOC dahil değil (skaler-rel.); SOC K-VBM'yi ~74 meV yükseltir, +0.477 eV'lik
Γ üstünlüğünü değiştirmez. Girdi/çıktılar ve ayrıştırıcı `kol_B/dft/` altındadır.

## Kanıt hattı 4 — Katman karakteri ek bastırma getirir

Γ değerlik durumu Mo d_z² + S p_z **katmanlar-arası hibrit** (bonding/antibonding;
katmanlara yayılmış); K iletim durumu Mo d-orbital ağırlıklı, **katman-lokalize**
(Splendiani 2010). Zarf örtüşmesi de zayıf: geçiş dipolü, momentum yasağına ek olarak
katman uzayında da bastırılır (nicel katsayı literatürde yok — nitel destek).

## Geçiş dipolü sıfırdan ayırt edilebilir mi?

Sıfırıncı mertebede hayır (momentum-yasaklı). Sonlu katkılar yalnızca fonon/kusur-destekli
yüksek mertebe süreçlerden gelir; deneysel görünürlük ancak PL-dışı kanallarla (fotoakım)
mümkün olmuştur (Westenberg 2025). Optik arayüz kurmaya yetecek osilatör şiddeti **yoktur**.

---

## Sonuç ve yönlendirme

| Soru | Cevap |
|---|---|
| Değerlik maksimumu | Γ (moiré düz bantları Γ-türevli, bal peteği) |
| Geçiş direkt mi | Hayır — Γ(v)→K/Q(c), momentum-indirekt, her açıda |
| Dipol | Sıfırdan ayırt edilemez (fonon-destekli artıklar hariç) |
| Etiket | **KARANLIK** |

**Mimari sonucu:** Bu kol optik qubit adayı değildir. Görev 4–6, "elektriksel qubit +
hyperfine-temiz uzun-ömürlü platform" çerçevesinde tamamlanmıştır (bkz. `rol_hukmu.md`).
Spin okuması için yol: kapı-tanımlı yük algılama / Pauli spin blokajı (TMD QD literatürü);
optik başlatma/okuma bu malzemede kapalıdır.

**Not — SOC ve spin yapısı (görev 2'ye girdi):** Γ noktasında değerlik SOC yarılması
**tam sıfır** (Kramers/zaman-tersinimi; Kormányos 2015, açıkça belirtilmiş). Γ deliğinin
vadi serbestlik derecesi yok, saf spin-1/2 Kramers çifti; tek-parçacık beklentisi g≈2
(d_z²/p_z, L_z=0). Karşılaştırma: K'da değerlik SOC 148 meV, iletim 3 meV, Q'da ~70 meV
(monolayer değerleri; Kormányos et al., 2D Mater. 2, 022001, 2015).

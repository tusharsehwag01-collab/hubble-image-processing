#!/usr/bin/env python3
"""
CAT'S EYE NEBULA (NGC 6543) — NASA Publication Grade Pipeline
=============================================================
Data: WFPC2/PC, Proposal 05403, HLA single-exposure DRZ files
Filters: F437N (H-gamma blue), F502N (OIII green), F658N (NII red)

Order of operations (NASA-style, lessons from Orion):
  1. Load + EXPTIME normalize
  2. Crop to 500x500 around bright center
  3. Cosmic ray removal (astroscrappy) per exposure
  4. Hot pixel cleanup per exposure
  5. Median stack per filter
  6. Background subtraction (photutils Background2D)
  7. Build RGB: R=F658N + F673N*0 (just NII), G=F502N, B=F437N*0.3 (suppressed)
  8. Pre-stretch percentile scaling
  9. make_lupton_rgb with Q/stretched params tuned for Cat's Eye
 10. HDR star core blend
 11. Background color noise reduction
 12. Final output
"""
import os, warnings, time, subprocess
import numpy as np
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
from astropy.stats import sigma_clipped_stats, SigmaClip
from scipy.ndimage import gaussian_filter, median_filter
from photutils.background import Background2D, MedianBackground
import astroscrappy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

DATA = "/Users/tusharsehwag/python projects/catseye_data"

FILTER_FILES = {
    'F437N': [
        'hst_05403_01_wfpc2_f437n_pc_01_drz.fits',
        'hst_05403_01_wfpc2_f437n_pc_02_drz.fits',
    ],
    'F502N': [
        'hst_05403_01_wfpc2_f502n_pc_01_drz.fits',
        'hst_05403_01_wfpc2_f502n_pc_02_drz.fits',
        'hst_05403_01_wfpc2_f502n_pc_03_drz.fits',
        'hst_05403_01_wfpc2_f502n_pc_04_drz.fits',
    ],
    'F658N': [
        'hst_05403_01_wfpc2_f658n_pc_01_drz.fits',
        'hst_05403_01_wfpc2_f658n_pc_02_drz.fits',
    ],
}

t_start = time.time()

print("=" * 70)
print("  CAT'S EYE NEBULA — NASA-GRADE PIPELINE")
print("  WFPC2/PC Proposal 05403 | HLA single-exposure DRZ")
print("=" * 70)

# ═══ STEP 1: LOAD + EXPTIME ═══
print(f"\n[1/12] Loading FITS + EXPTIME normalization...")
all_data = {}
for filt, flist in FILTER_FILES.items():
    exps = []
    for fname in flist:
        path = os.path.join(DATA, fname)
        hdul = fits.open(path)
        sci_ext = 0
        for i, hdu in enumerate(hdul):
            if hdu.data is not None and len(hdu.data.shape) >= 2:
                sci_ext = i; break
        d = hdul[sci_ext].data.astype(np.float64)
        hdr = hdul[sci_ext].header
        hdul.close()
        d[~np.isfinite(d)] = 0
        # NOTE: do NOT divide by EXPTIME — HLA DRZ values are already low (counts/sec)
        # Lupton stretch expects larger inputs (hundreds, not 0.0x)
        exps.append(d)
    all_data[filt] = exps
    print(f"  ✓ {filt}: {len(exps)} exps × {exps[0].shape}")

# ═══ STEP 2: LOCATE + CROP ═══
print(f"\n[2/12] Locating nebula and cropping...")
ref = all_data['F502N'][0]
bright = ref > np.percentile(ref, 99.5)
rows, cols = np.where(bright)
cy, cx = int(np.median(rows)), int(np.median(cols))
CROP = 500
y1, y2 = max(0, cy-CROP//2), min(ref.shape[0], cy+CROP//2)
x1, x2 = max(0, cx-CROP//2), min(ref.shape[1], cx+CROP//2)
for filt in all_data:
    all_data[filt] = [d[y1:y2, x1:x2] for d in all_data[filt]]
print(f"  Center ({cy},{cx}) → crop {y2-y1}×{x2-x1}")

# ═══ STEP 3: COSMIC RAYS ═══
print(f"\n[3/12] Cosmic ray removal (astroscrappy)...")
n_cr = 0
for filt in all_data:
    cl = []
    for d in all_data[filt]:
        mask, c = astroscrappy.detect_cosmics(
            d, sigclip=4.5, sigfrac=0.3, objlim=5.0,
            satlevel=np.inf, niter=3)
        cl.append(c); n_cr += int(np.sum(mask))
    all_data[filt] = cl
print(f"  ✓ Removed: {n_cr}")

# ═══ STEP 4: HOT PIXELS ═══
print(f"\n[4/12] Hot pixel removal (sigma=6, conservative)...")
n_hot = 0
for filt in all_data:
    cl = []
    for d in all_data[filt]:
        _, _, std = sigma_clipped_stats(d, sigma=3.0, maxiters=5)
        mi = median_filter(d, size=3)
        hm = (d - mi) > (6.0 * std)
        c = d.copy()
        c[hm] = mi[hm]
        cl.append(c); n_hot += int(np.sum(hm))
    all_data[filt] = cl
print(f"  ✓ Fixed: {n_hot}")

# ═══ STEP 5: MEDIAN STACK ═══
print(f"\n[5/12] Median stacking per filter...")
stacked = {}
for filt in all_data:
    stacked[filt] = np.median(all_data[filt], axis=0)
    d = stacked[filt]
    print(f"  ✓ {filt}: range [{d.min():.3f}, {d.max():.3f}]")

# ═══ STEP 6: BACKGROUND SUBTRACTION per filter ═══
print(f"\n[6/12] Background subtraction (Background2D per filter)...")
sc = SigmaClip(sigma=3.0, maxiters=10)
for filt in stacked:
    bkg = Background2D(
        stacked[filt], (50, 50),
        filter_size=(3, 3),
        sigma_clip=sc,
        bkg_estimator=MedianBackground(),
        exclude_percentile=10
    )
    sub = stacked[filt] - bkg.background
    sub[sub < 0] = 0
    stacked[filt] = sub
    print(f"  ✓ {filt}: median_bkg={bkg.background_median:.4f}")

# ═══ STEP 7: BUILD RGB ═══
# Palette 1 (Traditional/OIII→Green) — matches official Hubble:
#   R = F658N (NII) — red outer shells
#   G = F502N (OIII) — dominant green glow
#   B = F437N (H-gamma) — suppressed blue (very faint, noisy)
print(f"\n[7/12] Building RGB channels...")
print(f"  Palette: R=F658N  G=F502N  B=F437N×0.3")

h, w = stacked['F502N'].shape
red = stacked['F658N']
green = stacked['F502N']
blue = stacked['F437N'] * 0.3  # suppress noisy blue

# ═══ STEP 8: PER-CHANNEL SCALING ═══
# Equalize using 95th percentile of positive pixels (not median — that's background!)
print(f"\n[8/12] Per-channel scaling (p95 equalization for nebula, not background)...")
refs_95 = []
for ch in [red, green, blue]:
    pos = ch[ch > 0]
    if len(pos) > 0:
        refs_95.append(np.percentile(pos, 95))
    else:
        refs_95.append(1.0)
target = float(np.median(refs_95))
for i, r in enumerate(refs_95):
    if r > 0:
        if i == 0: red = red * (target / r)
        elif i == 1: green = green * (target / r)
        else: blue = blue * (target / r)
print(f"  p95 of positives: R={refs_95[0]:.4f}, G={refs_95[1]:.4f}, B={refs_95[2]:.4f}")
print(f"  All scaled to ~{target:.4f}")

# ═══ STEP 9: Lupton RGB ═══
# Cat's Eye is a small target, dominated by OIII — need stronger contrast
# NASA images show rich green with red outer shell and subtle blue hints
print(f"\n[9/12] make_lupton_rgb...")
MAIN_Q = 12
MAIN_STRETCH = 5.0

# Clip to avoid inf/negatives
red = np.clip(np.nan_to_num(red), 0, None)
green = np.clip(np.nan_to_num(green), 0, None)
blue = np.clip(np.nan_to_num(blue), 0, None)

main_rgb = make_lupton_rgb(
    red, green, blue,
    minimum=0.0,
    Q=MAIN_Q,
    stretch=MAIN_STRETCH,
    filename=None
).astype(np.float32) / 255.0
print(f"  Main: Q={MAIN_Q}, stretch={MAIN_STRETCH}")

# ═══ STEP 10: HDR CORE BLEND ═══
print(f"\n[10/12] HDR core blend...")
CORE_Q = 25
CORE_STRETCH = 1.5

core_rgb = make_lupton_rgb(
    red, green, blue,
    minimum=0.0,
    Q=CORE_Q,
    stretch=CORE_STRETCH,
    filename=None
).astype(np.float32) / 255.0

# Luminance mask — only very bright core
lum = 0.299 * main_rgb[:,:,0] + 0.587 * main_rgb[:,:,1] + 0.114 * main_rgb[:,:,2]
# Only top 2% brightest pixels get core blend
thresh = np.percentile(lum, 98)
mask = np.clip((lum - thresh) / max(1e-6, lum.max() - thresh), 0, 1)
mask = gaussian_filter(mask, sigma=6)
mask = np.clip(mask ** 1.7, 0, 1)

blended = main_rgb * (1 - mask[:,:,None]) + core_rgb * mask[:,:,None]
print(f"  ✓ Mask min={mask.min():.3f} max={mask.max():.3f}")

# ═══ STEP 11: BACKGROUND DENOISE ═══
print(f"\n[11/12] Background color noise reduction...")
dim = np.mean(blended, axis=2) < np.percentile(np.mean(blended, axis=2), 30)
for i in range(3):
    med = median_filter(blended[:,:,i], size=3)
    blended[:,:,i][dim] = 0.65 * blended[:,:,i][dim] + 0.35 * med[dim]
print(f"  ✓ SCNR on 30% dimmest pixels")

# ═══ STEP 12: SAVE ═══
print(f"\n[12/12] Saving...")
out_main = os.path.join(DATA, "catseye_final.png")
out_master = os.path.join(DATA, "catseye_final_master.png")

# Save 8-bit PNG
blended = np.clip(blended, 0, 1)
plt.imsave(out_main, blended, origin='lower', dpi=200)
plt.imsave(out_master, blended, origin='lower', vmin=0, vmax=1)

# Make 3-way comparison (different stretches)
fig, axes = plt.subplots(1, 3, figsize=(24, 8))
for ax, q, st, title in [(axes[0], 8, 0.10, "Gentle (Q=8, stretch=0.10)"),
                         (axes[1], MAIN_Q, MAIN_STRETCH, f"Standard (Q={MAIN_Q}, stretch={MAIN_STRETCH})"),
                         (axes[2], blended, None if True else None, "Final HDR blend")]:
    if isinstance(q, np.ndarray):
        img = q
    else:
        img = make_lupton_rgb(red, green, blue, minimum=0.0, Q=q, stretch=st, filename=None).astype(np.float32) / 255.0
    ax.imshow(np.clip(img, 0, 1), origin='lower')
    ax.set_title(title, color='white')
    ax.axis('off')

fig.patch.set_facecolor('black')
for ax in axes:
    ax.set_facecolor('black')
plt.tight_layout()
comp_path = os.path.join(DATA, "catseye_comparison_3ways.png")
plt.savefig(comp_path, facecolor='black', dpi=150, bbox_inches='tight')
plt.close()

t_elapsed = time.time() - t_start

# Quality diagnostics
print("\n" + "=" * 70)
print(f"  PIPELINE COMPLETE — {t_elapsed:.1f}s")
print(f"=" * 70)
print(f"  Target:    Cat's Eye Nebula (NGC 6543)")
print(f"  Telescope: HST / WFPC2 (PC chip)")
print(f"  Proposal:  05403")
print(f"  Filter map: R=NII(F658N) | G=OIII(F502N) | B=Hγ(F437N)×0.3")
print(f"  Palette:   Traditional (OIII→Green, matches official Hubble)")
print(f"  Lupton:    Q={MAIN_Q}, stretch={MAIN_STRETCH}")
print(f"  CR removed: {n_cr}")
print(f"  Hot px:     {n_hot}")
print(f"  Final size: {blended.shape}")
print(f"=" * 70)
print(f"  Channel balance:")
print(f"    R mean: {blended[:,:,0].mean():.3f}")
print(f"    G mean: {blended[:,:,1].mean():.3f}")
print(f"    B mean: {blended[:,:,2].mean():.3f}")
print(f"    Saturated (value=1): {int(np.sum(blended >= 1.0))}")
print(f"=" * 70)
print(f"  ✓ {out_main}")
print(f"  ✓ {out_master}")
print(f"  ✓ {comp_path}")
print(f"=" * 70)

subprocess.run(['open', out_main])

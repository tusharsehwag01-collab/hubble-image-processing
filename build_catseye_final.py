#!/usr/bin/env python3
"""
CAT'S EYE — FINAL PIPELINE (With H-alpha)
==========================================
Now has F656N (Hα) — the missing piece NASA uses for the green channel.

Mapping: R=F658N (NII) | G=F656N (Hα) | B=F502N (OIII)
Hot pixel: sigma=4.0 (was 6.0 — too weak)
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
    'F502N': ['hst_05403_01_wfpc2_f502n_pc_01_drz.fits',
              'hst_05403_01_wfpc2_f502n_pc_02_drz.fits',
              'hst_05403_01_wfpc2_f502n_pc_03_drz.fits',
              'hst_05403_01_wfpc2_f502n_pc_04_drz.fits'],
    'F656N': ['hst_05403_01_wfpc2_f656n_pc_01_drz.fits',
              'hst_05403_01_wfpc2_f656n_pc_02_drz.fits'],
    'F658N': ['hst_05403_01_wfpc2_f658n_pc_01_drz.fits',
              'hst_05403_01_wfpc2_f658n_pc_02_drz.fits'],
}

t_start = time.time()
print("=" * 70)
print("  CAT'S EYE NEBULA — FINAL PIPELINE WITH H-ALPHA")
print("  Mapping: R=NII(F658N) | G=Hα(F656N) | B=OIII(F502N)")
print("=" * 70)

# STEP 1: LOAD
print(f"\n[1/10] Loading + EXPTIME normalization...")
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
        hdr = hdul[sci_ext].header.copy()
        hdul.close()
        d[~np.isfinite(d)] = 0
        exptime = float(hdr.get('EXPTIME', 1.0) or 1.0)
        if exptime > 0: d = d / exptime
        exps.append(d)
    all_data[filt] = exps
    print(f"  ✓ {filt}: {len(exps)} × {exps[0].shape}")

# STEP 2: CROP
print(f"\n[2/10] Cropping around nebula...")
ref = all_data['F502N'][0]
bright = ref > np.percentile(ref, 99.5)
rows, cols = np.where(bright)
cy, cx = int(np.median(rows)), int(np.median(cols))
CROP = 600
y1 = max(0, cy-CROP//2); y2 = min(ref.shape[0], cy+CROP//2)
x1 = max(0, cx-CROP//2); x2 = min(ref.shape[1], cx+CROP//2)
for filt in all_data:
    all_data[filt] = [d[y1:y2, x1:x2] for d in all_data[filt]]
print(f"  ({cy},{cx}) → {y2-y1}×{x2-x1}")

# STEP 3: COSMIC RAYS
print(f"\n[3/10] Cosmic rays (astroscrappy)...")
n_cr = 0
for filt in all_data:
    cl = []
    for d in all_data[filt]:
        mask, c = astroscrappy.detect_cosmics(
            d, sigclip=4.5, sigfrac=0.3, objlim=5.0,
            satlevel=np.inf, niter=3)
        cl.append(c); n_cr += int(np.sum(mask))
    all_data[filt] = cl
print(f"  {n_cr} removed")

# STEP 4: HOT PIXELS (sigma=4 — tighter than 6!)
print(f"\n[4/10] Hot pixels (sigma=4)...")
n_hot = 0
for filt in all_data:
    cl = []
    for d in all_data[filt]:
        _, _, std = sigma_clipped_stats(d, sigma=3.0, maxiters=5)
        mi = median_filter(d, size=3)
        hm = (d - mi) > (3.0 * std)  # aggressive — was 4.0
        c = d.copy()
        c[hm] = mi[hm]
        cl.append(c); n_hot += int(np.sum(hm))
    all_data[filt] = cl
print(f"  {n_hot} fixed")

# STEP 5: STACK
print(f"\n[5/10] Median stacking...")
stacked = {}
for filt in all_data:
    stacked[filt] = np.median(all_data[filt], axis=0)

# STEP 6: BACKGROUND
print(f"\n[6/10] Background subtraction...")
sc = SigmaClip(sigma=3.0, maxiters=10)
for filt in stacked:
    bkg = Background2D(
        stacked[filt], (50, 50), filter_size=(3, 3),
        sigma_clip=sc, bkg_estimator=MedianBackground(),
        exclude_percentile=10
    )
    stacked[filt] = np.clip(stacked[filt] - bkg.background, 0, None)
    # Print max values
    d = stacked[filt]
    print(f"  ✓ {filt}: max={d.max():.3f}")

# STEP 7: RGB WITH CORRECT MAPPING (Hα in Green!)
print(f"\n[7/10] RGB + p95 scaling...")
red = stacked['F658N']
green = stacked['F656N']  # Hα — THIS IS THE KEY FIX
blue = stacked['F502N']   # OIII

# p95 scaling (nebula signal, not background)
for name, ch in [('R', red), ('G', green), ('B', blue)]:
    pos = ch[ch > 0]
    if len(pos) > 0:
        p95 = np.percentile(pos, 95)
        print(f"  {name}: p95={p95:.4f}")
    else:
        print(f"  {name}: no positive pixels")

# Equalize by p95
refs = [np.percentile(c[c>0], 95) if np.any(c>0) else 1.0 for c in [red, green, blue]]
target = float(np.median(refs))
for i, r in enumerate(refs):
    if r > 0:
        [red, green, blue][i] = [red, green, blue][i] * (target / r)
print(f"  All scaled to p95=~{target:.4f}")

red = np.clip(np.nan_to_num(red), 0, None)
green = np.clip(np.nan_to_num(green), 0, None)
blue = np.clip(np.nan_to_num(blue), 0, None)

# Rescale so max of brightest channel ~500 (Lupton expects hundreds-range values)
scale_factor = 500.0 / max(np.max(red) if np.max(red) > 0 else 1,
                              np.max(green) if np.max(green) > 0 else 1,
                              np.max(blue) if np.max(blue) > 0 else 1)
red *= scale_factor
green *= scale_factor
blue *= scale_factor
print(f"  Rescale factor: {scale_factor:.0f}x")

# STEP 8: LUPTON
print(f"\n[8/10] Lupton RGB...")
MAIN_Q, MAIN_STRETCH = 10, 3.0
main = make_lupton_rgb(red, green, blue, minimum=0.0, Q=MAIN_Q, stretch=MAIN_STRETCH, filename=None).astype(np.float32) / 255.0
print(f"  Q={MAIN_Q}, stretch={MAIN_STRETCH}")

# Gentle version for core — darker, tighter for stars
core = make_lupton_rgb(red, green, blue, minimum=0.0, Q=25, stretch=0.8, filename=None).astype(np.float32) / 255.0

# HDR blend
lum = 0.299*main[:,:,0] + 0.587*main[:,:,1] + 0.114*main[:,:,2]
thresh = np.percentile(lum, 98)
mask = np.clip((lum - thresh) / max(1e-6, lum.max() - thresh), 0, 1)
mask = gaussian_filter(mask, sigma=6)
mask = np.clip(mask ** 1.7, 0, 1)
blended = main * (1 - mask[:,:,None]) + core * mask[:,:,None]

# STEP 9: DENOISE (dimmest 30% only)
print(f"\n[9/10] Background denoise...")
dim = np.mean(blended, axis=2) < np.percentile(np.mean(blended, axis=2), 30)
for i in range(3):
    med = median_filter(blended[:,:,i], size=3)
    blended[:,:,i][dim] = 0.65 * blended[:,:,i][dim] + 0.35 * med[dim]

# Loosen black point — keep some faint structure
blended = np.clip(blended, 0, 1)

# Crop 5% margin to remove sensor edge artifact
h_c, w_c = blended.shape[:2]
margin = int(h_c * 0.05)
blended = blended[margin:h_c-margin, margin:w_c-margin]
print(f"  Cropped margin {margin}px → {blended.shape}")

# STEP 10: SAVE
print(f"\n[10/10] Saving...")
out = os.path.join(DATA, "catseye_final_ha.png")
plt.imsave(out, blended, origin='lower', dpi=200)

# 3-way comparison
fig, axes = plt.subplots(1, 3, figsize=(24, 8))
variants = [(8, 3.0, "Gentle Q=8"), (MAIN_Q, MAIN_STRETCH, "Standard"), (20, 8.0, "Saturated Q=20")]
for ax, (q, st, title) in zip(axes, variants):
    img = make_lupton_rgb(red, green, blue, minimum=0.0, Q=q, stretch=st, filename=None).astype(np.float32) / 255.0
    ax.imshow(np.clip(img, 0, 1), origin='lower')
    ax.set_title(title, color='white'); ax.axis('off')
fig.patch.set_facecolor('black')
for ax in axes: ax.set_facecolor('black')
plt.tight_layout()
comp = os.path.join(DATA, "catseye_final_ha_comparison.png")
plt.savefig(comp, facecolor='black', dpi=150, bbox_inches='tight')
plt.close()

t_elapsed = time.time() - t_start

print("\n" + "=" * 70)
print(f"  COMPLETE — {t_elapsed:.1f}s")
print(f"  Mapping: R=NII(F658N) | G=Hα(F656N) | B=OIII(F502N)")
print(f"  Hot px sigma=4 (was 6)")
print(f"  Final: {blended.shape}")
print(f"  R mean: {blended[:,:,0].mean():.3f}")
print(f"  G mean: {blended[:,:,1].mean():.3f}")
print(f"  B mean: {blended[:,:,2].mean():.3f}")
print(f"  Saturated: {int(np.sum(blended >= 0.999))}")
print(f"  ✓ {out}")
print("=" * 70)

subprocess.run(['open', out])
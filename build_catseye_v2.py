#!/usr/bin/env python3
"""
CAT'S EYE v2 — Individual channel tuning, no equalization
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
print("  CAT'S EYE v2 — Individual channel tuning")
print("=" * 70)

# Load (no EXPTIME divide — keep raw DRZ values)
print(f"\n[1/12] Loading FITS...")
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
        # EXPTIME normalization — puts all filters on consistent flux scale
        exptime = float(hdr.get('EXPTIME', 1.0) or 1.0)
        if exptime > 0:
            d = d / exptime
        exps.append(d)
    all_data[filt] = exps
    print(f"  ✓ {filt}: {len(exps)} exps × {exps[0].shape}")

# Crop
print(f"\n[2/12] Locating + cropping...")
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

# Cosmic rays
print(f"\n[3/12] Cosmic ray removal...")
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

# Hot pixels
print(f"\n[4/12] Hot pixel removal...")
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

# Stack
print(f"\n[5/12] Median stacking...")
stacked = {}
for filt in all_data:
    stacked[filt] = np.median(all_data[filt], axis=0)
print(f"  ✓ All filters stacked")

# Background
print(f"\n[6/12] Background subtraction...")
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
print(f"  ✓ All backgrounds subtracted")

# Hubble Palette for Cat's Eye:
# R = F658N [NII] — outer red shell
# G = F437N [Hγ] — very faint, must boost heavily
# B = F502N [OIII] — creates the famous blue-white inner glow
print(f"\n[7/12] Building RGB (NO equalization, CORRECTED MAPPING)...")
print(f"  Palette: R=F658N×1.0 | G=F437N×5.0 | B=F502N×1.0")

h, w = stacked['F502N'].shape
red = stacked['F658N'] * 1.0
green = stacked['F437N'] * 5.0  # Hγ is 10-20× fainter than OIII, must boost
blue = stacked['F502N'] * 1.0

# Clip
red = np.clip(np.nan_to_num(red), 0, None)
green = np.clip(np.nan_to_num(green), 0, None)
blue = np.clip(np.nan_to_num(blue), 0, None)

# Rescale back to "hundreds" range for Lupton stretch (EXPTIME-norm puts data in counts/sec, very small)
SCALE_FACTOR = 1000.0
red = red * SCALE_FACTOR
green = green * SCALE_FACTOR
blue = blue * SCALE_FACTOR
print(f"  Rescaled all channels ×{SCALE_FACTOR} for Lupton stretch")

print(f"\n[8/12] Lupton RGB (Q=15, stretch=8.0)...")
MAIN_Q = 15
MAIN_STRETCH = 8.0

main_rgb = make_lupton_rgb(
    red, green, blue,
    minimum=0.0,
    Q=MAIN_Q,
    stretch=MAIN_STRETCH,
    filename=None
).astype(np.float32) / 255.0

print(f"\n[9/12] HDR core blend...")
CORE_Q = 30
CORE_STRETCH = 2.0

core_rgb = make_lupton_rgb(
    red, green, blue,
    minimum=0.0,
    Q=CORE_Q,
    stretch=CORE_STRETCH,
    filename=None
).astype(np.float32) / 255.0

lum = 0.299 * main_rgb[:,:,0] + 0.587 * main_rgb[:,:,1] + 0.114 * main_rgb[:,:,2]
thresh = np.percentile(lum, 98)
mask = np.clip((lum - thresh) / max(1e-6, lum.max() - thresh), 0, 1)
mask = gaussian_filter(mask, sigma=6)
mask = np.clip(mask ** 1.7, 0, 1)

blended = main_rgb * (1 - mask[:,:,None]) + core_rgb * mask[:,:,None]

print(f"\n[10/12] Background denoise...")
dim = np.mean(blended, axis=2) < np.percentile(np.mean(blended, axis=2), 30)
for i in range(3):
    med = median_filter(blended[:,:,i], size=3)
    blended[:,:,i][dim] = 0.65 * blended[:,:,i][dim] + 0.35 * med[dim]

print(f"\n[11/12] Final crop...")
bright_final = np.mean(blended, axis=2)
thresh = np.percentile(bright_final, 55)
y_idx, x_idx = np.where(bright_final > thresh)
if len(y_idx) > 0:
    y1_f = max(0, y_idx.min()-20)
    y2_f = min(h, y_idx.max()+20)
    x1_f = max(0, x_idx.min()-20)
    x2_f = min(w, x_idx.max()+20)
    blended = blended[y1_f:y2_f, x1_f:x2_f]

blended = np.clip(blended, 0, 1)

print(f"\n[12/12] Saving...")
out_main = os.path.join(DATA, "catseye_v2.png")
out_master = os.path.join(DATA, "catseye_v2_master.png")
plt.imsave(out_main, blended, origin='lower', dpi=200)
plt.imsave(out_master, blended, origin='lower', vmin=0, vmax=1)

# 3-way comparison
fig, axes = plt.subplots(1, 3, figsize=(24, 8))
for ax, q, st, title in [(axes[0], 10, 5.0, "Gentle (Q=10, stretch=5)"),
                         (axes[1], MAIN_Q, MAIN_STRETCH, f"Standard (Q={MAIN_Q}, stretch={MAIN_STRETCH})"),
                         (axes[2], 20, 12.0, "Saturated (Q=20, stretch=12)")]:
    img = make_lupton_rgb(red, green, blue, minimum=0.0, Q=q, stretch=st, filename=None).astype(np.float32) / 255.0
    ax.imshow(np.clip(img, 0, 1), origin='lower')
    ax.set_title(title, color='white')
    ax.axis('off')

fig.patch.set_facecolor('black')
for ax in axes:
    ax.set_facecolor('black')
plt.tight_layout()
comp_path = os.path.join(DATA, "catseye_v2_comparison.png")
plt.savefig(comp_path, facecolor='black', dpi=150, bbox_inches='tight')
plt.close()

t_elapsed = time.time() - t_start

print("\n" + "=" * 70)
print(f"  PIPELINE COMPLETE — {t_elapsed:.1f}s")
print(f"=" * 70)
print(f"  Channel scaling: R×1.5, G×1.0 (reference), B×0.15")
print(f"  Lupton: Q={MAIN_Q}, stretch={MAIN_STRETCH}")
print(f"  Final: {blended.shape}")
print(f"=" * 70)
print(f"  R mean: {blended[:,:,0].mean():.3f}")
print(f"  G mean: {blended[:,:,1].mean():.3f}")
print(f"  B mean: {blended[:,:,2].mean():.3f}")
print(f"  Saturated: {int(np.sum(blended >= 1.0))}")
print(f"=" * 70)
print(f"  ✓ {out_main}")
print(f"  ✓ {comp_path}")
print(f"=" * 70)

subprocess.run(['open', out_main])

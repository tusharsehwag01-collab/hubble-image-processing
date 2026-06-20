# Hubble Image Processing

Processing real Hubble Space Telescope data using Python, astropy, and professional-grade pipelines.

## Targets

- **Cat's Eye Nebula (NGC 6543)** — Planetary nebula, 3 narrowband filters (F502N, F656N, F658N)
- **M51 Whirlpool Galaxy** — Face-on spiral, 3 broadband filters (F439W, F555W, F814W)
- **Orion Nebula (M42)** — Star-forming region (attempted)
- **Jupiter** — WFC3/UVIS narrowband (attempted)

## Pipeline Steps

1. MAST archive query via astroquery
2. FITS download and EXPTIME normalization
3. Cosmic ray rejection (astroscrappy)
4. Background subtraction (photutils Background2D)
5. WCS reprojection for alignment
6. Median stacking per filter
7. AsinhStretch per channel
8. RGB composite via make_lupton_rgb
9. HDR core blend for bright targets

## Tools

Python · astropy · astroquery · numpy · scipy · matplotlib · photutils · astroscrappy · reproject · FITS Liberator

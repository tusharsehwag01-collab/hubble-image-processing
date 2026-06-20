# Cat's Eye Nebula (NGC 6543) — Hubble Pipeline

Multi-filter composite of the Cat's Eye planetary nebula from real Hubble Space Telescope data.

## Filters Used

| Filter | Emission Line | Color | 
|--------|--------------|-------|
| F502N | [O III] (502nm) | Blue channel |
| F656N | Hα (656nm) | Green channel |
| F658N | [N II] (658nm) | Red channel |

## Pipeline Steps

1. Download from MAST archive via astroquery
2. EXPTIME normalization
3. Cosmic ray rejection (astroscrappy)
4. Background subtraction (photutils)
5. Median stacking
6. AsinhStretch per channel
7. RGB composite with make_lupton_rgb

## Data Source

Proposal 05403, WFPC2/PC instrument, Hubble Space Telescope.

# M51 · Whirlpool Galaxy research notes

These notes record the three-filter broadband composite explored in the Hubble-learning project. **This file is a research note, not an executable M51 pipeline.** No M51 processing script, FITS dataset or rendered composite is included in the current public repository.

## Recorded filter mapping

- **F439W** → blue: shorter-wavelength structure and young star clusters.
- **F555W** → green: the middle visible channel.
- **F814W** → red: the longer-wavelength / near-infrared channel.

These are display-channel choices, not a claim of exact human-eye colour or quantitative calibration.

## Recorded source context

The earlier project notes identify **Hubble proposal 07375, WFPC2/PC**, with the three filters from the same pointing. Those notes describe them as pre-aligned. Inspect the actual FITS headers and WCS before relying on that assumption in a new run; matching pointing alone does not verify pixel alignment.

No new archive download or scientific validation was performed in this documentation update.

## Continue from here

- [Repository overview and limitations](README.md)
- [Published Cat's Eye H-alpha experiment](build_catseye_final.py) — a different target, not a substitute M51 implementation.
- [Tushar's project profile](https://github.com/tusharsehwag01-collab)

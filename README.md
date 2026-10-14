# Hubble image processing

Python experiments with real Hubble Space Telescope FITS data. This repository preserves Cat's Eye Nebula processing variants and M51 source/filter notes as part of an astronomy-learning portfolio.

**Status:** experimental scripts and research notes. This is not a NASA-endorsed pipeline, a calibrated science product or a self-contained download-and-run package.

## What is actually in this repository

- [`build_catseye.py`](build_catseye.py) — an initial Cat's Eye variant using F437N, F502N and F658N, with stacking, cleanup, background subtraction and colour/contrast experiments.
- [`build_catseye_v2.py`](build_catseye_v2.py) — another colour-mapping/channel-tuning experiment using those filters.
- [`build_catseye_final.py`](build_catseye_final.py) — an H-alpha variant using F502N, F656N and F658N, mapped to blue, green and red respectively.
- [`M51-README.md`](M51-README.md) — Whirlpool Galaxy source/filter notes. There is no M51 processing script in this published repository.

The filenames record historical variants; **“final” does not certify scientific correctness**. FITS input files, rendered images, a dataset downloader, a portable configuration layer and an automated test suite are not included.

## Techniques explored

The included scripts use FITS loading, bright-region cropping, cosmic-ray rejection, median-filter hot-pixel cleanup, median stacking, background estimation/subtraction, Lupton RGB conversion, contrast experiments and bright-core blending.

Previous work also explored M51, Orion and Jupiter. Only the code and notes listed above are published here. A target being explored is not proof of a finished result, and WCS reprojection/archive-query features are not implemented by these included scripts.

## Setup for inspecting the scripts

A local Python environment and the following libraries are needed:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy astropy scipy photutils astroscrappy matplotlib
```

This is an unpinned dependency list derived from the actual imports, not a tested lockfile. The original scripts were written for a local macOS setup and attempt to open a rendered image with the macOS `open` command after processing.

Before running a variant:

1. Obtain the matching Hubble FITS inputs from an authoritative archive and review their product metadata.
2. Set that script's `DATA` value to your own local data directory. The historical scripts still contain an author-machine path.
3. Check every entry in `FILTER_FILES`; the chosen variant must match the available observations.
4. Verify FITS extensions, units, exposure-time treatment and alignment before processing. The scripts differ in their normalisation assumptions; a header/comment is not scientific validation.
5. Review clipping, cleanup thresholds, colour mapping and the effect on faint structure. Do not treat an attractive composite as a quantitative measurement.

After preparing and reviewing those inputs, run only the chosen variant, for example:

```sh
python build_catseye_final.py
```

The scripts write outputs into `DATA`. Keep raw data untouched and use a working copy. No data-processing run or dependency installation was performed as part of this documentation refresh.

## Scientific and provenance limits

- Phrases such as “NASA-grade” in historical source comments are the author's earlier wording, not NASA affiliation, endorsement or certification.
- Colour choices are visual mappings, not necessarily what a human eye would see.
- Exposure/unit handling and image alignment need dataset-specific review. These scripts should not be used as a general calibrated pipeline without that work.
- Links to target notes identify the recorded research context; they are not an independent verification of every dataset claim.
- No project-wide code licence has been added here. Refer to upstream data/library terms before reuse or publication.

## Related projects

[Project profile](https://github.com/tusharsehwag01-collab) · [Website design studies](https://github.com/tusharsehwag01-collab/local-biz-ai) · [The Living Mystics](https://www.youtube.com/@TheLivingMystics)

# ElectrophyAnalysis

Automated, reproducible analysis of whole-cell current-clamp recordings.

Three command-line pipelines and one point-and-click application that turn Axon
Binary Format (ABF) files into figure-ready Excel workbooks, with no cursor
placement and no manual parameter hunting.

| Pipeline | Protocol | What it returns |
|---|---|---|
| `intrinsic_firing.py` | 500 ms repetitive-firing steps | 63-column intrinsic phenotype: rheobase, F–I curve, biphasic adaptation fit, fast/slow AHP, ADP, passive properties (Rin, τm, Cm, sag) |
| `single_ap.py` | brief suprathreshold pulse | Per-spike kinetics via the Allen Institute IPFX detector, plus a full-width-at-half-rise metric (`Durée_PA`) |
| `slow_depol.py` | sustained recordings (≥ 3 min) | Detection and area-under-curve quantification of slow depolarizations against a drift-tracking local baseline |

All three are protocol-aware, batch-capable, deterministic, and fully logged.
Every numerical parameter is exposed and overridable; nothing is hard-coded.

---

## Two ways to run it

### Point and click — no Python required

Download `ElectrophyAnalysis-<version>-win64.zip` from the
[Releases](../../releases) page, unzip it anywhere, and run
`ElectrophyAnalysis.exe`. Everything is bundled; no installation, no Python, no
administrator rights.

The application wraps the three pipelines without duplicating any analysis
logic — it builds the command line for you and streams the pipeline's output
into a log panel. Results are identical to running the scripts directly.

### Command line

```bash
python intrinsic_firing.py --input  ./data/current_steps \
                           --output intrinsic \
                           --config examples/config_intrinsic.yaml

python single_ap.py        --input  ./data/single_ap \
                           --output SpikeFeatures.xlsx \
                           --mapping examples/group_mapping.csv

python slow_depol.py       --input  ./data/long_recordings \
                           --output ./results \
                           --bl-window 7 --threshold 10
```

Every option can also be supplied from a YAML file via `--config`, so a run is
fully described by a single file. Command-line flags take precedence over the
config file, which takes precedence over the built-in defaults.

Run any pipeline with `--help` for the complete option list.

---

## Installation

> **The install order matters.** `ipfx 2.0.0` breaks with `numpy >= 1.24`
> (type-promotion changes), and `efel 5.7.17` declares a `pynwb` requirement it
> does not actually use at runtime. Installing `requirements.txt` in one pass
> with a naive resolver produces a broken environment. Follow the sequence
> below, or use `build_exe.ps1`, which automates it.

Python 3.10.16 is the reference interpreter. 3.11 works — `numpy 1.23.5` and
`scipy 1.10.1` both ship cp310 and cp311 wheels — but 3.10.16 is what the
published results were produced with.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 1. numpy first, and pinned, so nothing can silently upgrade it
pip install numpy==1.23.5

# 2. everything else under the numpy constraint; ipfx before efel,
#    because ipfx pins pynwb==2.2.0 and efel would otherwise bump it
pip install -r requirements.txt -c constraints.txt --prefer-binary

# 3. efel's runtime deps and efel itself, without dependency resolution
pip install neo quantities --no-deps
pip install efel==5.7.17 --no-deps

# 4. verify
pip check
python -c "import numpy,scipy,pandas,matplotlib,pyabf,efel,ipfx; print('OK', numpy.__version__)"
```

A conda alternative is provided in `environment.yml`.

---

## Reproducibility

Every run writes a `Metadata` sheet into its output workbook and a timestamped
log into `runinfo/`, recording:

- the version of every scientific library actually loaded, not the pinned one
- the absolute input and output paths
- every parameter used, including the ones left at their default
- the SHA-256 of the group-mapping CSV and of the config file (`single_ap.py`)

Two runs that report the same checksums and the same parameter block were
driven by byte-identical inputs. A published result can therefore be tied to
the exact configuration that produced it, and re-executed from the stored
parameters alone.

---

## Validation

`validation/auc_validation_v3.py` checks the AUC estimator against analytically
known solutions:

- **24 baseline-stable cases** — six signal archetypes (three rectangles, a
  triangle, a trapezoid, a half-sine) at 5, 10, 15 and 20 kHz. Maximum error
  3 × 10⁻⁶ %.
- **5 drift cases** — linear baseline drift of 0, ±5, −10 mV over a 120 s
  recording. Maximum error 2.33 %, i.e. within tolerance at the most extreme
  drift rate simulated (≈ 0.08 mV·s⁻¹).

Run it with `python validation/auc_validation_v3.py`; it regenerates
`AUC_Validation_Synthetique_V3.xlsx` and three figures.

Validation against manual ClampFit analysis on real recordings is reported in
the accompanying manuscript.

---

## Running the tests

```bash
pip install pytest
pytest -v
```

The suite covers the AUC estimator against closed-form solutions, the
config/CLI precedence rules, the checksum helper, and a smoke test per pipeline
on the example dataset.

---

## Repository layout

```
intrinsic_firing.py      single_ap.py      slow_depol.py    the three pipelines
electrophy_gui.py                                           the GUI front-end
electrophy_gui.spec      build_exe.ps1     run_gui.bat      Windows build chain
examples/                                                   configs, mapping CSV, sample ABFs
validation/                                                 synthetic AUC validation
tests/                                                      pytest suite
requirements.txt         constraints.txt   environment.yml  pinned environment
```

---

## Citing

If this software contributed to your work, please cite both the software
release and the accompanying paper. See [`CITATION.cff`](CITATION.cff), or use
GitHub's "Cite this repository" button.

---

## Licence

MIT — see [`LICENSE`](LICENSE).

Built on [pyABF](https://github.com/swharden/pyABF),
[eFEL](https://github.com/BlueBrain/eFEL),
[IPFX](https://github.com/AllenInstitute/ipfx), NumPy, SciPy, pandas and
matplotlib. Please cite those projects as well.

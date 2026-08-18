"""The AUC estimator against closed-form solutions.

These are the assertions behind the synthetic-validation section of the paper.
`validation/auc_validation_v3.py` produces the figures and the workbook; this
file is the machine-checkable version of the same claims, so a regression fails
CI instead of quietly changing a published number.

The functions under test are imported from `slow_depol.py` itself — the actual
shipped pipeline, not a copy — so the test cannot drift away from the code it
is meant to protect.
"""
import numpy as np
import pytest

# Six archetypes x four sampling rates = the 24 baseline-stable cases.
SAMPLE_RATES = (5_000, 10_000, 15_000, 20_000)

# (label, shape, amplitude_mV, duration_s, analytic_area_mV_s)
ARCHETYPES = (
    ("rectangle_A", "rectangle", 10.0, 2.0, 10.0 * 2.0),
    ("rectangle_B", "rectangle", 25.0, 1.0, 25.0 * 1.0),
    ("rectangle_C", "rectangle", 5.0, 0.5, 5.0 * 0.5),
    ("triangle", "triangle", 20.0, 2.0, 0.5 * 20.0 * 2.0),
    ("trapezoid", "trapezoid", 15.0, 1.5, 15.0 * 1.5 * 0.75),
    ("half_sine", "half_sine", 10.0, 1.0, 10.0 * 1.0 * (2 / np.pi)),
)

PAD_S = 10.0        # > BL_WINDOW_S (7 s), so the estimator sees pure baseline
BL_WINDOW_S = 7.0
BL_STEP_S = 0.5
BL_PERCENTILE = 10

STABLE_TOL_PCT = 0.5    # baseline-stable cases
DRIFT_TOL_PCT = 3.0     # matches DRIFT_TOL_PCT in the validation script


def build_event(shape, amp, duration_s, fs):
    """A single event of known area, padded with silence on both sides."""
    n = int(duration_s * fs)
    x = np.linspace(0.0, 1.0, n, endpoint=False)
    if shape == "rectangle":
        body = np.full(n, amp)
    elif shape == "triangle":
        body = amp * (1.0 - np.abs(2.0 * x - 1.0))
    elif shape == "trapezoid":
        body = np.clip(amp * np.minimum(4.0 * x, 4.0 * (1.0 - x)), 0.0, amp)
    elif shape == "half_sine":
        body = amp * np.sin(np.pi * x)
    else:
        raise ValueError(shape)

    pad = int(PAD_S * fs)
    v = np.concatenate([np.zeros(pad), body, np.zeros(pad)])
    t = np.arange(v.size) / fs
    return t, v, pad - 1, pad + n


def measure(mod, t, v, i0, i1):
    baseline = mod.compute_local_baseline(v, t, BL_WINDOW_S, BL_STEP_S, BL_PERCENTILE)
    return mod.compute_global_area(v, t, baseline, i0, i1)


@pytest.mark.parametrize("label,shape,amp,dur,analytic", ARCHETYPES)
@pytest.mark.parametrize("fs", SAMPLE_RATES)
def test_area_matches_analytic_solution(slow_depol, label, shape, amp, dur, analytic, fs):
    """24 cases: every archetype, every sampling rate, stable baseline."""
    t, v, i0, i1 = build_event(shape, amp, dur, fs)
    measured = measure(slow_depol, t, v, i0, i1)
    err_pct = 100.0 * abs(measured - analytic) / analytic
    assert err_pct < STABLE_TOL_PCT, (
        f"{label} @ {fs // 1000} kHz: {measured:.6f} vs analytic {analytic:.6f} "
        f"({err_pct:.4f} % error)"
    )


@pytest.mark.parametrize("drift_mv", (0.0, -5.0, +5.0, -10.0))
def test_area_is_recovered_under_baseline_drift(slow_depol, drift_mv):
    """A drifting resting potential must not leak into the measured area."""
    fs, total_s, amp, dur = 10_000, 120.0, 10.0, 2.0
    analytic = amp * dur

    total = int(total_s * fs)
    baseline_true = np.linspace(0.0, drift_mv, total)
    v = baseline_true.copy()

    start = total // 2
    n = int(dur * fs)
    v[start:start + n] += amp
    t = np.arange(total) / fs

    measured = measure(slow_depol, t, v, start - 1, start + n)
    err_pct = 100.0 * abs(measured - analytic) / analytic
    assert err_pct < DRIFT_TOL_PCT, (
        f"drift {drift_mv:+.0f} mV over {total_s:.0f} s: {measured:.4f} vs "
        f"analytic {analytic:.4f} ({err_pct:.3f} % error)"
    )


def test_hyperpolarizations_do_not_subtract_from_the_area(slow_depol):
    """The AUC is one-sided: dips below baseline are clipped, not subtracted."""
    fs, amp, dur = 10_000, 10.0, 2.0
    t, v, i0, i1 = build_event("rectangle", amp, dur, fs)

    without_dip = measure(slow_depol, t, v, i0, i1)
    v_dip = v.copy()
    v_dip[i0 - int(0.2 * fs):i0] = -15.0        # hyperpolarizing artefact just before
    with_dip = measure(slow_depol, t, v_dip, i0, i1)

    assert with_dip == pytest.approx(without_dip, rel=1e-9)


def test_area_scales_linearly_with_amplitude(slow_depol):
    fs, dur = 10_000, 1.0
    areas = []
    for amp in (5.0, 10.0, 20.0):
        t, v, i0, i1 = build_event("rectangle", amp, dur, fs)
        areas.append(measure(slow_depol, t, v, i0, i1))
    assert areas[1] == pytest.approx(2 * areas[0], rel=1e-6)
    assert areas[2] == pytest.approx(4 * areas[0], rel=1e-6)

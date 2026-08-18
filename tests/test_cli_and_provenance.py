"""Command-line / config precedence, and the provenance metadata.

The manuscript claims a uniform `--input` / `--output` / `--config` interface
across the three pipelines, and a Metadata sheet that fingerprints the inputs.
These tests hold the code to both claims.
"""
import textwrap

import pytest


# --------------------------------------------------------------------------
# SHA-256 fingerprinting of the mapping CSV and the config file
# --------------------------------------------------------------------------

def test_checksum_is_stable_and_prefixed(single_ap, tmp_path):
    csv = tmp_path / "map.csv"
    csv.write_text("date,file_number_start,file_number_end,group\n"
                   "2024_01_15,0,19,Control\n", encoding="utf-8")

    digest = single_ap.file_sha256(str(csv))
    assert digest.startswith("sha256:")
    assert len(digest) == len("sha256:") + 64
    assert digest == single_ap.file_sha256(str(csv)), "must be deterministic"


def test_checksum_changes_when_the_mapping_changes(single_ap, tmp_path):
    csv = tmp_path / "map.csv"
    csv.write_text("date,file_number_start,file_number_end,group\n"
                   "2024_01_15,0,19,Control\n", encoding="utf-8")
    before = single_ap.file_sha256(str(csv))

    with csv.open("a", encoding="utf-8") as fh:
        fh.write("2024_01_15,20,39,Treatment\n")

    assert single_ap.file_sha256(str(csv)) != before, (
        "a result must not be able to claim the same mapping after an edit"
    )


def test_checksum_of_a_missing_file_is_none_not_a_crash(single_ap, tmp_path):
    """Provenance is best-effort: a missing file must not abort an analysis."""
    assert single_ap.file_sha256(str(tmp_path / "nope.csv")) is None


# --------------------------------------------------------------------------
# YAML config loading
# --------------------------------------------------------------------------

def test_config_is_read_as_a_mapping(single_ap, tmp_path):
    cfg = tmp_path / "params.yaml"
    cfg.write_text(textwrap.dedent("""
        input: /data/abf
        output: results.xlsx
        mapping: groups.csv
        overwrite: true
    """), encoding="utf-8")

    loaded = single_ap.load_config(str(cfg))
    assert loaded["input"] == "/data/abf"
    assert loaded["overwrite"] is True


def test_empty_config_is_an_empty_dict_not_none(single_ap, tmp_path):
    cfg = tmp_path / "empty.yaml"
    cfg.write_text("", encoding="utf-8")
    assert single_ap.load_config(str(cfg)) == {}


def test_missing_config_exits_with_a_clear_message(single_ap, tmp_path):
    with pytest.raises(SystemExit):
        single_ap.load_config(str(tmp_path / "absent.yaml"))


def test_non_mapping_config_is_rejected(single_ap, tmp_path):
    cfg = tmp_path / "list.yaml"
    cfg.write_text("- not\n- a\n- mapping\n", encoding="utf-8")
    with pytest.raises(SystemExit):
        single_ap.load_config(str(cfg))


# --------------------------------------------------------------------------
# slow_depol.py: CLI over YAML over defaults
# --------------------------------------------------------------------------

def test_defaults_apply_when_nothing_is_given(slow_depol, tmp_path):
    args = slow_depol.parse_args(["-i", str(tmp_path), "-o", str(tmp_path)])
    cfg = slow_depol.build_config(args)
    assert cfg["BL_WINDOW_S"] == slow_depol.DEFAULTS["BL_WINDOW_S"]
    assert cfg["DEPOL_THRESHOLD_MV"] == slow_depol.DEFAULTS["DEPOL_THRESHOLD_MV"]


def test_command_line_overrides_the_default(slow_depol, tmp_path):
    args = slow_depol.parse_args(
        ["-i", str(tmp_path), "-o", str(tmp_path), "--bl-window", "15", "--threshold", "8"])
    cfg = slow_depol.build_config(args)
    assert cfg["BL_WINDOW_S"] == 15
    assert cfg["DEPOL_THRESHOLD_MV"] == 8


def test_command_line_overrides_the_config_file(slow_depol, tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("BL_WINDOW_S: 20\nDEPOL_THRESHOLD_MV: 5\n", encoding="utf-8")

    args = slow_depol.parse_args(
        ["-i", str(tmp_path), "-o", str(tmp_path), "-c", str(cfg_file), "--bl-window", "9"])
    cfg = slow_depol.build_config(args)

    assert cfg["BL_WINDOW_S"] == 9, "explicit CLI flag must win over the config file"
    assert cfg["DEPOL_THRESHOLD_MV"] == 5, "config file must still supply the rest"


# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------

def test_run_leaves_a_log_on_disk(slow_depol, tmp_path):
    """The paper claims all three pipelines are fully logged."""
    log_path = slow_depol.setup_logging(str(tmp_path))
    slow_depol.log.info("Depolarizations detected: 28")
    slow_depol.log.warning("No baseline point found")

    for handler in slow_depol.log.handlers:
        handler.flush()

    contents = open(log_path, encoding="utf-8").read()
    assert "Depolarizations detected: 28" in contents
    assert "WARNING" in contents

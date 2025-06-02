# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Spike Feature Extraction Script
------------------------------

Author: Elias [@YourLab]
Affiliation: [Your Institution]
Version: 1.2.0
Date: 2024-06-12

Description:
Extracts spike features from ABF files (using IPFX), groups by experimental condition, saves to Excel,
logs all errors, checks file/data integrity, and generates QC plots (PDF/PNG) for every group.

Usage (CLI):
    python SingleSpikechat.py --data_dir "<path_to_abf_files>" --output "<output_excel_file.xlsx>"

Dependencies:
    - numpy==1.23.5
    - pandas==1.5.3
    - matplotlib==3.10.3
    - pyabf==2.3.8
    - ipfx==2.0.0
    - openpyxl==3.1.5

See README.md for details and test instructions.
"""
import getpass
import os
import sys
import shutil
import numpy as np
import pandas as pd
import pyabf
import argparse
import logging
import traceback
from datetime import datetime
import matplotlib.pyplot as plt

from ipfx.spike_detector import detect_putative_spikes
from ipfx.feature_extractor import SpikeFeatureExtractor
import yaml
# ===================== CONSTANTS AND CONFIGURATION ===========================

EXPECTED_COLS = [
    "Width (ms)", "Amplitude (mV)", "Threshold (mV)",
    "Duration (ms)", "Durée_PA", "Epoch_Level", "File_Info"
]

def setup_logging(log_file="analysis.log", error_log_file="error.log"):
    """Set up logging for info and error separately."""
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w', encoding='utf-8'),
            logging.StreamHandler()
        ])
    error_logger = logging.getLogger("errorLogger")
    error_handler = logging.FileHandler(error_log_file, mode='w', encoding='utf-8')
    error_handler.setLevel(logging.ERROR)
    error_logger.addHandler(error_handler)
    return error_logger
def export_run_params(params_dict, mapping_path, output_dir):
    """
    Sauvegarde les paramètres du run (dictionnaire params_dict)
    et le mapping utilisé dans un fichier run_params.txt ET run_params.yaml dans output_dir.
    """
    run_params_txt = os.path.join(output_dir, "run_params.txt")
    run_params_yaml = os.path.join(output_dir, "run_params.yaml")
    with open(run_params_txt, "w", encoding="utf-8") as f:
        f.write("### Spike Feature Extraction - Run Parameters ###\n\n")
        for key, val in params_dict.items():
            f.write(f"{key} : {val}\n")
        f.write("\n## Mapping file used :\n")
        with open(mapping_path, "r", encoding="utf-8") as m:
            f.write(m.read())
    # YAML export (plus lisible et exploitable par un autre script)
    params_dict['mapping_file_used'] = mapping_path
    try:
        with open(run_params_yaml, "w", encoding="utf-8") as f:
            yaml.dump(params_dict, f, allow_unicode=True, sort_keys=False)
    except Exception as e:
        print(f"Erreur export YAML: {e}")
    print(f"Run parameters exported to {run_params_txt} and {run_params_yaml}")

def parse_filename(filename):
    """Parse the ABF filename to extract date and file number."""
    try:
        date_str, number_str = filename.split("-")[0], filename.split("-")[-1].split(".")[0]
        date = date_str.split("_")
        file_number = int(number_str)
        return date, file_number
    except Exception as e:
        logging.warning(f"Failed to parse filename: {filename}, error: {e}")
        return None, None

def determine_group_from_csv(date, file_number, group_df):
    """
    Assign group from date (as list/str) and file number, based on CSV mapping.
    """
    if isinstance(date, list) or isinstance(date, tuple):
        date_str = "_".join(str(d) for d in date)
    else:
        date_str = str(date)
    sub = group_df[group_df['date'] == date_str]
    for _, row in sub.iterrows():
        if row['file_number_start'] <= file_number <= row['file_number_end']:
            return row['group']
    return "Other"

def detect_and_correct_scale(t, v):
    """Detects and corrects for voltage scale issues."""
    if np.any(v < -100):
        duration_over_threshold = t[np.where(v < -100)[0][-1]] - t[0]
        if duration_over_threshold > 0.05:
            logging.info("Scale correction applied.")
            v = v / 10
            return True, v
    return False, v

def extract_spike_features(t, v, abf, sweep_number, fin_epoch_time, error_logger, filename):
    """Extract spike features for one sweep. Returns dict of features."""
    try:
        if getattr(abf._adcSection, 'fTelegraphAdditGain', [1])[0] == 5:
            v = v * 5
        putative_spikes = detect_putative_spikes(v, t, filter=None, dv_cutoff=40)
        feature_extractor = SpikeFeatureExtractor(start=np.min(t), end=np.max(t), filter=None)
        features = feature_extractor.process(t, v, i=None)

        peak_t = features["peak_t"].values[0]
        threshold_t = features["threshold_t"].values[0]
        width_ms = features["width"].values[0] * 1000 if "width" in features else np.nan
        amplitude_pa = round(features["peak_v"].values[0] - features["threshold_v"].values[0], 3)
        threshold_pa = round(features["threshold_v"].values[0], 3)
        duration_pa_threshold_fast_trough_ms = round((features["fast_trough_t"].values[0] - features["threshold_t"].values[0]) * 1000, 3)
        duration_threshold_peak_ms = peak_t - threshold_t
        d1_t = threshold_t + (duration_threshold_peak_ms / 2)
        d1_t_index = np.argmin(np.abs(t - d1_t))
        d1_v = v[d1_t_index]
        index_peak = np.argmin(np.abs(t - peak_t))
        indices_after_peak = np.where(v[index_peak:] <= d1_v)[0]
        if len(indices_after_peak) > 0:
            first_index_below_d1_v = index_peak + indices_after_peak[0]
            time_below_d1_v = t[first_index_below_d1_v]
            duration_pa = time_below_d1_v - d1_t
        else:
            duration_pa = np.nan

        file_info_suffix = ""
        if (threshold_t is not np.nan) and (threshold_t < fin_epoch_time):
            file_info_suffix = " threshold < de2"

        return {
            "Width (ms)": width_ms,
            "Amplitude (mV)": amplitude_pa,
            "Threshold (mV)": threshold_pa,
            "Duration (ms)": duration_pa_threshold_fast_trough_ms,
            "Durée_PA": duration_pa,
            "File_Info_Suffix": file_info_suffix
        }
    except Exception as e:
        error_logger.error(f"Error extracting features: file {filename} sweep {sweep_number} -- {str(e)}\n{traceback.format_exc()}")
        return {
            "Width (ms)": np.nan,
            "Amplitude (mV)": np.nan,
            "Threshold (mV)": np.nan,
            "Duration (ms)": np.nan,
            "Durée_PA": np.nan,
            "File_Info_Suffix": ""
        }

def process_abf_file(filepath, error_logger, group_df):
    filename = os.path.basename(filepath)
    try:
        date, file_number = parse_filename(filename)
        if not date or file_number is None:
            raise ValueError("Parsing error")
        group_name = determine_group_from_csv(date, file_number, group_df)
        if not group_name:
            error_logger.error(f"File {filename} -- Group not found, assigned to 'Other'")
            group_name = "Other"
        abf = pyabf.ABF(filepath)
        # Integrity checks
        if abf.sweepCount == 0:
            error_logger.error(f"File {filename} -- No sweeps found (corrupt or empty file)")
            return None, []
        try:
            epoch_level = round(abf.sweepEpochs.levels[2], 3)
        except Exception:
            epoch_level = np.nan
        try:
            fin_epoch_index_2 = abf.sweepEpochs.p2s[2]
            fin_epoch_time = fin_epoch_index_2 / abf.dataRate
        except Exception:
            fin_epoch_time = np.nan

        sweep_results = []
        for sweep_number in abf.sweepList:
            try:
                abf.setSweep(sweep_number)
                t = abf.sweepX
                v = abf.sweepY
                if np.max(v) <= 0:
                    error_logger.error(f"File {filename} sweep {sweep_number} -- No positive voltage detected (empty or corrupt sweep)")
                    continue
                scale_correction, v_corr = detect_and_correct_scale(t, v)
                features = extract_spike_features(
                    t, v_corr, abf, sweep_number, fin_epoch_time, error_logger, filename)
                file_info = f"{filename}_{sweep_number}{features['File_Info_Suffix']}"
                if scale_correction:
                    file_info += "/10"
                sweep_results.append({
                    "Width (ms)": features["Width (ms)"],
                    "Amplitude (mV)": features["Amplitude (mV)"],
                    "Threshold (mV)": features["Threshold (mV)"],
                    "Duration (ms)": features["Duration (ms)"],
                    "Durée_PA": features["Durée_PA"],
                    "Epoch_Level": epoch_level,
                    "File_Info": file_info
                })
            except Exception as e:
                error_logger.error(f"Error in file {filename} sweep {sweep_number}: {str(e)}\n{traceback.format_exc()}")
                continue
        return group_name, sweep_results
    except Exception as e:
        error_logger.error(f"Critical error for file {filename}: {str(e)}\n{traceback.format_exc()}")
        return None, []

def process_all_files(data_dir, error_logger, group_df, group_list):
    """Process all ABF files in a directory with error resilience."""
    grouped_results = {g: [] for g in group_list}
    all_files = [f for f in os.listdir(data_dir) if f.endswith(".abf")]
    if len(all_files) == 0:
        error_logger.error(f"No ABF files found in directory: {data_dir}")
    logging.info(f"Found {len(all_files)} ABF files in {data_dir}")
    for filename in all_files:
        filepath = os.path.join(data_dir, filename)
        group_name, results = process_abf_file(filepath, error_logger, group_df)
        if group_name and results:
            grouped_results[group_name].extend(results)
    return grouped_results

def save_results_to_excel(grouped_results, excel_file, meta_dict):
    """Saves all group DataFrames to an Excel file."""
    with pd.ExcelWriter(excel_file) as writer:
        summary_rows = []
        for group, results in grouped_results.items():
            df = pd.DataFrame(results)
            for col in EXPECTED_COLS:
                if col not in df.columns:
                    df[col] = np.nan
            df = df[EXPECTED_COLS]
            if not df.empty:
                df_ok = df[~df["File_Info"].astype(str).str.contains("threshold < de2", na=False)]
            else:
                df_ok = df
            df_ok.to_excel(writer, sheet_name=group, index=False)
            summary_rows.append({
                "Group": group,
                "N_sweeps": len(df_ok),
                "Width mean": df_ok["Width (ms)"].mean() if not df_ok.empty else np.nan,
                "Amplitude mean": df_ok["Amplitude (mV)"].mean() if not df_ok.empty else np.nan
            })
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name="Summary", index=False)
        # Metadata
        pd.DataFrame(meta_dict).to_excel(writer, sheet_name="Metadata", index=False)
    logging.info(f"Excel file '{excel_file}' created successfully.")

def export_script_copy(script_path, output_dir):
    """Export the current script file to the results directory for reproducibility."""
    try:
        shutil.copy(script_path, os.path.join(output_dir, os.path.basename(script_path)))
        logging.info(f"Script copied to {output_dir} for reproducibility.")
    except Exception as e:
        logging.warning(f"Could not copy script to results folder: {e}")

def generate_qc_plots(grouped_results, output_dir):
    """Generate QC plots (hist, boxplots) for each group."""
    for group, results in grouped_results.items():
        df = pd.DataFrame(results)
        if not df.empty:
            for feature in ["Width (ms)", "Amplitude (mV)", "Threshold (mV)"]:
                plt.figure(figsize=(8, 4))
                plt.hist(df[feature].dropna(), bins=20, alpha=0.7, edgecolor='k')
                plt.title(f"{group}: Histogram of {feature}")
                plt.xlabel(feature)
                plt.ylabel("Count")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{group}_{feature.replace(' ','_')}_hist.png"))
                plt.close()

                plt.figure(figsize=(4, 6))
                plt.boxplot(df[feature].dropna(), vert=True, patch_artist=True)
                plt.title(f"{group}: Boxplot of {feature}")
                plt.ylabel(feature)
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{group}_{feature.replace(' ','_')}_boxplot.png"))
                plt.close()

def save_full_excel(grouped_results, full_excel_file):
    """Save all data (sans filtre) pour archivage complet, avec nom personnalisé s'il existe déjà."""
    print(">>> Appel à save_full_excel")
    dfs = {}
    for group, results in grouped_results.items():
        df = pd.DataFrame(results)
        for col in EXPECTED_COLS:
            if col not in df.columns:
                df[col] = np.nan
        dfs[group] = df[EXPECTED_COLS]
    # Gestion du nom de fichier
    file_exists = os.path.isfile(full_excel_file)
    if file_exists:
        print(f"Le fichier Excel '{full_excel_file}' existe déjà.")
        replace_file = input("Voulez-vous le remplacer? (O/N): ").strip().upper()
        if replace_file == "N":
            new_file_name = input("Veuillez entrer un nouveau nom de fichier Excel: ").strip()
            full_excel_file = new_file_name + ".xlsx"
    # Sauvegarde
    with pd.ExcelWriter(full_excel_file) as writer:
        for group, df in dfs.items():
            print(f"Group: {group}, shape: {df.shape}")
            df.to_excel(writer, sheet_name=group, index=False)
    print(f"Le fichier Excel '{full_excel_file}' a été créé avec succès.")

def validate_group_mapping(group_df):
    """Valide l'absence de range qui se chevauchent et de noms de groupes dupliqués, affiche un warning si mapping ambigu."""
    warnings = []
    # Vérif groupes uniques
    group_names = group_df["group"].tolist()
    if len(set(group_names)) < len(group_names):
        warnings.append("WARNING: Some group names are duplicated in the mapping.")
    # Vérif overlap de ranges pour chaque date
    for date, subdf in group_df.groupby("date"):
        ranges = []
        for _, row in subdf.iterrows():
            start, end = row["file_number_start"], row["file_number_end"]
            ranges.append((start, end, row["group"]))
        # Test overlap
        ranges = sorted(ranges)
        for i in range(len(ranges)-1):
            _, end1, g1 = ranges[i]
            start2, _, g2 = ranges[i+1]
            if start2 <= end1:
                warnings.append(f"WARNING: Overlapping file number ranges for date {date}: groups '{g1}' and '{g2}' overlap.")
    # Affiche warnings
    if warnings:
        for w in warnings:
            print(w)
        print("Please check 'group_mapping.csv' for ambiguities.\n")
    else:
        print("Mapping validation: OK (no overlap, all group names unique)")

def generate_summary_report(grouped_results, output_dir):
    """Génère un rapport texte et markdown récapitulatif dans le dossier résultats."""
    report_lines = []
    report_lines.append("# Spike Feature Extraction – Summary Report\n")
    report_lines.append(f"Date : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    report_lines.append("## Number of files and sweeps per group\n")
    report_lines.append("| Group | N files | N sweeps | Width mean | Amplitude mean | Amplitude min | Amplitude max |")
    report_lines.append("|-------|---------|----------|------------|---------------|---------------|---------------|")

    for group, results in grouped_results.items():
        n_files = len(set([r["File_Info"].split("_")[0] for r in results])) if results else 0
        n_sweeps = len(results)
        if n_sweeps > 0:
            amplitudes = [r["Amplitude (mV)"] for r in results if r["Amplitude (mV)"] is not None]
            width_mean = np.nanmean([r["Width (ms)"] for r in results])
            amp_mean = np.nanmean(amplitudes)
            amp_min = np.nanmin(amplitudes)
            amp_max = np.nanmax(amplitudes)
        else:
            width_mean = amp_mean = amp_min = amp_max = float('nan')
        report_lines.append(f"| {group} | {n_files} | {n_sweeps} | {width_mean:.2f} | {amp_mean:.2f} | {amp_min:.2f} | {amp_max:.2f} |")

    # Groupes non utilisés
    unused_groups = [g for g, results in grouped_results.items() if len(results) == 0]
    if unused_groups:
        report_lines.append("\n## Groups with no files or sweeps found\n")
        report_lines.append(", ".join(unused_groups) + "\n")

    # Fichiers "Other"
    if "Other" in grouped_results and grouped_results["Other"]:
        files_other = sorted(set([r["File_Info"].split("_")[0] for r in grouped_results["Other"]]))
        report_lines.append("\n## Files classified as 'Other' (no group found in mapping)\n")
        for f in files_other:
            report_lines.append(f"- {f}")
    else:
        report_lines.append("\n## No files were classified as 'Other'.")

    # Sauvegarde du rapport
    summary_file = os.path.join(output_dir, "summary_report.md")
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print(f"Summary report written to {summary_file}")

def main():
    import platform
    import matplotlib
    import ipfx
    from datetime import datetime

    # Arguments
    parser = argparse.ArgumentParser(description="Spike Feature Extraction from ABF files (see --help).")
    parser.add_argument('--data_dir', type=str, required=True, help='Path to directory with ABF files')
    parser.add_argument('--output', type=str, default="SpikeFeatureResults.xlsx", help='Output Excel filename (will be inside results folder)')
    parser.add_argument('--script_copy', action='store_true', help='Copy script to results directory')
    args = parser.parse_args()

    # 1. Création du dossier horodaté
    now = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    results_dir = f"results_{now}"
    os.makedirs(results_dir, exist_ok=True)
    qc_dir = os.path.join(results_dir, "QC_Plots")
    os.makedirs(qc_dir, exist_ok=True)
    runinfo_dir = os.path.join(results_dir, "runinfo")
    os.makedirs(runinfo_dir, exist_ok=True)
    # 2. Chemins redirigés
    output_excel = os.path.join(results_dir, os.path.basename(args.output))
    log_file = os.path.join(runinfo_dir, "analysis.log")
    error_log_file = os.path.join(runinfo_dir, "error.log")
    full_excel_file = os.path.join(results_dir, "AnalyseSpikeCC-C4PicroCPG.xlsx")

    # 3. Logging (dans le bon dossier)
    error_logger = setup_logging(log_file, error_log_file)
    logging.info("=== Spike Feature Extraction Script Started ===")
    logging.info(f"Input directory: {args.data_dir}")

    # 4. Lecture du mapping groupes depuis CSV + Validation mapping
    group_mapping_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "group_mapping.csv")
    try:
        group_df = pd.read_csv(group_mapping_path)
        validate_group_mapping(group_df)
    except Exception as e:
        logging.error(f"Impossible de lire le fichier group_mapping.csv : {e}")
        sys.exit(1)

    # 5. Génération dynamique de la liste des groupes
    group_list = sorted(group_df["group"].unique().tolist()) + ["Other"]

    # 6. Reproducibilité: export du script si demandé
    if args.script_copy:
        try:
            export_script_copy(__file__, results_dir)
        except Exception as e:
            logging.warning(f"Script copy failed: {e}")

    # 7. Meta-data pour l'onglet Excel
    meta = {
        "Run timestamp": [datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        "Script version": ["1.2.0"],
        "Input folder": [os.path.abspath(args.data_dir)],
        "Python version": [platform.python_version()],
        "Pandas version": [pd.__version__],
        "Numpy version": [np.__version__],
        "pyabf version": [pyabf.__version__],
        "ipfx version": [getattr(ipfx, '__version__', 'unknown')],
        "matplotlib version": [matplotlib.__version__],
        "Log file": [os.path.abspath(log_file)],
        "Error log file": [os.path.abspath(error_log_file)]
    }

    # 8. Traitement principal
    grouped_results = process_all_files(args.data_dir, error_logger, group_df, group_list)
    save_results_to_excel(grouped_results, output_excel, meta)
    try:
        save_full_excel(grouped_results, full_excel_file)
    except Exception as e:
        print(f"Erreur pendant save_full_excel : {e}")

    # 9. QC plots
    generate_qc_plots(grouped_results, qc_dir)

    # 10. Rapport récapitulatif automatique
    generate_summary_report(grouped_results, runinfo_dir)

    # 11. (Optionnel) Copie du README dans le dossier results pour traçabilité
    try:
        readme_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")
        if os.path.isfile(readme_src):
            shutil.copy(readme_src, os.path.join(results_dir, "README.md"))
    except Exception as e:
        logging.warning(f"README copy failed: {e}")
      # 12. Export run parameter (TXT et YAML)
    run_params = {
        "Run timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Script name": os.path.basename(__file__),
        "Script version": "1.2.0",
        "Input folder": os.path.abspath(args.data_dir),
        "Output Excel": output_excel,
        "QC folder": qc_dir,
        "Log file": os.path.abspath(log_file),
        "Error log file": os.path.abspath(error_log_file),
        "Full Excel export": os.path.abspath(full_excel_file),
        "Python version": sys.version,
        "Pandas version": pd.__version__,
        "Numpy version": np.__version__,
        "pyabf version": pyabf.__version__,
        "ipfx version": getattr(__import__('ipfx'), '__version__', 'unknown'),
        "matplotlib version": plt.matplotlib.__version__,
        "Mapping used": os.path.abspath(group_mapping_path),
    }
    export_run_params(run_params, group_mapping_path, runinfo_dir)

    logging.info("=== Script completed successfully ===")
    
if __name__ == "__main__":
    main()


# pour lancer le code suivre la synthaxe suivante :  python SingleSpikechat.py --data_dir "<chemin_vers_tes_fichiers_ABF>" --output "<nom_de_ton_fichier_excel_de_sortie>.xlsx"
 
# python SingleSpikechat.py --data_dir "E:/Electrophy/Data culture/Data culture 4 PicroCP +G/Spike CC/" --output "SpikeFeatureResults.xlsx"

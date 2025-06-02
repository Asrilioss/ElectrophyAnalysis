import os
import efel
import pyabf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
import logging
import sys
import yaml

# ---------------------------
# Logger configuration
# ---------------------------
def get_logger(runinfo_dir):
    logfile_path = os.path.join(runinfo_dir, "pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[
            logging.FileHandler(logfile_path, mode='w', encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)
    return logger

# ---------------------------
# Argument parsing
# ---------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Pipeline analyse ABF Depol 500ms - extraction, QC, export")
    parser.add_argument('--folder', type=str, required=True,
                        help='Chemin du dossier contenant les fichiers ABF')
    parser.add_argument('--export_name', type=str, default="Depol500_CC",
                        help='Nom racine pour le dossier d’export')
    parser.add_argument('--stim_start', type=float, default=499,
                        help='Début de la stimulation (ms, par défaut 499)')
    parser.add_argument('--stim_end', type=float, default=1040,
                        help='Fin de la stimulation (ms, par défaut 1040)')
    parser.add_argument('--features', type=str, nargs="+",
                        default=['Spikecount', 'adaptation_index', 'time_to_first_spike', 'ISI_CV', 'ISI_values'],
                        help="Liste des features efel à extraire")
    parser.add_argument('--retry_errors', action='store_true',
                        help="Si présent, relance seulement les fichiers listés dans ABF_errors.csv du dernier run")
    return parser.parse_args()

# ---------------------------
# Extraction de features
# ---------------------------
def extract_sweeps_from_abf(file_path, file_name, stim_start_time, stim_end_time, features_efel, error_list=None, logger=None):
    abf = pyabf.ABF(file_path)
    sweeps = []
    for sweep in abf.sweepList:
        try:
            abf.setSweep(sweep)
            if abf._adcSection.fTelegraphAdditGain[0] == 5:
                abf.setSweep(sweep, channel=0)
                v = abf.sweepY
                v *= 5

            stim_start_index = np.where(abf.sweepX*1000 >= stim_start_time)[0][0]
            stim_end_index = np.where(abf.sweepX*1000 <= stim_end_time)[-1][-1]
            trace = {'T': abf.sweepX[stim_start_index:stim_end_index+1]*1000,
                     'V': abf.sweepY[stim_start_index:stim_end_index+1],
                     'stim_start': [stim_start_time],
                     'stim_end': [stim_end_time]}
            traces = [trace]

            efel.api.set_setting('Threshold', 0)
            efel.api.set_setting('DerivativeThreshold', 20)

            feature_values = efel.get_feature_values(
                traces,
                features_efel,
                raise_warnings=None
            )[0]
            current_mean = np.average(
                abf.sweepC[int(500 * abf.dataPointsPerMs):int(1000 * abf.dataPointsPerMs)]
            )
            
            sweep_dict = {
                'File_name': file_name,
                'Sweep': sweep,
                'Current_step': current_mean,
                'Spikecount': feature_values.get('Spikecount', [None])[0] if feature_values.get('Spikecount') else 0,
                'adaptation': None,
                'Latency_ms': feature_values.get('time_to_first_spike', [None])[0] if feature_values.get('time_to_first_spike') else None,
                'ISI_CV': None,
                'ISI_mean_ms': None
            }
            spikecount = sweep_dict['Spikecount']
            if spikecount is not None and spikecount > 4:
                adaptation_index = feature_values.get('adaptation_index')
                if adaptation_index is not None and len(adaptation_index) > 0 and adaptation_index[0] is not None:
                    sweep_dict['adaptation'] = adaptation_index[0]
                ISI_CV = feature_values.get('ISI_CV')
                if ISI_CV is not None and len(ISI_CV) > 0 and ISI_CV[0] is not None:
                    sweep_dict['ISI_CV'] = ISI_CV[0]
                ISI_values = feature_values.get('ISI_values')
                if ISI_values is not None and len(ISI_values) > 0 and ISI_values[0] is not None:
                    sweep_dict['ISI_mean_ms'] = ISI_values[0] / 1000

            sweeps.append(sweep_dict)
        except Exception as e:
            if logger:
                logger.error(f"Erreur sur {file_name} sweep {sweep}: {e}")
            if error_list is not None:
                error_list.append({'File_name': file_name, 'Sweep': sweep, 'Error': str(e)})
    return sweeps

# ---------------------------
# Extraction tous sweeps
# ---------------------------
def extract_all_sweeps(folder_path, stim_start_time, stim_end_time, features_efel, error_list=None, logger=None, files_to_analyze=None):
    if files_to_analyze is not None:
        abf_files = files_to_analyze
    else:
        abf_files = [f for f in os.listdir(folder_path) if f.endswith('.abf')]
    all_sweeps = []
    for file in abf_files:
        file_path = os.path.join(folder_path, file)
        if not os.path.exists(file_path):
            if logger:
                logger.warning(f"Fichier {file_path} introuvable, ignoré.")
            continue
        sweeps = extract_sweeps_from_abf(file_path, file, stim_start_time, stim_end_time, features_efel, error_list=error_list, logger=logger)
        all_sweeps.extend(sweeps)
    return pd.DataFrame(all_sweeps), abf_files

# ---------------------------
# Moyenne spikecount
# ---------------------------
def compute_spikecount_means(final_data, abf_files):
    abf_files_sorted = sorted(abf_files, key=lambda x: int(x.split("-")[-1].split(".")[0]))
    average_spike_counts = []
    moyenne_files = []

    previous_file_number = None
    previous_group_files = []

    for file in abf_files_sorted:
        file_number = int(file.split("-")[-1].split(".")[0])
        if previous_file_number is not None and file_number == previous_file_number + 1:
            previous_group_files.append(file)
        else:
            if previous_group_files:
                moyenne_files.append(previous_group_files)
                group_data = final_data[final_data['File_name'].isin(previous_group_files)]
                max_sweeps = int(group_data['Sweep'].max()) + 1
                all_sweeps = np.arange(max_sweeps)
                all_files = previous_group_files
                idx = pd.MultiIndex.from_product([all_files, all_sweeps], names=["File_name", "Sweep"])
                df_full = pd.DataFrame(index=idx).reset_index()
                df_merged = pd.merge(
                    df_full,
                    group_data[["File_name", "Sweep", "Spikecount"]],
                    on=["File_name", "Sweep"],
                    how="left"
                )
                df_merged["Spikecount"] = df_merged["Spikecount"].fillna(0)
                group_average_spike_counts = group_data.groupby('Sweep')['Spikecount'].mean().reset_index()
                group_average_spike_counts['Files_Moyennes'] = [previous_group_files] * len(group_average_spike_counts)
                average_spike_counts.append(group_average_spike_counts)
            previous_group_files = [file]
        previous_file_number = file_number

    if previous_group_files:
        moyenne_files.append(previous_group_files)
        group_data = final_data[final_data['File_name'].isin(previous_group_files)]
        max_sweeps = int(group_data['Sweep'].max()) + 1
        all_sweeps = np.arange(max_sweeps)
        all_files = previous_group_files
        idx = pd.MultiIndex.from_product([all_files, all_sweeps], names=["File_name", "Sweep"])
        df_full = pd.DataFrame(index=idx).reset_index()
        df_merged = pd.merge(
            df_full,
            group_data[["File_name", "Sweep", "Spikecount"]],
            on=["File_name", "Sweep"],
            how="left"
        )
        df_merged["Spikecount"] = df_merged["Spikecount"].fillna(0)
        group_average_spike_counts = df_merged.groupby("Sweep")["Spikecount"].mean().reset_index()
        group_average_spike_counts['Files_Moyennes'] = [previous_group_files] * len(group_average_spike_counts)
        average_spike_counts.append(group_average_spike_counts)

    average_spike_count_df = pd.concat(average_spike_counts, ignore_index=True)
    average_spike_count_df = average_spike_count_df.reindex(sorted(average_spike_count_df.columns), axis=1)
    return average_spike_count_df, average_spike_counts

# ---------------------------
# Export paramètres de run et summary
# ---------------------------
def export_run_params(params_dict, output_dir):
    params_txt = os.path.join(output_dir, "run_params.txt")
    params_yaml = os.path.join(output_dir, "run_params.yaml")
    with open(params_txt, "w", encoding="utf-8") as f:
        for k, v in params_dict.items():
            f.write(f"{k}: {v}\n")
    with open(params_yaml, "w", encoding="utf-8") as f:
        yaml.dump(params_dict, f, allow_unicode=True)

def generate_summary_report(final_data, abf_files, error_list, output_dir):
    lines = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines.append(f"# ABF Depol 500ms Summary Report\n")
    lines.append(f"- Run date: {now}")
    lines.append(f"- Number of files analyzed: {len(abf_files)}")
    lines.append(f"- Number of sweeps: {len(final_data)}")
    lines.append(f"- Number of error files: {len(error_list)}\n")
    lines.append("## Spikecount Statistics (global)")
    lines.append(str(final_data['Spikecount'].describe()))
    lines.append("\n## Files with 0 spike sweeps:")
    for file in abf_files:
        nb0 = len(final_data[(final_data['File_name'] == file) & (final_data['Spikecount'] == 0)])
        if nb0 > 0:
            lines.append(f"- {file}: {nb0} sweeps à 0 spike")
    if error_list:
        lines.append("\n## Files in error:")
        for err in error_list:
            lines.append(f"- {err['File_name']} sweep {err['Sweep']}: {err['Error']}")
    summary_path = os.path.join(output_dir, "summary_report.md")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

# ---------------------------
# MAIN
# ---------------------------
if __name__ == "__main__":
    args = parse_args()

    folder_path = args.folder
    export_name = args.export_name
    stim_start_time = args.stim_start
    stim_end_time = args.stim_end
    features_efel = args.features
    retry_errors = args.retry_errors

    dt = datetime.now().strftime("%Y%m%d_%H%M%S")
    export_dir = f"export_{export_name}_{dt}"
    os.makedirs(export_dir, exist_ok=True)
    stats_dir = os.path.join(export_dir, "stats_groupes")
    os.makedirs(stats_dir, exist_ok=True)
    runinfo_dir = os.path.join(export_dir, "runinfo")
    os.makedirs(runinfo_dir, exist_ok=True)

    logger = get_logger(runinfo_dir)
    OUTLIER_SEM_THRESHOLD = 5

    # ----------- 1. Relance sur erreurs si demandé -----------
    error_list = []
    files_to_analyze = None
    if retry_errors:
        error_file = os.path.join(export_dir, "runinfo", "ABF_errors.csv")
        prev_export_dirs = sorted([d for d in os.listdir(".") if d.startswith(f"export_{export_name}_") and d != os.path.basename(export_dir)])
        if prev_export_dirs:
            prev_dir = prev_export_dirs[-1]
            prev_error_file = os.path.join(prev_dir, "runinfo", "ABF_errors.csv")
            if not os.path.exists(prev_error_file):
                logger.error("Aucun ABF_errors.csv trouvé dans le dossier d'export précédent. Impossible de relancer uniquement les erreurs.")
                sys.exit(1)
            else:
                df_errors = pd.read_csv(prev_error_file)
                files_to_analyze = df_errors['File_name'].unique().tolist()
                logger.info(f"Relance uniquement sur les fichiers ayant échoué : {files_to_analyze}")
        else:
            logger.error("Aucun dossier d'export précédent trouvé.")
            sys.exit(1)

    final_data, abf_files = extract_all_sweeps(
        folder_path, stim_start_time, stim_end_time, features_efel,
        error_list=error_list, logger=logger, files_to_analyze=files_to_analyze
    )

    # ----------- 2. Export Excel et CSV des résultats bruts -----------
    excel_file_path = os.path.join(export_dir, 'Depol500_defchat.xlsx')
    final_data.to_excel(excel_file_path, index=False)
    logger.info(f"Les données ont été enregistrées dans {excel_file_path}")

    csv_file_path = os.path.join(export_dir, 'Depol500_defchat.csv')
    final_data.to_csv(csv_file_path, index=False)
    logger.info(f"Les données ont aussi été enregistrées en CSV sous {csv_file_path}")

    # ----------- 3. Calcul et export des moyennes par groupe -----------
    average_spike_count_df, average_spike_counts = compute_spikecount_means(final_data, abf_files)
    average_excel_file_path_by_sweep = os.path.join(export_dir, 'testfulldefchat.xlsx')
    average_spike_count_df.to_excel(average_excel_file_path_by_sweep, index=False)
    logger.info(f"Les moyennes Spike Counts par sweep pour les groupes de fichiers consécutifs ont été enregistrées dans {average_excel_file_path_by_sweep}")

    csv_avg_path = os.path.join(export_dir, 'testfulldefchat.csv')
    average_spike_count_df.to_csv(csv_avg_path, index=False)
    logger.info(f"Les moyennes Spike Counts par sweep ont aussi été enregistrées en CSV sous {csv_avg_path}")

    # ----------- 4. Export détaillé des erreurs (CSV) -----------
    error_report_path = os.path.join(runinfo_dir, "ABF_errors.csv")
    if error_list:
        pd.DataFrame(error_list).to_csv(error_report_path, index=False)
        logger.warning(f"Export des erreurs détaillées sous {error_report_path}")
    else:
        logger.info("Aucune erreur critique à exporter (ABF_errors.csv non créé)")

    # ----------- 5. CONTROLE QUALITE + LOGGING -----------
    log_path = os.path.join(runinfo_dir, "QC_report.txt")
    img_path = os.path.join(export_dir, "QC_SpikecountMoyen.png")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("========== Rapport Qualité Analyse ABF ==========\n")
        f.write(f"Date/heure analyse : {now}\n")
        f.write(f"Dossier de travail : {folder_path}\n")
        f.write(f"Nombre de fichiers ABF analysés : {len(abf_files)}\n")
        f.write(f"Nombre total de sweeps extraits : {len(final_data)}\n")
        f.write(f"Nombre de groupes : {len(average_spike_counts)}\n\n")
        f.write("---- Nombre de sweeps par fichier (QC) ----\n")
        sweeps_per_file = final_data.groupby("File_name")["Sweep"].nunique()
        for file, n in sweeps_per_file.items():
            f.write(f"{file}: {n} sweeps\n")
            if n < final_data['Sweep'].max() + 1:
                f.write(f"   >> WARNING: Sweep manquant (max={final_data['Sweep'].max() + 1})\n")
        f.write("\n")
        # QC : Nb de sweeps à 0 spike
        f.write("---- Sweeps à 0 spike (par fichier) ----\n")
        for file in abf_files:
            nb0 = len(final_data[(final_data['File_name'] == file) & (final_data['Spikecount'] == 0)])
            f.write(f"{file}: {nb0} sweeps à 0 spike\n")
        f.write("\n")
        # Stats de base
        f.write("---- Statistiques spikecount (global) ----\n")
        desc = final_data['Spikecount'].describe()
        f.write(str(desc) + "\n\n")
        f.write("---- Distribution spikecount par sweep ----\n")
        by_sweep = final_data.groupby("Sweep")["Spikecount"].describe()
        f.write(str(by_sweep) + "\n\n")
        # Résumé erreurs
        f.write("---- Récapitulatif erreurs critiques (voir ABF_errors.csv) ----\n")
        f.write(f"Nombre total d'erreurs lors de l'analyse : {len(error_list)}\n")
        if error_list:
            for err in error_list:
                f.write(f"{err['File_name']} sweep {err['Sweep']}: {err['Error']}\n")
        else:
            f.write("Aucune erreur critique signalée.\n")
        f.write("\n")

    logger.info(f"Rapport qualité automatique sauvegardé sous {log_path}")

    # ----------- 6. GRAPHIQUE QC -----------
    plt.figure(figsize=(10,5))
    final_data.groupby('Sweep')['Spikecount'].mean().plot(kind='bar')
    plt.title("Spikecount moyen par sweep (tous groupes)")
    plt.xlabel("Sweep")
    plt.ylabel("Spikecount moyen")
    plt.tight_layout()
    plt.savefig(img_path, dpi=300)
    plt.close()
    logger.info(f"Graphique QC sauvegardé sous {img_path}")
    logger.info("Analyse terminée. Vérifie les fichiers de log et graphique pour le contrôle qualité.")

    # ----------- 7. README exporté dans le sous-dossier runinfo -----------
    readme_path = os.path.join(runinfo_dir, "README_ANALYSE.txt")
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write("=================== README - Analyse ABF Depol 500ms ===================\n\n")
        f.write(f"Date de l'analyse : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Dossier de provenance : {folder_path}\n\n")
        f.write("Paramètres de l’analyse :\n")
        f.write("- Extraction des fichiers ABF par pyabf\n")
        f.write("- Extraction des features électrophysiologiques par efel\n")
        f.write("- Durée de la fenêtre d’analyse : {} ms à {} ms\n".format(stim_start_time, stim_end_time))
        f.write("- Features extraites : {}\n".format(", ".join(features_efel)))
        f.write("- Seuil efel (Threshold) : 0\n")
        f.write("- Seuil efel (DerivativeThreshold) : 20\n")
        f.write("- Fichiers manquants/erreurs listés dans QC_report.txt\n\n")
        f.write(f"Version du script : Depol500final_chat.py\n")
        f.write(f"Date du script : {datetime.now().strftime('%Y-%m-%d')}\n\n")
        f.write("Hash du code : [À compléter manuellement si tu utilises git ou md5sum]\n")
        f.write("------------------------------------------------------------------------\n")
        f.write("Ce README décrit les paramètres utilisés pour cette analyse. Toute modification du script ou du pipeline devra être documentée ici pour assurer la traçabilité scientifique.\n")
        f.write("------------------------------------------------------------------------\n")
    logger.info(f"README généré sous {readme_path}")

    # ----------- 8. Export run params et résumé (runinfo) -----------
    params_dict = {
        "Run timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Script name": os.path.basename(__file__),
        "Input folder": os.path.abspath(folder_path),
        "Output folder": os.path.abspath(export_dir),
        "Stimulation window (ms)": [stim_start_time, stim_end_time],
        "Features": features_efel,
        "Python version": sys.version,
        "Numpy version": np.__version__,
        "Pandas version": pd.__version__,
        "pyabf version": pyabf.__version__,
        "efel version": efel.__version__ if hasattr(efel, "__version__") else "unknown",
    }
    export_run_params(params_dict, runinfo_dir)
    generate_summary_report(final_data, abf_files, error_list, runinfo_dir)

    # ----------- 9. Copie du script dans runinfo (optionnel) -----------
    try:
        script_path = os.path.abspath(__file__)
        shutil.copy(script_path, os.path.join(runinfo_dir, os.path.basename(script_path)))
        logger.info("Copie du script pour traçabilité effectuée dans runinfo.")
    except Exception as e:
        logger.warning(f"Impossible de copier le script : {e}")

    # ----------- 10. Génération stats & graphiques par groupe + QC_flag -----------
    logger.info("Génération des statistiques et graphiques par groupe...")

    for i, group_df in enumerate(average_spike_counts, 1):
        group_name = f"Groupe_{i}"
        files = group_df['Files_Moyennes'].iloc[0]
        if isinstance(files, list):
            files_str = ", ".join(files)
        else:
            files_str = str(files)

        group_files = files if isinstance(files, list) else eval(files)
        group_data = final_data[final_data['File_name'].isin(group_files)]
        max_sweeps = int(group_data['Sweep'].max()) + 1
        all_sweeps = np.arange(max_sweeps)
        idx = pd.MultiIndex.from_product([group_files, all_sweeps], names=["File_name", "Sweep"])
        df_full = pd.DataFrame(index=idx).reset_index()
        df_merged = pd.merge(
            df_full,
            group_data[["File_name", "Sweep", "Spikecount"]],
            on=["File_name", "Sweep"],
            how="left"
        )
        df_merged["Spikecount"] = df_merged["Spikecount"].fillna(0)

        stats = df_merged.groupby('Sweep')['Spikecount'].agg(['mean','std','sem']).reset_index()
        stats.columns = ['Sweep', 'Spikecount_mean', 'Spikecount_std', 'Spikecount_sem']
        stats["QC_flag"] = np.where(stats["Spikecount_sem"] > OUTLIER_SEM_THRESHOLD, "outlier", "ok")
        stats["is_valid"] = stats["Spikecount_sem"] <= OUTLIER_SEM_THRESHOLD

        stats_path = os.path.join(stats_dir, f"Stats_{group_name}.csv")
        stats.to_csv(stats_path, index=False)
        logger.info(f"Table stats exportée pour {group_name} sous {stats_path}")

        plt.figure(figsize=(8,4))
        plt.errorbar(stats['Sweep'], stats['Spikecount_mean'], yerr=stats['Spikecount_sem'], fmt='-o', capsize=3, label='Spikecount moyen ± SEM')
        plt.title(f"{group_name} ({files_str})")
        plt.xlabel("Sweep")
        plt.ylabel("Spikecount moyen")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(stats_dir, f"{group_name}_mean_sem.png"), dpi=200)
        plt.close()
        logger.info(f"Graphique exporté sous {os.path.join(stats_dir, f'{group_name}_mean_sem.png')}")

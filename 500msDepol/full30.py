import os
import efel
import pyabf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import logging
import sys
import yaml
import shutil
import warnings
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
warnings.filterwarnings('ignore', category=OptimizeWarning)

# ignore les RuntimeWarning spécifiques à EFEL (pic non détecté, sag_amplitude invalide, etc.)
warnings.filterwarnings(
    'ignore',
    message=r'Error while calculating .*',
    module=r'efel\.pyfeatures\.cppfeature_access'
)
# --------------------------- Logger configuration ---------------------------
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
# ── Helpers ──────────────────────────────────────────────────────────────
def get_current_steps_info(file_path):
    """Retourne les paliers de courant et leur incrément pour un fichier ABF."""
    abf = pyabf.ABF(file_path)
    all_steps = set()
    for sw in [0, 1, abf.sweepCount - 1]:
        abf.setSweep(sw)
        all_steps.update(np.unique(np.round(abf.sweepC, 2)))
    steps = np.array(sorted(all_steps))
    incr  = steps[1] - steps[0] if len(steps) > 1 else np.nan
    return {'steps': steps, 'incr': incr}

def safe_get(d, key, idx=0):
    """
    Renvoie d[key][idx] que ce soit une liste Python ou un numpy array,
    ou NaN si la clé n'existe pas ou si l'index est hors-limites.
    """
    arr = d.get(key)
    if arr is None:
        return np.nan
    try:
        # fonctionne pour list ou numpy.ndarray
        return arr[idx]
    except Exception:
        return np.nan

def safe_list(d, key):
    """
    Renvoie toujours une Python list, même si EFEL retourne un ndarray,
    et [] si la clé n'existe pas.
    """
    arr = d.get(key)
    if arr is None:
        return []
    try:
        return list(arr)
    except Exception:
        return []


def extract_sweeps_from_abf(file_path, file_name, stim_start_time, stim_end_time,
                             features_efel, error_list=None, logger=None):
    abf = pyabf.ABF(file_path)
    sweeps = []
    for sweep in abf.sweepList:
        try:
            abf.setSweep(sweep)
            # amplify if needed
            if abf._adcSection.fTelegraphAdditGain[0] == 5:
                v = abf.sweepY * 5
            else:
                v = abf.sweepY
            t = abf.sweepX * 1000

            # frame indices
            stim_i = np.where(t >= stim_start_time)[0][0]
            stim_f = np.where(t <= stim_end_time)[-1][-1]

            trace = {'T': t.tolist(),
                     'V': v.tolist(),
                     'stim_start': [stim_start_time],
                     'stim_end': [stim_end_time]}
            I = abf.sweepC
            midx1 = int(stim_start_time * abf.dataPointsPerMs)
            midx2 = int(stim_end_time * abf.dataPointsPerMs)
            mean_I_pA = np.mean(I[midx1:midx2])
            trace['stimulus_current'] = [mean_I_pA / 1000.0]  # en nA
            # add stimulus_current if needed
            if 'ohmic_input_resistance' in features_efel:
                I = abf.sweepC
                midx1 = int(stim_start_time * abf.dataPointsPerMs)
                midx2 = int(stim_end_time * abf.dataPointsPerMs)
                mean_I_pA = np.mean(I[midx1:midx2])
                trace['stimulus_current'] = [mean_I_pA / 1000.0]

            efel.api.set_setting('Threshold', 0)
            efel.api.set_setting('DerivativeThreshold', 20)
            res = efel.get_feature_values([trace], features_efel)[0]

            # assemble result dict
            sweep_dict = {'File_name': file_name, 'Sweep': sweep}
            # current step mean
            sweep_dict['Injected_current_pA'] = mean_I_pA if 'ohmic_input_resistance' in features_efel else np.nan
            for feat in features_efel:
                val = res.get(feat)
                if isinstance(val, list):
                    sweep_dict[feat] = val[0] if val else None
                else:
                    sweep_dict[feat] = val
            sweeps.append(sweep_dict)
        except Exception as e:
            if logger:
                logger.error(f"Erreur sur {file_name} sweep {sweep}: {e}")
            if error_list is not None:
                error_list.append({'File_name': file_name, 'Sweep': sweep, 'Error': str(e)})
    return sweeps

# --------------------------- Extraction tous sweeps ---------------------------
def extract_all_sweeps(folder_path, stim_start_time, stim_end_time,
                       error_list=None, logger=None, files_to_analyze=None):
    if files_to_analyze:
        abf_files = files_to_analyze
    else:
        abf_files = [f for f in os.listdir(folder_path) if f.endswith('.abf')]
    all_sweeps = []
    for file in abf_files:
        file_path = os.path.join(folder_path, file)
        if not os.path.exists(file_path):
            logger.warning(f"Fichier {file_path} introuvable, ignoré.")
            continue
        proto = detect_protocol(file_path, stim_start_time, stim_end_time)
        features = FEATURE_MAP.get(proto, [])
        logger.info(f"Detected protocol '{proto}' for file {file}")
        sweeps = extract_sweeps_from_abf(file_path, file, stim_start_time,
                                         stim_end_time, features,
                                         error_list=error_list, logger=logger)
        all_sweeps.extend(sweeps)
    df = pd.DataFrame(all_sweeps)
    # ensure all possible feature columns exist
    all_features = set().union(*FEATURE_MAP.values())
    for feat in all_features.union({'Injected_current_pA'}):
        if feat not in df.columns:
            df[feat] = np.nan
    return df, abf_files

# (the rest of the pipeline: compute_spikecount_means, exports, summary, QC, etc.)
def compute_spikecount_means(final_data, abf_files):
    # Tri des fichiers selon leur numéro
    abf_files_sorted = sorted(
        abf_files,
        key=lambda x: int(x.split("-")[-1].split(".")[0])
    )
    average_spike_counts = []
    current_group = []
    previous_num = None

    def process_group(group_files):
        # Construit une table complète (toutes files × tous sweeps)
        sub = final_data[final_data['File_name'].isin(group_files)]
        max_sweep = int(sub['Sweep'].max()) + 1
        all_idx = pd.MultiIndex.from_product(
            [group_files, np.arange(max_sweep)],
            names=["File_name", "Sweep"]
        )
        full = pd.DataFrame(index=all_idx).reset_index()
        merged = pd.merge(
            full,
            sub[["File_name", "Sweep", "Spikecount"]],
            on=["File_name", "Sweep"],
            how="left"
        )
        merged["Spikecount"] = merged["Spikecount"].fillna(0)
        # Moyenne et SEM par sweep
        stats = merged.groupby("Sweep")["Spikecount"] \
                      .agg(['mean', 'sem']) \
                      .reset_index() \
                      .rename(columns={
                          'mean': 'Spikecount_mean',
                          'sem':  'Spikecount_sem'
                      })
        stats['Files_Moyennes'] = [group_files] * len(stats)
        return stats

    # Parcours et regroupement par files consécutifs
    for file in abf_files_sorted:
        num = int(file.split("-")[-1].split(".")[0])
        if previous_num is not None and num == previous_num + 1:
            current_group.append(file)
        else:
            if current_group:
                average_spike_counts.append(process_group(current_group))
            current_group = [file]
        previous_num = num

    # dernier groupe
    if current_group:
        average_spike_counts.append(process_group(current_group))

    # DataFrame globale (si besoin)
    average_spike_count_df = pd.concat(average_spike_counts, ignore_index=True)
    average_spike_count_df = average_spike_count_df.reindex(
        sorted(average_spike_count_df.columns),
        axis=1
    )
    return average_spike_count_df, average_spike_counts


# ── Extraction complète des 32 paramètres Petilla ────────────────────────
def extract_all_parameters(folder_path, stim_start=220, stim_end=1000):
    rows = []
    abf_files = [f for f in os.listdir(folder_path) if f.endswith('.abf')]
    for file in abf_files:
        path = os.path.join(folder_path, file)
        # calcul incrément global pour R_in
        info   = get_current_steps_info(path)
        incr   = info['incr']
        abf    = pyabf.ABF(path)
        for sw in abf.sweepList:
            abf.setSweep(sw)
            T = abf.sweepX * 1000  # ms
            V = abf.sweepY         # mV
            I = abf.sweepC         # pA

            # indices de la fenêtre stimulus
            i0 = np.argmin(np.abs(T - stim_start))
            i1 = np.argmin(np.abs(T - stim_end))

            # 1) Resting, v_peak, v_ss, ΔV, current
            v_rest  = np.mean(V[:i0])
            v_peak  = np.min(V[i0:i1])
            v_ss    = np.mean(V[i1-5:i1])
            v_delta = v_ss - v_rest
            current = float(I[(i0 + i1)//2])

            # 2) Input resistance R = ΔV / ΔI
            R_in = v_delta / incr if not np.isnan(incr) and incr!=0 else np.nan
            if R_in <= 0 or np.isclose(R_in,0): R_in = np.nan

            # 3) tau & 4) C_m
            decay = efel.get_feature_values([{
                'T': T.tolist(), 'V': V.tolist(),
                'stim_start': [stim_start], 'stim_end': [stim_end],
                'decay_start_after_stim': [1.0], 'decay_end_after_stim': [30.0]
            }], ['decay_time_constant_after_stim'])[0]
            tau_decay = safe_get(decay, 'decay_time_constant_after_stim')
            C_m = tau_decay / R_in if not np.isnan(R_in) and not np.isnan(tau_decay) else np.nan

            # 5) Ghyp, Gsag & sag_index
            if current < 0 and not np.isclose(v_peak, v_rest):
                ghyp = current / (v_peak - v_rest)
            else:
                ghyp = np.nan

            if current < 0 and not np.isclose(v_ss, v_rest):
                gsag = current / (v_ss - v_rest)
            else:
                gsag = np.nan

            # sag_index selon Halabisky et al. (2006) = (Gsag - Ghyp) / Gsag
            if current < 0 and not np.isclose(gsag, 0):
                sag_index = (gsag - ghyp) / gsag
            else:
                sag_index = np.nan

            # eFEL features de base + AHP/ADP
            feats = [
                'Spikecount','time_to_first_spike',
                'AP_amplitude','AP_width','peak_time',
                'AHP_depth_from_peak','AHP_time_from_peak','min_AHP_indices',
                'ADP_peak_amplitude','ADP_peak_indices'
            ]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                vals = efel.get_feature_values([{
                    'T':T.tolist(),'V':V.tolist(),'I':I.tolist(),
                    'stim_start':[stim_start],'stim_end':[stim_end]
                }], feats)[0]

            # 6) Rheobase (sera calculé après regroupement par fichier)

            # 7) First‐spike latency
            first_lat = safe_get(vals, 'time_to_first_spike')

            # 8–9–10) Adaptation près du seuil, F_min & class_firing
            peaks = np.array(safe_list(vals, 'peak_time'))
            if len(peaks) >= 3:
                ISIs = np.diff(peaks[:4])
                inst = 1000.0/ISIs
                tfit = peaks[1:4] - stim_start
                m_thr, F_min = np.polyfit(tfit, inst, 1)
                if m_thr < 0:
                    cls = 'accelerating'
                elif F_min > 10 and m_thr > 0.5:
                    cls = 'burst'
                else:
                    cls = 'continuous'
            else:
                m_thr = F_min = np.nan
                cls   = np.nan

            # 11–14) Biphasic adaptation fit
            if len(peaks) > 1:
                ISI2 = np.diff(peaks)
                finst = 1000.0/ISI2
                tinst = ((peaks[1:]+peaks[:-1])/2) - stim_start
                def fsat(t, A, tau, m, F): return A*np.exp(-t/tau) + m*t + F
                p0 = [finst[0]-finst[-1], 200.0, 0.0, finst[-1]]
                bnds = ([0,1e-3,-np.inf,0], [p0[0]*1.5,1000,np.inf,np.max(finst)])
                Asat, tau_sat, m_sat, Fmax = curve_fit(fsat, tinst, finst, p0=p0, bounds=bnds, maxfev=5000)[0]
            else:
                Asat = tau_sat = m_sat = Fmax = np.nan

            # 15–18) A1, A2, D1, D2
            amps = safe_list(vals, 'AP_amplitude')
            wds  = safe_list(vals, 'AP_width')
            A1, A2 = (amps + [np.nan, np.nan])[:2]
            D1, D2 = (wds  + [np.nan, np.nan])[:2]

            # 19) amp_reduction & 20) dur_increase
            amp_red = (A1 - A2)/A1 if not np.isnan(A1) and not np.isnan(A2) else np.nan
            dur_inc = (D2 - D1)/D1 if not np.isnan(D1) and not np.isnan(D2) else np.nan

            # 21–24) AHP amplitudes & latencies
            ahpd = safe_list(vals, 'AHP_depth_from_peak')
            ahpt = safe_list(vals, 'AHP_time_from_peak')
            AHP1_amp, AHP2_amp     = (ahpd + [np.nan, np.nan])[:2]
            AHP1_lat, AHP2_lat     = (ahpt + [np.nan, np.nan])[:2]

            # 25–28) AHP indices
            ahpi = safe_list(vals, 'min_AHP_indices')
            AHP1_idx, AHP2_idx     = (ahpi + [np.nan, np.nan])[:2]

            # 29–30) ADP amplitudes
            adpa = safe_list(vals, 'ADP_peak_amplitude')
            ADP1_amp, ADP2_amp     = (adpa + [np.nan, np.nan])[:2]

            # 31–32) ADP latencies & indices
            adpi = safe_list(vals, 'ADP_peak_indices')
            adpt = T[[int(ii) for ii in adpi]] if adpi else []
            ADP1_lat, ADP2_lat     = (list(adpt) + [np.nan, np.nan])[:2]
            ADP1_idx, ADP2_idx     = (adpi + [np.nan, np.nan])[:2]

            rows.append({
                'File_name': file, 'Sweep': sw,
                'Resting_mV': v_rest, 'V_peak_mV': v_peak, 'V_ss_mV': v_ss,
                'Delta_V_mV': v_delta, 'Injected_current_pA': current,
                'R_in_Gohm': R_in, 'decay_tau_ms': tau_decay, 'C_m_pF': C_m,
                'Ghyp_nS': ghyp, 'Gsag_nS': gsag, 'sag_index': sag_index,
                'Spikecount': safe_get(vals,'Spikecount'),
                'first_spike_lat_ms': first_lat,
                'm_threshold': m_thr, 'F_min': F_min, 'class_firing': cls,
                'Asat_Hz': Asat, 'tau_sat_ms': tau_sat,
                'm_sat_per_ms': m_sat, 'Fmax_Hz': Fmax,
                'AP1_amp_mV': A1, 'AP2_amp_mV': A2,
                'AP1_width_ms': D1, 'AP2_width_ms': D2,
                'amp_reduction': amp_red, 'dur_increase': dur_inc,
                'AHP1_amp': AHP1_amp, 'AHP2_amp': AHP2_amp,
                'AHP1_latency_ms': AHP1_lat, 'AHP2_latency_ms': AHP2_lat,
                'AHP1_idx': AHP1_idx, 'AHP2_idx': AHP2_idx,
                'ADP1_amp': ADP1_amp, 'ADP2_amp': ADP2_amp,
                'ADP1_latency_ms': ADP1_lat, 'ADP2_latency_ms': ADP2_lat,
                'ADP1_idx': ADP1_idx, 'ADP2_idx': ADP2_idx
            })

    # ── Calcul de la rhéobase (par fichier) ────────────────────────────────
    df = pd.DataFrame(rows)
    df['Rheobase_pA'] = np.nan
    for file in df['File_name'].unique():
        sub = df[df['File_name']==file]
        # premier courant positif avec Spikecount>0
        positives = sub[(sub['Injected_current_pA']>0) & (sub['Spikecount']>0)]
        if not positives.empty:
            rheo = positives['Injected_current_pA'].min()
            df.loc[df['File_name']==file, 'Rheobase_pA'] = rheo

    return df
# --------------------------- Export paramètres de run et summary ---------------------------
# --------------------------- Export paramètres de run et summary ---------------------------
def export_run_params(params_dict, output_dir):
    params_txt  = os.path.join(output_dir, "run_params.txt")
    params_yaml = os.path.join(output_dir, "run_params.yaml")
    with open(params_txt,  "w", encoding="utf-8") as f:
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
    return summary_path
# --------------------------- MAIN ---------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import shutil
    import sys

    # 0) Paramètres à ajuster
    folder_path     = r"F:\Electrophy\Data culture\Data culture 4 PicroCP +G\testrename\IC Firing High Sampling 20pA"
    export_name     = "Depol500_CC20pa"
    stim_start_time = 220    # ms
    stim_end_time   = 1000   # ms

    # 1) Préparation des dossiers
    dt = datetime.now().strftime("%Y%m%d_%H%M%S")
    export_dir  = f"export_{export_name}_{dt}"
    os.makedirs(export_dir,  exist_ok=True)
    stats_dir   = os.path.join(export_dir, "stats_groupes"); os.makedirs(stats_dir,   exist_ok=True)
    runinfo_dir = os.path.join(export_dir, "runinfo");     os.makedirs(runinfo_dir, exist_ok=True)

    logger = get_logger(runinfo_dir)
    OUTLIER_SEM_THRESHOLD = 5

    # 2) Extraction Petilla (32+ paramètres) pour chaque sweep
    logger.info("Extraction Petilla (rhéobase_full) pour chaque sweep…")
    petilla_df = extract_all_parameters(folder_path, stim_start_time, stim_end_time)
    abf_files  = sorted(petilla_df['File_name'].unique().tolist())

    # 3) Extraction des features « FI » / rheobase précédentes
    logger.info("Extraction des features FI / rheobase d'origine…")
    BASE_FEATS = [
    'adaptation_index','mean_frequency','ISI_values',
    'AP_amplitude','AP_width','voltage_deflection',
    'sag_ratio1','sag_amplitude','ohmic_input_resistance',
    'time_to_first_spike',
    'voltage_base',  # nécessaire pour sag_index_corrected
    'steady_state_voltage_stimend'  # idem
    ]

    raw_rows = []
    for file in abf_files:
        file_path = os.path.join(folder_path, file)
        raw_rows.extend(
            extract_sweeps_from_abf(
                file_path, file,
                stim_start_time, stim_end_time,
                BASE_FEATS,
                error_list=None, logger=logger
            )
        )
    raw_df = pd.DataFrame(raw_rows)

    # 4) On retire les colonnes dupliquées et on prefixe le reste
    cols_to_drop = ['Injected_current_pA','Spikecount']
    keep = [c for c in raw_df.columns if c not in cols_to_drop]
    raw_df = raw_df[keep]
    # prefixer toutes les features d'origine hors File_name/Sweep
    rename_map = {c: 'orig_'+c for c in raw_df.columns if c not in ['File_name','Sweep']}
    raw_df = raw_df.rename(columns=rename_map)

    # 5) Merge Petilla + features d'origine
    # Fusion des données Petilla avec les features d'origine
    final_data = petilla_df.merge(raw_df, on=['File_name', 'Sweep'], how='left')

    

    # 6) Export Excel & CSV des résultats bruts
    excel_path = os.path.join(export_dir, f"{export_name}.xlsx")
    final_data.to_excel(excel_path, index=False)
    logger.info(f"Données brutes exportées vers {excel_path}")

    csv_path = os.path.join(export_dir, f"{export_name}.csv")
    final_data.to_csv(csv_path, index=False)
    logger.info(f"Données brutes exportées vers {csv_path}")

    # 4) Moyennes de Spikecount par groupe
    avg_df, average_spike_counts = compute_spikecount_means(final_data, abf_files)
    avg_excel = os.path.join(export_dir, f"{export_name}_mean_spikecount.xlsx")
    avg_df.to_excel(avg_excel, index=False)
    logger.info(f"Moyennes spikecount exportées vers {avg_excel}")

    avg_csv = os.path.join(export_dir, f"{export_name}_mean_spikecount.csv")
    avg_df.to_csv(avg_csv, index=False)
    logger.info(f"Moyennes spikecount exportées vers {avg_csv}")

    # 5) QC report & logging
    log_path = os.path.join(runinfo_dir, "QC_report.txt")
    img_qc   = os.path.join(export_dir, "QC_SpikecountMoyen.png")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("====== QC Analyse ABF ======\n")
        f.write(f"Date       : {datetime.now()}\n")
        f.write(f"Dossier    : {folder_path}\n")
        f.write(f"Fichiers   : {len(abf_files)}\n")
        f.write(f"Sweeps     : {len(final_data)}\n\n")
        f.write("-- Sweeps par fichier --\n")
        for fname, n in final_data.groupby("File_name")["Sweep"].nunique().items():
            f.write(f"{fname}: {n} sweeps\n")
        f.write("\n-- Sweeps à 0 spike --\n")
        zeros = final_data[final_data["Spikecount"] == 0].groupby("File_name").size()
        for fname in abf_files:
            f.write(f"{fname}: {zeros.get(fname,0)} sweeps à 0 spike\n")
    logger.info(f"QC report généré : {log_path}")

    # 6) Graphique QC
    plot_data = final_data.groupby('Sweep')['Spikecount'].mean()
    plt.figure(figsize=(10,5))
    plot_data.plot(kind='bar')
    plt.title("Spikecount moyen par sweep")
    plt.xlabel("Sweep"); plt.ylabel("Spikecount moyen")
    plt.tight_layout()
    plt.savefig(img_qc, dpi=300)
    plt.close()
    logger.info(f"Graphique QC sauvegardé : {img_qc}")

    # 7) README ANALYSE
    readme = os.path.join(runinfo_dir, "README_ANALYSE.txt")
    with open(readme, "w", encoding="utf-8") as f:
        f.write(f"Analyse ABF - {datetime.now()}\n")
        f.write(f"Stim window: {stim_start_time}-{stim_end_time} ms\n")
        f.write("Extraction de 32+ paramètres Petilla pour chaque sweep.\n")
    logger.info(f"README exporté : {readme}")

    # 8) Export run params & summary
    params  = {
        "run_timestamp": datetime.now().isoformat(),
        "input_folder":  os.path.abspath(folder_path),
        "stim_start_ms": stim_start_time,
        "stim_end_ms":   stim_end_time,
        "export_dir":    os.path.abspath(export_dir)
    }
    export_run_params(params, runinfo_dir)
    summary = generate_summary_report(final_data, abf_files, [], runinfo_dir)
    logger.info(f"Summary report : {summary}")

    # 9) Copie du script pour traçabilité
    try:
        shutil.copy(os.path.abspath(__file__),
                    os.path.join(runinfo_dir, os.path.basename(__file__)))
        logger.info("Script copié dans runinfo pour traçabilité")
    except Exception as e:
        logger.warning(f"Impossible de copier le script : {e}")

    # 10) Stats & graphiques par groupe (+ QC_flag)
    for idx, grp in enumerate(average_spike_counts, 1):
        name      = f"Groupe_{idx}"
        stats_df  = grp.rename(columns={'Spikecount':'Spikecount_mean'})
        stats_csv = os.path.join(stats_dir, f"{name}_stats.csv")
        stats_df.to_csv(stats_csv, index=False)
        logger.info(f"Stats exportées : {stats_csv}")
        fig, ax = plt.subplots(figsize=(8,4))
        ax.errorbar(stats_df['Sweep'], stats_df['Spikecount_mean'],
                    yerr=stats_df['Spikecount_sem'], fmt='-o', capsize=3)
        ax.set(title=name, xlabel='Sweep', ylabel='Mean Spikecount')
        fig_path = os.path.join(stats_dir, f"{name}_mean_sem.png")
        fig.tight_layout(); fig.savefig(fig_path, dpi=200); plt.close(fig)
        logger.info(f"Graphique {name} : {fig_path}")

    logger.info("Pipeline terminé, toutes les colonnes Petilla sont dans l'Excel.")

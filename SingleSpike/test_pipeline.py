import os
import pandas as pd
import numpy as np
from SingleSpikechat import process_all_files, setup_logging

def test_process_abf_files():
    # Chemin vers ton dossier test
    data_dir = r"E:\Electrophy\Code_excel\SinglePA Spyke CC Code\test"
    log_file = "test_analysis.log"
    error_log_file = "test_error.log"
    group_mapping_path = os.path.join(os.path.dirname(__file__), "group_mapping.csv")

    error_logger = setup_logging(log_file, error_log_file)
    group_df = pd.read_csv(group_mapping_path)

    # Exécution du pipeline sur le dossier de test
    grouped_results = process_all_files(data_dir, error_logger, group_df)
    
    # --------- Tests de base ----------
    assert isinstance(grouped_results, dict), "Le résultat doit être un dictionnaire de groupes"
    n_total = sum([len(val) for val in grouped_results.values()])
    print(f"Nombre total de sweeps trouvés : {n_total}")
    assert n_total > 0, "Aucun sweep détecté : il y a un souci avec les fichiers ou le parsing."
    for group, res in grouped_results.items():
        if res:  # Si non vide
            for entry in res:
                for col in ["Width (ms)", "Amplitude (mV)", "Threshold (mV)"]:
                    assert col in entry, f"Colonne manquante : {col} pour groupe {group}"

if __name__ == "__main__":
    test_process_abf_files()
    print("✅ Test basique du pipeline terminé (aucune erreur détectée).")

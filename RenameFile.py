import os
import datetime

def rename_abf_by_creation_date(folder_path):
    # Lister tous les fichiers ABF (quel que soit leur nom d'origine)
    abf_files = [f for f in os.listdir(folder_path) if f.lower().endswith('.abf')]
    # Associer chaque fichier à sa date de création (timestamp)
    abf_files_with_ctime = [
        (f, os.path.getctime(os.path.join(folder_path, f))) for f in abf_files
    ]
    # Trier par date de création (croissante = du plus ancien au plus récent)
    abf_files_sorted = sorted(abf_files_with_ctime, key=lambda x: x[1])

    for idx, (file_name, ctime) in enumerate(abf_files_sorted, 1):
        # Conversion de la date de création en jour/mois/année
        dt = datetime.datetime.fromtimestamp(ctime)
        jour = dt.day
        mois = dt.month
        annee = dt.year % 100  # pour avoir 2 chiffres seulement
        # Format du nouveau nom : 17_5_24-001.abf
        new_name = f"{jour}_{mois}_{annee}-{idx:03d}.abf"
        old_file_path = os.path.join(folder_path, file_name)
        new_file_path = os.path.join(folder_path, new_name)
        if not os.path.exists(new_file_path):
            os.rename(old_file_path, new_file_path)
            print(f"{file_name} → {new_name}")
        else:
            print(f"ATTENTION : {new_name} existe déjà, {file_name} ignoré !")

folder_path = r"F:\Electrophy\Data culture\Data culture 4 PicroCP +G\testrename"
rename_abf_by_creation_date(folder_path)

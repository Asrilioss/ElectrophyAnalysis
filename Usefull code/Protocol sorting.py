# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:42:55 2024

@author: Elias
"""

import os
import shutil
from pyabf import ABF

# Chemin du dossier contenant les fichiers ABF
folder_path = r"E:/Electrophy/Data culture 4 CP +G"

# Obtenir la liste des fichiers dans le dossier
file_list = os.listdir(folder_path)

# Créer un dictionnaire pour stocker les protocoles associés à chaque fichier
protocols = {}

# Parcourir tous les fichiers ABF et stocker les protocoles associés
for file in file_list:
    # Chemin complet du fichier ABF
    abf_file_path = os.path.join(folder_path, file)

    # Charger le fichier ABF
    abf = ABF(abf_file_path)

    # Accéder au protocole
    protocol = abf.protocol

    # Stocker le protocole associé au fichier
    protocols[file] = protocol

# Créer un dossier pour chaque protocole différent
for protocol in set(protocols.values()):
    # Chemin du nouveau dossier pour le protocole
    protocol_folder_path = os.path.join(folder_path, protocol)

    # Vérifier si le dossier existe déjà, sinon le créer
    if not os.path.exists(protocol_folder_path):
        os.makedirs(protocol_folder_path)

# Déplacer les fichiers ABF dans les dossiers correspondants aux protocoles
for file, protocol in protocols.items():
    # Chemin complet du fichier ABF
    abf_file_path = os.path.join(folder_path, file)

    # Chemin du dossier pour le protocole
    protocol_folder_path = os.path.join(folder_path, protocol)

    # Déplacer le fichier ABF dans le dossier correspondant au protocole
    shutil.move(abf_file_path, protocol_folder_path)

print("Tous les fichiers ont été déplacés dans les dossiers correspondants aux protocoles.")

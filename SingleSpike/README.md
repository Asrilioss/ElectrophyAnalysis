Spike Feature Extraction Script

Author : Elias [@YourLab]
Affiliation : [Your Institution]
Version : 1.2.0
Date : 2025-06-02

Overview

This Python script automates the extraction of spike features from Axon Binary Files (ABF, electrophysiology recordings), using the IPFX library.

Classifies recordings by experimental group (date + file index)

Extracts key spike features for each sweep, using robust feature detection

Saves results to Excel, by group, with detailed metadata

Generates quality control plots (histograms, boxplots) for each group

Logs all errors and info for transparent troubleshooting

Traceable outputs : All results, logs, and metadata are stored in a uniquely timestamped output folder for full reproducibility

Installation

We recommend using a virtual environment (conda or venv).

Clone this repository (or download the script)

Change directory to the project folder

Create and activate a new environment (for example with venv)

Install dependencies with pip install -r requirements.txt

Required packages : numpy, pandas, matplotlib, pyabf, ipfx

Note : When you install ipfx, all secondary dependencies (such as h5py, pynwb, etc.) are installed automatically via pip.

Usage

Run from the command line :

python SingleSpikechat.py --data_dir "<path_to_abf_files>" --output "<output_excel_file.xlsx>"

Arguments :
--data_dir : Path to the folder containing your ABF files
--output : Name for the Excel output file (will be placed inside the results folder)

Example :

python SingleSpikechat.py --data_dir "E:/Electrophy/Data culture/Data culture 4 PicroCP +G/Spike CC/" --output "SpikeFeatureResults.xlsx"

All results and plots will be saved in a new folder named like results_YYYY-MM-DD_HH-MM-SS.

Outputs

Excel file : One sheet per group, with all extracted features

Summary sheet : Quick stats per group (number of sweeps, means)

Metadata sheet : Script version, parameters, software versions, run timestamp

QC plots : Group-wise histograms and boxplots (PNG images)

Logs : Detailed info and error logs for traceability

Features Extracted

Spike width (ms)

Amplitude (mV)

Threshold (mV)

Duration (ms)

“Durée_PA” (see code for definition)

Epoch level, file info

For a full explanation of features, see the code (SingleSpikechat.py) or contact the author.

Reproducibility

All outputs, logs, and metadata are saved in a timestamped folder.

Optionally, the script and README are copied to the results folder for full traceability.

The script is modular and easily adaptable for new group assignments or feature sets.


How to add or modify experimental groups

This script assigns each ABF file to an experimental group using the mapping table in group_mapping.csv.
If you want to adapt the script for your own experiments, simply edit this file.

Open group_mapping.csv

This CSV file is located in the main project folder. Each row defines one group.
The file has four columns:

date: The date of the experiment as it appears in your ABF filenames (format: dd_mm_yy, for example 19_9_23)

file_number_start: The first file number in the range (inclusive)

file_number_end: The last file number in the range (inclusive)

group: The group name to assign

Add a new group

To add a group, add a new line to the CSV file.
For example, to assign files with date 01_7_25 and file numbers from 10 to 20 to a new group called MyNewGroup, add:

01_7_25,10,20,MyNewGroup

Modify an existing group

To change an existing group, find the corresponding line and update the file number range or group name as needed.
For example, to rename group NT to NewName for date 19_9_23 and files 0–39, change the line to:

19_9_23,0,39,NewName

Save the file

After making your changes, save group_mapping.csv.
You do not need to edit the Python script. The next time you run the analysis, your new groups will be used automatically.

Tips

Each row must correspond to a unique combination of date and file number range.

If a file does not match any line, it will be assigned to the group “Other”.

You can add as many groups as you need.

To cover all files for a date, use a wide range, for example: 0,9999,AllSamples

Example of a complete mapping file:

date,file_number_start,file_number_end,group
19_9_23,0,39,NT
19_9_23,40,70,C99
20_9_23,0,29,Poly A




How to cite

If you use this script in your research, please cite :

Elias, et al. (2025). Spike Feature Extraction Script (v1.2.0). [Software]. https://github.com/YourLab/spike-feature-extractor

License

[Specify your license here: MIT, GPL, etc.]

Questions or feedback ?

Open an issue on GitHub or contact Elias at [your.email@yourinstitution.edu].
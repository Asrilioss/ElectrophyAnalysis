==============================
ABF AUTOMATED ANALYSIS PIPELINE - USER GUIDE
==============================

This Python script enables **automated, reproducible, and parameterizable analysis of ABF electrophysiology files**: spike detection, statistics, group averages, and full quality control.

----------------------------------------------------
1. PREREQUISITES
----------------------------------------------------
- **Python 3.10
- Required packages: `pyabf`, `numpy`, `pandas`, `matplotlib`, `efel`
  - You can install them with:  
    `pip install pyabf numpy pandas matplotlib efel`
- The script file: `500mschatlast.py`
- Your `.abf` files ready in a folder

----------------------------------------------------
2. LAUNCHING THE SCRIPT (Command Line)
----------------------------------------------------
- Open a terminal:
  - **Windows**: Start > “Command Prompt” (or type “cmd” and hit Enter)
  - **Mac/Linux**: Terminal app
- Navigate to the folder containing your script (using `cd ...`)
- Basic usage example:
  
      python 500mschatlast.py --folder "E:/YourFolder/ABF" --export_name "Test2025"

----------------------------------------------------
3. OPTIONS & PARAMETERS
----------------------------------------------------
| Option         | Required? | Default           | Description                                             |
|:---------------|:---------:|:------------------|:--------------------------------------------------------|
| --folder       | Yes       | -                 | Path to the folder containing `.abf` files              |
| --export_name  | No        | Depol500_CC       | Name of the export/results folder                       |
| --stim_start   | No        | 499               | Start of analysis window (ms)                           |
| --stim_end     | No        | 1040              | End of analysis window (ms)                             |
| --features     | No        | Spikecount adaptation_index time_to_first_spike ISI_CV ISI_values | List of EFEL features to extract, space-separated |
| --retry_errors | No        | off               | (Advanced) Rerun only on files with previous errors     |

**See all options and their descriptions with:**
    
    python Depol500final_chat.py --help

**Example run with all parameters:**

    python Depol500final_chat.py --folder "E:/Data/ABF" --export_name "MyExperiment" --stim_start 500 --stim_end 1100 --features Spikecount ISI_CV

----------------------------------------------------
4. RESULTS & OUTPUTS
----------------------------------------------------
After each run, a results folder is created, e.g.:
  
    export_<export_name>_<date-time>/

Inside, you’ll find:
- **Excel and CSV files**: all sweep results, and group means
- **README_ANALYSE.txt**: all run parameters and context (in `runinfo/`)
- **QC_report.txt**: detailed QC/summary report (in `runinfo/`)
- **pipeline.log**: complete log of all steps and issues (in `runinfo/`)
- **ABF_errors.csv**: list of any sweeps with errors (in `runinfo/`)
- **Stats_groupes/**: group-level stats and plots
- **Graphs**: global and per-group quality plots

**All “meta” and logs are organized in a `runinfo/` subfolder for clarity.**

----------------------------------------------------
5. REPRODUCIBILITY & TRACEABILITY
----------------------------------------------------
- Each analysis is timestamped and saved in its own folder
- All settings, logs, and even the script version are archived with each run
- You can always revisit or share the exact context of any result

----------------------------------------------------
6. TROUBLESHOOTING & CONTACT
----------------------------------------------------
- For error messages or issues, check the `pipeline.log` and `ABF_errors.csv` in the runinfo folder.
- If you encounter missing files, permissions errors, or strange results, check the README and logs first.
- For help, contact the script maintainer or your lab’s e-phys expert.

----------------------------------------------------
7. BEST PRACTICES
----------------------------------------------------
- Change `--export_name` for each experiment/run to avoid overwriting and keep your workflow organized.
- Keep your data and export folders well-structured:  
  `E:/ABF_Data/` and `E:/ABF_Exports/`
- For collaborative projects, consider using version control (Git) and archiving your requirements.txt

----------------------------------------------------
8. CREDITS
----------------------------------------------------
Developed by [Your Name/Lab]  
Version: [update date/version]  
Feel free to adapt/extend this pipeline for your own research!

----------------------------------------------------

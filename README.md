PyADM1ODE-based Pohang Food Waste Anaerobic Digestion Simulator
This project models and simulates the anaerobic digestion process of food waste in Pohang City using the PyADM1ODE solver.

The pohang_adm1_colab_v6_en_2.py code included in this repository is the English edition of the v6 logic. It retains identical simulation logic but all comments, messages, and labels have been translated to English for broader accessibility.
Key Features
High-Resolution 1-Second Time Series
The hourly ODE results calculated by PyADM1ODE are linearly interpolated to generate 1-second resolution time series data.
Data is saved monthly in parquet/zstd format for optimal performance (automatically falls back to csv.gz if the pyarrow package is missing).
Seamless Google Colab Integration & Automated Pathing
Designed to run seamlessly in Google Colab without manual local environment setup.
The google.colab.drive mounting feature is built into the script. The code automatically reads inputs and writes outputs to your Colab local disk (/content/ADM1_QRA_results) and backs them up to your mounted Google Drive (/content/drive/MyDrive/ADM1_QRA_results).
Dynamic Feedstock Input via Excel (Automated YAML Generation)
No need to write complex YAML files by hand.
The script parses the /content/pohang_fw_yaml_monthly.xlsx file. It reads 17 variables across 12 months from the sheet exactly named 표2_YAML변수 and automatically generates 12 monthly feedstock files (pohang_waste_MM.yaml).
H2S & NH3 Concentration Estimation
While ADM1 lacks a sulphur cycle, this code estimates the upper-bound equilibrium concentrations of Hydrogen Sulphide (H2S) and Ammonia (NH3) in the off-gas (ppmv). This is calculated using the model's protein degradation extent, pH/temperature partitioning (dissociation), and Henry's law. Dedicated plots for these concentrations are automatically generated.
Repository Structure
The project folder structure is organized as follows:

.

├── Hourly meteorological observ.../ # Meteorological Dataset

├── src/                             # Simulation source code (including pohang_adm1_colab_v6_en_2.py)

├── LICENSE                          # License information (Initial commit)

├── README.md                        # Project documentation (this file)

└── Requirements.txt                 # List of package dependencies
Getting Started
1. Environment Setup
This code is highly recommended to be run in a Google Colab environment.

Open the Python script (pohang_adm1_colab_v6_en_2.py) located in the src/ folder using Google Colab.
When prompted during execution, grant permission to mount your Google Drive.
2. Prepare the Data File
Upload the required Excel file using the Colab file panel on the left, or place it in the root of your mounted Google Drive.

Required File: pohang_fw_yaml_monthly.xlsx
Target Sheet Name: 표2_YAML변수 (This exact sheet name must exist in the workbook).
3. Run the Simulation
You can control the simulation duration using the RUN_MODE setting at the top of the script:

RUN_MODE = "test" : Runs a short 2-day simulation (recommended for initial setup and debugging).
RUN_MODE = "full" : Runs the complete 12-month simulation.
Assumptions & Variables
Digester Volume: Working volume 6000 m³, Headspace 800 m³
Temperature: 35 °C (308.15 K)
Feed Flow (Q_FEED): 200 m³/d (YAML concentrations represent 1:9 diluted values)
H2S Estimation Parameters:
F_S_PROTEIN: 0.010 (Sulphur mass fraction of protein)
SO4_IN_mgS_L: 0.0 (Sulphate in the feed slurry)
FE_DOSE_mol_m3: 0.0 (Iron salt dose)
References
Feedstock Data: Park, Sang-Hyuk et al. (2022). Monthly characteristics of food waste in Pohang City (Raw material, 2020.11~2021.10). The simulation inputs for the monthly biochemical compositions and parameters used in the pohang_fw_yaml_monthly.xlsx workbook are digitized and derived from this study.

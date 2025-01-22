# Comparative analysis of Im7 dynamics with polarizable and non-polarizable CHARMM family of force fields
This repository contains input files for simulation of Im7 protein with CHARMM36m and DRUDE2019 force field along with the analysis scripts for comparing the protein dynamics, salt bridge dynamics and cation-protein interactions between the CHARMM36m and DRUDE2019 force fields. The trajectory files used in the analysis can be found at https://doi.org/10.5281/zenodo.14715013

## Contents
- **`simulation_setup`**: Starter files to setup the simulation with CHARMM36m and DRUDE2019 force field.
- **`codes`**: important python codes used in the paper which are not directly available in any MD simulation analysis package.
- **`plots`**: Contains all figures from the main text.
- **`plots_main_figure.ipynb`**: Jupyter notebook containing codes for plotting of the figures in main text.

## Requirements

We used the following Python packages for the analysis:

- `mdanalysis` (version 2.7.0): For trajectory analysis.
- `mdtraj` (version 1.9.9): For dssp secondary structure prediction.
- `openmm` (version 7.7.0): For simulations with DRUDE force field.
- `numpy` (version 1.23.5): For array operations and data handling.
- `scipy` (version 1.9.3): statistical analysis.
- `matplotlib` (version 3.8.4): For data visualization and plotting.
- `pandas` (version 2.2.1): For data handling.
- `seaborn` (version 0.12.2): For data visualization and plotting.
- `tqdm` (version 4.66.4): Progress bar.

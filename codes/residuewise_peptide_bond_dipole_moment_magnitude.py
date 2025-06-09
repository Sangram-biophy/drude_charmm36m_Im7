import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from tqdm import tqdm
from MDAnalysis.lib.log import ProgressBar

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='psffile', help='Input psf file', required=True)
parser.add_argument('-f', dest='traj', help='concatenated trajectory', required=True)
parser.add_argument('-fftype', dest='type', help='CHARMM/DRUDE', required=True)
parser.add_argument('-skip', dest='skip', help='CHARMM/DRUDE', required=True,type=int)
args = parser.parse_args()

top = args.psffile
traj = args.traj
ff=args.type
skip=args.skip
eA2D=4.803205467200591

u=mda.Universe(top,traj)
protein=u.select_atoms("protein")
# Precompute selections for all residues
if ff=="CHARMM":
    selections = [
        u.select_atoms(f"protein and (resid {i} and name C O H NH CA HA)")
        for i in range(np.unique(protein.resids)[1], np.unique(protein.resids)[-1])
    ]
elif ff=="DRUDE":
    selections = [
        u.select_atoms(f"protein and (resid {i} and name LPOA LPOB C DC O DO N DN HN DCA CA HA)",updating=True)
        for i in range(np.unique(protein.resids)[1], np.unique(protein.resids)[-1])
    ]
    
# Preallocate the array

n_frames = u.trajectory.n_frames//skip
dipole_moment_peptide = np.zeros((n_frames, protein.n_residues - 1))  # Store x, y, z dipole components

for frame_idx, ts in enumerate(ProgressBar(u.trajectory[::skip])):
    for i, peptide_bond in enumerate(selections):
        dipole_moment_peptide[frame_idx, i] = peptide_bond.dipole_moment() * eA2D

if ff=="CHARMM":
    np.save("residuewise_peptide_bond_dipole_moment_time_evol.npy",dipole_moment_peptide,allow_pickle=True)
if ff=="DRUDE":
    np.save(f"residuewise_DRUDE_LP_include_peptide_bond_dipole_moment_time_evol_skip{skip}.npy",dipole_moment_peptide,allow_pickle=True)

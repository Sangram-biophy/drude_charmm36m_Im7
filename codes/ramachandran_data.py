import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import argparse
from MDAnalysis.analysis.dihedrals import Ramachandran

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='psffile', help='Input psf file', required=True)
parser.add_argument('-f', dest='traj', help='concatenated trajectory', required=True)
parser.add_argument('-res', dest='res', help='protein residues', default=None)
args = parser.parse_args()

refpsf = args.psffile
traj = args.traj
prot_residue = args.res

u=mda.Universe(refpsf,traj)
atom_selection = u.select_atoms(f"protein and resid {prot_residue}")
R = Ramachandran(atom_selection,check_protein=True).run(verbose=True)
all_frame_phi_psi=np.concatenate(R.results.angles)
np.save(f"all_frame_phi_psi_resid_{prot_residue}.npy",all_frame_phi_psi)

fig, ax = plt.subplots(figsize=(5,3),dpi=1000,layout="tight")
R.plot(ax=ax, color='k',s=0.2,ref=True,marker=".",alpha=0.5)
ax.set_aspect('equal')
fig.savefig(f"ramachandran_residue{prot_residue}.png",dpi=1000)


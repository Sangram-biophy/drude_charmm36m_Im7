import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import argparse
from MDAnalysis.analysis.dihedrals import Dihedral

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='psffile', help='Input psf file', required=True)
parser.add_argument('-r', dest='refcrd', help='Reference coordinates file', required=True)
parser.add_argument('-f', dest='traj', help='concatenated trajectory', required=True)
parser.add_argument('-res', dest='res', help='protein residues', default=None)
args = parser.parse_args()

refpsf = args.psffile
refcrd = args.refcrd
traj = args.traj
prot_residue = args.res

u=mda.Universe(refpsf,traj)
ref=mda.Universe(refpsf,refcrd)

def dihedral_calc(u, ref, residues=None):
    if residues is None:
        selection = u.select_atoms("protein")
        selection_ref = ref.select_atoms("protein")
    else:
        selection = u.select_atoms(f"protein and resid {residues}")
        selection_ref = ref.select_atoms(f"protein and resid {residues}")
    
    phi = [res.phi_selection() for res in selection.residues if res.phi_selection() is not None]
    psi = [res.psi_selection() for res in selection.residues if res.psi_selection() is not None]
    
    phi_ref = [res.phi_selection() for res in selection_ref.residues if res.phi_selection() is not None]
    psi_ref = [res.psi_selection() for res in selection_ref.residues if res.psi_selection() is not None]
    
    phi_val = Dihedral(phi).run(verbose=True, step=1)
    psi_val = Dihedral(psi).run(verbose=True, step=1)
    
    phi_val_ref = Dihedral(phi_ref).run(verbose=True, step=1)
    psi_val_ref = Dihedral(psi_ref).run(verbose=True, step=1)
    
    return phi_val, psi_val, phi_val_ref, psi_val_ref

def dihedral_rmsd(phi_val,psi_val,phi_val_ref,psi_val_ref):
    #diff_phi = np.abs(phi_val.results["angles"]-phi_val_ref.results["angles"])
    #diff_phi = np.where(diff_phi > 180, 360 - diff_phi, diff_phi)
    #np.savetxt(f"phi_diff_resid{}.txt",diff_phi)
    diff_phi = phi_val.results["angles"]-phi_val_ref.results["angles"]
    diff_phi = np.where(diff_phi > 180, 360 - diff_phi, np.where(diff_phi<-180,-360-diff_phi,diff_phi))
    diff_phi_sq=np.square(diff_phi)
    
    # diff_psi = np.abs(psi_val.results["angles"]-psi_val_ref.results["angles"])
    # diff_psi = np.where(diff_psi > 180, 360 - diff_psi, diff_psi)
    #np.savetxt(f"psi_diff_resid{}.txt",diff_psi)
    diff_psi = psi_val.results["angles"]-psi_val_ref.results["angles"]
    diff_psi = np.where(diff_psi > 180, 360 - diff_psi, np.where(diff_psi<-180,-360-diff_psi,diff_psi))
    diff_psi_sq=np.square(diff_psi)
    
    
    RMSD_phi=np.sqrt(np.mean(diff_phi_sq,axis=1))
    RMSD_psi=np.sqrt(np.mean(diff_psi_sq,axis=1))
    
    t_phi=np.arange(0,len(RMSD_phi))*0.02/1e3 #in microsec
    t_psi=np.arange(0,len(RMSD_psi))*0.02/1e3
    
    t_RMSD_phi=np.column_stack((t_phi,RMSD_phi))
    t_RMSD_psi=np.column_stack((t_psi,RMSD_psi))
    
    return t_RMSD_phi,t_RMSD_psi

phi_val,psi_val,phi_val_ref,psi_val_ref=dihedral_calc(u,ref,prot_residue)
RMSD_phi,RMSD_psi=dihedral_rmsd(phi_val,psi_val,phi_val_ref,psi_val_ref)

if prot_residue is None:
    np.save(f"rmsd_phi_prot_all.npy",RMSD_phi)
    np.save(f"rmsd_psi_prot_all.npy",RMSD_psi)
else:
    np.save(f"rmsd_phi_prot_resid{prot_residue}.npy",RMSD_phi)
    np.save(f"rmsd_psi_prot_resid{prot_residue}.npy",RMSD_psi)

import os
import numpy as np
import argparse
import MDAnalysis as mda
import pandas as pd
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='refpsffile', help='Reference psf file', default="step3_charmm2omm.psf",type=str)
parser.add_argument('-r', dest='refcrd', help='Reference coordinates file', default="step3_charmm2omm.crd",type=str)
parser.add_argument('-o', dest='out', help='Output data file', default="iRED_S2.dat",type=str)
parser.add_argument('-xtc', dest='xtc', help='trajectory in xtc format', required=True,type=str)
parser.add_argument('-dt', dest='dt', help='time step while saving', default=20,type=int)

args = parser.parse_args()

refpsf = args.refpsffile
refcrd = args.refcrd
traj_xtc=args.xtc
out = args.out
dt = args.dt

u=mda.Universe(refpsf,traj_xtc)
ref=mda.Universe(refpsf,refcrd)
stop_time=u.trajectory.n_frames*dt/1000

def S2_iRED(traj,sel1,sel2):
    M=np.zeros((sel1.n_atoms,sel1.n_atoms))
    for ts in traj:
        NHvec=(sel1.positions-sel2.positions)/(np.linalg.norm(sel1.positions-sel2.positions,axis=1)[:,None])
        corrmatt=np.dot(NHvec,NHvec.T)
        M+=(3*np.square(corrmatt)-1)
    M=0.5*(M/len(traj))
    eigval,eigvec=np.linalg.eigh(M)   #linalg.eigh is used as M matrix is real symmetric 
    S2=np.zeros(sel1.n_atoms)
    for i in range(len(eigval[:-5])):
        S2+=eigval[i]*np.square(eigvec[:,i])
    return (1-S2)

HN=u.select_atoms("name HN")
N=u.select_atoms("name N")
NHvecresid=[]
for i in HN.resids:
    for j in N.resids:
        if i==j:
            NHvecresid.append(i)

NHbondlist=" ".join(str(x) for x in NHvecresid)
HN_sel=u.select_atoms(f"resid {NHbondlist} and name HN")
N_sel=u.select_atoms(f"resid {NHbondlist} and name N")

Tired=[1,5,10,25,50,100,250,500,1000]  #it is the block size in ns
avg_S2_iRED=np.zeros((len(NHvecresid),len(Tired)+1))
avg_S2_iRED[:,0]=np.array(NHvecresid)
for n,tblock in enumerate(tqdm(Tired)):
    blocks=np.arange(start=0,stop=stop_time,step=tblock)*1000/dt #converting ns to ps and then to frames
    avg_S2=np.zeros(len(NHvecresid))
    for j in range(len(blocks)-1):
        traj=u.trajectory[int(blocks[j]):int(blocks[j+1]):]
        avg_S2+=S2_iRED(traj,sel1=HN_sel,sel2=N_sel)
    avg_S2=avg_S2/(len(blocks)-1)
    avg_S2_iRED[:,n+1]=avg_S2
delimiter="        "
header_2=delimiter.join(str(x)+"ns" for x in Tired)
np.savetxt(f"{out}",avg_S2_iRED,delimiter="\t",header="residue"+"\t"+header_2)
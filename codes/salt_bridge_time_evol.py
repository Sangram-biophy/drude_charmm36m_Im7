import os
import numpy as np
import argparse
import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis import contacts
from MDAnalysis.lib.log import ProgressBar

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='refpsffile', help='Reference psf file', default="step3_charmm2omm.psf",type=str)
parser.add_argument('-o', dest='out', help='Output data file', default="framewise_salt_bridge_pair.npy",type=str)
parser.add_argument('-xtc', dest='xtc', help='trajectory in xtc format', required=True,type=str)
parser.add_argument('-dt', dest='dt', help='time step while saving', default=20,type=int)
parser.add_argument('-basic', dest='sel_basic', help='basic group selection', default="(resname ARG LYS HIS) and (name NH* NZ)",type=str)
parser.add_argument('-acidic', dest='sel_acidic', help='acidic group selection', default="(resname ASP GLU) and (name OE* OD*)",type=str)

args = parser.parse_args()

refpsf = args.refpsffile
traj_xtc=args.xtc
out = args.out
dt = args.dt
sel_basic = args.sel_basic
sel_acidic = args.sel_acidic

u=mda.Universe(refpsf,traj_xtc)
acidic = u.select_atoms(sel_acidic)
basic = u.select_atoms(sel_basic)

def distinct_contact_pairs(group_a,group_b,dimension,cutoff):
    dist=contacts.distance_array(group_a.positions, group_b.positions,box=dimension)
    indices=np.array(np.nonzero(dist<4.0))
    contact_pairs=set()
    for i in range(len(indices.T)):
        contact_pairs.add((group_a.resids[indices[:,i][0]],group_b.resids[indices[:,i][1]]))
    return contact_pairs

contact_pair_timeseries=[]
count_contact_timeseries=[]
for ts in ProgressBar(u.trajectory):
    box_d=ts.dimensions
    contact_pairs_framewise=distinct_contact_pairs(group_a=acidic,group_b=basic,dimension=box_d,cutoff=4.0)
    contact_pair_timeseries.append([ts.frame*dt/1e6,contact_pairs_framewise])
    count_contact_timeseries.append([ts.frame*dt/1e6,len(contact_pairs_framewise)])

np.save(f"{out}",np.array(contact_pair_timeseries))
np.save(f"count_{out}",np.array(count_contact_timeseries))

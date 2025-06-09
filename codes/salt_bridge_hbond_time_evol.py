import os
import numpy as np
import argparse
import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis import contacts
from MDAnalysis.lib.log import ProgressBar
from collections import defaultdict
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='refpsffile', help='Reference psf file', default="step3_charmm2omm.psf",type=str)
parser.add_argument('-o', dest='out', help='Output data file', default="framewise_salt_bridge_pair.npy",type=str)
parser.add_argument('-xtc', dest='xtc', help='trajectory in xtc format', required=True,type=str)
parser.add_argument('-dt', dest='dt', help='time step while saving', default=20,type=int)
parser.add_argument('-basic', dest='sel_basic', help='basic group selection', default="(resname ARG LYS HIS) and (name NH* NZ)",type=str)
parser.add_argument('-acidic', dest='sel_acidic', help='acidic group selection', default="(resname ASP GLU) and (name OE* OD*)",type=str)
parser.add_argument('-hsel', dest='sel_hyd', help='hydrogen group selection', default="(resname ARG LYS HIS) and (name HH* HZ*)",type=str)

args = parser.parse_args()

refpsf = args.refpsffile
traj_xtc=args.xtc
out = args.out
dt = args.dt
sel_basic = args.sel_basic
sel_acidic = args.sel_acidic
sel_hyd = args.sel_hyd

u=mda.Universe(refpsf,traj_xtc)
hbonds = HydrogenBondAnalysis(
    universe=u,
    donors_sel=sel_basic,
    hydrogens_sel=sel_hyd,
    acceptors_sel=sel_acidic,
    d_a_cutoff=4.0,
    d_h_a_angle_cutoff=150,
    update_selections=False
)

hbonds.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
)

hbonds1 = hbonds.results.hbonds

frame_to_unique_pairs = defaultdict(set)

for frame, donor_idx, hydrogen_idx, acceptor_idx, dist, angle in hbonds1:
    donor_res = u.atoms[int(donor_idx)].residue
    acceptor_res = u.atoms[int(acceptor_idx)].residue

    pair = tuple(sorted((donor_res.resid, acceptor_res.resid)))

    frame_to_unique_pairs[int(frame)].add(pair)

dt_ps = 20
num_frames = len(u.trajectory)
framewise_count = np.zeros((num_frames, 2))

for frame in tqdm(range(num_frames)):
    time = frame * dt_ps
    count = len(frame_to_unique_pairs.get(frame, []))
    framewise_count[frame] = [time, count]

np.save(f"hcriteria_{out}",framewise_count)

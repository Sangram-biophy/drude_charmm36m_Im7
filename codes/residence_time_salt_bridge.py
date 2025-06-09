import MDAnalysis as mda
from MDAnalysis.analysis import distances
from collections import Counter
import multiprocessing as mp
import numpy as np
import sys, os
from itertools import product
from tqdm import tqdm
import argparse

def product_of_groups(lp1,lp2):
    return [(x,y) for x,y in list(product(lp1,lp2)) if x!=y]

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='refpsffile', help='Reference psf file', default="step3_charmm2omm.psf",type=str)
#parser.add_argument('-o', dest='out', help='Output data file', default="sod_res_time",type=str)
parser.add_argument('-xtc', dest='xtc', help='trajectory in xtc format', required=True,type=str)
#parser.add_argument('-sel', dest='sel', help='which part of protein selection', default="protein",type=str)
parser.add_argument('-cutoff', dest='cutoff', help='cutoff', default=4.0,type=float)

args = parser.parse_args()

refpsf = args.refpsffile#"../step3_input.psf"#args.refpsffile
traj = args.xtc #"../prod_fit.xtc"#args.xtc
#out = args.out
#sel = args.sel
cutoff = args.cutoff

u=mda.Universe(refpsf,traj)

basic_resid=u.select_atoms("(resname ARG LYS HIS) and (name CA)")
acid_resid=u.select_atoms("(resname ASP GLU) and (name CA)")

allresiduepairlist = list(product_of_groups(acid_resid.resids, basic_resid.resids))
residence_time_counters = Counter(allresiduepairlist)
for pair in allresiduepairlist:
    residence_time_counters[pair] = 0

basic_res=u.select_atoms("(resname ARG LYS HIS) and (name NH* NZ)")
acid_res=u.select_atoms("(resname ASP GLU) and (name OE* OD*)")

fstart=0
fend=u.trajectory.n_frames
frames=np.arange(fstart,fend)
#cutoff=5.0
residence_times = []
residence_times_pair = []
for frame in tqdm(frames):
    #print(f"current frame number: {frame}/{fend}")
    u.trajectory[frame]

    distarr = distances.distance_array(acid_res.positions, basic_res.positions)
    paircontactindices = np.where(distarr < cutoff)
    paircontactindices = np.array(paircontactindices).swapaxes(0,1)

    paircontacts = []
    for pairindices in paircontactindices:
        pair = (acid_res.resids[pairindices[0]], basic_res.resids[pairindices[1]])
        if pair not in paircontacts:
            paircontacts.append((acid_res.resids[pairindices[0]], basic_res.resids[pairindices[1]]))

    if (frame > fstart) and (frame < fend):
        for pair in allresiduepairlist:
            if pair in paircontacts:
                residence_time_counters[pair] += 1

            elif (pair in prev_paircontacts) and (pair not in paircontacts):
                residence_times.append(residence_time_counters[pair])
                residence_times_pair.append(pair)
                residence_time_counters[pair] = 0

    if frame == fend:
        for pair in allresiduepairlist:
            if residence_time_counters[pair] > 0:
                residence_times.append(residence_time_counters[pair])
                residence_times_pair.append(pair)

    prev_paircontacts = paircontacts
    
np.savetxt(f'residence_time_salt_bridge.dat', residence_times, fmt='%i')
np.save(f'residence_time_salt_bridge_pairs.npy',np.array(residence_times_pair))


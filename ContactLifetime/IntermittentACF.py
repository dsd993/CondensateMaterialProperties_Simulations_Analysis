#Eg usage: python IntermittentACF.py 0 1 1000 10 ../../EKV_2.dat EKV_2

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 
import MDAnalysis.analysis as ana 
import MDAnalysis.analysis.distances
import sys 
np.set_printoptions(threshold=sys.maxsize)

start = int(sys.argv[1])
stride = int(sys.argv[2])
stop = int(sys.argv[3])
interval_time = int(sys.argv[4])
seq = open(str(sys.argv[5]),'r').read().strip()
fout = str(sys.argv[6])
array = np.load("{}_CLifetimeACFArray.txt.npy".format(fout))

nsteps = int((stop-start)/stride)

#Preparing a matrix to average over only E-K residue pair autocorrelation values
temp_seq_array_opppair = np.zeros((len(seq),len(seq)))
for idx, iseq in enumerate(seq):
    for jdx, jseq in enumerate(seq):
        if ((iseq == "E" and jseq == "K") or (iseq == "K" and jseq == "E")):
            temp_seq_array_opppair[idx,jdx]=1.0
npair_opppair = np.sum(temp_seq_array_opppair)
print("npair_opppair:", npair_opppair)

#Preparing a matrix to average over only E-E residue pair/K-K residue pair autocorrelation values
temp_seq_array_samepair = np.zeros((len(seq),len(seq)))
for idx, iseq in enumerate(seq):
    for jdx, jseq in enumerate(seq):
        if ((iseq == "E" and jseq == "E") or (iseq == "K" and jseq == "K")):
            temp_seq_array_samepair[idx,jdx]=1.0
npair_samepair = np.sum(temp_seq_array_samepair)
print("npair_samepair:", npair_samepair)

acf_opppair = np.zeros(nsteps)
for i in range(0,nsteps):
    acf_opppair[i] = np.sum(array[i,:,:]*temp_seq_array_opppair[:,:])/npair_opppair

acf_samepair = np.zeros(nsteps)
for i in range(0,nsteps):
    acf_samepair[i] = np.sum(array[i,:,:]*temp_seq_array_samepair[:,:])/npair_samepair

acf_combinedall = np.zeros(nsteps)
for i in range(0,nsteps):
    acf_combinedall[i] = np.sum(array[i,:,:])/(len(seq)*len(seq))

normacf_opppair = acf_opppair/acf_opppair[0]
normacf_samepair = acf_samepair/acf_samepair[0]
normacf_combinedall = acf_combinedall/acf_combinedall[0]

time_ps = np.arange(0,nsteps)*interval_time*1e-2 #1e-2 considering time units used is 10 fs
time_ns = np.arange(0,nsteps)*interval_time*1e-5 #1e-5 considering time units used is 10 fs

final_data = np.transpose(np.vstack((time_ps,time_ns,acf_opppair,normacf_opppair,acf_samepair,normacf_samepair,acf_combinedall,normacf_combinedall)))

np.savetxt('{}_ContactLifetimeACF.txt'.format(fout), final_data, delimiter=' ', header = "Time_ps Time_ns ACF_OppPair NormACF_OppPair ACF_SamePair NormACF_SamePair ACF_AllCombined NormACF_AllCombined")
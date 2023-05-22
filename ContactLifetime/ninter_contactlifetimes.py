#Eg usage: python ninter_contactlifetimes.py ../../veryshort_dumps.gsd ../../veryshort_dumps.gsd 0 1 1000 10 500 500 ../../EKV_2.dat EKV_2

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 
import MDAnalysis.analysis as ana 
import MDAnalysis.analysis.distances
import sys 
np.set_printoptions(threshold=sys.maxsize)

traj = sys.argv[1]
top = sys.argv[2]
start = int(sys.argv[3])
stride = int(sys.argv[4])
stop = int(sys.argv[5])
interval_time = int(sys.argv[6])
nchains = int(sys.argv[7])
nchains_ofinterest = int(sys.argv[8]) #No. of chains as reference I will use to compute rather than all chains in the system
seq = open(str(sys.argv[9]),'r').read().strip()
fout = str(sys.argv[10])

def CalcCutoffMatrix(seq1,seq2):
    cutoff = np.zeros((len(seq1),len(seq2)))
    for i in range(len(seq1)):
        if seq1==seq2:
            for j in range(i,len(seq1)):
                cutoff[i,j] = (siglist[AAlist.index(seq1[i])] + siglist[AAlist.index(seq1[j])])/2.0
                cutoff[j,i] = (siglist[AAlist.index(seq1[i])] + siglist[AAlist.index(seq1[j])])/2.0

        else:
            for j in range(len(seq2)):
                cutoff[i,j] = (siglist[AAlist.index(seq1[i])] + siglist[AAlist.index(seq2[j])])/2.0
    return cutoff 

univ = mda.Universe(top, traj)
global siglist, AAlist 
siglist = np.loadtxt('aminoacids_vdwdiameter.dat',usecols=1)
AAlist = list('ARNDCQEGHILKMFPSTWYV') #Must be of same order as that appears in the 'aminoacids_vdwdiameter.dat' file

natoms = int(len(univ.atoms))
nres = int(natoms/nchains)
nsteps = int((stop-start)/stride)
cutoff = CalcCutoffMatrix(seq,seq)*1.5

print("natoms:", natoms)
print("nchains:", nchains)
print("nres:", nres)
print("nsteps:", nsteps)
print("seq:", seq)
print("cutoff_shape:", np.shape(cutoff))

pair_count_array = np.zeros((nchains_ofinterest,nsteps,nres,nres))
chain_count_array = np.zeros((nsteps,nres,nres))
ncc = nchains 

for k in range(0,nchains_ofinterest):
    ncc=ncc-1 #Reducing the no. of pairs as we keep moving from chain 0 to chain 1 to chain 2 as reference chains
    npair=0
    count_array = np.zeros((ncc,nres,nres,nsteps))
    for l in range(k+1,nchains):
        print("pair_ids:", k,l)
        for ts in univ.trajectory[start:stop:stride]:
            moli = univ.atoms.positions[k*nres:(k+1)*nres,:]
            molj = univ.atoms.positions[l*nres:(l+1)*nres,:]
            dist = ana.distances.distance_array(moli, molj, box=univ.dimensions, backend="serial")
            count_array[npair,:,:,ts.frame] = (dist <= cutoff)

        npair += 1 #Counting the number of pairs for a given reference chain
    
    for framestride in range(0,nsteps):
        for npair_loop in range(0,npair):
            #Below, 0 added to framestride indicates only 0th frame is taken as reference, no moving avg. is done for a specific reason
            pair_count_array[k,framestride,:,:] += (np.where(count_array[npair_loop,:,:,0] > 0, count_array[npair_loop,:,:,0], -1) == count_array[npair_loop,:,:,0+framestride]).astype(int)
        
        if (npair > 0):   
            pair_count_array[k,framestride,:,:] = pair_count_array[k,framestride,:,:]/npair #Averaging over the no. of pairs for a given reference chain
        
        chain_count_array[framestride,:,:] += pair_count_array[k,framestride,:,:]

chain_count_array[:,:,:] = chain_count_array[:,:,:]/nchains_ofinterest #This array has for each residue the non-normalized (raw) autocorrelation function values after averaging over no. of pairs and no. of chains

np.save('{}_CLifetimeACFArray.txt'.format(fout), chain_count_array)

#Preparing a matrix to average over only E-K residue pair autocorrelation values
#temp_seq_array_opppair = np.zeros((len(seq),len(seq)))
#for idx, iseq in enumerate(seq):
    #for jdx, jseq in enumerate(seq):
        #if ((iseq == "E" and jseq == "K") or (iseq == "K" and jseq == "E")):
            #temp_seq_array_opppair[idx,jdx]=1.0
#npair_opppair = np.sum(temp_seq_array_opppair)
#print("npair_opppair:", npair_opppair)

#Preparing a matrix to average over only E-E residue pair/K-K residue pair autocorrelation values
#temp_seq_array_samepair = np.zeros((len(seq),len(seq)))
#for idx, iseq in enumerate(seq):
    #for jdx, jseq in enumerate(seq):
        #if ((iseq == "E" and jseq == "E") or (iseq == "K" and jseq == "K")):
            #temp_seq_array_samepair[idx,jdx]=1.0
#npair_samepair = np.sum(temp_seq_array_samepair)
#print("npair_samepair:", npair_samepair)

#acf_opppair = np.zeros(nsteps)
#for i in range(0,nsteps):
    #acf_opppair[i] = np.sum(chain_count_array[i,:,:]*temp_seq_array_opppair[:,:])/npair_opppair

#acf_samepair = np.zeros(nsteps)
#for i in range(0,nsteps):
    #acf_samepair[i] = np.sum(chain_count_array[i,:,:]*temp_seq_array_samepair[:,:])/npair_samepair

#acf_combinedall = np.zeros(nsteps)
#for i in range(0,nsteps):
    #acf_combinedall[i] = np.sum(chain_count_array[i,:,:])/(nres*nres)

#normacf_opppair = acf_opppair/acf_opppair[0]
#normacf_samepair = acf_samepair/acf_samepair[0]
#normacf_combinedall = acf_combinedall/acf_combinedall[0]

#time_ps = np.arange(0,nsteps)*interval_time*1e-2 #1e-2 considering time units used is 10 fs
#time_ns = np.arange(0,nsteps)*interval_time*1e-5 #1e-5 considering time units used is 10 fs

#final_data = np.transpose(np.vstack((time_ps,time_ns,acf_opppair,normacf_opppair,acf_samepair,normacf_samepair,acf_combinedall,normacf_combinedall)))

#np.savetxt('{}_ContactLifetimeACF.txt'.format(fout), final_data, delimiter=' ', header = "Time_ps Time_ns ACF_OppPair NormACF_OppPair ACF_SamePair NormACF_SamePair ACF_AllCombined NormACF_AllCombined")
import msprime as msp
import numpy as np
import random
import os
import sys
import scipy.special as sp
from numpy import log
from scipy.special import betaln
import argparse
from joblib import Parallel, delayed
import scipy.stats

#Sim parameters from Moorjani et al 2016
#Ne0 Neanderthal Ne 2500
#Ne1 Europe Ne 10000
#Ne2 East Asia Ne 10000
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#t1 split time_1 12000 gen
#t2 split time_2 2300 gen
#t3 split time 3 1500 gen
#f1 0.022 - original neanderthal pulse
#f2 0.01 - second pulse to east asia
#f3 0.01 - second pulse to europe
#f4 0.20 - dilution pulse to europe
#m1 f1 time 2000 gen
#m2 f2 time 1000 gen
#m3 f3 time 1000 gen
#m4 f4 time 1000 gen
#eu=european pop 0, as=asian pop 1, ba=basaleur pop 2, nean pop 3		

#TODO: Add Ne to the random parameters, figure out overlapping dates

def sim_pipeline(ID,m1,m2,m3,m4,t1,t2,t3,f1,f2,f3,f4,w,n):
	print(ID)
	outfile = open('outfile_sim%s.bed' %(ID), 'w+')
	outfile.close()

	N_admix = neanderthal_admixture_model(ID,seed=ID,mix_time1=m1,mix_time2=m2,mix_time3=m3,mix_time4=m4,split_time_1=t1,split_time_2=t2,split_time_3=t3,f1=f1,f2=f2,f3=f3,f4=f4,window_size =w,num_rep=n)

	#bedops
	B_ops = bedops(ID)
	
	
	#sys_stat
	S_stat = sys_stat(ID)

	#outfile reference and matrix
	O_file = ofile(ID,m1,m2,m3,m4,t1,t2,t3,f1,f2,f3,f4)
	
	#O_matrix = outmatrix(EU_AS)


def neanderthal_admixture_model(ID=1,seed=1,num_eu=170,num_as=394,num_nean = 1,anc_time=900,mix_time1=2000,mix_time2=1000,mix_time3=1000,mix_time4=1000,split_time_1=120000,split_time_2=2300,split_time_3=1500,f1=0.022,f2=0.00,f3=0.00,f4=0.20,Ne0=10000,Ne1=2500,Ne2=10000,mu=1.5e-8,window_size = 100000,num_SNP = 1,num_rep=1,coverage=False):
	for chr in range(1,23):
		infile = "/mnt/md0/villanea/MSprime/chr%s_map" %(chr)
		rho_map = msp.RecombinationMap.read_hapmap(infile)
		samples = [msp.Sample(population=0,time=0)]*num_eu
		samples.extend([msp.Sample(population=1,time=0)]*num_as) #no sampling of Basal Eurasian pop
		samples.extend([msp.Sample(population=3,time=anc_time)]*(num_nean)) #sample 1 Neanderthal for comparison
		pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1)]
		divergence = [msp.MassMigration(time=mix_time4,source=0,destination=2,proportion = f4), #BE dilution into EU
				msp.MassMigration(time=mix_time3,source=0,destination=3,proportion = f3), #second pulse EU
				msp.MassMigration(time=mix_time2,source=1,destination=3,proportion = f2), #second pulse AS
				msp.MassMigration(time=split_time_3,source=0,destination=1,proportion=1.0), #EU AS split
				msp.MassMigration(time=mix_time1,source=1,destination=3,proportion = f1), #first pulse
				msp.MassMigration(time=split_time_2,source=1,destination=2,proportion=1.0), #BE AS split
				msp.MassMigration(time=split_time_1,source=3,destination=2,proportion=1.0)] # Neand AS split
		divergence = sorted(divergence, key = lambda x: x.time)
		sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,recombination_map=rho_map,num_replicates=num_rep)
		chrom = "chr%s" %(chr)
		pos = []
		pos1 = []
		freq_EU = []
		freq_AS = []
		last = np.array(rho_map.get_positions()[-1])
		cur_sim = 0
		for sim in sims:
			cur_win = 1
			cur_start = 0
			cur_end = window_size-1
			cur_site = int(((cur_start+cur_end)+1)/2.0)
			cur_sim += 1
			for tree in sim.trees():
				F_int = tree.get_interval()
				while cur_site >= F_int[0] and cur_site < F_int[1]:
					cur_node = len(samples)-1  #the very last leaf, when adding more modern pops make sure Neanderthal is still last
					while tree.get_time(tree.get_parent(cur_node)) < split_time_1:
						cur_node = tree.get_parent(cur_node)
					N_freq_EU = 0
					N_freq_AS = 0
					for leaf in tree.leaves(cur_node):
						if tree.get_population(leaf)== 0:
							N_freq_EU += 1
						elif tree.get_population(leaf) == 1:
							N_freq_AS += 1
					pos.append(cur_site)
					pos1.append(cur_site+1)
					freq_EU.append(N_freq_EU)
					freq_AS.append(N_freq_AS)
					cur_start += window_size
					cur_end += window_size
					if cur_end > last:
						break
					cur_win += 1
					print cur_win
					cur_site = int(((cur_start+cur_end)+1)/2.0)
		outfile = open('outfile_sim%s.bed' %(ID), 'a')
		for line in range(0,len(freq_AS)):
			outfile.write(chrom)
			outfile.write('\t')
			outfile.write(str(pos[line]))
			outfile.write('\t')
			outfile.write(str(pos1[line]))
			outfile.write('\t')
			outfile.write(str(freq_EU[line]))
			outfile.write('\t')
			outfile.write(str(freq_AS[line]))
			outfile.write('\n')
		outfile.close()
	return np.array(pos), np.array(pos1), np.array(freq_EU), np.array(freq_AS)


def bedops(ID):
	os.system("sort-bed outfile_sim%s.bed > outfile_sim%s_sorted.bed" %(ID,ID))
	os.system("rm outfile_sim%s.bed" %(ID))
	os.system("bedops --element-of 1 outfile_sim%s_sorted.bed human_genome_mask_sorted.bed > outfile_sim%s_masked.bed" %(ID,ID))
	os.system("rm outfile_sim%s_sorted.bed" %(ID))

def sys_stat(ID):
		EU = np.genfromtxt('outfile_sim%s_masked.bed' %(ID), usecols=3)
		AS = np.genfromtxt('outfile_sim%s_masked.bed' %(ID), usecols=4)

		#delete sim file
		os.system("rm outfile_sim%s_masked.bed" %(ID))

		#initialize and fill the matrix
		EU_AS = np.zeros((171, 395)) #170+1, 394+1: +1 to include fixed alleles
		for i in range(0,len(AS)):
			EU_freq = EU[i]	
			AS_freq = AS[i]
			EU_AS[(EU_freq), (AS_freq)] = EU_AS[(EU_freq),(AS_freq)]+1
		np.savetxt('symmetry_matrix_%s.txt' %(ID), EU_AS, delimiter='\t')

def ofile(ID,m1,m2,m3,m4,t1,t2,t3,f1,f2,f3,f4):	
	outfile = open('symmetry_stat_%s.txt' %(ID), 'w+')
	outfile.write(str(ID))
	outfile.write('\t')
	outfile.write(str(t1))
	outfile.write('\t')
	outfile.write(str(t2))
	outfile.write('\t')
	outfile.write(str(t3))
	outfile.write('\t')
	outfile.write(str(f1))
	outfile.write('\t')
	outfile.write(str(f2))
	outfile.write('\t')
	outfile.write(str(f3))
	outfile.write('\t')
	outfile.write(str(f4))
	outfile.write('\t')
	outfile.write(str(m1))
	outfile.write('\t')
	outfile.write(str(m2))
	outfile.write('\t')
	outfile.write(str(m3))
	outfile.write('\t')
	outfile.write(str(m4))
	outfile.write('\t')


#Ne0 Neanderthal Ne 2500
#Ne1 Europe Ne 10000
#Ne2 East Asia Ne 10000
#t1 split time_1 12000 gen 
#t2 split time_2 2300 gen
#t3 split time 3 1500 gen
#f1 0.022 - original neanderthal pulse
#f2 0.01 - second pulse to east asia
#f3 0.01 - second pulse to europe
#f4 0.20 - dilution pulse to europe
#m1 f1 time 2000 gen
#m2 f2 time 1000 gen
#m3 f3 time 1000 gen
#m4 f4 time 1000 gen
#eu=european pop 0, as=asian pop 1, ba=basaleur pop 2, nean pop 3    

ID = np.random.randint(1,100000000,size=2)
m1 = scipy.stats.uniform.rvs(loc=1701, scale=2500, size=2)
m2 = scipy.stats.uniform.rvs(loc=500, scale=1300, size=2)
m3 = scipy.stats.uniform.rvs(loc=500, scale=1300, size=2)
m4 = scipy.stats.uniform.rvs(loc=500, scale=1300, size=2)
t1 = scipy.stats.uniform.rvs(loc=10000, scale=14000, size=2)
t2 = scipy.stats.uniform.rvs(loc=2100, scale=2501, size=2)
t3 = scipy.stats.uniform.rvs(loc=1301, scale=1700, size=2)
f1 = scipy.stats.uniform.rvs(loc=0, scale=0.1, size=2)
f2 = scipy.stats.uniform.rvs(loc=0, scale=0.1, size=2)
f3 = scipy.stats.uniform.rvs(loc=0, scale=0.1, size=2)
f4 = scipy.stats.uniform.rvs(loc=0, scale=0.5, size=2)

Sim = Parallel(n_jobs=2)(delayed(sim_pipeline)(ID[i],m1=m1[i],m2=m2[i],m3=m3[i],m4=m4[i],t1=t1[i],t2=t2[i],t3=t3[i],f1=f1[i],f2=f2[i],f3=f3[i],f4=f4[i],w=100000,n=1) for i in range(2))

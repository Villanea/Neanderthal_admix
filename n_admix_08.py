import msprime as msp
import numpy as np
import random

import sys

#One Model to rule them all, One Model to find them, One Model to bring them all and in the darkness bind them

#Sim parameters from Moorjani et al 2016
#Ne0 Neanderthal Ne 2500
#Ne1 Europe Ne 10000
#Ne2 East Asia Ne 10000
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#split time_1 12000 gen
#split time_2 2300 gen
#split time 3 1500 gen
#f1 0.02 - original neanderthal pulse
#f2 0.01 - second pulse to east asia
#f3 0.01 - second pulse to europe
#f4 0.20 - dilution pulse to europe
#f1 time 2000 gen
#f2 time 1000 gen
#f3 time 1000 gen
#f4 time 1000 gen
#eu=european pop 0, as=asian pop 1, ba=basaleur pop 2, neand pop 3		

#TODO: output in bed format

outfile = open('outfile_map_wholegen_dil.bed', 'w+')
#outfile.write("window\tfrequency_EU\tfrequency_AS")
#outfile.write('\n')
outfile.close()
def neanderthal_admixture_model(num_eu=170,num_as=394,num_nean = 1,anc_time=900,mix_time1=2000,mix_time2=1000,mix_time3=1000,mix_time4=1000,split_time_1=120000,split_time_2=2300,split_time_3=1500,f1=0.022,f2=0.00,f3=0.00,f4=0.25,Ne0=10000,Ne1=2500,Ne2=10000,mu=1.5e-8,window_size = 100000,num_SNP = 1,num_rep=1,coverage=False):
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
		sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,recombination_map=rho_map,num_replicates=num_rep)
		chrom = "chr%s" %(chr)
		#win = []
		pos = []
		pos1 = []
		freq_EU = []
		freq_AS = []
		last = np.array(rho_map.get_positions()[-1])
		#leng = []
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
					#F_length = tree.get_length()
					N_freq_EU = 0
					N_freq_AS = 0
					for leaf in tree.leaves(cur_node):
						if tree.get_population(leaf)== 0:
							N_freq_EU += 1
						elif tree.get_population(leaf) == 1:
							N_freq_AS += 1
					#win.append(cur_win)
					pos.append(cur_site)
					pos1.append(cur_site+1)
					freq_EU.append(N_freq_EU)
					freq_AS.append(N_freq_AS)
					#leng.append(F_length)
					cur_start += window_size
					cur_end += window_size
					print cur_end
					print last
					if cur_end > last:
						break
					cur_win += 1
					print cur_win
					cur_site = int(((cur_start+cur_end)+1)/2.0) #random.randint(cur_start,cur_end)
					print cur_site
		outfile = open('outfile_map_wholegen_dil.bed', 'a')
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
			#outfile.write(str(leng[line]))
			#outfile.write('\n')
		outfile.close()
	return np.array(pos), np.array(pos1), np.array(freq_EU), np.array(freq_AS)#, np.array(leng)

num_rep = 1
window_size = 100000
if len(sys.argv) > 1: 
	num_rep = int(sys.argv[1]) # take some command line arguments
if len(sys.argv) > 2:
	window_size = int(sys.argv[2]) # take some command line arguments
N_admix = neanderthal_admixture_model(window_size = window_size, num_rep = num_rep)

import msprime as msp
import numpy as np
import random

import sys

#Sim parameters from Moorjani et al 2016
#Modern humans Ne 10000, samples 10
#Neanderthal Ne 2500, samples 0
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#split time_1 12000 gen
#split time_2 1500 gen
#f1 0.02
#f2 0.01
#f1 time 2000 gen
#f2 time 1000 gen
#eu=european pop 0, as=asian pop 1		

		
def neanderthal_admixture_model(num_eu=100,num_as=100,anc_num = 1, anc_time=900,mix_time1=2000,mix_time2=1000,split_time_1=120000,split_time_2=1500,f1=0.02,f2=0.01,Ne0=10000,Ne1=2500,mu=1.5e-8,rho=1.0e-8,length=10000000,window_size = 1000000,num_SNP = 1,num_rep=100,coverage=False): #when is best time to sample Neanderthal? 100 gen before f?
	if f1 < 0 or f1 > 1 or f2 < 0 or f2 > 1: #error catching, leave there for now
		print "Admixture fraction is not in [0,1]"
		return None
	samples = [msp.Sample(population=0,time=0)]*num_eu
	samples.extend([msp.Sample(population=1,time=0)]*num_as)
	samples.extend([msp.Sample(population=anc_pop,time=anc_time)]*(anc_num)) #sample 1 Neanderthal for comparison
	pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1)]
	divergence = [msp.MassMigration(time=mix_time2,source=1,destination=2,proportion = f2), #second pulse
			msp.MassMigration(time=split_time_2,source=0,destination=1,proportion=1.0), #EU AS split
			msp.MassMigration(time=mix_time1,source=1,destination=2,proportion = f1), #first pulse
			msp.MassMigration(time=split_time_1,source=1,destination=2,proportion=1.0)] # Neand EU split. Do Neanderthals exchange into pop 1 then pop 1 into 2? or should Neanderthal exchange with both pops?
	sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,recombination_rate=rho,length=length,num_replicates=num_rep)
	print "done simulating"
	win = []
	freq_EU = []
	freq_AS = []
	leng = []
	cur_sim = 0
	for sim in sims:
		cur_win = 1
		cur_start = 0
		cur_end = window_size-1
		cur_site = (cur_start+cur_end)/2.0 #random.randint(cur_start,cur_end)
		cur_sim += 1
		print "current simulation"
		print cur_sim
		for tree in sim.trees():
			F_int = tree.get_interval()
			if cur_site >= F_int[0] and cur_site < F_int[1]:
				#raw_input()
				cur_node = len(samples)-1  #the very last leaf, when adding more modern pops make sure Neanderthal is still last
				while tree.get_time(tree.get_parent(cur_node)) < split_time:
					cur_node = tree.get_parent(cur_node)
				F_length = tree.get_length()
				N_freq_EU = 0
				N_freq_AS = 0
				for leaf in tree.leaves(cur_node):
					if tree.get_population(leaf)== 0:
						N_freq_EU += 1
					elif tree.get_population(leaf) == 1:
						N_freq_AS += 1
				win.append(cur_win)
				freq_EU.append(N_freq_EU)
				freq_AS.append(N_freq_AS)
				leng.append(F_length)
				cur_start += window_size
				cur_end += window_size
				if cur_end > length:
					break
				cur_win += 1
				print cur_win
				cur_site = (cur_start+cur_end)/2.0 #random.randint(cur_start,cur_end)
	outfile = open('outfile_ea.txt', 'w')
	outfile.write("window\tfrequency_EU\tfrequency_AS\tlength")
	outfile.write('\n')
	for line in range(0,len(leng)):
		outfile.write(str(win[line]))
		outfile.write('\t')
		outfile.write(str(freq_EU[line]))
		outfile.write('\t')
		outfile.write(str(freq_AS[line]))
		outfile.write('\t')
		outfile.write(str(leng[line]))
		outfile.write('\n')
	outfile.close()
	return np.array(win), np.array(freq_EU), np.array(freq_AS), np.array(leng)

num_rep = 1
window_size = 100000
if len(sys.argv) > 1: 
	num_rep = int(sys.argv[1]) # take some command line arguments
if len(sys.argv) > 2:
	window_size = int(sys.argv[2]) # take some command line arguments
N_admix = neanderthal_admixture_model(window_size = window_size, num_rep = num_rep)

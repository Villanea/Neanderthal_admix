import msprime as msp
import numpy as np
import random

import sys

#Sim parameters from Moorjani et al 2016
#Modern humans Ne 10000, samples 10
#Neanderthal Ne 2500, samples 0
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#split time 12000 gen
#f 0.03
#f time 100-2500 gen
		
#TODO: loop to automate between windows and replicates
		
def neanderthal_admixture_model(num_modern=1000,anc_pop = 1, anc_num = 1, anc_time=900,mix_time=2000,split_time=120000,f=0.03,Ne0=10000,Ne1=2500,mu=1.5e-8,rho=1.0e-8,length=10000000,window_size = 1000000,num_SNP = 1,num_rep=100,coverage=False):
	#when is best time to sample Neanderthal? 100 gen before f?
	#error catching, leave there for now
	if f < 0 or f > 1:
		print "Admixture fraction is not in [0,1]"
		return None
	samples = [msp.Sample(population=0,time=0)]*num_modern #sample 1 Neanderthal for comparison
	samples.extend([msp.Sample(population=anc_pop,time=anc_time)]*(anc_num))
	pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1)]
	divergence = [msp.MassMigration(time=mix_time,source=0,destination=1,proportion = f),
			msp.MassMigration(time=split_time,source=1,destination=0,proportion=1.0)]
	sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,recombination_rate=rho,length=length,num_replicates=num_rep)
	win = []
	freq = []
	leng = []
	#FYI mean fragment length from test_2 model ~6000 bp
	for sim in sims:
		cur_win = 1
		cur_start = 0
		cur_end = window_size-1
		cur_site = (cur_start+cur_end)/2.0 #random.randint(cur_start,cur_end)
		#print cur_start, cur_end, cur_site
		for tree in sim.trees():
			F_int = tree.get_interval()
			if cur_site >= F_int[0] and cur_site < F_int[1]:
				#print cur_site, F_int
				#raw_input()
				cur_node = len(samples)-1  #the very last leaf, when adding more modern pops make sure Neanderthal is still last
				while tree.get_time(tree.get_parent(cur_node)) < split_time:
					cur_node = tree.get_parent(cur_node)
				F_length = tree.get_length()
				N_freq = (tree.get_num_leaves(cur_node) - 1) #minus our lone Neanderthal
				win.append(cur_win)
				freq.append(N_freq)
				leng.append(F_length)
				cur_start += window_size
				cur_end += window_size
				if cur_end > length:
					break
				cur_win += 1
				cur_site = (cur_start+cur_end)/2.0 #random.randint(cur_start,cur_end)
				#print cur_start, cur_end, cur_site
	outfile = open('outfile_s.txt', 'w')
	outfile.write("window\tfrequency\tlength")
	outfile.write('\n')
	for line in range(0,len(leng)):
		outfile.write(str(win[line]))
		outfile.write('\t')
		outfile.write(str(freq[line]))
		outfile.write('\t')
		outfile.write(str(leng[line]))
		outfile.write('\n')
	outfile.close()
	return np.array(win), np.array(freq), np.array(leng)

num_rep = 1000
window_size = 100000
if len(sys.argv) > 1: 
	num_rep = int(sys.argv[1]) # take some command line arguments
if len(sys.argv) > 2:
	window_size = int(sys.argv[2]) # take some command line arguments
N_admix = neanderthal_admixture_model(window_size = window_size, num_rep = num_rep, num_modern = 100)

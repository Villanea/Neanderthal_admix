#!/usr/bin/env python

#test_2 is the simple model including recombination
#To do: concatenate "invisible" fragments - sequential fragments with the same frequency 
import msprime as msp
import numpy as np
import re
#based on Schraiber Admixture model: https://github.com/Schraiber/continuity/blob/master/ancient_genotypes.py
class FreqError(Exception):
	pass
#Sim parameters from Moorjani et al 2016
#Modern humans Ne 10000, samples 10
#Neanderthal Ne 2500, samples 0
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#split time 12000 gen
#f 0.03
#f time 100-2500 gen
		
def neanderthal_admixture_model(num_modern=1000,anc_pop = 1, anc_num = 1, anc_time=900,mix_time=1000,split_time=12000,f=0.03,Ne0=10000,Ne1=2500,mu=1.5e-8,rho=1.0e-8,length=10000000,num_rep=100,coverage=False):
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
	outfile = open('outfile.csv', 'w')
	outfile.write("frequency,length")
	outfile.write('\n')
	freq = [0]
	length = []
	sim_num = 0	
	for sim in sims:
		for tree in sim.trees():
			F_length = tree.get_length()
			cur_node = len(samples)-1  # the very last leaf, when adding more modern pops make sure Neanderthal is still last
			while tree.get_time(tree.get_parent(cur_node)) < split_time:
				cur_node = tree.get_parent(cur_node)
			N_freq = (tree.get_num_leaves(cur_node) - 1) #minus our lone Neanderthal
			if N_freq == freq[-1]: #hope the very first loop does not crash
				F_length = F_length+length[-1]
				length[-1] = F_length
			else:			
				length.append(F_length)
				freq.append(N_freq)
			outfile.write(str(N_freq))
			outfile.write(",")
			outfile.write(str(F_length))
			outfile.write("\n")
	del freq[0]
	outfile.close()
	return np.array(freq), np.array(length)

N_admix = neanderthal_admixture_model()

#!/usr/bin/env python

import msprime as msp
import numpy as np

#based on Schraiber Admixture model: https://github.com/Schraiber/continuity/blob/master/ancient_genotypes.py
#class FreqError(Exception):
#	pass
#Sim parameters from Moorjani et al 2016
#Modern humans Ne 10000, samples 10
#Neanderthal Ne 2500, samples 0
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#split time 12000 gen
#f 0.03
#f time 100-2500 gen

def neanderthal_admixture_model(num_modern=10,anc_pop = 1, anc_num = 1, anc_time=900,mix_time=1000,split_time=12000,f=0.03,Ne0=10000,Ne1=2500,mu=1.5e-8,length=1000,num_rep=1,coverage=False):
	#when is best time to sample Neanderthal? 100 gen before f?
	outFile = open('outfile.csv', 'w')
	outFile.write("Frequency")
	outFile.write('\n')
	#error catching, leave there for now
	if f < 0 or f > 1:
		print "Admixture fraction is not in [0,1]"
		return None
	samples = [msp.Sample(population=0,time=0)]*num_modern 
	#sample 1 neanderthal for comparison 
	samples.extend([msp.Sample(population=anc_pop,time=anc_time)]*(anc_num))
	pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1)]
	divergence = [msp.MassMigration(time=mix_time,source=0,destination=1,proportion = f),
			msp.MassMigration(time=split_time,source=1,destination=0,proportion=1.0)]
	sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,length=length,num_replicates=num_rep)
	print "msp is done simulating"
	freq = []
	sim_num = 0	
	for sim in sims:
	#what is the correct syntax to get a specific fragment? is it for variant in sim.variants?
		for variant in sim.variants():
			cur_node = get_samples(population_id=1)
			lenght = tree.get_branch_length(cur_node)
			#check that is returning a length in gens
			print lenght
			while lenght < mix_time:
				cur_node = tree.get_parent(cur_node)
				lenght = lenght + tree.get_branch_length(cur_node)
			N_freq = (get_num_leaves(cur_node) - 1)/num_modern
			freq.append(N_freq)
			#check that is printing a fraction
			print N_freq
			#write inside loop more inefficient, but ok for now
			outFile.write(N_freq)
			outFile.write('\n')
	#To do: We need mean frequency across all replicates
	return np.array(freq)
	outFile.close()
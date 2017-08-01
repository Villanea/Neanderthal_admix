import msprime as msp
import numpy as np
import random

import sys
import argparse

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
#f4 0.01 - dilution pulse to europe
#f1 time 2000 gen
#f2 time 1000 gen
#f3 time 1000 gen
#f4 time 1000 gen
#eu=european pop 0, as=asian pop 1, ba=basaleur pop 2, neand pop 3		

#TOD: simulate len chromosomes in genome
#TODO: fit admixture maps
#Do we need to record fragment length anymore?

parser = argparse.ArgumentParser("Simulate introgression with various parameters in a simple model of recombination")
parser.add_argument("-t1", default=12000, type = float, help = "split time of humans and Neandertals, in generations")
parser.add_argument("-t2", default=2300, type = float, help = "split time of basal Eurasian population")
parser.add_argument("-t3", default=1500, type = float, help = "split time of East Asians and Europeans")
parser.add_argument("-f1", default = 0.022, type = float, help = "introgression from Neandertals into ancestor of Europeans and Asiasn")
parser.add_argument("-f2", default = 0.003, type = float, help = "introgression from Neandertals into just East Asians")
parser.add_argument("-f3", default = 0.0, type = float, help = "introgression from Neandertals into just Europeans")
parser.add_argument("-f4", default = 0.0, type = float, help = "dilution from Basal Eurasians into Europeans")
parser.add_argument("-m1", default = 2000, type = float, help = "time of Neandertal to ancestor of European and Asian admxiture")
parser.add_argument("-m2", default = 1000, type = float, help = "time of admixture from Neandertal into East Asian")
parser.add_argument("-m3", default = 1000, type = float, help = "time of admixture from Neandertal into European")
parser.add_argument("-m4", default = 1000, type = float, help = "time of dilution from Basal Eurasian into European")
parser.add_argument("-w", default = 100000, type = int, help = "window size for pulling out admixed bases")
parser.add_argument("-n", default = 1000, type = int, help = "number of replicate simulations to run")
parser.add_argument("-l", default = 10000000, type = int, help = "length of fragment to simulate")

args = parser.parse_args() 

		
def neanderthal_admixture_model(num_eu=170,num_as=394,num_nean = 1,anc_time=900,mix_time1=2000,mix_time2=1000,mix_time3=1000,mix_time4=1000,split_time_1=120000,split_time_2=2300,split_time_3=1500,f1=0.022,f2=0.003,f3=0.00,f4=0.00,Ne0=10000,Ne1=2500,Ne2=10000,mu=1.5e-8,rho=1.0e-8,length=10000000,window_size = 100000,num_SNP = 1,num_rep=1000,coverage=False):
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
				while tree.get_time(tree.get_parent(cur_node)) < split_time_1:
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
				#print cur_win
				cur_site = (cur_start+cur_end)/2.0 #random.randint(cur_start,cur_end)
	outfile = open('outfile_all_2f.txt', 'w')
	outfile.write("window\tfrequency_EU\tfrequency_AS")
	outfile.write('\n')
	for line in range(0,len(leng)):
		outfile.write(str(win[line]))
		outfile.write('\t')
		outfile.write(str(freq_EU[line]))
		outfile.write('\t')
		outfile.write(str(freq_AS[line]))
		#outfile.write('\t')
		#outfile.write(str(leng[line]))
		outfile.write('\n')
	outfile.close()
	return np.array(win), np.array(freq_EU), np.array(freq_AS), np.array(leng)

N_admix = neanderthal_admixture_model(mix_time1=args.m1,mix_time2=args.m2,mix_time3=args.m3,mix_time4=args.m4,split_time_1=args.t1,split_time_2=args.t2,split_time_3=args.t3,f1=args.f1,f2=args.f2,f3=args.f3,f4=args.f4,length=args.l,window_size = args.w,num_rep=args.n)

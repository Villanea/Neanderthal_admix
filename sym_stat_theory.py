import numpy as np
from scipy.sparse.linalg import expm_multiply as expma
import scipy.special as sp
import scipy.optimize as opt

def generate_Q_het(n):
	Q = np.zeros((n+1,n+1))
	for k in range(n+1):
		Q[k,k] = -k*(n-k)
		if k > 0: Q[k,k-1] = .5*k*(k-1)
		if k < n: Q[k,k+1] = .5*(n-k-1)*(n-k)
	return Q

def generate_Qm_het(n):
	Qm = np.zeros((n+1,n+1))
	for k in range(n+1):
		if k > 0: Qm[k,k-1] = -k*(n-k)
		if k > 1: Qm[k,k-2] = .5*k*(k-1)
		Qm[k,k] = .5*(n-k-1)*(n-k)
	return Qm

def generate_het(x, n):
	pows = np.arange(0,n+1)
	het = x**pows*(1-x)**(n-pows)
	return het

def generate_het_one(f1, n):
	return generate_het(f1,n)

def generate_het_dilution(f1,f2,n):
	return generate_het((1-f2)*f1,n)
	
def generate_het_two(f1,f2,n):
	return generate_het(f2+(1-f2)*f1,n)

def one_pulse(f,t,n):
	h = generate_het_one(f,n)
	Q = generate_Q_het(n)
	sample = expma(Q*t,h)
	comb = sp.binom(n,np.arange(0,n+1))
	return comb*sample

def two_pulse(f1,t1,f2,t2,n):
	h = generate_het_two(f1,f2,n)
	Q = generate_Q_het(n)
	Qm = generate_Qm_het(n)
	sample = expma((Q-f2*Qm)*t1,h)
	sample = expma(Q*t2,sample)
	comb = sp.binom(n,np.arange(0,n+1))
	return comb*sample

def dilution(f1,t1,f2,t2,n):
	h = generate_het_dilution(f1,f2,n)
	Q = generate_Q_het(n)
	sample = expma(Q*((1-f2)*t1+t2),h)
	comb = sp.binom(n,np.arange(0,n+1))
	return comb*sample

def sym_stat_two_pulse(f1,f2,tp,ta1,ta2,te,n):
	t1AS = tp+ta1
	t2AS = ta2
	tEU = tp+te
	freqEU = one_pulse(f1,tEU,n)
	freqAS = two_pulse(f1,t1AS,f2,t2AS,n)
	sym = (freqEU-freqAS)/(freqEU+freqAS+1)
	print f1,f2,tp,ta1,ta2,te
	print sym
	return sym

def sym_stat_dilution(f1,f2,tp,ta,te1,te2,n):
	tAS = tp+ta
	t1EU = tp+te1
	t2EU = te2
	freqEU = dilution(f1,t1EU,f2,t2EU,n)
	freqAS = one_pulse(f1,tAS,n)	
	sym = (freqEU-freqAS)/(freqEU+freqAS+1)
	return sym

#NB: DON'T INCLUDE THE ZERO CATEGORY!!!!
def fit_two_pulse(dat):
	n = len(dat)
	return opt.curve_fit(lambda x, f1,f2,tp,ta1,ta2,te: sym_stat_two_pulse(f1,f2,tp,ta1,ta2,te,n)[1:], 1, dat,bounds=((0,0,0,0,0,0),(1,1,1,1,1,1)))
	

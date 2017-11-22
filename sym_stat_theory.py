import numpy as np
from scipy.sparse.linalg import expm_multiply as expma
import scipy.special as sp
import scipy.optimize as opt
import scipy.stats as st

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

def generate_Qsel_het(n):
	Qsel = np.zeros((n+1,n+1))
	for k in range(n+1):
		Qsel[k,k] = k
		if k < n: Qsel[k,k+1] = -(n-k)
	return Qsel

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

def beta_binom(k,n,alpha,beta, log=False):
	logP = sp.gammaln(n+1) - sp.gammaln(k+1) - sp.gammaln(n-k+1)
	logP += sp.gammaln(alpha + beta) - sp.gammaln(alpha) - sp.gammaln(beta)
	logP += sp.gammaln(k + alpha) + sp.gammaln(n - k + beta) - sp.gammaln(n + alpha + beta)
	if log:
		return logP
	else:
		return np.exp(logP)

def error_model(p,neg,pos,alpha,beta):
	n = len(p)-1
	res = np.zeros(len(p))
	res[0] = (1-pos)*p[0] + neg*(1-p[0])
	f = beta_binom(np.arange(n+1),n,alpha,beta)
	f = f[1:]/(1-f[0])
	res[1:] = (1-neg)*p[1:] + pos*p[0]*f	
	return res

def one_pulse(f,t,n, error = None):
	h = generate_het_one(f,n)
	Q = generate_Q_het(n)
	sample = expma(Q*t,h)
	comb = sp.binom(n,np.arange(0,n+1))
	sample *= comb
	if error is not None:
		neg, pos, alpha, beta = error
		sample = error_model(sample,neg,pos,alpha,beta)
	return sample

def two_pulse(f1,t1,f2,t2,n, error = None):
	h = generate_het_two(f1,f2,n)
	Q = generate_Q_het(n)
	Qm = generate_Qm_het(n)
	sample = expma((Q-f2*Qm)*t1,h)
	sample = expma(Q*t2,sample)
	comb = sp.binom(n,np.arange(0,n+1))
	sample *= comb
	if error is not None:
		neg, pos, alpha, beta = error
		sample = error_model(sample, neg, pos, alpha, beta)
	return sample

def dilution(f1,t1,f2,t2,n):
	h = generate_het_dilution(f1,f2,n)
	Q = generate_Q_het(n)
	sample = expma(Q*((1-f2)*t1+t2),h)
	comb = sp.binom(n,np.arange(0,n+1))
	return comb*sample

def one_pulse_sel(f,t,gamma,n,error = None):
	h = generate_het_one(f,n)
	Q = generate_Q_het(n)
	Qsel = generate_Qsel_het(n)
	sample = expma((Q+gamma*Qsel)*t,h)
	comb = sp.binom(n,np.arange(0,n+1))
	sample *= comb
	if error is not None:
		print error
		neg, pos, alpha, beta = error
		sample = error_model(sample,neg,pos,alpha,beta)
	return sample

def sym_stat_two_pulse(f1,f2,tp,ta1,ta2,te,n,norm=False):
	t1AS = tp+ta1
	t2AS = ta2
	tEU = tp+te
	freqEU = one_pulse(f1,tEU,n)
	freqAS = two_pulse(f1,t1AS,f2,t2AS,n)
	if norm:
		freqEU = freqEU[1:]/sum(freqEU[1:])
		freqAS = freqAS[1:]/sum(freqAS[1:])
	sym = (freqEU-freqAS)/(freqEU+freqAS+1)
	print f1,f2,tp,ta1,ta2,te
	return sym

def sym_stat_dilution(f1,f2,tp,ta,te1,te2,n,norm=False):
	tAS = tp+ta
	t1EU = tp+te1
	t2EU = te2
	freqEU = dilution(f1,t1EU,f2,t2EU,n)
	freqAS = one_pulse(f1,tAS,n)	
	if norm:
		freqEU = freqEU[1:]/sum(freqEU[1:])
		freqAS = freqAS[1:]/sum(freqAS[1:])
	sym = (freqEU-freqAS)/(freqEU+freqAS+1)
	print f1,f2,tp,ta,te1,te2
	return sym

def sym_stat_two_two(f1,fa,fe,tp,ta1,ta2,te1,te2,n,norm=False):
	t1AS = tp+ta1
	t2AS = ta2
	t1EU = tp+te1
	t2EU = te2
	freqEU = two_pulse(f1,t1EU,fe,t2EU,n)
	freqAS = two_pulse(f1,t1AS,fa,t2AS,n)
	if norm:
		freqEU = freqEU[1:]/sum(freqEU[1:])
		freqAS = freqAS[1:]/sum(freqAS[1:])
	sym = (freqEU-freqAS)/(freqEU+freqAS+1)
	print f1,fa,fe,tp,ta1,ta2,te1,te2
	return sym

#NB: DON'T INCLUDE THE ZERO CATEGORY!!!!
def fit_two_pulse(dat):
	n = len(dat)
	return opt.curve_fit(lambda x, f1,f2,tp,ta1,ta2,te: sym_stat_two_pulse(f1,f2,tp,ta1,ta2,te,n,True), 1, dat, p0 = st.uniform.rvs(size=6,scale=.5), bounds=((0,0,0,0,0,0),(1,1,.5,.5,.5,.5)))

def fit_dilution(dat):
	n = len(dat)
	return opt.curve_fit(lambda x, f1,f2,tp,ta,te1,te2: sym_stat_dilution(f1,f2,tp,ta,te1,te2,n,True), 1, dat, p0 = st.uniform.rvs(size=6,scale=.5), bounds=((0,0,0,0,0,0),(.5,.8,.5,.5,.5,.5)))#,ftol=1e-10,gtol=1e-10)

def fit_two_two(dat):
	n = len(dat)
	return opt.curve_fit(lambda x, f1,fa,fe,tp,ta1,ta2,te1,te2: sym_stat_two_two(f1,fa,fe,tp,ta1,ta2,te1,te2,n,True), 1, dat, p0 = st.uniform.rvs(size=8,scale=.5), bounds=((0,0,0,0,0,0,0,0),(.5,.5,.5,.5,.5,.5,.5,.5)))#,ftol=1e-10,gtol=1e-10)

def one_like(dat, f, t,start=1, end = None,error = 0):
	n = len(dat)-1
	if end is None:
		end = len(dat)
	dat = dat[start:end]
	theory = one_pulse(f,t,n,error=error)[start:end]
	theory = theory/np.sum(theory)
	lnL = np.sum(dat*np.log(theory))
	print f, t, error,  lnL
	return lnL

def two_like(dat, f1, t1, f2, t2, start = 1, end = None, error = 0):
	n = len(dat) - 1
	if end is None:
		end = len(dat)
	dat = dat[start:end]
	theory = two_pulse(f1,t1,f2,t2,n, error=error)[start:end]
	theory = theory/np.sum(theory)
	lnL = np.sum(dat*np.log(theory))
	print f1, t1, f2, t2, error, lnL
	return lnL

def sel_like(dat, f, t, gamma, start = 1, end = None, error = 0):
	n = len(dat) - 1
	if end is None:
		end = len(dat)
	dat = dat[start:end]
	theory = one_pulse_sel(f,t,gamma,n,error=error)[start:end]
	theory = theory/np.sum(theory)
	lnL = np.sum(dat*np.log(theory))
	print f, t, gamma, error, lnL
	return lnL


#project down to 100 by 100 matrix
def lchoose(N,k):
	#return -betaln(1 + int(N) - k, 1 + k) - log(int(N) + 1)
	return sp.gammaln(N+1) - sp.gammaln(N-k+1) - sp.gammaln(k+1)

def project_down(d,m):
	n = len(d)-1 #check if -1 because matrix dimensions are 170+1, 394+1
	l = np.arange(0,n+1)
	res = np.zeros(m+1)#initializes res array? check:numeric(m+1), is +1 bc R is 1 offset?
	for i in np.arange(0,m+1):
		res[i] = np.sum(d*np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))) #check this line: res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
	return res

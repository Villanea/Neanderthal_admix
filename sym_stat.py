import numpy as np
from scipy import special
from numpy import log
from scipy.special import betaln

#def symmetry_stat():
EU = np.genfromtxt('outfile_map_wholegen_dil_masked.bed', usecols=3)
EUav = np.mean(EU)/170
print EUav
AS = np.genfromtxt('outfile_map_wholegen_dil_masked.bed', usecols=4)
ASav = np.mean(AS)/394
print ASav

#initialize and fill the matrix
EU_AS = np.zeros((171, 395)) #170+1, 394+1: +1 to include fixed alleles
for i in range(0,len(AS)):
	EU_freq = EU[i]	
	AS_freq = AS[i]
	EU_AS[(EU_freq), (AS_freq)] = EU_AS[(EU_freq),(AS_freq)]+1

#project down to 100 by 100 matrix
def lchoose(N, k):
    # return -betaln(1 + int(N) - k, 1 + k) - log(int(N) + 1)
	return special.gammaln(N+1) - special.gammaln(N-k+1) - special.gammaln(k+1)

def project_down(d,m):
	n = len(d)-1 #check if because matrix dimensions are 170+1, 394+1
	l = range(0,n)
	res = np.zeros(m)#initializes res array? check:numeric(m+1), is +1 bc R is 1 offset?
	for i in range(0,m):
		res[i] = np.sum(d*np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))) #check this line: res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
	return res

EU_AS_d = np.zeros((101, 395))
for i in range(0,393):
	EU_AS_d[i,:] = project_down(EU_AS[:,i],100)

EU_AS_pd = np.zeros((101, 395))
for i in range(0,100):
	EU_AS_pd[i,:] = project_down(EU_AS_d[i,:],100)

EU_AS_pd[0,0] = 0

#calculate and write symmetry stat
sym_stat = []
for i in range(0,100):
	sym_stat = c(sym_stat, np.sum((EU_AS_pd[i,:] - EU_AS_pd[:,i]))/np.sum((EU_AS_pd[i,:] + EU_AS_pd[:,i]+1)))
outfile = open('symmetry_stat', 'a')
outfile.write(sym_stat)
outfile.write('\n')
outfile.close()
#return


#Symm = symmetry_stat

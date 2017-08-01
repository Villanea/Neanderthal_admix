import numpy as np
import scipy.special as sp
from numpy import log
from scipy.special import betaln

#def symmetry_stat():
EU = np.genfromtxt('outfile_map_wholegen_dil_masked.bed', usecols=3)
EUav = np.mean(EU)/170 #remove in final version
print EUav #remove in final version
AS = np.genfromtxt('outfile_map_wholegen_dil_masked.bed', usecols=4)
ASav = np.mean(AS)/394 #remove in final version
print ASav #remove in final version

#initialize and fill the matrix
EU_AS = np.zeros((171, 395)) #170+1, 394+1: +1 to include fixed alleles
for i in range(0,len(AS)):
	EU_freq = EU[i]	
	AS_freq = AS[i]
	EU_AS[(EU_freq), (AS_freq)] = EU_AS[(EU_freq),(AS_freq)]+1

#project down to 100 by 100 matrix
def lchoose(N,k):
	#return -betaln(1 + int(N) - k, 1 + k) - log(int(N) + 1)
	return sp.gammaln(N+1) - sp.gammaln(N-k+1) - sp.gammaln(k+1)

#test = np.exp(lchoose(np.arange(1,100),100))
#print test
def project_down(d,m):
	n = len(d)-1 #check if -1 because matrix dimensions are 170+1, 394+1
	l = np.arange(0,n+1)
	res = np.zeros(m+1)#initializes res array? check:numeric(m+1), is +1 bc R is 1 offset?
	for i in np.arange(0,m+1):
		res[i] = np.sum(d*np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))) #check this line: res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
	return res

EU_AS_d = np.zeros((101, 394))
for i in range(0,394):
	EU_AS_d[:,i] = project_down(EU_AS[:,i],100)

EU_AS_pd = np.zeros((101, 101))
for i in range(0,101):
	EU_AS_pd[i,:] = project_down(EU_AS_d[i,:],100)

EU_AS_pd[0,0] = 0

#TODO:calculate and write symmetry stat
sym_stat = []
for i in range(0,101):
	stat =  np.sum((EU_AS_pd[i,:] - EU_AS_pd[:,i]))/np.sum((EU_AS_pd[i,:] + EU_AS_pd[:,i]+1))
	sym_stat.append(stat)
outfile = open('symmetry_stat', 'a')
outfile.write(str(sym_stat))
outfile.write('\n')
outfile.close()
#return


#Symm = symmetry_stat

import numpy as np
from scipy import special

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
      return special.gammaln(N+1) - special.gammaln(N-k+1) - special.gammaln(k+1)

def project_down(d,m):
	n = len(d)-1 #check if because matrix dimensions are 170+1, 394+1
	l = range(0,n)
	res = np.zeros(m)#initializes res array? check:numeric(m+1), is +1 bc R is 1 offset?
	for i in range(0,m):
		res[i] = np.sum(d*np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))) #check this line: res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
	return n,l,res

test = project_down(EU_AS[3],100)
print test
#print EU_AS
	#return EU, AS

#Symm = symmetry_stat

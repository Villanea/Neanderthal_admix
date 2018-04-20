import numpy as np

def Pk_site(site_probs,recal = lambda x: x):
	n = len(site_probs)
	z = np.zeros(n+1)
	site_probs = recal(site_probs)
	z[0] = 1-site_probs[0]
	z[1] = site_probs[0]
	for i in range(1,n):
		for j in reversed(range(i+2)):
			if j == 0: z[j] = z[j]*(1-site_probs[i])
			elif j == i: z[j] = z[j-1]*site_probs[i]
			else:  z[j] = z[j]*(1-site_probs[i])+z[j-1]*site_probs[i]

	return z

def Pk_genome(fn,step=100000,useMax=False, recal = lambda x: x):
	chrom = "chr0"
	genomeProbs = np.zeros(0)
	i = 0
	for line in open(fn):
		splitLine = line.split()
		curChrom = splitLine[0]
		pos = int(splitLine[1])
		probs = np.array(splitLine[3:],dtype=np.float64)
		if len(genomeProbs) == 0: genomeProbs = np.zeros(len(probs)+1)
		if curChrom != chrom:
			chrom = curChrom
			end = pos
		if pos < end:
			continue
		curPk = Pk_site(probs,recal)
		if useMax:
			whichMax = np.argmax(curPk)
			genomeProbs[whichMax] += 1
		else: 
			genomeProbs += curPk
		end += step +1
		if i % 1000 == 0: print i, curChrom, pos
		i += 1
	return genomeProbs

	
		

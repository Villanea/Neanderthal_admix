import os
import re
import numpy as np
import glob

path = '/mnt/md0/villanea/neanderthal_data/CEU_lax/CEU_lax_chr22/'

chr=22
pop='CEU_lax'

#commence the header
head = []
head.append('chromosome')
head.append('position')

#make the position list
for filename in sorted(glob.glob('%s%s' %(path,'*.pos'))):
        data_file = open('%s' %(filename),"r")
cur_pos = re.split('\n', data_file.read())
cur_pos = np.array(cur_pos[:-1], dtype=np.float64)
cur_pos = np.reshape(cur_pos,(len(cur_pos),1))
#print cur_pos

end_pos = cur_pos
end_pos = end_pos + 1

#populate the chromosome list
data = np.full((len(cur_pos),1),'chr%s' %(chr),dtype=np.chararray)
#print len(data)
#print len(cur_pos)
#print len(end_pos)
#append position column
data = np.hstack((data,cur_pos))

#append second position column
data = np.hstack((data,end_pos))

#loop through all IND files in a CHR folder
for filename in sorted(glob.glob('%s%s' %(path,'*.filtered'))):
        name = os.path.basename(filename)
        name = name.strip('.filtered')
        head.append(name)
        data_file = open('%s' %(filename),"r")
        ind = re.split('\n', data_file.read())
        ind = np.array(ind[:-1])
        ind = np.reshape(ind,(len(ind),1))
	#print len(ind)
        data = np.hstack((data,ind))
	print data

#Write the array line by line, header first
outfile = open('%s_chr%s_out_bed' %(pop,chr), 'w')
outfile.write('\t'.join(head))
outfile.write('\n')
for line in range(0,(len(data)-1)):
        outfile.write('\t'.join(map(str,data[line])))
        outfile.write('\n')
outfile.close()








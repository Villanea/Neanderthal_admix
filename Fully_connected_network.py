#fully connected network
import numpy as np
import tensorflow as tf
import scipy
import glob
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten, Masking
from keras.layers import Convolution2D, MaxPooling2D
from keras.utils import np_utils
from keras.datasets import mnist
import scipy.special as sp
from keras.callbacks import EarlyStopping  

#project down
def make_PD_operator(n,m):
	PD = np.zeros((m+1,n+1))
	l = np.arange(n+1)
	for i in range(m+1):
		PD[i,:] = np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))
	return PD
	
def project_down_matrix(d,m1,m2):
	n1,n2=np.array(d.shape)-1
	PD1 = make_PD_operator(n1,m1)
	PD2 = make_PD_operator(n2,m2)
	res = np.dot(PD1,d)
	res = np.dot(res,PD2.transpose())
	return res

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

#set up number of sims to load
sims = 10000
num_cat = 5

#load up data into 4 dimensional array
data = np.empty([(sims*num_cat),64,64])

for i, file in enumerate(glob.glob('/mnt/md0/villanea/MSprime/1_pulse/symmetry_matrix_*.txt')):
	EU_AS = np.loadtxt(file)
	EU_AS_pd = project_down_matrix(EU_AS,63,63)
	EU_AS_pd[0,0] = 0
	EU_AS_pd = EU_AS_pd/np.sum(EU_AS_pd)
	data[i] = EU_AS_pd
	del EU_AS
	del EU_AS_pd
	if i >= (sims-1):
		break

for i, file in enumerate(glob.glob('/mnt/md0/villanea/MSprime/2_pulse/symmetry_matrix_*.txt')):
	EU_AS = np.loadtxt(file)
	EU_AS_pd = project_down_matrix(EU_AS,63,63)
	EU_AS_pd[0,0] = 0
	EU_AS_pd = EU_AS_pd/np.sum(EU_AS_pd)
	data[i+sims] = EU_AS_pd
	del EU_AS
	del EU_AS_pd
	if i >= (sims-1):
		break

for i, file in enumerate(glob.glob('/mnt/md0/villanea/MSprime/3_pulse/symmetry_matrix_*.txt')):
	EU_AS = np.loadtxt(file)
	EU_AS_pd = project_down_matrix(EU_AS,63,63)
	EU_AS_pd[0,0] = 0
	EU_AS_pd = EU_AS_pd/np.sum(EU_AS_pd)
	data[i+(sims*2)] = EU_AS_pd
	del EU_AS
	del EU_AS_pd
	if i >= (sims-1):
		break

for i, file in enumerate(glob.glob('/mnt/md0/villanea/MSprime/4_pulse/symmetry_matrix_*.txt')):
	EU_AS = np.loadtxt(file)
	EU_AS_pd = project_down_matrix(EU_AS,63,63)
	EU_AS_pd[0,0] = 0
	EU_AS_pd = EU_AS_pd/np.sum(EU_AS_pd)
	data[i+(sims*3)] = EU_AS_pd
	del EU_AS
	del EU_AS_pd
	if i >= (sims-1):
		break

for i, file in enumerate(glob.glob('/mnt/md0/villanea/MSprime/5_pulse/symmetry_matrix_*.txt')):
	EU_AS = np.loadtxt(file)
	EU_AS_pd = project_down_matrix(EU_AS,63,63)
	EU_AS_pd[0,0] = 0
	EU_AS_pd = EU_AS_pd/np.sum(EU_AS_pd)
	data[i+(sims*4)] = EU_AS_pd
	del EU_AS
	del EU_AS_pd
	if i >= (sims-1):
		break

data = data.reshape(data.shape[0],data.shape[1],data.shape[2],1)

#set up labels and transform into keras categorical
labels = np.concatenate((np.zeros(sims), np.ones(sims),np.full(sims,2),np.full(sims,3),np.full(sims,4)))
indices = range(sims*num_cat)
np.random.shuffle(indices)
data_random = data[indices]
labels_random = labels[indices]
labels_random = keras.utils.to_categorical(labels_random, num_cat)

data_flat = data_random.reshape((num_cat*sims), 4096)

model = Sequential()
model.add(Dense(1024, activation='relu', input_shape=(4096,)))
model.add(Dropout(0.20))
model.add(Dense(512, activation='relu'))
model.add(Dropout(0.20))
model.add(Dense(64, activation='relu'))
model.add(Dropout(0.20))
model.add(Dense(num_cat, activation='softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])

dense = model.fit(x=data_flat, y=labels_random, epochs=150, verbose=1, batch_size=64, validation_split=0.25)

yun_data = np.loadtxt('/mnt/md0/villanea/MSprime/Yun_SNP_matrix_2.txt')
yun_data_pd = project_down_matrix(yun_data,63,63)
yun_data_pd[0,0] = 0
yun_data_pd = yun_data_pd.reshape(1,yun_data_pd.shape[0],yun_data_pd.shape[1],1)
yun_data_pd = yun_data_pd/np.sum(yun_data_pd)
yun_data_flat = yun_data_pd.reshape(1, 4096)

sriram_data = np.loadtxt('/mnt/md0/schraiber/Neanderthal_admix/Srirams_SNP_matrix.trimmed.txt')
sriram_data_pd = project_down_matrix(sriram_data,63,63)
sriram_data_pd[0,0] = 0
sriram_data_pd = sriram_data_pd.reshape(1,sriram_data_pd.shape[0],sriram_data_pd.shape[1],1)
sriram_data_pd = sriram_data_pd/np.sum(sriram_data_pd)
sriram_data_flat = sriram_data_pd.reshape(1, 4096)

model.predict(yun_data_flat, batch_size=None, verbose=0,)
model.predict(sriram_data_flat, batch_size=None, verbose=0,)

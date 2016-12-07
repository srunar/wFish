import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

#Make bins of redshift
zmin = 0.2
zmax = 0.46
bin_num = 20

#Load the Fisher matrix
Fisher = np.matrix(np.load('/Users/EPICAC/Dropbox/Project/Fisher20_new.npy'))

#Calculate the inverse of the Fisher matrix
fishinv = linalg.inv(Fisher)

#Get w(z) submatrix
w_marg = fishinv[3:,3:]

#Invert w(z) submatrix to marginalize over other parameters
winv = linalg.inv(w_marg)
fw = linalg.eig(winv)[1]

#Make array of midpoints of redshift bins
zbins = np.linspace(zmin,zmax, bin_num+1)
zbin_new = []

for i in range(0,bin_num):
	zbin_new.append(0.5*(zbins[i+1]+zbins[i]))

#Plot the PCs
for i in range(0, bin_num):
	plt.plot(zbin_new, float(2*i)+np.squeeze(np.asarray(fw[:,i])))
	plt.axhline(float(2*i), 0, 0.85, color='k', linestyle='--', linewidth = 0.3)
	plt.text(0.46, float(2*i), 'i = '+str(i+1))

plt.xlabel('Redshift $(z)$', fontsize=14)
plt.ylabel('$e_i(z)$', fontsize=20)
plt.tick_params(axis='x', labelsize=18, length=10)
plt.tick_params(axis='y',left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off', labelright='off', labelbottom='off')#, labelsize=18, length=10)
plt.tick_params('both', length=5, which='minor')
plt.rc('font', family = 'serif')
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True
plt.savefig('zbins_pc.png')
plt.show()
import numpy as np 
from classint import *
from compute import *

####################################
############USER INPUT##############
####################################

#Input parameters
zmin = 0.2
zmax = 0.46
bin_num = 20

#K values should be appropriate for the choice of redshift bin
kmin = 0.0071
kmax = 0.08

#Make sure delta is chosen such that none of the parameters will be too large or too small
delta = -0.01

fsky = 0.40
sigmaov = 0.0019
sigmagz = 0.04

#for the moment, only these parameters can be varied and they can all be varied the same amount
params_name = ['H0','Omega_b', 'Omega_cdm', 'w0_fld']
params = [69.0,0.0222,0.1197,-1.0]

#feel free to change this, it is the array of output k values
karr = np.linspace(0.001, 0.1, 500)

#these should all be tailored to the individual computer
param_filename = "param.ini"
path = "/Users/EPICAC/Documents/class_public/"
output = 'sruthi_test_'

#######CREATING ORIGINAL PARAM FILE########

file_exists(path+param_filename)

write_paramfile(params,params_name , zmax,path+param_filename, output)

#######RUNNING CLASS########

run_class(path,param_filename)

#######SETTING FIDUCIAL ARRAYS########

k, P_fid = class_out(path,output,karr)

#######CREATING EMPTY FISHER MATRIX########

Fisher = np.zeros((len(params)+bin_num-1,len(params)+bin_num-1))

#######GENERATING NECESSARY INPUTS FOR FISHER########

print 'Making array of redshifts...'
zbins = np.linspace(zmin, zmax, bin_num+1)
print 'Redshift bin endpoints are', zbins
print 'Generating array of fiducial P(k) in redshift bins...'
P_total_fid = all_fiducial_zbins(zmin, zmax, bin_num, P_fid,params,k, path, param_filename, output)
print 'Generating array of varied P(k) in redshift bins...'
P_total_vary = all_varied_zbins(zmin, zmax, bin_num, delta, P_fid,params,k, path, param_filename, output)
print 'Generating array of perturbed P(k) in redshift bins...'
Pk_new = pk_wz(zmin,zmax, bin_num,P_total_fid, P_total_vary)
print 'Generating total array of perturbed P(k) for all parameters'
Pk_varied = pk_all(Pk_new, params, params_name, delta,zmax,k, path, param_filename, output)
mu_arr = np.linspace(-3,3,50)
volume = volumelement(fsky,zmax,zmin,params)
a = 0.5*(2*np.pi)**-2

#######WRITING FISHER MATRIX########
print 'Generating Fisher Matrix...'

for i in range(0,len(params)+bin_num-1):
	for j in range(0,len(params)+bin_num-1):
		Fisher[i,j] = fishij(zbins,P_total_fid, P_total_vary, Pk_new, Pk_varied, mu_arr, volume, a, P_fid, params,i,j,k, sigmaov, sigmagz,kmin,kmax,delta,fsky)
		Fisher[j,i] = Fisher[i,j]
		

print 'Fisher matrix is:' ,Fisher

#######SAVING TO OUTPUT FILE########
#should also be personalized
print 'Writing Fisher matrix to output file...'
np.save("/Users/EPICAC/Dropbox/Project/Fisher20_new.npy", Fisher)

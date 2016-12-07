import os
import numpy as np 
import subprocess
from subprocess import call
import shutil
from scipy.interpolate import UnivariateSpline


###########################################
#############CLASS INTERFACE###############
###########################################

#This module serves to connect python with CLASS. It writes the parameter files, runs CLASS and parses the output into a format that python can read and work with.

#Determine if a param file already exists, if so remove it.
def file_exists(path_to_param_file):
	output_exists = os.path.isfile(path_to_param_file)
	print output_exists
	if output_exists == True:
		subprocess.Popen("rm "+path_to_param_file, shell=True).wait()

#Write specific parameter file
#param_names=(H0, Omega_b, Omega_cdm, w0_fld)
def write_paramfile(params, param_names, redshift, param_filename, output_file):
	print "Writing param file...."
	f=open(param_filename, 'w')
	f.write('\n')
	f.write('--------------')
	f.write('\n')
	f.write('--->parameters')
	f.write('\n')
	f.write('--------------')
	f.write('\n')
	for i in range(0, len(params)):
		f.write(param_names[i] + ' ' + '=' + ' ' + str(params[i]))
		f.write('\n')
	f.write('Omega_fld = 0.')
	f.write('\n')
	f.write('wa_fld = 0.')
	f.write('\n')
	f.write('cs2_fld = 1')
	f.write('\n')
	f.write('YHe = 0.25')
	f.write('\n')
	f.write('z_reio = 10.0')
	f.write('\n')
	f.write('------------------------')
	f.write('\n')
	f.write('--->output perturbations')
	f.write('\n')
	f.write('------------------------')
	f.write('\n')
	f.write('output = tCl,pCl,lCl,mPk')
	f.write('\n')
	f.write('P_k_ini type = analytic_Pk')
	f.write('\n')
	f.write('lensing = yes')
	f.write('\n')
	f.write('k_pivot = 0.05')
	f.write('\n')
	f.write('A_s = 2.3e-9')
	f.write('\n')
	f.write('n_s = 1.0')
	f.write('\n')
	f.write('alpha_s = 0')
	f.write('\n')
	f.write('l_max_scalars = 3000')
	f.write('\n')
	f.write('l_max_tensors = 500')
	f.write('\n')
	f.write('P_k_max_h/Mpc = 1.')
	f.write('\n')
	f.write('z_pk = '+str(redshift))
	f.write('\n')
	f.write('-------------------')
	f.write('\n')
	f.write('--->redshift of dCl')
	f.write('\n')
	f.write('-------------------')
	f.write('\n')
	f.write('selection = gaussian')
	f.write('\n')
	f.write('selection_mean = 1.0')
	f.write('\n')
	f.write('selection_width = 0.5')
	f.write('\n')
	f.write('-------------------')
	f.write('\n')
	f.write('--->where to output')
	f.write('\n')
	f.write('-------------------')
	f.write('\n')
	f.write('root = output/'+output_file)
	f.write('\n')
	f.write('headers = yes')
	f.write('\n')
	f.write('format = class')
	f.write('\n')
	f.write('bessel file = yes')
	f.write('\n')
	f.write('write parameters = yeap')
	f.write('\n')
	f.write('background_verbose = 1')
	f.write('\n')
	f.write('thermodynamics_verbose = 1')
	f.write('\n')
	f.write('perturbations_verbose = 1')
	f.write('\n')
	f.write('bessels_verbose = 1')
	f.write('\n')
	f.write('transfer_verbose = 1')
	f.write('\n')
	f.write('primordial_verbose = 1')
	f.write('\n')
	f.write('spectra_verbose = 1')
	f.write('\n')
	f.write('nonlinear_verbose = 1')
	f.write('\n')
	f.write('lensing_verbose = 1')
	f.write('\n')
	f.write('output_verbose = 1')
	f.close()
	print "Finished writing param file"

#Runs class with given parameter file and waits until it finishes
def run_class(path, param_filename):
	os.chdir(path)
	print 'Currently in directory:'
	call("pwd")
	print 'Running CLASS...'
	subprocess.Popen(["./class",param_filename]).wait()
	print 'Finished running CLASS'

#Obtains the array of scales and matter power spectra. Interpolates to an input array of scales (k).
def class_out(path, output_file,k_fid):
	f = open(path+'output/'+output_file+'pk.dat','r')
	lines = f.readlines()
	k_arr = []
	Pk_arr = []
	for i in range(4, len(lines)):
		k_arr.append(float(lines[i][:25]))
		Pk_arr.append(float(lines[i][26:]))
	spl = UnivariateSpline(k_arr, Pk_arr)
	return k_fid, spl(k_fid)



	
		





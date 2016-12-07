import math
import numpy as np 
from classint import *
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline
import time


###########################################
#############COMPUTE SPECTRA###############
###########################################

#This module serves to compute the elements of the Fisher matrix

#Changing the parameter file to incorporate a perturbed parameter other than w
def param_vary(param_name, delta, params,zmax, path, param_filename):
	with open(path+param_filename,'r') as file:
		data = file.readlines()
	for i in range(0, len(data)):
		if "H0" in data[i]:
			hval = i
		if "Omega_b" in data[i]:
			obval = i
		if "Omega_cdm" in data[i]:
			ocval = i
		if "z_pk" in data[i]:
			zval = i
		if "w0_fld" in data[i]:
			wval = i
		if "Omega_fld" in data[i]:
			fval = i
	data[wval] = 'w0_fld = '+str(params[3])+'\n'
	data[zval] = 'z_pk = '+str(zmax)+'\n'
	data[fval] = 'Omega_fld = '+str(0.0)+'\n'
	if param_name == 'H0':
		data[hval] = 'H0 = '+str(params[0]+delta)+'\n'
		data[obval] = 'Omega_b = '+str(params[1])+'\n'
		data[ocval] = 'Omega_cdm = '+str(params[2])+'\n'
	if param_name == 'Omega_b':
		data[hval] = 'H0 = '+str(params[0])+'\n'
		data[obval] = 'Omega_b = '+str(params[1]+delta)+'\n'
		data[ocval] = 'Omega_cdm = '+str(params[2])+'\n'
	if param_name == 'Omega_cdm':
		data[hval] = 'H0 = '+str(params[0])+'\n'
		data[obval] = 'Omega_b = '+str(params[1])+'\n'
		data[ocval] = 'Omega_cdm = '+str(params[2]+delta)+'\n'
	with open(path+param_filename,'w') as file:
		file.writelines(data)

#Changing the parameter file to incorporate a perturbed w and a modified upper bound for the integration of P(k)
def wz_vary(z, delta,params, path, param_filename):
	with open(path+param_filename,'r') as file:
		data = file.readlines()
	for i in range(0, len(data)):
		if "H0" in data[i]:
			hval = i
		if "Omega_b" in data[i]:
			obval = i
		if "Omega_cdm" in data[i]:
			ocval = i
		if "z_pk" in data[i]:
			zval = i
		if "w0_fld" in data[i]:
			wval = i
		if "Omega_fld" in data[i]:
			fval = i
	data[hval] = 'H0 = '+str(params[0])+'\n'
	data[obval] = 'Omega_b = '+str(params[1])+'\n'
	data[ocval] = 'Omega_cdm = '+str(params[2])+'\n'
	data[wval] = 'w0_fld = '+str(params[3]+delta)+'\n'
	data[zval] = 'z_pk = '+str(z)+'\n'
	if delta == 0:
		data[fval] = 'Omega_fld = '+str(0.)+'\n'
	else:
		data[fval] = 'Omega_fld = '+str(0.1)+'\n'
	with open(path+param_filename,'w') as file:
		file.writelines(data)

#Calculates P(k) for each redshift upper bound for no modified parameters
def all_fiducial_zbins(zmin, zmax, bin_num, P_fid,params,k_fid, path, param_filename,output):
	zbins = np.linspace(zmin, zmax, bin_num+1)
	P_total_fid = np.zeros((bin_num, len(P_fid)))
	for i in range(1, len(zbins)):
		wz_vary(zbins[i], 0,params, path, param_filename)
		run_class(path,param_filename)
		k_arr, Pk_arr = class_out(path,output,k_fid)
		P_total_fid[i-1,:] = Pk_arr
	return P_total_fid

#Calculates P(k) for each redshift upper bound for modified w and redshift upper bound for P(k)
def all_varied_zbins(zmin, zmax, bin_num, delta, P_fid,params,k_fid, path,param_filename,output):
	zbins = np.linspace(zmin, zmax, bin_num+1)
	P_total_vary = np.zeros((bin_num, len(P_fid)))
	for i in range(1, len(zbins)):
		wz_vary(zbins[i], delta,params, path, param_filename)
		run_class(path,param_filename)
		k_arr, Pk_arr = class_out(path,output,k_fid)
		P_total_vary[i-1,:] = Pk_arr
	return P_total_vary

#Calculates P(k) for w(z) modified in each bin
def pk_wz(zmin ,zmax,bin_num, P_total_fid, P_total_vary):
	zbins = np.linspace(zmin, zmax, bin_num+1)
	Pk_new = 0.0*P_total_fid
	for i in range(0, bin_num):
		if i == 0:
			Pk_new[i,:] = P_total_vary[i,:]+P_total_fid[bin_num-1,:]-P_total_fid[i,:]
		elif i == bin_num-1:
			Pk_new[i,:] = P_total_fid[i-1,:]+P_total_vary[i,:]-P_total_vary[i-1,:]
		else:
			Pk_new[i,:] = P_total_fid[i-1,:]+P_total_vary[i,:]-P_total_vary[i-1,:]+P_total_fid[bin_num-1,:]-P_total_fid[i,:]
	return Pk_new

#Calculates P(k) with variations in parameters other than w(z)
def pk_all(Pk_new,params, params_name,delta,zmax,k_fid,path,param_filename, output):
	Pk_varied = np.zeros((len(Pk_new[:,0])+len(params)-1, len(Pk_new[0,:])))
	for i in range(0, len(params)):
		param_vary(params_name[i], delta,params,zmax, path, param_filename)
		run_class(path,param_filename)
		k_arr, Pk_arr = class_out(path,output,k_fid)
		Pk_varied[i,:] = Pk_arr
	for i in range(0, len(Pk_new[:,0])):
		Pk_varied[i+len(params)-1,:] = Pk_new[i,:]
	return Pk_varied

#Takes in the matter power spectrum and calculates the galaxy power spectrum for a given k, mu
def galpower(Pk, k, zmax, z,params,delta, sigmaov, sigmagz,mu):
	omega_m = params[1]+params[2]
	b = 1.0+0.84*z
	beta = (omega_m)**(0.56)/b
	if z == zmax:
		Hsquare = (params[0]**2)*((omega_m)*(1.0+z)**3 + (1.0-omega_m)*(1.0+z)**(-3.0*(1.0+params[3]+delta)))
	else:
		Hsquare = (params[0]**2)*((omega_m)*(1.0+z)**3 + (1.0-omega_m))
	sigmazsquare = ((1.0+z)**2)*(sigmaov**2 + sigmagz**2)
	return np.log(((1.0+beta*mu**2)**2)*(b**2)*Pk)-((299792.0**2)*sigmazsquare*(k**2)*(mu**2)/Hsquare)

#Calculates the invariant volume element
def volumelement(fsky, zmax, zmin, params):
	c = 299792.0
	omega_m = params[1]+params[2]
	dcmax = integrate.quad(lambda x: c/np.sqrt((params[0]**2)*((omega_m)*(1+x)**3 + (1-omega_m))), 0, zmax)[0]
	dcmin = integrate.quad(lambda x: c/np.sqrt((params[0]**2)*((omega_m)*(1+x)**3 + (1-omega_m))), 0, zmin)[0]
	return ((4.0*np.pi)/3.0)*fsky*(dcmax**3 - dcmin**3)

#Calculates the element of the Fisher matrix for parameters theta_i, theta_j
def fishij(zbins,P_total_fid, P_total_vary, Pk_new, Pk_varied, mu_arr, volume, a, P_fid, params,i,j,k, sigmaov, sigmagz,kmin,kmax,delta,fsky):
	zmax = np.max(zbins)
	zmin = np.min(zbins)
	bin_num = len(zbins)-1
	nbar = (4.0*np.pi*fsky/volume)*integrate.quad(lambda x: 640.0*(x**2)*np.exp(-x/0.35), zmin, zmax)[0]
	#print nbar
	Pgi = np.zeros((len(k), len(mu_arr)))
	Pgj = np.zeros((len(k), len(mu_arr)))
	for n in range(0,len(k)):
		for m in range(0,len(mu_arr)):
			if i == len(params)+bin_num-2:
				Pg_f = galpower(P_fid[n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pg_v = galpower(Pk_varied[i,n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pgi[n,m] = (nbar*Pg_f/(nbar*Pg_f + 1.0))*(k[n]/delta)*(Pg_v-Pg_f)
			elif i <= len(params)-2:
				Pg_f = galpower(P_fid[n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pg_v = galpower(Pk_varied[i,n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m]) 
				Pgi[n,m] = (nbar*Pg_f/(nbar*Pg_f + 1.0))*(k[n]/delta)*(Pg_v-Pg_f)
			else:
				Pg_f = galpower(P_fid[n], k[n], zmax, zbins[i-len(params)-1], params, delta, sigmaov, sigmagz, mu_arr[m])
				Pg_v = galpower(Pk_varied[i,n], k[n], zmax, zbins[i-len(params)-1], params, delta, sigmaov, sigmagz, mu_arr[m])
				Pgi[n,m] = (nbar*Pg_f/(nbar*Pg_f + 1.0))*(k[n]/delta)*(Pg_v-Pg_f)
			if j == len(params)+bin_num-2:
				Pg_f = galpower(P_fid[n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pg_v = galpower(Pk_varied[j,n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pgj[n,m] = (nbar*Pg_f/(nbar*Pg_f + 1.0))*(k[n]/delta)*(Pg_v-Pg_f)
			if j <= len(params)-2:
				Pg_f = galpower(P_fid[n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pg_v = galpower(Pk_varied[j,n], k[n], zmax, zmax, params, delta, sigmaov, sigmagz, mu_arr[m])
				Pgj[n,m] = (nbar*Pg_f/(nbar*Pg_f + 1.0))*(k[n]/delta)*(Pg_v-Pg_f)
			else:
				Pg_f = galpower(P_fid[n], k[n], zmax, zbins[j-len(params)-1], params, delta, sigmaov, sigmagz, mu_arr[m])
				Pg_v = galpower(Pk_varied[j,n], k[n], zmax, zbins[j-len(params)-1], params, delta, sigmaov, sigmagz, mu_arr[m])
				Pgj[n,m] = (nbar*Pg_f/(nbar*Pg_f + 1.0))*(k[n]/delta)*(Pg_v-Pg_f)
	pmu = 0.0*mu_arr
	for n in range(0,len(mu_arr)):
		s = UnivariateSpline(k,a*volume*Pgi[:,n]*Pgj[:,n], k=3)
		pmu[n] = s.integral(kmin,kmax)
	s = UnivariateSpline(mu_arr,pmu,k=3)
	return -1.0*s.integral(-1,1)






	

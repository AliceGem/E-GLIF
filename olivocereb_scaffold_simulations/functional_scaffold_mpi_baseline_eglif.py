# -*- coding: utf-8 -*-
"""
Created on Tue May 15 17:53:43 2018

@author: Alice Geminiani and Claudia Casellato 

"""

from __future__ import print_function  # python 2 & 3 compatible
import os
import sys

import time
import glob
import logging
import numpy as np
import h5py
import auxiliary_functions
import matplotlib.pyplot as plt
from scaffold_paramsIO import *
# import pyNN.nest as sim # only nest as simulator
import nest
import nest.raster_plot
import random
from operator import itemgetter

nest.Install("eglif_module")		# To check/change based on the PC where we run the code

# CORES = int(sys.argv[1])
CORES = 6
NODES = 1
NEU_MOD = 'eglif_cond_alpha_multisyn'            # Neuron model name


# Control variables
RECORDING_CELLS = True
MULTI_CORE = True

from mpi4py import MPI

COMM = MPI.COMM_WORLD
using_mpi4py = True


# PARAMETERS
# Single neurons
if NEU_MOD=='eglif_cond_alpha_multisyn':
	param_goc = {'t_ref': 2.0, 'C_m': 145.0,'tau_m': 44.0,'V_th': -55.0,'V_reset': -75.0,'Vinit': -62.0,'E_L': -62.0,'lambda':1.0, 'deltaV':0.4,'Ie': 16.214,'k_adap': 0.217,'k1': 0.031, 'k2': 0.023,'A1': 259.988,'A2':178.01}
	param_grc = {'t_ref': 1.5, 'C_m': 7.0,'tau_m': 24.15,'V_th': -41.0,'V_reset': -70.0,'Vinit': -62.0,'E_L': -62.0,'lambda':1.0, 'deltaV':0.3,'Ie': -0.888,'k_adap': 0.022,'k1': 0.311, 'k2': 0.041,'A1': 0.01,'A2':-0.94}
	param_pc = {'t_ref': 0.5, 'C_m': 334.0,'tau_m': 47.0,'V_th': -43.0,'V_reset': -69.0,'Vinit': -59.0,'E_L': -59.0,'lambda':4.0, 'deltaV':3.5,'Ie': 742.54,'k_adap': 1.492,'k1': 0.1950, 'k2': 0.041,'A1': 157.622,'A2':172.622}
	param_mli = {'t_ref': 1.59, 'C_m': 14.6,'tau_m': 9.125,'V_th': -53.0,'V_reset': -78.0,'Vinit': -68.0,'E_L': -68.0,'lambda':1.8, 'deltaV':1.1,'Ie': 3.711,'k_adap': 2.025,'k1': 1.887, 'k2': 1.096,'A1': 5.953,'A2':5.863}
	param_dcn = {'t_ref': 1.5, 'C_m': 142.0,'tau_m': 33.0,'V_th': -36.0,'V_reset': -55.0,'Vinit': -45.0,'E_L': -45.0,'lambda':3.5, 'deltaV':3.0,'Ie': 75.385,'k_adap': 0.408,'k1': 0.697, 'k2': 0.047,'A1': 13.857,'A2':3.477}
	param_dcnp = {'t_ref': 3.0, 'C_m': 56.0,'tau_m': 56.0,'V_th': -39.0,'V_reset': -55.0,'Vinit': -40.0,'E_L': -40.0,'lambda':0.9, 'deltaV':1.0,'Ie': 2.384,'k_adap': 0.079,'k1': 0.041, 'k2': 0.044,'A1': 176.358,'A2':176.358}
	param_io = {'t_ref': 1.0, 'C_m': 189.0,'tau_m': 11.0,'V_th': -35.0,'V_reset': -45.0,'Vinit': -45.0,'E_L': -45.0,'lambda':1.2, 'deltaV':0.8,'Ie': -18.101,'k_adap': 1.928,'k1': 0.191, 'k2': 0.091,'A1': 1810.93,'A2':1358.197}

# Parameters to avoid threshold variation
a = 50.0
b = 0.0
c = 0.0

# Synapse parameters (alpha conductances in NEST have the same rise and decay time; we use ref decay time from literature)
# In E-GLIF, 3 synaptic receptos per neuron: the first is always associated to exc, the second to inh, the third to remaining synapse type
Erev_exc = 0.0		# [mV]	#[Cavallari et al, 2014]
Erev_inh = -80.0		# [mV]
tau_exc = {'goc': 0.23, 'grc': 5.8, 'pc': 1.1, 'mli': 0.64, 'dcn': 1.0, 'dcnp': 3.64, 'io': 1.0}		#tau_exc for pc is for pf input; tau_exc for goc is for mf input; tau_exc for mli is for pf input
tau_inh = {'goc': 10.0, 'grc': 13.61, 'pc': 2.8, 'mli': 2.0, 'dcn': 0.7, 'dcnp': 1.14, 'io': 60.0}
tau_exc_cfpc = 0.4		# decay time of the CF-mono-EPSC (so tau_exc_cf for pc is for cf input)
tau_exc_pfgoc = 0.5
tau_exc_cfmli = 1.2

# Connection delays
conn_delays = {'aa_goc': 2.0, 'aa_pc': 2.0, 'bc_pc': 4.0, 'dcnp_io': 20.0, 'gj_bc': 1.0, 'gj_sc': 1.0, 'glom_dcn': 4.0,
               'glom_goc': 4.0, 'glom_grc': 4.0, 'goc_glom': 0.5, 'gj_goc': 1.0, 'goc_grc': 2.0, 'io_dcn': 4.0, 'io_dcnp': 5.0,
               'io_bc': 70.0,'io_sc': 70.0, 'io_pc': 4.0, 'pc_dcn': 4.0, 'pc_dcnp': 4.0, 'pf_bc': 5.0, 'pf_goc': 5.0,'pf_pc': 5.0,
               'pf_sc': 5.0, 'sc_pc':5.0}

sd_iomli = 10.0				# IO-MLI weights set as normal distribution to reproduce the effect of spillover-based transmission
min_iomli = 40.0

# Connection weights
# Tuned weights with EGLIF
conn_weights = {'aa_goc': 1.2, 'aa_pc': 0.7, 'bc_pc': 0.3, 'dcnp_io': 3.0, 'gj_bc': 0.2, 'gj_sc': 0.2, 'glom_dcn': 0.05,
				'glom_goc': 1.5, 'glom_grc': 0.15, 'goc_glom': 0.0, 'gj_goc': 0.3,'goc_grc': 0.6, 'io_dcn': 0.1, 'io_dcnp': 0.2,			# -0.48 gjgoc    goc-grc = 2.0
				'io_bc': 1.0,'io_sc': 1.0, 'io_pc': 350.0, 'pc_dcn': 0.4, 'pc_dcnp': 0.12, 'pf_bc': 0.015, 'pf_goc': 0.05,'pf_pc': 0.007,
				'pf_sc': 0.015, 'sc_pc': 0.3}


''' VARIABLES INITIALIZATION '''
nest.set_verbosity('M_ERROR')
nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads': CORES,
                      'total_num_virtual_procs': CORES*NODES,
                      'resolution': 0.01,
                      'overwrite_files': True})


msd = 1000  # master seed
msdrange1 = range(msd, msd + CORES)
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2 = range(msd + CORES + 1, msd + 1 + 2 * CORES)

# NEST random
nest.SetKernelStatus({'grng_seed': msd + CORES,
                      'rng_seeds': msdrange2})

# Python random (taken from NEST doc on "random_numbers")
N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
pyrngs = [np.random.RandomState(s) for s in range(msd,msd+N_vp)]
#random.seed()
#nest.SetKernelStatus({'grng_seed' : random.randint(10,10000)})


start_netw = time.time()
print(start_netw)


################## CEREBELLAR NEURONS #######################################################À

filename = 'scaffold_full_IO_400.0x400.0_v3_microzone.hdf5'
f = h5py.File(filename, 'r+')
positions = np.array(f['positions'])

if using_mpi4py:
	COMM.Barrier()

sorted_nrn_types = sorted(list(cell_type_ID.values()))

# Create a 'inverse map' of cell_type_IDaa_pc
id_2_cell_type = {val: key for key, val in cell_type_ID.items()}

# Create a dictionary with all cell names (keys)
# and lists that will contain nest models (values)
neuron_models = {key: [] for key in cell_type_ID.keys()}
for cell_id in sorted_nrn_types:
	cell_name = id_2_cell_type[cell_id]
	if cell_name != 'glomerulus':
		if cell_name not in nest.Models():
			if cell_name != 'golgi' and cell_name != 'granule' and cell_name != 'io':
				nest.CopyModel(NEU_MOD, cell_name)
			else:
				nest.CopyModel('eglif_cond_alpha_multisyn_vmin', cell_name)				# GoC, GR and IO neurons have lower limit on Vm
		if cell_name == 'golgi':
			nest.SetDefaults(cell_name, {'t_ref': param_goc['t_ref'],  # ms
                                         'C_m': param_goc['C_m'],  # pF
										 'tau_m': param_goc['tau_m'],  # ms
                                         'Vth_init': param_goc['V_th'],  # mV
										 'Vth_inf' : param_goc['V_th'],
										 'Vth_reset' : param_goc['V_th'],
										 'a': a,							# These needs to be specified to avoid threshold variation effects that are not in the final model having Vth constant
										 'b': b,
										 'c': c,
                                         'V_reset': param_goc['V_reset'],  # mV
										 'Vinit': param_goc['Vinit'],  # mV
										 'Vmin': -150.0,
                                         'E_L': param_goc['E_L'],  # mV
										 'lambda_0' : param_goc['lambda'],
										 'delta_V' : param_goc['deltaV'],
                                         'Ie_const': param_goc['Ie'],  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
										 'adaptC' :  param_goc['k_adap'],
										 'k1' : param_goc['k1'],
										 'k2' : param_goc['k2'],
										 'A1' : param_goc['A1'],
										 'A2' : param_goc['A2'],
                                         'tau_syn1': tau_exc['goc'],		# glom-goc (mf-goc)
							 			 'tau_syn2': tau_inh['goc'],		# goc-goc
                                         'tau_syn3': tau_exc_pfgoc,			# pf-goc (grc-goc)
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh,
										 'E_rev3': Erev_exc})
		elif cell_name == 'granule':
			nest.SetDefaults(cell_name, {'t_ref': param_grc['t_ref'],  # ms
										 'C_m': param_grc['C_m'],  # pF
										 'tau_m': param_grc['tau_m'],  # ms
										 'Vth_init': param_grc['V_th'],  # mV
										 'Vth_inf' : param_grc['V_th'],
										 'Vth_reset' : param_grc['V_th'],
										 'a': a,
										 'b': b,
										 'c': c,
										 'V_reset': param_grc['V_reset'],  # mV
										 'Vinit': param_grc['Vinit'],  # mV
										 'Vmin': -150.0,
										 'E_L': param_grc['E_L'],  # mV
										 'lambda_0' : param_grc['lambda'],
										 'delta_V' : param_grc['deltaV'],
										 'Ie_const': param_grc['Ie'],  # pA  
										 'adaptC' :  param_grc['k_adap'],
										 'k1' : param_grc['k1'],
										 'k2' : param_grc['k2'],
										 'A1' : param_grc['A1'],
										 'A2' : param_grc['A2'],
										 'tau_syn1': tau_exc['grc'],		# glom-grc (mf-grc)
										 'tau_syn2': tau_inh['grc'],		# goc-grc
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh})
		elif cell_name == 'purkinje':
			nest.SetDefaults(cell_name, {'t_ref': param_pc['t_ref'],  # ms
										 'C_m': param_pc['C_m'],  # pF
										 'tau_m': param_pc['tau_m'],  # ms
										 'Vth_init': param_pc['V_th'],  # mV
										 'Vth_inf' : param_pc['V_th'],
										 'Vth_reset' : param_pc['V_th'],
										 'a': a,
										 'b': b,
										 'c': c,
										 'V_reset': param_pc['V_reset'],  # mV
										 'Vinit': param_pc['Vinit'],  # mV
										 'E_L': param_pc['E_L'],  # mV
										 'lambda_0' : param_pc['lambda'],
										 'delta_V' : param_pc['deltaV'],
										 'Ie_const': param_pc['Ie'],  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
										 'adaptC' :  param_pc['k_adap'],
										 'k1' : param_pc['k1'],
										 'k2' : param_pc['k2'],
										 'A1' : param_pc['A1'],
										 'A2' : param_pc['A2'],
										 'tau_syn1': tau_exc['pc'],		# grc-pc (pf/aa-pc)
										 'tau_syn2': tau_inh['pc'],		# mli-pc
										 'tau_syn3': tau_exc_cfpc,			# io-pc (cf-pc)
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh,
										 'E_rev3': Erev_exc})
		elif cell_name == 'stellate' or cell_name == 'basket':
			nest.SetDefaults(cell_name, {'t_ref': param_mli['t_ref'],  # ms
										 'C_m': param_mli['C_m'],  # pF
										 'tau_m': param_mli['tau_m'],  # ms
										 'Vth_init': param_mli['V_th'],  # mV
										 'Vth_inf' : param_mli['V_th'],
										 'Vth_reset' : param_mli['V_th'],
										 'a': a,
										 'b': b,
										 'c': c,
										 'V_reset': param_mli['V_reset'],  # mV
										 'Vinit': param_mli['Vinit'],  # mV
										 'E_L': param_mli['E_L'],  # mV
										 'lambda_0' : param_mli['lambda'],
										 'delta_V' : param_mli['deltaV'],
										 'Ie_const': param_mli['Ie'],  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
										 'adaptC' :  param_mli['k_adap'],
										 'k1' : param_mli['k1'],
										 'k2' : param_mli['k2'],
										 'A1' : param_mli['A1'],
										 'A2' : param_mli['A2'],
										 'tau_syn1': tau_exc['mli'],		# grc-mli
										 'tau_syn2': tau_inh['mli'],		# mli-mli
										 'tau_syn3' : tau_exc_cfmli,		# io-mli (cf-mli)
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh,
										 'E_rev3': Erev_exc})
		elif cell_name == 'dcn':
			nest.SetDefaults(cell_name, {'t_ref': param_dcn['t_ref'],  # ms
										 'C_m': param_dcn['C_m'],  # pF
										 'tau_m': param_dcn['tau_m'],  # ms
										 'Vth_init': param_dcn['V_th'],  # mV
										 'Vth_inf' : param_dcn['V_th'],
										 'Vth_reset' : param_dcn['V_th'],
										 'a': a,
										 'b': b,
										 'c': c,
										 'V_reset': param_dcn['V_reset'],  # mV
										 'Vinit': param_dcn['Vinit'],  # mV
										 'E_L': param_dcn['E_L'],  # mV
										 'lambda_0' : param_dcn['lambda'],
										 'delta_V' : param_dcn['deltaV'],
										 'Ie_const': param_dcn['Ie'],  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
										 'adaptC' :  param_dcn['k_adap'],
										 'k1' : param_dcn['k1'],
										 'k2' : param_dcn['k2'],
										 'A1' : param_dcn['A1'],
										 'A2' : param_dcn['A2'],
										 'tau_syn1': tau_exc['dcn'],
										 'tau_syn2': tau_inh['dcn'],
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh})
		elif cell_name == 'dcnp':
			nest.SetDefaults(cell_name, {'t_ref': param_dcnp['t_ref'],  # ms
										 'C_m': param_dcnp['C_m'],  # pF
										 'tau_m': param_dcnp['tau_m'],  # ms
										 'Vth_init': param_dcnp['V_th'],  # mV
										 'Vth_inf' : param_dcnp['V_th'],
										 'Vth_reset' : param_dcnp['V_th'],
										 'a': a,
										 'b': b,
										 'c': c,
										 'V_reset': param_dcnp['V_reset'],  # mV
										 'Vinit': param_dcnp['Vinit'],  # mV
										 'E_L': param_dcnp['E_L'],  # mV
										 'lambda_0' : param_dcnp['lambda'],
										 'delta_V' : param_dcnp['deltaV'],
										 'Ie_const': param_dcnp['Ie'],  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
										 'adaptC' :  param_dcnp['k_adap'],
										 'k1' : param_dcnp['k1'],
										 'k2' : param_dcnp['k2'],
										 'A1' : param_dcnp['A1'],
										 'A2' : param_dcnp['A2'],
										 'tau_syn1': tau_exc['dcnp'],
										 'tau_syn2': tau_inh['dcnp'],
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh})
		else:		# Case: IO neurons
			nest.SetDefaults(cell_name, {'t_ref': param_io['t_ref'],  # ms
										 'C_m': param_io['C_m'],  # pF
										 'tau_m': param_io['tau_m'],  # ms
										 'Vth_init': param_io['V_th'],  # mV
										 'Vth_inf' : param_io['V_th'],
										 'Vth_reset' : param_io['V_th'],
										 'V_reset': param_io['V_reset'],  # mV
										 'Vinit': param_io['Vinit'],  # mV
										 'E_L': param_io['E_L'],  # mV
										 'lambda_0' : param_io['lambda'],
										 'delta_V' : param_io['deltaV'],
										 'Ie_const': param_io['Ie'],  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
										 'adaptC' :  param_io['k_adap'],
										 'k1' : param_io['k1'],
										 'k2' : param_io['k2'],
										 'A1' : param_io['A1'],
										 'A2' : param_io['A2'],
										 'a':50.0,
										 'b':0.0,
										 'c': 0.0,
										 'Vmin': -60.0,
										 'tau_syn1': tau_exc['io'],		# glom-goc (mf-goc)
										 'tau_syn2': tau_inh['io'],		# goc-goc
										 'tau_syn3': tau_exc['io'],			# pf-goc (grc-goc)
										 'E_rev1': Erev_exc,
										 'E_rev2' : Erev_inh,
										 'E_rev3': Erev_exc})
	else:
		if cell_name not in nest.Models():
			nest.CopyModel('parrot_neuron', cell_name)				# glomeruli are parrot neurons
	cell_pos = positions[positions[:, 1] == cell_id, :]
	neuron_models[cell_name] = nest.Create(cell_name, cell_pos.shape[0])


# Random initialization  (between EL-half1 and EL+half2; being half1 the half of the range between Vreset and EL, half2 the half of the range between EL and Vth)
for x in range(1,len(neuron_models['golgi']),2):
	nest.SetStatus(neuron_models['golgi'][x-1:x],{'Vinit':param_goc['Vinit']+random.randint(-13,3)})

for x in range(1,len(neuron_models['granule']),2):
	nest.SetStatus(neuron_models['granule'][x-1:x],{'Vinit':param_grc['Vinit']+random.randint(-13,10)})

for x in range(1,len(neuron_models['purkinje']),2):
	nest.SetStatus(neuron_models['purkinje'][x-1:x],{'Vinit':param_pc['Vinit']+random.randint(-10,8)})

for x in range(1,len(neuron_models['stellate']),2):
	nest.SetStatus(neuron_models['stellate'][x-1:x],{'Vinit':param_mli['Vinit']+random.randint(-10,7)})

for x in range(1,len(neuron_models['basket']),2):
	nest.SetStatus(neuron_models['basket'][x-1:x],{'Vinit':param_mli['Vinit']+random.randint(-10,7)})

for x in range(1,len(neuron_models['dcn']),2):
	nest.SetStatus(neuron_models['dcn'][x-1:x],{'Vinit':param_dcn['Vinit']+random.randint(-10,4)})

for x in range(1,len(neuron_models['dcnp']),2):
	nest.SetStatus(neuron_models['dcnp'][x-1:x],{'Vinit':param_dcnp['Vinit']+random.randint(-15,0)})

for x in range(1,len(neuron_models['io']),2):
	nest.SetStatus(neuron_models['io'][x-1:x],{'Vinit':param_io['Vinit']+random.randint(0,5)})

############################ CEREBELLAR CONNECTIVITY #########################################################

# connections in parallel on #cores (reading only pre-prepared matrices "connMatrixForEachCore.py")
def connect_neuron(conn_mat, pre, post, syn_param, conn_param='one_to_one'):
    pre_idx = [np.where(pre == x)[0][0] for x in conn_mat[:, 0] + 1]
    post_idx = [np.where(post == x)[0][0] for x in conn_mat[:, 1] + 1]
    prenn = [pre[idx] for idx in pre_idx]
    postnn = [post[idx] for idx in post_idx]
    nest.Connect(prenn, postnn, conn_spec=conn_param, syn_spec=syn_param)
    check = nest.GetConnections(pre, post)
    return check


rank = COMM.Get_rank()
print(rank)

f=h5py.File(filename, 'r+')


if MULTI_CORE:
    nameMatrix = "aa_goc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "aa_goc"

conn_aa_goc = np.array(f['connections/aa_goc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['aa_goc'], "delay": conn_delays['aa_goc'],"receptor_type":3}
aa_goc_conn = connect_neuron(conn_aa_goc, neuron_models['granule'], neuron_models['golgi'], syn_param)



if MULTI_CORE:
    nameMatrix = "aa_pc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "aa_pc"
conn_aa_pc = np.array(f['connections/aa_pc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['aa_pc'], "delay": conn_delays['aa_pc'],"receptor_type":1}
a_pc_conn = connect_neuron(conn_aa_pc, neuron_models['granule'], neuron_models['purkinje'],syn_param)



if MULTI_CORE:
    nameMatrix = "bc_pc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "bc_pc"
conn_bc_pc = np.array(f['connections/bc_pc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['bc_pc'], "delay": conn_delays['bc_pc'],"receptor_type":2}  # -5.0
bc_pc_conn = connect_neuron(conn_bc_pc, neuron_models['basket'], neuron_models['purkinje'], syn_param)



if MULTI_CORE:
    nameMatrix = "gj_bc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "gj_bc"

conn_gj_bc = np.array(f['connections/gj_bc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['gj_bc'], "delay": conn_delays['gj_bc'],"receptor_type":2}
gj_bc_conn = connect_neuron(conn_gj_bc, neuron_models['basket'], neuron_models['basket'], syn_param)



if MULTI_CORE:
    nameMatrix = "gj_sc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "gj_sc"

conn_gj_sc = np.array(f['connections/gj_sc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['gj_sc'], "delay": conn_delays['gj_sc'],"receptor_type":2}
gj_sc_conn = connect_neuron(conn_gj_sc, neuron_models['stellate'], neuron_models['stellate'], syn_param)



if MULTI_CORE:
    nameMatrix = "glom_goc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "glom_goc"

conn_glom_goc = np.array(f['connections/glom_goc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['glom_goc'], "delay": conn_delays['glom_goc'],"receptor_type":1}
glom_goc_conn = connect_neuron(conn_glom_goc, neuron_models['glomerulus'], neuron_models['golgi'], syn_param)



if MULTI_CORE:
    nameMatrix = "glom_grc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "glom_grc"

conn_glom_grc = np.array(f['connections/glom_grc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['glom_grc'], "delay": conn_delays['glom_grc'],"receptor_type":1}
glom_grc_conn = connect_neuron(conn_glom_grc, neuron_models['glomerulus'], neuron_models['granule'], syn_param)



if MULTI_CORE:
    nameMatrix = "goc_glom"+str(rank)
    exec("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "goc_glom"

conn_goc_glom = np.array(f['connections/goc_glom'])
syn_param = {"model" : "static_synapse", "weight" : conn_weights['goc_glom'], "delay": conn_delays['goc_glom'],"receptor_type":1}
goc_glom_conn = connect_neuron(conn_goc_glom, neuron_models['golgi'], neuron_models['glomerulus'],syn_param)


if MULTI_CORE:
	nameMatrix = "gj_goc" + str(rank)
	exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
	nameMatrix = "gj_goc"

conn_gj_goc = np.array(f['connections/gj_goc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['gj_goc'], "delay": conn_delays['gj_goc'],"receptor_type":2}
gj_goc_conn = connect_neuron(conn_gj_goc, neuron_models['golgi'], neuron_models['golgi'], syn_param)



if MULTI_CORE:
	nameMatrix = "goc_grc" + str(rank)
	exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "goc_grc"

conn_goc_grc = np.array(f['connections/goc_grc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['goc_grc'], "delay": conn_delays['goc_grc'],"receptor_type":2}
goc_grc_conn = connect_neuron(conn_goc_grc, neuron_models['golgi'], neuron_models['granule'], syn_param)



if MULTI_CORE:
    nameMatrix = "pc_dcn" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "pc_dcn"

conn_pc_dcn = np.array(f['connections/pc_dcn'])
syn_param = {"model": "static_synapse", "weight": conn_weights['pc_dcn'], "delay": conn_delays['pc_dcn'],"receptor_type":2}
pc_dcn_conn = connect_neuron(conn_pc_dcn, neuron_models['purkinje'], neuron_models['dcn'], syn_param)



if MULTI_CORE:
    nameMatrix = "pc_dcnp" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "pc_dcnp"

conn_pc_dcnp = np.array(f['connections/pc_dcnp'])
syn_param = {"model": "static_synapse", "weight": conn_weights['pc_dcnp'], "delay": conn_delays['pc_dcnp'],"receptor_type":2}
pc_dcnp_conn = connect_neuron(conn_pc_dcnp, neuron_models['purkinje'], neuron_models['dcnp'], syn_param)



if MULTI_CORE:
    nameMatrix = "pf_bc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "pf_bc"

conn_pf_bc = np.array(f['connections/pf_bc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['pf_bc'], "delay": conn_delays['pf_bc'],"receptor_type":1}
pf_bc_conn = connect_neuron(conn_pf_bc, neuron_models['granule'], neuron_models['basket'], syn_param)



if MULTI_CORE:
    nameMatrix = "pf_goc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "pf_goc"

conn_pf_goc = np.array(f['connections/pf_goc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['pf_goc'], "delay": conn_delays['pf_goc'],"receptor_type":3}  # 10.0
pf_goc_conn = connect_neuron(conn_pf_goc, neuron_models['granule'], neuron_models['golgi'], syn_param)



if MULTI_CORE:
    nameMatrix = "pf_pc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "pf_pc"

conn_pf_pc = np.array(f['connections/pf_pc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['pf_pc'], "delay": conn_delays['pf_pc'],"receptor_type":1}  # 10.0
pf_pc_conn = connect_neuron(conn_pf_pc, neuron_models['granule'], neuron_models['purkinje'], syn_param)



if MULTI_CORE:
    nameMatrix = "pf_sc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "pf_sc"

conn_pf_sc = np.array(f['connections/pf_sc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['pf_sc'], "delay": conn_delays['pf_sc'],"receptor_type":1}
pf_sc_conn = connect_neuron(conn_pf_sc, neuron_models['granule'], neuron_models['stellate'], syn_param)



if MULTI_CORE:
    nameMatrix = "sc_pc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "sc_pc"

conn_sc_pc = np.array(f['connections/sc_pc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['sc_pc'], "delay": conn_delays['sc_pc'],"receptor_type":2}  # -5.0
sc_pc_conn = connect_neuron(conn_sc_pc, neuron_models['stellate'], neuron_models['purkinje'], syn_param)



if MULTI_CORE:
    nameMatrix = "glom_dcn" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "glom_dcn"

conn_glom_dcn = np.array(f['connections/glom_dcn'])
syn_param = {"model": "static_synapse", "weight": conn_weights['glom_dcn'], "delay": conn_delays['glom_dcn'],"receptor_type":1}
glom_dcn_conn = connect_neuron(conn_glom_dcn, neuron_models['glomerulus'], neuron_models['dcn'], syn_param)


############################ OLIVARY CONNECTIVITY #############################################

if MULTI_CORE:
    nameMatrix = "io_pc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "io_pc"

conn_io_pc = np.array(f['connections/io_pc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['io_pc'], "delay": conn_delays['io_pc'],"receptor_type":3}
io_pc_conn = connect_neuron(conn_io_pc, neuron_models['io'], neuron_models['purkinje'], syn_param)



if MULTI_CORE:
    nameMatrix = "io_dcn" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "io_dcn"

conn_io_dcn = np.array(f['connections/io_dcn'])
syn_param = {"model": "static_synapse", "weight": conn_weights['io_dcn'], "delay": conn_delays['io_dcn'],"receptor_type":1}
io_dcn_conn = connect_neuron(conn_io_dcn, neuron_models['io'], neuron_models['dcn'], syn_param)



if MULTI_CORE:
    nameMatrix = "io_dcnp" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "io_dcnp"

conn_io_dcnp = np.array(f['connections/io_dcnp'])
syn_param = {"model": "static_synapse", "weight": conn_weights['io_dcnp'], "delay": conn_delays['io_dcnp'],"receptor_type":1}
io_dcnp_conn = connect_neuron(conn_io_dcnp, neuron_models['io'], neuron_models['dcnp'], syn_param)



if MULTI_CORE:
    nameMatrix = "dcnp_io" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "dcnp_io"

conn_dcnp_io = np.array(f['connections/dcnp_io'])
syn_param = {"model": "static_synapse", "weight":conn_weights['dcnp_io'], "delay": conn_delays['dcnp_io'],"receptor_type":2}
dcnp_io_conn = connect_neuron(conn_dcnp_io, neuron_models['dcnp'], neuron_models['io'], syn_param)



if MULTI_CORE:
    nameMatrix = "io_bc" + str(rank)
    exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
    nameMatrix = "io_bc"

conn_io_bc = np.array(f['connections/io_bc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['io_bc'], "delay": {'distribution': 'normal_clipped', 'low': min_iomli, 'mu': conn_delays['io_bc'], 'sigma': sd_iomli},"receptor_type":3}
io_bc_conn = connect_neuron(conn_io_bc, neuron_models['io'], neuron_models['basket'], syn_param)


if MULTI_CORE:
	nameMatrix = "io_sc" + str(rank)
	exec ("f=h5py.File('%s.hdf5', 'r+')" % (nameMatrix))
else:
	nameMatrix = "io_sc"

conn_io_sc = np.array(f['connections/io_sc'])
syn_param = {"model": "static_synapse", "weight": conn_weights['io_sc'], "delay": {'distribution': 'normal_clipped', 'low': min_iomli, 'mu': conn_delays['io_sc'], 'sigma': sd_iomli},"receptor_type":3}
io_sc_conn = connect_neuron(conn_io_sc, neuron_models['io'], neuron_models['stellate'], syn_param)


end_netw = time.time()  # end of network creation (placement and connectome)
print(end_netw-start_netw)


###################################################  Stim and simulation for TUNING ######################################################################

TOT_DURATION = 1760.  # mseconds
CS_START = 1000.  # beginning of stimulation 200.0 #
CS_END = 1260.  # end of stimulation 460.0 #
CS_FREQ = 36.  # Frequency in Hz - considering the background at 4 Hz (sum of Poisson processes = Poisson proc with the sum of rates)
US_START = 1250.  # beginning of stimulation 450.#
US_END = 1260.  # end of stimulation 460.#
US_FREQ = 500.  # Frequency in Hz
glom_num = len(neuron_models['glomerulus'])
io_num = len(neuron_models['io'])


CS = nest.Create('poisson_generator',params={'rate':CS_FREQ, 'start': CS_START, 'stop': CS_END})


# Localized CS

RADIUS = 150
gloms_pos = positions[positions[:,1]==cell_type_ID['glomerulus'], :]
x_c, z_c = 200., 200.

# Find glomeruli falling into the selected volume
target_gloms_idx = np.sum((gloms_pos[:,[2,4]] - np.array([x_c, z_c]))**2, axis=1).__lt__(RADIUS**2)
target_gloms = gloms_pos[target_gloms_idx,0]+1
id_stim = [glom for glom in neuron_models['glomerulus'] if glom in target_gloms]
id_stim = sorted(list(set(id_stim)))
n = len(id_stim)

nest.Connect(CS, id_stim)

# US not as Poisson to avoid that some IO do not fire:
spike_nums = np.int(np.round((US_FREQ * (US_END - US_START)) / 1000.))
US_array = np.round(np.linspace(US_START, US_END, spike_nums))				# US_array = np.random.sample(range(US_START, US_END), spike_nums) 		#

US = nest.Create("spike_generator", io_num/2,
					   params = {'spike_times': US_array})

background = nest.Create('poisson_generator',params={'rate':4.0, 'start': 0.0, 'stop': TOT_DURATION})		# Rancz

#nest.Connect(CS,neuron_models['glomerulus'])			# Input to parrot neurons
nest.Connect(background,neuron_models['glomerulus'])			# Input to parrot neurons
syn_param = {"model": "static_synapse", "weight":55.0, "delay": 0.1,"receptor_type":1}
nest.Connect(US,neuron_models['io'][:io_num/2],{'rule':'one_to_one'},syn_param)			# Input to first half of IO, corresponding to first microzone

print('Create recording devices')
## Record spikes from all cells
glom_spikes = nest.Create("spike_detector",
                          params={"withgid": True, "withtime": True, "to_file": True, "label": "glom_spikes"})
grc_spikes = nest.Create("spike_detector",
                         params={"withgid": True, "withtime": True, "to_file": True, "label": "grc_spikes"})
goc_spikes = nest.Create("spike_detector",
                         params={"withgid": True, "withtime": True, "to_file": True, "label": "goc_spikes"})
pc_spikes = nest.Create("spike_detector",
                        params={"withgid": True, "withtime": True, "to_file": True, "label": "pc_spikes"})
sc_spikes = nest.Create("spike_detector",
                        params={"withgid": True, "withtime": True, "to_file": True, "label": "sc_spikes"})
bc_spikes = nest.Create("spike_detector",
                        params={"withgid": True, "withtime": True, "to_file": True, "label": "bc_spikes"})
dcn_spikes = nest.Create("spike_detector",
                         params={"withgid": True, "withtime": True, "to_file": True, "label": "dcn_spikes"})
dcnp_spikes = nest.Create("spike_detector",
						 params={"withgid": True, "withtime": True, "to_file": True, "label": "dcnp_spikes"})
io_spikes = nest.Create("spike_detector",
						 params={"withgid": True, "withtime": True, "to_file": True, "label": "io_spikes"})

device_list = [grc_spikes, goc_spikes, glom_spikes, pc_spikes, sc_spikes, bc_spikes, dcn_spikes, dcnp_spikes, io_spikes]

# Connect all the devices to cell models
campione = random.sample(neuron_models['granule'],int(round(0.01*len(neuron_models['granule']))))
nest.Connect(campione, grc_spikes)
nest.Connect(neuron_models['golgi'], goc_spikes)
nest.Connect(neuron_models['glomerulus'], glom_spikes)
nest.Connect(neuron_models['purkinje'], pc_spikes)
nest.Connect(neuron_models['stellate'], sc_spikes)
nest.Connect(neuron_models['basket'], bc_spikes)
nest.Connect(random.sample(neuron_models['basket'],int(round(0.1*len(neuron_models['basket'])))), bc_spikes)
nest.Connect(neuron_models['dcn'], dcn_spikes)
nest.Connect(neuron_models['dcnp'], dcnp_spikes)
nest.Connect(neuron_models['io'], io_spikes)

'''
# Multimeter for DCN and IO
dcn_vm = nest.Create('multimeter')
nest.SetStatus(dcn_vm, {'withtime': True, 'record_from': ['V_m','I_stc1','I_stc2','G1','G2','G3'], 'to_file': True, 'label': 'dcn_vm'})
nest.Connect(dcn_vm, neuron_models['dcn'])

io_vm = nest.Create('multimeter')
nest.SetStatus(io_vm, {'withtime': True, 'record_from': ['V_m','I_stc2','G1','G2'], 'to_file': True, 'label': 'io_vm'})
nest.Connect(io_vm, neuron_models['io'])
'''

RECORD_VM = False
if RECORD_VM:
	grc_vm = nest.Create('multimeter')
	goc_vm = nest.Create('multimeter')
	pc_vm = nest.Create('multimeter')
	bc_vm = nest.Create('multimeter')
	sc_vm = nest.Create('multimeter')
	dcn_vm = nest.Create('multimeter')
	dcnp_vm = nest.Create('multimeter')
	io_vm = nest.Create('multimeter')

	nest.SetStatus(grc_vm, {'withtime': True, 'record_from': ['V_m','I_stc1','I_stc2','G1','G2'], 'to_file': True, 'label': 'granule_vm'})
	nest.SetStatus(goc_vm, {'withtime': True, 'record_from': ['V_m','I_stc1','I_stc2','G1','G2','G3'], 'to_file': True, 'label': 'golgi_vm'})
	nest.SetStatus(pc_vm, {'withtime': True, 'record_from': ['V_m'], 'to_file': True, 'label': 'purkinje_vm'})
	nest.SetStatus(bc_vm, {'withtime': True, 'record_from': ['V_m'], 'to_file': True, 'label': 'basket_vm'})
	nest.SetStatus(sc_vm, {'withtime': True, 'record_from': ['V_m'], 'to_file': True, 'label': 'stellate_vm'})
	nest.SetStatus(dcn_vm, {'withtime': True, 'record_from': ['V_m','I_stc1','I_stc2'], 'to_file': True, 'label': 'dcn_vm'})
	nest.SetStatus(dcnp_vm, {'withtime': True, 'record_from': ['V_m','I_stc1','I_stc2'], 'to_file': True, 'label': 'dcnp_vm'})
	nest.SetStatus(io_vm, {'withtime': True, 'record_from': ['V_m'], 'to_file': True, 'label': 'io_vm'})

	nest.Connect(grc_vm, campione)
	nest.Connect(goc_vm, neuron_models['golgi'])
	nest.Connect(pc_vm, neuron_models['purkinje'])
	nest.Connect(bc_vm, neuron_models['basket'])
	nest.Connect(sc_vm, neuron_models['stellate'])
	nest.Connect(dcn_vm, neuron_models['dcn'])
	nest.Connect(dcnp_vm, neuron_models['dcnp'])
	nest.Connect(io_vm, neuron_models['io'])

print('Simulation start')
nest.Simulate(TOT_DURATION)
print('Simulation finished')



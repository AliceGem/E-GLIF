# -*- coding: utf-8 -*-
#
# one_neuron.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

'''
Single neuron example
------------------

This script simulates a neuron driven by a constant external current
and records its membrane potential.
'''
# First, we import all necessary modules for simulation, analysis and
# plotting. Additionally, we set the verbosity to suppress info
# messages and reset the kernel.
# Resetting the kernel allows you to execute the script several
# times in a Python shell without interferences from previous NEST
# simulations. Thus, without resetting the kernel the network status
# including connections between nodes, status of neurons, devices and
# intrinsic time clocks, is kept and influences the next simulations.
import os

import nest
import nest.voltage_trace
import pylab as pl
import numpy as np
from matplotlib.pylab import *
import random


nest.Install("cerebmodule")
nest.set_verbosity("M_WARNING")
nest.ResetKernel()
nest.SetKernelStatus({"overwrite_files": True,				# Parameters for writing on files
                      "data_path": "/home/alice/workspace/E-GLIF/single_neu_simulations",
                      "data_prefix": "eglif_GoC"})

random.seed()
nest.SetKernelStatus({'grng_seed' : random.randint(10,10000)})
# Second, the nodes (neurons and devices) are created using `Create()`.
# We store the returned handles in variables for later reference.
# The `Create` function also allow you to create multiple nodes
# e.g. nest.Create('iaf_neuron',5)
# Also default parameters of the model can be configured using 'Create'
# by including a list of parameter dictionaries
# e.g. `nest.Create("iaf_neuron", params=[{'I_e':376.0}])`
# or `nest.Create("voltmeter", [{"withgid": True, "withtime": True}])`.
# In this example we will configure these parameters in an additional
# step, which is explained in the third section.


# Neurons

param_goc = {'t_ref': 2.0, 'C_m': 145.0, 'tau_m': 44.0, 'V_th': -55.0, 'V_reset': -75.0,'E_L': -62.0}
eglif_cond_alpha_multisyn = {           
            "t_ref": 2.0,
            "C_m": 145.0,
            "V_th": -55.0,
            "V_reset": -75.0,
            "E_L": -62.0,
            "Vinit": -62.0,
            "lambda_0": 1.0,             
            "tau_V": 0.4,
            "tau_m": 44.0,
            "I_e": 16.21,
            "kadap": 0.22,
            "k1": 0.03,
            "k2": 0.02,
            "A1": 259.99,
            "A2": 178.01,
            "tau_syn1": 0.1,
            "tau_syn2": 0.1,
            "E_rev1": 0.0,
            "E_rev2": -80.0,
            "E_rev3": 0.0
        }
iaf_cond_alpha = {
            "t_ref": 1.5,
            "C_m": 142.0,
            "V_th": -36.0,
            "V_reset": -55.0,
            "E_L": -45.0,
            "I_e": 180.0,           # Ref: paper "Response" - tuned with trial and error only for network conditions
            "tau_syn_ex": 1.0,
            "tau_syn_in": 0.7,
            "g_L": 4.3
        }

single_neuron_GoC = nest.Create("eglif_cond_alpha_multisyn")		# Create E-GLIF neuron (then the parameters will define the neuron type)

print(nest.GetDefaults("eglif_cond_alpha_multisyn"))
print("eglif_cond_alpha_multisyn: {0}".format( 					# Print the recordable variables in the E-GLIF model
       nest.GetDefaults("eglif_cond_alpha_multisyn")["recordables"]))

# External input current
num_dc = 9
num_freq = 6		# Number of different frequency considered - literature protocol
num_step = 10    	# Number of current step in each square wave period - literature protocol
num_tot = num_freq*num_step

current_dc_GoC = []


for i in range(num_dc):
     current_dc_GoC.append(1)


for i in range(0,num_dc):
    current_dc_GoC[i] = nest.Create("dc_generator")	# Provide a positive input current (DC)


# Square wave input
current_sw_GoC = []


for i in range(num_tot):
     current_sw_GoC.append(1)

for i in range(0,num_tot):
    current_sw_GoC[i] = nest.Create("dc_generator")	# Provide a positive input current (DC)

print('massimo:',len(current_sw_GoC))

# Measurement devices
# voltmeter = nest.Create("voltmeter",params={"interval": 0.1})		# Measure membrane voltage -> not necessary 'cos we have the multimeter
sd = nest.Create('spike_detector',
			params = {"to_file": True,
			"withgid": True,
			"label": "spikes"})

m = nest.Create("multimeter",
		        params = {"interval": 1.0,
		                 "record_from": ["V_m", "V_th", "I_dep", "I_adap", "I_gen"],
		                  "withgid": True,
		                  "to_file": True,
		                  "label": "multimeter"})

# Third, the neuron and the voltmeter are configured using
# `SetStatus()`, which expects a list of node handles and a list of
# parameter dictionaries.
# In this example we use `SetStatus()` to configure the constant
# current input to the neuron. We also want to record the global id of
# the observed nodes and set the withgid flag of the voltmeter to
# True.


# Experiments linear adaptive neuron model:

Cm_GoC = 200.0

ratio_current_Cm_GoC = [0.0,1.0,0.0,2.0,0.0,3.0,0.0,-1,0.0]

amp_res = 2.47*Cm_GoC		# [pA] amplitude sw current steps


dur_res = 0.03     # [s]

sw_period = [round(1/0.3,2),round(1/3.0,2),round(1/6.0,2),round(1/9.0,2),round(1/12.0,2),round(1/15.0,2)]		# [s] - da protocollo GR layer ([D'Angelo et al., 2001])

#sw_period = [round(1/0.5,2),round(1/2.0,2),round(1/3.5,2),round(1/5.0,2),round(1/6.25,2),round(1/7.7,2),round(1/10.0,2)]	# protocollo sperimentale
print(sw_period)


# Durate intervalli dell'intero protocollo
durate = [10, 1, 1, 1, 1, 1, 1, 1, 1]


print(np.sum(durate)*1000)

# Current for GoC
nest.SetStatus(sd,[{"withgid": True, "withtime": True}])


start_res = sum(durate[:])

cont = 1
for i in durate[:num_dc]:
    print(durate[:cont])
    nest.SetStatus(current_dc_GoC[cont-1],{'amplitude' :ratio_current_Cm_GoC[cont-1]*Cm_GoC,'start' : (np.sum(durate[:cont-1]))*1000.0, 'stop' : (np.sum(durate[:cont]))*1000.0})
    cont += 1

csw = 0
for i in range(0,num_freq):
    for k in range(0,num_step):
        if i==0:
            starting = start_res + (k*dur_res+k*sw_period[i])
            ending = starting + dur_res
            print(starting)
            print(ending)
        else:
            starting = start_res+(i)*(num_step*dur_res)+num_step*sum(sw_period[:i])+(k*dur_res+k*sw_period[i])           # [ms]
            ending = starting + dur_res
            print(starting)
            print(ending)
        nest.SetStatus(current_sw_GoC[csw],{'amplitude' :amp_res,'start' : starting*1000, 'stop' : ending*1000})
        csw += 1

print('number of csw:',csw)

#GoC neuron: instrinsic tonic activity;
nest.SetStatus(single_neuron_GoC, eglif_cond_alpha_multisyn)


# Fourth, the neuron is connected to the devices. The command
# `Connect()` has different variants. Plain `Connect()` just takes the
# handles of pre- and post-synaptic nodes and uses the default values
# for weight and delay. Note that the connection direction for the voltmeter is
# reversed compared to the spike detector, because it observes the
# neuron instead of receiving events from it. Thus, `Connect()`
# reflects the direction of signal flow in the simulation kernel
# rather than the physical process of inserting an electrode into the
# neuron. The latter semantics is presently not available in NEST.

for i in range(0,num_dc):
    nest.Connect(current_dc_GoC[i], single_neuron_GoC)


nest.Connect(m, single_neuron_GoC)
nest.Connect(single_neuron_GoC, sd)


for i in range(0,num_tot):
   nest.Connect(current_sw_GoC[i], single_neuron_GoC)


print(starting-start_res)

# Now we simulate the network using `Simulate()`, which takes the
# desired simulation time in milliseconds.

nest.Simulate(18000.0)

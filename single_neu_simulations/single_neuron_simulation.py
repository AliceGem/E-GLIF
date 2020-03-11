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

import auxiliary_functions as aux
import nest
import nest.voltage_trace
import pylab as pl
import numpy as np
from matplotlib.pylab import *
import random

nest.Install("models")
nest.set_verbosity("M_WARNING")
nest.ResetKernel()
nest.SetKernelStatus({"overwrite_files": True,				# Parameters for writing on files
                      "data_path": "/home/workspace/E-GLIF/Simulation",
                      "data_prefix": "eglif_"})

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
single_neuron_MLI = nest.Create("eglif_cond_alpha_multisyn")		# Create E-GLIF neuron (then the parameters will define the neuron type)	
 
print(nest.GetDefaults("eglif_cond_alpha_multisyn"))	
print("eglif_cond_alpha_multisyn: {0}".format( 					# Print the recordable variables in the E-GLIF model
       nest.GetDefaults("eglif_cond_alpha_multisyn")["recordables"]))

# External input current
num_dc = 9
num_freq = 6		# Number of different frequency considered - literature protocol
num_step = 10    	# Number of current step in each square wave period - literature protocol
num_tot = num_freq*num_step

current_dc_MLI = []
       

for i in range(num_dc):
     current_dc_MLI.append(1)
        
 
for i in range(0,num_dc):
    current_dc_MLI[i] = nest.Create("dc_generator")	# Provide a positive input current (DC)
      

# Square wave input
current_sw_MLI = []
      

for i in range(num_tot):
     current_sw_MLI.append(1)
        
for i in range(0,num_tot):
    current_sw_MLI[i] = nest.Create("dc_generator")	# Provide a positive input current (DC)
  


# Measurement devices
# voltmeter = nest.Create("voltmeter",params={"interval": 0.1})		# Measure membrane voltage -> not necessary 'cos we have the multimeter
sd = nest.Create('spike_detector',
			params = {"to_file": True,
			"withgid": True,
			"label": "spikes"})

m = nest.Create("multimeter",
		        params = {"interval": 1.0,
		                 "record_from": ["V_m", "V_th", "I_stc1", "I_stc2", "sum_buffer","I_gen"],
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

Cm_MLI = 14.6;
     
ratio_current_Cm_MLI = [0.0,0.82,0.0,1.64,0.0,2.47,0.0,-1.64,0.0]	
  
amp_res = 2.47*Cm_MLI		# [pA] amplitude sw current steps
        

dur_res = 0.03     # [s]

sw_period = [round(1/0.3,2),round(1/3.0,2),round(1/6.0,2),round(1/9.0,2),round(1/12.0,2),round(1/15.0,2)]		# [s] - da protocollo GR layer ([D'Angelo et al., 2001])

#sw_period = [round(1/0.5,2),round(1/2.0,2),round(1/3.5,2),round(1/5.0,2),round(1/6.25,2),round(1/7.7,2),round(1/10.0,2)]	# protocollo sperimentale
print(sw_period)


# Durate intervalli dell'intero protocollo
durate = [10, 1, 1, 1, 1, 1, 1, 1, 1]		


print(np.sum(durate)*1000)

# Current for MLI
nest.SetStatus(sd,[{"withgid": True, "withtime": True}])


start_res = sum(durate[:])

cont = 1
for i in durate[:num_dc]:
    print(durate[:cont])
    nest.SetStatus(current_dc_MLI[cont-1],{'amplitude' :ratio_current_Cm_MLI[cont-1]*Cm_MLI,'start' : (np.sum(durate[:cont-1]))*1000.0, 'stop' : (np.sum(durate[:cont]))*1000.0})
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
        nest.SetStatus(current_sw_MLI[csw],{'amplitude' :amp_res,'start' : starting*1000, 'stop' : ending*1000})
	csw += 1
       
        
        
print(csw)


#Da ottimizzazioni MATLAB
param_all = [2.0246, 1.0959, 5.8634, 1.8868, 5.9532, 3.7113]				

#MLI neuron: instrinsic tonic activity; 
nest.SetStatus(single_neuron_MLI, {'bT_size':185,
			           't_ref' : 1.59,
		          	   'C_m' : Cm_MLI,				
				   'tau_m' : 9.125,
				   'E_L' : -68.0,			
				   'Ie_const': param_all[5],
				   'Ie_sinA' : 0.0*Cm_MLI,
				   'Ie_sinF': 0.0,	
				   'adaptC' : param_all[0],
				   'Vinit':-68.0,
				   'V_reset' : -78.0,		
				   'Vth_init' : -53.0,
				   'Vth_inf' : -53.0,
				   'Vth_reset' : -53.0,
				   'lambda_0' : 1.8,	
				   'delta_V' : 1.1, 	
			 	   'k1' : param_all[3],	
				   'k2': param_all[1],
				   'R1' : 0.0,				
				   'R2init' : 1.0,	
				   'A1' : param_all[4],
			           'A2' : param_all[2],
				   })


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
    nest.Connect(current_dc_MLI[i], single_neuron_MLI)
	    

nest.Connect(m, single_neuron_MLI)		
nest.Connect(single_neuron_MLI, sd)
   
for i in range(0,num_tot):
   nest.Connect(current_sw_MLI[i], single_neuron_MLI)
 

print(starting-start_res)

# Now we simulate the network using `Simulate()`, which takes the
# desired simulation time in milliseconds.

nest.Simulate(18000.0)



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
import time

def remove_files():
    for f in os.listdir('.'):
        if '.gdf' in f or '.dat' in f:
            os.remove(f)

remove_files()

# os.environ["NEST_MODULE_PATH"] = "/home/massimo/nest-simulator-2.18.0-build/lib/nest"
# os.environ["SLI_PATH"] = "/home/massimo/nest-simulator-2.18.0-build/share/nest/sli"
# os.environ["LD_LIBRARY_PATH"] = "/home/massimo/nest-simulator-2.18.0-build/lib/nest:/home/massimo/bin/lib"
# os.environ["PATH"] = "/home/massimo/bin/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
# os.environ["SPATIALINDEX_C_LIBRARY"] = "/home/massimo/bin/lib/libspatialindex.so"
# os.environ["PYTHONPATH"] = "/home/massimo/extra-cereb-nest/Tests:/opt/amber18/lib/python3.6/site-packages/:/home/massimo/.local/nrn/lib/python:"
# import sys
# sys.path[0] = '/home/massimo/nest-simulator-2.18.0-build/lib/python3.6/site-packages/'

import nest
import numpy as np
import random

nest.Install("cerebmodule")
nest.set_verbosity("M_WARNING")

n_simulation = 0
while(n_simulation < 10):
    n_simulation += 1
    nest.ResetKernel()
    nest.SetKernelStatus({"overwrite_files": True,				# Parameters for writing on files
                        "data_path": "/home/nrp/workspace/E-GLIF/single_neu_simulations/CA1PC_sim",
                        "data_prefix": "eglif_CA1PC_"+str(n_simulation)+'_'})

    random.seed()
    seed = random.randint(10, 10000)
    print(seed)
    nest.SetKernelStatus({'grng_seed' : seed})
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

    # param_all =  K_adap, K2, A2, K1, A1, I_e
    #optim1
    #param_all = [3.3313, 0.3556, 351.3812, 0.6634, 168.2647, -56.6884 ] #  last nopt (lambda 1.5 e tau 1.0)
    param_all = [3.4178, 0.3265, 283.0757, 1.2063, 62.5273, 70.3318] # median
    #param_all = [2.0178, 0.2265, 13.0757, 1.2063, 32.5273, -30.0] # median with trial and error
    # param_all = [4.3029, 0.4314, 496.0874, 0.0007, 102.9319, 55.1788]         # optim7 with min cost and all constr negative
    #param_all = [4.3029, 0.4314, 96.0874, 0.0007, 102.9319, -10.1788]         # trial and error from optim7 with min cost and all constr negative


    param_CA1PC = {'t_ref': 2.15, 'C_m': 90.0, 'tau_m': 15.0, 'V_th': -48.0, 'V_reset': -55.0,'E_L': -65.0}

    eglif_cond_alpha_multisyn = {
                    "t_ref": param_CA1PC['t_ref'],
                    "C_m": param_CA1PC['C_m'],
                    "V_th": param_CA1PC['V_th'],
                    "V_reset": param_CA1PC['V_reset'],
                    "E_L": param_CA1PC['E_L'],
                    "Vinit": param_CA1PC['E_L'],
                    "lambda_0": 1.1,
                    "tau_V": 0.7,
                    "tau_m": param_CA1PC['tau_m'],
                    "I_e": param_all[5],
                    "kadap": param_all[0],
                    "k1": param_all[3],
                    "k2": param_all[1],
                    "A1": param_all[4],
                    "A2":param_all[2],
                    "tau_syn1": 1.0,
                    "tau_syn2": 0.7,
                    "E_rev1": 0.0,
                    "E_rev2": -80.0,
                    "E_rev3": 0.0
    }


    single_neuron_CA1PC = nest.Create("eglif_cond_alpha_multisyn")		# Create E-GLIF neuron (then the parameters will define the neuron type)

    # External input current
    num_dc = 9
    num_freq = 6		# Number of different frequency considered - literature protocol
    num_step = 10    	# Number of current step in each square wave period - literature protocol
    num_tot = num_freq*num_step

    current_dc_CA1PC = []


    for i in range(num_dc):
        current_dc_CA1PC.append(1)


    for i in range(0,num_dc):
        current_dc_CA1PC[i] = nest.Create("dc_generator")	# Provide a positive input current (DC)


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
    current_amplitude = [0.0,100.0,0.0,200.0,0.0,300.0,0.0,-100.0,0.0]

    # Durate intervalli dell'intero protocollo
    durate = [10, 0.4, 1, 0.4, 1, 0.4, 1, 1, 1]

    # Current for CA1PC
    nest.SetStatus(sd,[{"withgid": True, "withtime": True}])

    cont = 1
    for i in durate[:num_dc]:
        nest.SetStatus(current_dc_CA1PC[cont-1],{'amplitude' :current_amplitude[cont-1],'start' : (np.sum(durate[:cont-1]))*1000.0, 'stop' : (np.sum(durate[:cont]))*1000.0})
        cont += 1

    #CA1PC neuron: instrinsic tonic activity;
    nest.SetStatus(single_neuron_CA1PC, eglif_cond_alpha_multisyn)


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
        nest.Connect(current_dc_CA1PC[i], single_neuron_CA1PC)


    nest.Connect(m, single_neuron_CA1PC)
    nest.Connect(single_neuron_CA1PC, sd)

    # Now we simulate the network using `Simulate()`, which takes the
    # desired simulation time in milliseconds.

    nest.Simulate(18000.0)

    time.sleep(0.5)

print(n_simulation)

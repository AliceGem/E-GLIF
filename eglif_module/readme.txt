# Sostituire cereb_module con eglif_module
The module "cereb_module" contains neuron and synapse models that are required to simulate the cerebellar SNN 2.0:
- eglif neuron model (eglif_cond_alpha_multisyn, previously called "gif_psc_exp_multisynapse_mihalas")
- connection with SinExp kernel
- connection with CosExp kernel
- connection with Alpha kernel
- volume transmitter modified for cerebellar plasticity


TO COMPILE and INSTALL the module
$ cd <path_cereb_module>
$ cd ..
$ mkdir cereb_module_build
$ cd cereb_module_build

$ cmake -Dwith-nest=<nest_install_dir>/bin/nest-config ../cereb_module
		cmake -Dwith-nest=/home/claudia/workspace/nest-simulator/b/bin/nest-config ../cereb_module
$ make
$ make install



TO USE the module in PyNEST
import nest

nest.Install("cereb_module")

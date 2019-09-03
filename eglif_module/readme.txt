The module "eglif_module" contains neuron and synapse models that are required to simulate the cerebellar SNN 2.0:
- eglif neuron model (eglif_cond_alpha_multisyn)



TO COMPILE and INSTALL the module
$ cd <path_eglif_module>
$ cd ..
$ mkdir ceglif_module_build
$ cd eglif_module_build

$ cmake -Dwith-nest=<nest_install_dir>/bin/nest-config ../eglif_module
		
$ make
$ make install



TO USE the module in PyNEST
import nest

nest.Install("eglif_module")

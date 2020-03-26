The module "eglif_module" contains neuron and synapse models that are required to simulate cerebellar single neuron and Spiking Neural Networks with advanced point neuron dynamics:
- eglif neuron model (eglif_cond_alpha_multisyn)

To use the latest version of custom NEST models for cerebellar simulations, refer to https://github.com/dbbs-lab/cereb-nest, compatible with NEST 2.18.0

The version in the current folder is compatible with NEST 2.14 and is now deprecated.

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

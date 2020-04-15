# E-GLIF module

Source C++ NEST code for the E-GLIF module installation and usage in PyNEST (compatible with NEST 2.18.0) - cfr. with the cerebNEST module containing also other cerebellum-specific NEST models (https://github.com/dbbs-lab/cereb-nest)

To install the module (https://nest.github.io/nest-simulator/extension_modules):
* install NEST 2.18.0 following the instructions at https://nest-simulator.readthedocs.io/en/latest/installation/
* export an environment variable containing the installation directory of NEST. E.g.:
```
export NEST_INSTALL_DIR=$HOME/nest-simulator-install
```
* download the E-GLIF module source code files and save them in the folder `eglif_module`
* create a new folder to build the module and make it the current directory:
```
mkdir $HOME/eglif_module-build
cd $HOME/eglif_module-build
```

* run cmake (tested with CMake 3.10.2):
```
cmake -Dwith-nest=${NEST_INSTALL_DIR}/bin/nest-config $HOME/eglif_module
```

* make and install the module
```
make
make install
```


* set environment variable:
```
export LD_LIBRARY_PATH=${NEST_INSTALL_DIR}/share/nest/sli:$LD_LIBRARY_PATH

```

* Every time you need the module, you can install it in this way (inside your python script):
```
import nest
nest.Install("eglifmodule")
```

* You can now use the E-GLIF neuron model inside your PyNEST script, e.g.:
```
nest.Create("eglif_cond_alpha_multisyn")
```

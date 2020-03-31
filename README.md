# E-GLIF repository

Code for simulations of cerebellar Spiking Neural Networks using the NEST simulator and the E-GLIF custom neuron model, as described in [Geminiani et al., Front Neuroinform, 2018],[Geminiani et al., Front Comput Neurosci, 2019a, b]

Specifically:
* folder `eglif_module` contains source C++ NEST code for the module installation and usage in PyNEST - deprecated (refer to https://github.com/dbbs-lab/cereb-nest for updated version compatible with NEST 2.18.0).
* folder `optimMATLAB` contains MATLAB (version R2015b) code for optimization of E-GLIF neuron parameters.
* folder `single_neu_simulations` contains PyNEST code for simulations of a single neuron activity and the corresponding MATLAB code to extract quantitative information and plots.
* folder `olivocereb_scaffold_simulations` contains PyNEST code for simulations of an olivocerebellar spiking neural networks with E-GLIF point neurons, providing sensory stimuli typical of eyeblink classical conditioning, a cerebellum-driven task.


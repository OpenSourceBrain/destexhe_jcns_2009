=============================================================================================================
Network simulations of self-sustained activity in networks of adaptive exponential integrate and fire neurons
=============================================================================================================

Demo files implemented using both NEURON and PyNN.

demo_cx-lts
-----------

Simulations of self-sustained AI states in a small N=500 network of
excitatory and inhibitory neurons, described by Adaptive
Exponential (Brette-Gerstner-Izhikevich) type neurons with
exponential approach to threshold.  The connectivity is random and
there is a small proportion (5%) of LTS cells among the excitatory
neurons.  This simulation reproduces Fig. 7 of the paper below.

demo_cx_up-down
---------------

Simulations of Up-Down states in a two-layer cortical network, with
one N=2000 network and a smaller N=500 network.  Both networks have
excitatory and inhibitory neurons described by Adaptative
Exponential (Brette-Gerstner-Izhikevich) type neurons with
exponential approach to threshold.  The connectivity is random
within each network as well as between them.  In the N=500 network,
there is a small proportion (5%) of LTS cells among the excitatory
neurons.  This simulation reproduces Fig. 13 of the paper below.

See details in the following article:

Destexhe, A. Self-sustained asynchronous irregular states and
Up/Down states in thalamic, cortical and thalamocortical networks
of nonlinear integrate-and-fire neurons.  Journal of Computational
Neuroscience 27: 493-506, 2009. 

arXiv preprint: http://arxiv.org/abs/0809.0654

Original NEURON implementation by Alain Destexhe
    destexhe@unic.cnrs-gif.fr
    http://cns.iaf.cnrs-gif.fr
    
Converted to PyNN by Andrew Davison
    davison@unic.cnrs-gif.fr
and Lyle Muller
    lyle.e.muller@gmail.com


Usage (NEURON version):

    nrnivmodl
    nrngui <file.oc>

Usage (Python version):

    python <file.py> <simulator>
    
where <file.py> is one of the demo files, and <simulator>
is one of neuron, nest, pcsim, brian, facets_hardware2, etc...


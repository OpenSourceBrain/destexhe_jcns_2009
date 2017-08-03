## Network simulations of self-sustained activity in networks of adaptive exponential integrate and fire neurons

Demo files implemented using both NEURON and PyNN.

### demo_cx-lts

Simulations of self-sustained AI states in a small N=500 network of
excitatory and inhibitory neurons, described by Adaptive
Exponential (Brette-Gerstner-Izhikevich) type neurons with
exponential approach to threshold.  The connectivity is random and
there is a small proportion (5%) of LTS cells among the excitatory
neurons.  This simulation reproduces Fig. 7 of the paper below.

**Original NEURON version**

This is the code originally written for the model by Alain Destexhe. To run:

    cd demo_cx-lts
    nrnivmodl ../demo_cx_up-down
    nrngui demo_cx05_N=500b_LTS.oc

**PyNN version**

This is a conversion of the model to PyNN by Andrew Davison (davison@unic.cnrs-gif.fr). See http://andrewdavison.info/notes/porting-NEURON-PyNN/

This can be run on NEURON:

    python demo_cx05_N=500b_LTS.py neuron 

but other PyNN backends should work also, e.g. NEST, Brian. 

### demo_cx_up-down

Simulations of Up-Down states in a two-layer cortical network, with
one N=2000 network and a smaller N=500 network.  Both networks have
excitatory and inhibitory neurons described by Adaptative
Exponential (Brette-Gerstner-Izhikevich) type neurons with
exponential approach to threshold.  The connectivity is random
within each network as well as between them.  In the N=500 network,
there is a small proportion (5%) of LTS cells among the excitatory
neurons.  This simulation reproduces Fig. 13 of the paper below.

**Original NEURON version**

This is the code originally written for the model by Alain Destexhe. To run:

    cd demo_cx_up-down
    nrnivmodl
    nrngui demo_cxcx01b_N=2500_LTS.oc

**PyNN version**

This is a conversion of the model to PyNN by Lyle Muller (lyle.e.muller@gmail.com). 

This can be run on NEURON:

    python demo_cx_Up-Down.py neuron 

### Publication

Destexhe, A. Self-sustained asynchronous irregular states and
Up/Down states in thalamic, cortical and thalamocortical networks
of nonlinear integrate-and-fire neurons.  [Journal of Computational
Neuroscience 27: 493-506, 2009](https://link.springer.com/article/10.1007%2Fs10827-009-0164-4). 

arXiv preprint: http://arxiv.org/abs/0809.0654

    


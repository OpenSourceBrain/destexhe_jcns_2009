"""
  2Layer Cortical Network Simulation
  Destexhe Lab, Lyle Muller
  LTS 2-Layer Demo
  18 Mar 2010

  This file implements a PyNN version of the model detailed in

  Destexhe, A. Self-sustained asynchronous irregular states and
  Up/Down states in thalamic, cortical and thalamocortical
  networks of nonlinear integrate-and-fire neurons.
  Journal of Computational Neuroscience 27: 493-506, 2009. 

  arXiv preprint: http://arxiv.org/abs/0809.0654

  This file reproduces Fig 13 of that paper


Updated to work with PyNN 0.8, August 2016 by Andrew Davison

"""

import numpy.random, os
from pyNN.nest import *
from numpy import *
import matplotlib.pyplot as plot
from pyNN.utility import Timer, plotting

setup()

### Switch Statements
DistanceDep = True
run_time = 5000.0                 # ms
b = .005                         # b = .05 SA, .005 WA
#

# Population Numbers
py_n = 1600                   
inh_n = 400                   
pyB_n = 400
inhB_n = 100
#

# Synaptic Conductances / External Stimulation
g_e = 6e-3        # nS
g_i = 67e-3        # nS       
g_ext = 6e-3
stim_dur = 50.0
p_c = .02
scale_factor = round((py_n+inh_n) / (pyB_n+inhB_n)) 
inter_p_c = .01
#

# Parameters
py_params = {'tau_m'      : 20.0,             # ms
               'tau_syn_E'  : 5.0,
               'tau_syn_I'  : 10.0,
               'tau_refrac' : 2.5,
               'v_rest'     : -60.0,
               'v_reset'    : -60.0,
               'v_thresh'   : -50.0,
               'delta_T'    : 2.5,
               'tau_w'      : 600.0,
               'cm'         : 0.200,
               'a'          : 0.001e3,
               'b'          : b   }         

pyB_params = {'tau_m'      : 20.0,             # ms
               'tau_syn_E'  : 5.0,
               'tau_syn_I'  : 10.0,
               'tau_refrac' : 2.5,
               'v_rest'     : -60.0,
               'v_reset'    : -60.0,
               'v_thresh'   : -50.0,
               'delta_T'    : 2.5,
               'tau_w'      : 600.0,
               'cm'         : 0.200,
               'a'          : 0.001e3,
               'b'          : b   }         

inh_params = {'tau_m'      : 20.0,             # ms
               'tau_syn_E'  : 5.0,
               'tau_syn_I'  : 10.0,
               'tau_refrac' : 2.5,
               'v_rest'     : -60.0,
               'v_reset'    : -60.0,
               'v_thresh'   : -50.0,
               'delta_T'    : 2.5,
               'tau_w'      : 600.0,
               'cm'         : 0.200,
               'a'          : 0.001e3,
               'b'          : 0.0    }

V_INIT = -60

print("Building Network")

# Create Populations
py = Population(py_n, EIF_cond_alpha_isfa_ista, cellparams=py_params, initial_values={'v': V_INIT})
inh = Population(inh_n, EIF_cond_alpha_isfa_ista, cellparams=inh_params, initial_values={'v': V_INIT})
pyB = Population(pyB_n, EIF_cond_alpha_isfa_ista, cellparams=pyB_params, initial_values={'v': V_INIT})
inhB = Population(inhB_n, EIF_cond_alpha_isfa_ista, cellparams=inh_params, initial_values={'v': V_INIT})
#

# External Stimulation - Start of Simulation
ext_py = Population(1, SpikeSourcePoisson,
                    cellparams={'start':0.0, 'rate':50., 'duration':stim_dur})
ext_prj_py = Projection(ext_py, py, FixedProbabilityConnector(.02),
                        synapse_type=StaticSynapse(weight=g_ext, delay=None),
                        receptor_type='excitatory')
ext_pyB = Population(1, SpikeSourcePoisson,
                     cellparams={'start':0.0, 'rate':50., 'duration':stim_dur})
ext_prj_pyB = Projection(ext_pyB, pyB, FixedProbabilityConnector(.02),
                         synapse_type=StaticSynapse(weight=g_ext, delay=None),
                         receptor_type='excitatory')
#
                  
# LTS Subgroup - 10% of Layer B
cells = pyB.local_cells
#cells = numpy.random.permutation(cells)[0:int(0.1*len(cells))]
cells = cells[0:int(0.05*len(cells))]

for cell in cells:
    cell.a = 0.02e3
    cell.b = 0.0

# Connect Groups - Random Connect 
print("Random Connect")
rng = NumpyRNG(1235342134, parallel_safe=False)
py_py = Projection(py, py,
                   FixedProbabilityConnector(p_c, allow_self_connections=False, rng=rng),
                   synapse_type=StaticSynapse(weight=g_e, delay=None),
                   receptor_type='excitatory')
print("Number of Synapses (py_py):", len(py_py))
py_inh = Projection(py, inh,
                    FixedProbabilityConnector(p_c, allow_self_connections=False, rng=rng),
                    synapse_type=StaticSynapse(weight=g_e, delay=None),
                    receptor_type='excitatory')
print("Number of Synapses (py_inh):", len(py_inh))
inh_py = Projection(inh, py,
                    FixedProbabilityConnector(p_c, allow_self_connections=False, rng=rng),
                    synapse_type=StaticSynapse(weight=g_i, delay=None),
                    receptor_type='inhibitory')
print("Number of Synapses (inh_py):", len(inh_py))
inh_inh = Projection(inh, inh,
                     FixedProbabilityConnector(p_c, allow_self_connections=False, rng=rng),
                     synapse_type=StaticSynapse(weight=g_i, delay=None),
                     receptor_type='inhibitory')
print("Number of Synapses (inh_inh):", len(inh_inh))

pyB_pyB = Projection(pyB, pyB,
                     FixedProbabilityConnector(p_c*scale_factor, allow_self_connections=False, rng=rng),
                     synapse_type=StaticSynapse(weight=g_e, delay=None),
                     receptor_type='excitatory')
print("Number of Synapses (pyB_pyB):", len(pyB_pyB))
pyB_inhB = Projection(pyB, inhB,
                      FixedProbabilityConnector(p_c*scale_factor, allow_self_connections=False, rng=rng),
                      synapse_type=StaticSynapse(weight=g_e, delay=None),
                      receptor_type='excitatory')
print("Number of Synapses (pyB_inhB):", len(pyB_inhB))
inhB_pyB = Projection(inhB, pyB,
                      FixedProbabilityConnector(p_c*scale_factor, allow_self_connections=False, rng=rng),
                      synapse_type=StaticSynapse(weight=g_i, delay=None),
                      receptor_type='inhibitory')
print("Number of Synapses (inhB_pyB):", len(inhB_pyB))
inhB_inhB = Projection(inhB, inhB,
                       FixedProbabilityConnector(p_c*scale_factor, allow_self_connections=False, rng=rng),
                       synapse_type=StaticSynapse(weight=g_i, delay=None),
                       receptor_type='inhibitory')
print("Number of Synapses (inhB_inhB):", len(inhB_inhB))

### Excitatory Connection Between the Layers
py_pyB = Projection(py, pyB,
                    FixedProbabilityConnector(inter_p_c, allow_self_connections=False, rng=rng),
                    synapse_type=StaticSynapse(weight=g_e, delay=None),
                    receptor_type='excitatory')
pyB_py = Projection(pyB, py,
                    FixedProbabilityConnector(inter_p_c, allow_self_connections=False, rng=rng),
                    synapse_type=StaticSynapse(weight=g_e, delay=None),
                    receptor_type='excitatory')
py_inhB = Projection(py, inhB,
                     FixedProbabilityConnector(inter_p_c, allow_self_connections=False, rng=rng),
                     synapse_type=StaticSynapse(weight=g_e, delay=None),
                     receptor_type='excitatory')
pyB_inh = Projection(pyB, inh,
                     FixedProbabilityConnector(inter_p_c, allow_self_connections=False, rng=rng),
                     synapse_type=StaticSynapse(weight=g_e, delay=None),
                     receptor_type='excitatory')

# Recording
#pyB.record_v(10)
#inhB.record_v(10)
py.record('spikes')
inh.record('spikes')
pyB.record('spikes')
inhB.record('spikes')
#

print("Running Network")
timer = Timer()
timer.reset()
run(run_time)
simCPUtime = timer.elapsedTime()

print("Simulation Time: %s" % str(simCPUtime))



py_sp = py.get_data().segments[0].spiketrains
# print("Layer A Pyramidal Mean Rate (initial stimulation): %s" % str(py_sp.mean_rate(t_start=0,
#                                                               t_stop=stim_dur)))
# print("Layer A Pyramidal Mean Rate: %s" % str(py_sp.mean_rate(t_start=stim_dur,
#                                                               t_stop=run_time)))
# print("Layer A Pyramidal Mean CV: %s" % str(mean(py_sp.cv_isi(float_only=True))))
plotting.plot_spiketrains(plot.subplot(221), py_sp, label='PY Layer A', ms=0.1)

inh_sp = inh.get_data().segments[0].spiketrains
# print("Layer A Interneuron Mean Rate (initial stimulation): %s" % str(inh_sp.mean_rate(t_start=0,
#                                                               t_stop=stim_dur)))
# print("Layer A Interneuron Mean Rate: %s" % str(inh_sp.mean_rate(t_start=stim_dur,
#                                                               t_stop=run_time)))
# print("Layer A Interneuron Mean CV: %s" % str(mean(inh_sp.cv_isi(float_only=True))))
plotting.plot_spiketrains(plot.subplot(222), inh_sp, label='INH Layer A', ms=0.1)

pyB_sp = pyB.get_data().segments[0].spiketrains
# print("Layer B Pyramidal Mean Rate (initial stimulation): %s" % str(pyB_sp.mean_rate(t_start=0,
#                                                               t_stop=stim_dur)))
# print("Layer B Pyramidal Mean Rate: %s" % str(pyB_sp.mean_rate(t_start=stim_dur,
#                                                               t_stop=run_time)))
# print("Layer B Pyramidal Mean CV: %s" % str(mean(pyB_sp.cv_isi(float_only=True))))
plotting.plot_spiketrains(plot.subplot(223), pyB_sp, label='PY Layer B', ms=0.1)
plot.axhline(y=.05*pyB_n, linewidth=2, color='r')

inhB_sp = inhB.get_data().segments[0].spiketrains
# print("Layer B Interneuron Mean Rate (initial stimulation): %s" % str(inhB_sp.mean_rate(t_start=0,
#                                                               t_stop=stim_dur)))
# print("Layer B Interneuron Mean Rate: %s" % str(inhB_sp.mean_rate(t_start=stim_dur,
#                                                               t_stop=run_time)))
# print("Layer B Interneuron Mean CV: %s" % str(mean(inhB_sp.cv_isi(float_only=True))))
plotting.plot_spiketrains(plot.subplot(224), inhB_sp, label='INH Layer B', ms=0.1, xlabel="Time(ms)", xticks=True)

plot.savefig("demo_cx_Up-Down.png")

# Cleanup
end()

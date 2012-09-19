"""
 Network of IF cells - Sustained Activity

 Cortical network consisting of N=500 EX and IN cells, with
 80-20% proportion, and random connectivity
 The excitatory cells also include a proportion of LTS cells

 calculate the nb of spikes for each cell -> "numspikes_cx05_LTS500b.dat"
 calculate spiketimes -> "spiketimes_cx05_LTS500b.dat"
 write the Vm of one cell to "Vm170_cx05_LTS500b.dat" for control
 print the values of the connectivity

 cortical excitatory inputs 61.03425  -> from 1.9 % of exc cells
 cortical inhibitory inputs 14.465    -> from 1.8 % of inh cells

 This file: all cells described by IF-BG4 mechanism (RS + FS cells
 for cortex), with correct storage of spike times (also faster)
 This file: interneurons are FS, adaptation of 0.005 for RS.
 Proportion of LTS cells: 5%
  => sustained activity (AI state)
  
 Refactoring of original Python conversion of Alain's Hoc file, putting more
 of the code into the cell classes.
 Replaced multiStimexp, multiAMPAexp and multiGABAAexp with ExpSyn
"""

from neuron import h, nrn, gui, load_mechanisms
from math import sqrt, pi
from pyNN import __path__ as pyNN_path

load_mechanisms(pyNN_path[0]+"/neuron/nmodl")

h.load_file("nrngui.hoc")
    
#-----------------------------------------------------------------
#  Parameters
#-----------------------------------------------------------------

# General parameters

SEED_LTS = 428577
SEED_CONN = 193566
SEED_GEN = 983649

DT = 0.1                                        # (ms) Time step
TSTART  = 0                                     
TSTOP   = 5000
V_INIT  = -60
h.celsius = 36

# Cell parameters

LENGTH          = sqrt(20000/pi)                # in um
DIAMETER        = sqrt(20000/pi)                # in um
AREA            = 1e-8 * pi * LENGTH * DIAMETER # membrane area in cm2
TAU             = 20                            # time constant in ms
CAPACITANCE     = 1                             # capacitance in muF/cm2
G_L             = 1e-3 * CAPACITANCE / TAU      # leak conductance in S/cm2
V_REST          = -60                           # resting potential     

a_RS            = 0.001 
b_RS            = 0.1   # full adaptation
b_RS            = 0.005 # weaker adaptation
a_LTS           = 0.02
b_LTS           = 0.0
a_FS            = 0.001
b_FS            = 0.0

TAU_W           = 600
DELTA           = 2.5


# Spike parameters

VTR             = -50           # threshold in mV
VTOP            = 40            # top voltage during spike in mV
VBOT            = -60           # reset voltage in mV
REFRACTORY      = 5.0/2         # refractory period in ms (correction for a bug in IF_CG4)

# Synapse parameters

scale = 1               # scaling factor (>1 = more synapses)

TAU_E           = 5             # excitation time constant in ms
TAU_I           = 10            # inhibition time constant in ms
V_E             = 0             # excitatory reversal potential
V_I             = -80           # inhibitory reversal potential
AMPA_GMAX       = 0.006/scale
GABA_GMAX       = 0.067/scale

# Network parameters
# Cortex
N_CX = 500              # Number of cortical cells
N_I = int(N_CX/5.0)     # Number of Cx inibitory cells
N_E = N_CX - N_I        # Number of excitatory cells
PROB_CONNECT = 0.02*scale       # Connection probability in cortex
PROB_CONNECT = 0.02*2000/N_CX   # prob renormalized to size

C_I = int(N_I*PROB_CONNECT)     # nb inh synapses per neuron
C_E = int(N_E*PROB_CONNECT)     # nb exc synapses per neuron
N_GEN = N_CX            # total number of cells

PROP = 0.05             # proportion of cortical LTS cells

# Stimulation parameters

N_STIM          = N_CX/5        # number of neurons stimulated
STOPSTIM        = 50            # duration of stimulation (ms)
NSYN_STIM       = 20            # nb of stim (exc) synapses per neuron
STIM_INTERVAL   = 70            # mean interval between stims (ms)

MODEL_ID        = "cx05_LTS500b"

NEURONS_TO_PLOT = [0, 10, 20, 30, N_E, N_E+10]
NEURONS_TO_RECORD = [170, 0, N_STIM-1]

#-----------------------------------------------------------------
#  Create cells
#-----------------------------------------------------------------

class AdExpNeuron(object):
  
    def __init__(self, length=100, diameter=100/pi, g_pas=1e-5, e_pas=-60,
                 v_thresh=-50, t_refrac=2.0, v_peak=40, v_reset=-70, v_spike=-40,
                 spikewidth=1e-12, a=0.001, b=0, tau_w=600, delta=2.5,
                 tau_esyn=2.0, e_e=0.0, tau_isyn=2.0, e_i=-70):
        self.soma = nrn.Section()
        self.soma.insert('pas')
        self.adexp = h.AdExpIF(0.5, sec=self.soma)
    
        self.soma.L         = length
        self.soma.diam      = diameter
        self.soma.e_pas     = e_pas
        self.soma.g_pas     = g_pas
        self.adexp.vthresh = v_thresh # spike threshold for exponential calculation purposes
        self.adexp.trefrac = t_refrac # add the DT for compatibility with IF_CG4
        self.adexp.vpeak   = v_peak
        self.adexp.vreset  = v_reset
        self.adexp.vspike  = v_spike # spike-detection threshold
        self.adexp.spikewidth = spikewidth # for display purposes
        self.adexp.a        = a
        self.adexp.b        = b
        self.adexp.tauw    = tau_w
        self.adexp.EL       = self.soma.e_pas
        self.adexp.GL       = self.soma.g_pas*self.area()*1e-2
        self.adexp.delta    = delta
    
        self.esyn = h.ExpSyn(0.5, sec=self.soma)
        self.esyn.tau = tau_esyn
        self.esyn.e = e_e
        
        self.isyn = h.ExpSyn(0.5, sec=self.soma)
        self.isyn.tau = tau_isyn
        self.isyn.e = e_i
    
        self.stimsyn_list = []
        self.ampa_list = []
        self.gabaa_list = []
        self.spike_times = h.Vector()
        self.rec = h.NetCon(self.soma(0.5)._ref_v, None, v_spike,
                            0.0, 0.0, sec=self.soma)
        self.rec.record(self.spike_times)

    def area(self):
        return self.soma.L * self.soma.diam * pi

    def record_v(self):
        """record the Vm"""
        npt = int(float(h.tstop)/h.dt)
        self.Vm = h.Vector(npt)
        self.Vm.record(self.soma(0.5)._ref_v, h.dt)    

    def write_v(self, filename):
        f = open(filename, 'w')
        tt=0
        npt = int(float(h.tstop)/h.dt)
        f.write("%g %g\n" % (npt, h.dt))
        for i in range(0, npt):
            f.write("%g %g\n" % (tt, self.Vm.get(i)))
            tt = tt + h.dt
        f.close()                     # close file


class CXcell(AdExpNeuron):
    pass

class THcell(AdExpNeuron):
    pass

class SpikeGen(object):
    
    def __init__(self, section, latency=0, shutoff=1e6, invl=10, min_val=0, max_val=50):
        self.g = h.gen(0.5, sec=section)
        self.g.latency = latency
        self.g.shutoff = shutoff
        self.g.invl = invl
        self.g.noise = 1           # noisy stimulus        
        self.g.min_val = min_val
        self.g.max_val = max_val


#-----------------------------------------------------------------
#  Create Network
#-----------------------------------------------------------------

neuron = []
rLTS = h.Random(SEED_LTS)
nLTS = 0

def netCreate ():
    global nLTS
    RS_parameters = {
        'length': LENGTH, 'diameter': DIAMETER, 'g_pas': G_L, 'e_pas': V_REST,
        'v_thresh': VTR, 't_refrac': REFRACTORY+DT, 'v_peak': VTOP,
        'v_reset': VBOT, 'v_spike': VTR, 'spikewidth': DT, 'a': a_RS, 'b': b_RS,
        'tau_w': TAU_W, 'delta': DELTA, 'tau_esyn': TAU_E, 'e_e': V_E,
        'tau_isyn': TAU_I, 'e_i': V_I
    }
    LTS_parameters = RS_parameters.copy()
    LTS_parameters.update({'a': a_LTS, 'b': b_LTS})
    FS_parameters = RS_parameters.copy()
    FS_parameters.update({'a': a_FS, 'b': b_FS})

    for nbactual in range(0, N_E):      # create cortical cells (excitatory)
        # check if LTS cell
        if rLTS.uniform(0,1) < PROP:
            print "Cell ", nbactual, " is LTS"
            neuron.append(CXcell(**LTS_parameters))
            nLTS = nLTS + 1
        else:
            neuron.append(CXcell(**RS_parameters))

    for nbactual in range(N_E, N_CX):     # create cortical cells (inhibitory)
        neuron.append(CXcell(**FS_parameters))

#  Connect cells

rCon = h.Random(SEED_CONN)

PRINT = 2        # flag to print; 0=minimal, 1=verbose, 2=summary

def netConnect(): # local i, j, rand, distvert, nbconn
    ne = 0
    ni = 0
    ie = 0
    ii = 0
    print "Calculate connectivity of cortical cells..."
    # scan cortical cells
    for i in range(0, N_CX):
        if PRINT==1:
           if i<N_E:
                print "Cortical EX cell ",i
           else:
                print "Cortical IN cell ",i
        nbconex = 0
        nbconin = 0

        # Insert excitatory inputs
        j = 0
        while (nbconex < C_E) and (j < N_E):
            rand = rCon.uniform(0.0, 1.0)
            if (i != j) and (rand <= PROB_CONNECT):
                nc = h.NetCon(neuron[j].soma(0.5)._ref_v, neuron[i].esyn,
                          neuron[j].adexp.v_spike, 0.0, AMPA_GMAX,
                          sec=neuron[j].soma)
                neuron[i].ampa_list.append(nc)
                nbconex = nbconex + 1    
            j = j + 1
        if PRINT==1:
            print " - exc inputs from CX:", nbconex
        ne = ne + nbconex
        ie = ie + 1

        # Insert inhibitory inputs
        j =  N_E
        while (nbconin < C_I) and (j < N_CX):
            rand = rCon.uniform(0.0, 1.0)
            if (i != j) and (rand <= PROB_CONNECT):
                nc = h.NetCon(neuron[j].soma(0.5)._ref_v, neuron[i].isyn,
                          neuron[j].adexp.v_spike, 0.0, GABA_GMAX,
                          sec=neuron[j].soma)
                neuron[i].gabaa_list.append(nc)
                nbconin = nbconin + 1
            j = j + 1
        if PRINT==1:
            print " - inh inputs from CX:", nbconin
        ni= ni + nbconin
        ii = ii + 1

    if PRINT==2:
        print "MEAN SYNAPSES PER NEURON:"
        print "cortical excitatory inputs ", float(ne)/ie
        print "cortical inhibitory inputs ", float(ni)/ii

#-----------------------------------------------------------------
#  External Input
#-----------------------------------------------------------------

nstim = NSYN_STIM
stim = []

def insertStimulation():
    print "Add stimulation of cortical neurons..."
    spike_gen_parameters = {'latency': TSTART, 'shutoff': STOPSTIM,
                            'invl': STIM_INTERVAL, 'min_val': VBOT,
                            'max_val': VTOP}
    for i in range(0, N_STIM):
        for j in range(0, nstim):
            g = SpikeGen(neuron[i].soma, **spike_gen_parameters)
            stim.append(g)
            nc = h.NetCon(g.g._ref_x, neuron[i].esyn,
                          neuron[i].adexp.v_spike, 0.0, AMPA_GMAX*scale,
                          sec=neuron[i].soma)
            neuron[i].stimsyn_list.append(nc)
    g.g.seed(SEED_GEN)

#-----------------------------------------------------------------
# Simulation settings
#-----------------------------------------------------------------

h.dt = DT
h.steps_per_ms = 1.0/DT
tstart = TSTART
h.tstop = TSTOP
h.v_init = V_INIT

#-----------------------------------------------------------------
#  Add graphs
#-----------------------------------------------------------------

g = [None]*20
ngraph = 0

def addgraph(v_min, v_max, label, colour):
    global ngraph
    ngraph = ngraph+1
    ii = ngraph-1
    g[ii] = h.Graph()
    g[ii].size(tstart, h.tstop, v_min, v_max)
    g[ii].xaxis()
    g[ii].yaxis()
    g[ii].addexpr(label, colour, 0)
    g[ii].save_name("graphList[0].")
    h.graphList[0].append(g[ii])
   

print ""
print "======================================================================="
print "            Network of ",N_GEN,"IF neurons in an active state"
print "======================================================================="
print ""

#------------------------------------------------------------------------------
#  creating cells
#------------------------------------------------------------------------------
print "----[ CREATING CELLS ]----"
netCreate()

#------------------------------------------------------------------------------
#  creating network
#------------------------------------------------------------------------------
print "----[ CREATING NETWORK ]----"
netConnect()

#------------------------------------------------------------------------------
#  adding network input
#------------------------------------------------------------------------------
print "----[ ADDING NETWORK INPUT ]----"
insertStimulation()

#------------------------------------------------------------------------------
#  procedures to write spike times
#------------------------------------------------------------------------------

nspikes = []

def write_spikes():
    f = open("spiketimes_%s.dat" % MODEL_ID, 'w')
    for i in range(0, N_GEN):
        nspikes.append(neuron[i].spike_times.size())
        for j in range(0, int(nspikes[i])):
            f.write("%g %g\n" % (i, neuron[i].spike_times.x[j]))
    f.close()
  
    f = open("numspikes_%s.dat" % MODEL_ID, 'w')
    f.write("%g %g\n" % (N_GEN, h.t))             # write nb of cells and time
    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    for i in range(0, N_GEN):
        f.write("%g\n" % nspikes[i])       # write tot number of spikes
        rate = nspikes[i] * 1000.0 / TSTOP
        if i<N_E:
            sum1 = sum1 + rate
            sum2 = sum2 + rate**2
        else:
            sum3 = sum3 + rate
            sum4 = sum4 + rate**2
    f.close()

    sum1 = float(sum1) / N_E
    sum2 = sqrt( float(sum2)/N_E - sum1**2 )
    sum3 = float(sum3) / N_I
    sum4 = sqrt( float(sum4)/N_I - sum3**2 )
    return sum1, sum2, sum3, sum4



#-----------------------------------------------------------------
#  Graphs
#-----------------------------------------------------------------

h('objref py')
h.py = h.PythonObject() # lets Hoc access Python
h.nrnmainmenu()
h.nrncontrolmenu()

# adding graphs
for id in NEURONS_TO_PLOT:
    addgraph(-80, 40, "py.neuron[%d].soma(0.5).v" % id, 4)

# record the Vm
for id in NEURONS_TO_RECORD:
    neuron[id].record_v()          

#-----------------------------------------------------------------
# Procedure to run simulation and menu
#-----------------------------------------------------------------

def run_sim():
    h.init()
    h.run()

    print "Writing spikes to file..."
    rate_RS, std_RS, rate_FS, std_FS = write_spikes()
    
    print "Mean rate per RS cell (Hz) = ", rate_RS
    print " standard deviation = ", std_RS
    print "Mean rate per FS cell (Hz) = ", rate_FS
    print " standard deviation = ", std_FS
    
    for id in NEURONS_TO_RECORD:
        neuron[id].write_v("Vm%d_%s.dat" % (id, MODEL_ID))

def make_Vpanel():                    # make panel
    h.xpanel("Brette-Gerstner network")
    h.xbutton("Run simulation", "py.run_sim()")
    h.xpanel()
    

make_Vpanel()



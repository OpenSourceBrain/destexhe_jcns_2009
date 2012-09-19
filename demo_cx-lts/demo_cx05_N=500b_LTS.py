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
  
 Direct conversion of Alain's Hoc file to Python with minimal changes.
"""

from neuron import h, nrn, gui
from math import sqrt, pi

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

# Spike parameters

VTR             = -50           # threshold in mV
VTOP            = 40            # top voltage during spike in mV
VBOT            = -60           # reset voltage in mV
REFRACTORY      = 5             # refractory period in ms

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


#-----------------------------------------------------------------
#  Create cells
#-----------------------------------------------------------------

class CXcell(object):
  
  def __init__(self):
    self.soma = nrn.Section()
    self.soma.insert('pas')
    self.soma.insert('IF_BG4')
    self.nclist = []


class THcell(object):
  
  def __init__(self):
    self.soma = nrn.Section()
    self.soma.insert('pas')
    self.soma.insert('IF_BG4')
    self.nclist = []

#-----------------------------------------------------------------
#  Create Network
#-----------------------------------------------------------------

neuron = []
rLTS = h.Random(SEED_LTS)
nLTS = 0

def area(sec):
    return sec.L * sec.diam * pi

def netCreate ():
    global nLTS
    for nbactual in range(0, N_E):      # create cortical cells (excitatory)
        neuron.append(CXcell())
        assert len(neuron) == nbactual+1
        soma = neuron[nbactual].soma
        soma.L = LENGTH
        soma.diam = DIAMETER
        soma.e_pas = V_REST
        soma.g_pas = G_L
        soma.Vtr_IF_BG4 = VTR
        soma.Ref_IF_BG4 = REFRACTORY
        soma.Vtop_IF_BG4 = VTOP
        soma.Vbot_IF_BG4 = VBOT

        # Alain parameters (RS cell)
        soma.a_IF_BG4        = .001
        soma.b_IF_BG4        = 0.1           # full adaptation
        soma.b_IF_BG4        = 0.005         # weaker adaptation

        # check if LTS cell
        if rLTS.uniform(0,1) < PROP:
            print "Cell ",nbactual," is LTS"
            soma.a_IF_BG4     = .02           # LTS cell
            soma.b_IF_BG4     = 0             # LTS cell
            nLTS = nLTS + 1

        soma.tau_w_IF_BG4    = 600
        soma.EL_IF_BG4       = soma.e_pas
        soma.GL_IF_BG4       = soma.g_pas
        soma.delta_IF_BG4    = 2.5
        soma.surf_IF_BG4     = area(soma)
           
        setExpAMPA(nbactual)
        setExpGABA(nbactual)
        setExpStim(nbactual)

    for nbactual in range(N_E, N_CX):     # create cortical cells (inhibitory)
        neuron.append(CXcell())
        assert len(neuron) == nbactual+1
        soma = neuron[nbactual].soma
        soma.L = LENGTH
        soma.diam = DIAMETER
        soma.e_pas = V_REST
        soma.g_pas = G_L
        soma.Vtr_IF_BG4 = VTR
        soma.Ref_IF_BG4 = REFRACTORY
        soma.Vtop_IF_BG4 = VTOP
        soma.Vbot_IF_BG4 = VBOT

        # Alain parameters (FS cell)
        soma.a_IF_BG4        = .001
        soma.b_IF_BG4        = 0             # no adaptation
        soma.tau_w_IF_BG4    = 600
        soma.EL_IF_BG4       = soma.e_pas
        soma.GL_IF_BG4       = soma.g_pas
        soma.delta_IF_BG4    = 2.5
        soma.surf_IF_BG4     = area(soma)
           
        setExpAMPA(nbactual)
        setExpGABA(nbactual)
        setExpStim(nbactual)

def setExpAMPA(id):
    neuron[id].ampa = h.multiAMPAexp(0.5, sec=neuron[id].soma)
    neuron[id].ampa.allocate(C_E)           # allocate space for synapse
    neuron[id].ampa.q = 1
    neuron[id].ampa.gmax = AMPA_GMAX        # max conductance
    neuron[id].ampa.id = id                 # id of cell
    h.Erev_multiAMPAexp = V_E         # excitatory reversal (mV)
    h.Prethresh_multiAMPAexp = VTR    # voltage treshold for release (mV)
    h.Deadtime_multiAMPAexp = 2 * DT  # synapse "refractory"
    h.Beta_multiAMPAexp = 1.0 / TAU_E   # inhibition time constant


def setExpGABA(id):
    neuron[id].gaba = h.multiGABAAexp(0.5, sec=neuron[id].soma)
    neuron[id].gaba.allocate(C_I)           # allocate space for synapse
    neuron[id].gaba.q = 1
    neuron[id].gaba.gmax = GABA_GMAX        # max conductance
    neuron[id].gaba.id = -id                # id of cell
    h.Erev_multiGABAAexp = V_I        # inhibitory reversal (mV)
    h.Prethresh_multiGABAAexp = VTR   # voltage treshold for release (mV)
    h.Deadtime_multiGABAAexp = 2 * DT # synapse "refractory"
    h.Beta_multiGABAAexp = 1.0 / TAU_I  # inhibition time constant

def setExpStim(id):
    neuron[id].stimsyn = h.multiStimexp(0.5, sec=neuron[id].soma)
    neuron[id].stimsyn.allocate(NSYN_STIM)
    neuron[id].stimsyn.q = 1
    neuron[id].stimsyn.gmax = AMPA_GMAX*scale
    neuron[id].stimsyn.id = 15000 + id
    h.Erev_multiStimexp = V_E         # excitatory reversal potential
    h.Prethresh_multiStimexp = VTR    # voltage treshold for release (mV)
    h.Deadtime_multiStimexp = 0       # no synapse response to input
    h.Beta_multiStimexp = 1.0 / TAU_E

#  Connect cells

rCon = h.Random(SEED_CONN)

PRINT = 2        # flag to print; 0=minimal, 1=verbose, 2=summary

ne = 0
ni = 0
ie = 0
ii = 0

def netConnect(): # local i, j, rand, distvert, nbconn
    global ne, ni, ie, ii   
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
                neuron[i].ampa.addlink(neuron[j].soma(0.5)._ref_v)
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
                neuron[i].gaba.addlink(neuron[j].soma(0.5)._ref_v)
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
    for i in range(0, N_STIM):
        #access neuron[i].soma
        for j in range(0, nstim):
            g = h.gen(0.5, sec=neuron[i].soma)
            stim.append(g)
            g.latency = TSTART
            g.shutoff = STOPSTIM
            g.invl = STIM_INTERVAL
            g.noise = 1           # noisy stimulus        
            g.min_val = VBOT
            g.max_val = VTOP
            neuron[i].stimsyn.addlink(g._ref_x)
            neuron[i].nclist.extend(stim)
    neuron[0].nclist[0].seed(SEED_GEN)

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
  f = open("spiketimes_cx05_LTS500b.dat", 'w')
  for i in range(0, N_GEN):
    nspikes.append(neuron[i].soma.nspike_IF_BG4)
    for j in range(0, int(nspikes[i])):
        f.write("%g %g\n" % (i, neuron[i].soma.spiketimes_IF_BG4[j]))
  f.close()

#-----------------------------------------------------------------
#  Graphs
#-----------------------------------------------------------------

h('objref py')
h.py = h.PythonObject() # lets Hoc access Python
h.nrnmainmenu()
h.nrncontrolmenu()
  
# access origin of network
##access neuron[0].soma

# adding graphs
addgraph(-80, 40, "py.neuron[0].soma(0.5).v", 4)   # excitatory CX
addgraph(-80, 40, "py.neuron[10].soma(0.5).v", 4)
addgraph(-80, 40, "py.neuron[20].soma(0.5).v", 4)
addgraph(-80, 40, "py.neuron[30].soma(0.5).v", 4)

addgraph(-80, 40, "py.neuron[%d].soma(0.5).v" % N_E, 4)  # inhibitory CX
addgraph(-80, 40, "py.neuron[%d].soma(0.5).v" % (N_E+10,), 4)

npt = int(float(h.tstop)/h.dt)
Vm = h.Vector(npt)
Vm.record(neuron[170].soma(0.5)._ref_v, h.dt)          # record the Vm


#-----------------------------------------------------------------
# Procedure to run simulation and menu
#-----------------------------------------------------------------

def run_sim():
    h.init()
    h.run()

    print "Writing spikes on file..."
    write_spikes()
    
    f = open("numspikes_cx05_LTS500b.dat", 'w')
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
    print "Mean rate per RS cell (Hz) = ", sum1
    print " standard deviation = ", sum2
    sum3 = float(sum3) / N_I
    sum4 = sqrt( float(sum4)/N_I - sum3**2 )
    print "Mean rate per FS cell (Hz) = ", sum3
    print " standard deviation = ", sum4

    f = open("Vm170_cx05_LTS500b.dat", 'w')
    tt=0
    f.write("%g %g\n" % (npt, h.dt))
    for i in range(0, npt):
        f.write("%g %g\n" % (tt, Vm.get(i)))
        tt = tt + h.dt
  
    f.close()                     # close file

def make_Vpanel():                    # make panel
    h.xpanel("Brette-Gerstner network")
    h.xbutton("Run simulation", "py.run_sim()")
    h.xpanel()
    

make_Vpanel()



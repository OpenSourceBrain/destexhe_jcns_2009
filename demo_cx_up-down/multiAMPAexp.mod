TITLE multiple AMPA receptors

COMMENT
-----------------------------------------------------------------------------
Mechanism for handling multiple presynaptic entries to the same compartment;
up to 1000 synapses can be handled using a different pointer that must be set
for each presynaptic variable using addlink.  

Model of AMPA receptors: "exponential" model modified from a first-order kinetic
model with pulse of transmitter (see Destexhe, A., Mainen, Z. and Sejnowski, T.J.
Neural Computation, 6: 14-18, 1994).  The parameters were obtained from fitting
the model to whole-cell recorded AMPA postsynaptic currents (Xiang et al., Soc
Neurosci. Abstracts 18: 1350, 1992).  The fit was performed using a simplex
algorithm using short pulses of transmitter (0.5 mM during 0.3 ms).

The kinetic model was modified to yield a simple "exponential" model, with
instantaneous rise and exponential decay (the decay phase is the same as in the
kinetic model).  Another simplification is that there is no "saturation"
mechanism of presynaptic inputs arriving at the same terminals, so the
exponentials just summate linearly.  This mechanism was implemented in hardware
(ASIC circuits) by the group of S. Renaud-LeMasson (University of Bordeaux,
France).

-----------------------------------------------------------------------------
EXAMPLE OF HOW TO USE:

create POST,PRE[10]		// create compartments
objectvar c			// create an object
c = new multiAMPAexp()		// create multiple AMPA kinetic synapses
POST c.loc(0.5)			// localize synapse on postsyn compartment
c.gmax = 0.001			// assign max conductance of each syn (mu S)
c.allocate(10)			// allocate space for 10 presyn variables
for i=0,9 { 			// link presynaptic variables
   c.addlink(&PRE[i].v)
}  
-----------------------------------------------------------------------------
WARNINGS:

  This mechanism only works for cases where all weights are equal

  Alain Destexhe, CNRS, 2002

-----------------------------------------------------------------------------
ENDCOMMENT

: defines maximal number of possible links to presynaptic variables
: this number should correpond to the number of pointers pre00, pre01, ...
: defined in the NEURON block

DEFINE MAXSYNAMPA 1000
VERBATIM
static int MAXSYNAMPA = 1000;
ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS multiAMPAexp
	NONSPECIFIC_CURRENT i
	RANGE R, ri, nsyn, g, id, gmax, q
	GLOBAL Beta, Erev, Prethresh, Deadtime
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt		(ms)

	Beta	= 0.18	(/ms)		: backward (unbinding) rate
	Erev	= 60	(mV)		: reversal potential
        gmax            (umho)          : maximum conductance of each synapse
	q	= 0.01			: quantal increase of R at each release
					: (quantal conductance is q*gmax)
	id

	Prethresh = 0 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
}


ASSIGNED {
	lastrelease[MAXSYNAMPA] (ms)	: last release for each synapse
	R				: sum of all synapses
	nsyn				: number of synapses

	v		(mV)		: postsynaptic voltage
	i 		(nA)		: total current = g*(v - Erev)
	g 		(umho)		: total conductance

	ptr_array_ampa			: pointer array
}

INITIAL { LOCAL j
	FROM j=0 TO nsyn-1 {
		lastrelease[j] = -9e9
	}
	R = 0
}

BREAKPOINT {
   if(gmax > 0) {
	SOLVE release
	g = gmax * R
	i = g*(v - Erev)
   } else {
	i = 0
   }
:VERBATIM
:    if(id == 228) {
:	printf("AMPA: id, t, i = %2.2f\t%2.2f\t%2.4f\n", id, t, i);
:    }
:ENDVERBATIM
}

VERBATIM
#define ppampa ((double***)(&(ptr_array_ampa)))
extern double* hoc_pgetarg();
ENDVERBATIM


PROCEDURE release() { LOCAL j,trel

  FROM j=0 TO nsyn-1 {

    trel = (t - lastrelease[j])		: time since last release

    if (trel > Deadtime) {			: ready for another release?
	VERBATIM			
	if (*((*ppampa)[_lj]) > Prethresh) {	// spike occured?
	  lastrelease[_lj] = t;		// memorize release time
	  R = R + q;			// increase state variable by quantum
	}
	ENDVERBATIM					
    }
  }

  if(R > 0) {				: exp decay of R
     R = R * exptable(- Beta * dt)
     if(R < 1e-9) { R = 0 }		: prevent roundoff errors
  }

}


FUNCTION exptable(x) { 
	TABLE  FROM -25 TO 25 WITH 10000

	if ((x > -25) && (x < 25)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}



:FUNCTION exptable(x) { 
:	TABLE  FROM -10 TO 10 WITH 10000
:
:	if ((x > -10) && (x < 10)) {
:		exptable = exp(x)
:	} else {
:		exptable = 0.
:	}
:}



:FUNCTION exptable(x) {
:	if(x > -50) {
:		exptable = exp(x)
:	} else {
:		exptable = 0.
:	}
:}



:-------------------------------------------------------------------
:  Procedures for pointer arrays in nmodl 
:  create a pointer array and link its pointers to variables passed
:  from hoc (adapted from Mike Hines)
:-------------------------------------------------------------------


VERBATIM
#define ppampa ((double***)(&(ptr_array_ampa)))
extern double* hoc_pgetarg();
ENDVERBATIM


:
: Procedure to allocate space for n pointers
:
PROCEDURE allocate(n) {
  VERBATIM
	if (*ppampa) {
	   free(*ppampa);
	}
	*ppampa = ((double**) hoc_Ecalloc((int)_ln, sizeof(double *))), hoc_malchk();
  ENDVERBATIM
}

:
: procedure to get the value of a presynaptic variable
: index is the number of the presynaptic var
:
FUNCTION presynaptic(index) {
  VERBATIM
	if(_lindex >= nsyn) {
	   printf("Warning: attempt to use pointer outside range\n");
	   printf(" trying to use pointer number %d\n",(int)_lindex);
	   printf(" but number of defined pointers was nsyn=%d.\n",(int) nsyn);
	}
	_lpresynaptic = *((*ppampa)[(int)_lindex]);
  ENDVERBATIM
}


:
: procedure to add a new presynaptic variable
: the address of the variable is passed as argument (from hoc)
: a new pointer is then linked to that variable
:
PROCEDURE addlink() {
  VERBATIM
	if(++nsyn > MAXSYNAMPA) {
	  printf("Exceeding maximum of allowed links MAXSYNAMPA=%d\n",MAXSYNAMPA);
	  printf("  edit the nmodl code to increase the maximum allowed.\n");
	  exit(-1);
	}
	(*ppampa)[(int)(nsyn-1)] = hoc_pgetarg(1);
  ENDVERBATIM
}

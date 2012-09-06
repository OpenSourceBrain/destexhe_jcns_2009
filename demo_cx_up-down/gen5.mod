COMMENT

Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a fast spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed), and 
either tonic or bursting.

Parameters;
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   fast_invl: 	fast interval, mean time between spikes (ms)
   slow_invl:	slow interval, mean burst silent period (ms), 0=tonic train
   burst_len: 	mean burst length (nb. spikes)
   latency:	latency at which spikes begin
   shutoff:	time at which the generator is shutoff

Control Variable:
   on:		on=1 the generator is on, on=0 it is interrupted

Note: the "slow" interval is actually equal to slow_invl+fast_invl


Written by A. Destexhe and Z. Mainen, The Salk Institute, 1994

ENDCOMMENT




INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}


NEURON	{ POINT_PROCESS gen
	  RANGE x, latency, shutoff, on
	  RANGE invl
	  RANGE noise, min_val, max_val
	  RANGE lcntr, fac

}

PARAMETER {
	invl	= 10		: burst period (msec)
	noise		= 1		: amount of randomness (0.0 - 1.0)
	dt		(ms)
	min_val		= 0		: min value of presynaptic variable
	max_val		= 50		: max value of presynaptic variable
	latency		= 0 (ms)	: time at which spikes begin
	on		= 1		: logical on/off variable	
	shutoff		= 1e6 (ms)	: shutoff time
	fac 		= 1000		: 
}

ASSIGNED {
	lcntr
	x
}

INITIAL {
	lcntr	= fpoisrand(invl, noise)
	x 	= min_val
}


BREAKPOINT {
	SOLVE generate
}
 

PROCEDURE generate() {	

    x = min_val

    if( (on>0) && (t>=latency) && (t<=shutoff) )  {
        
	lcntr = lcntr - dt	
	if (lcntr <= dt + 1e-6) {	:+1e-6 for numerical stability
	    x = max_val			
	    lcntr = fpoisrand(invl, noise)
:    VERBATIM
:        printf("t, x = %2.2f\t%2.2f\n", t, x);
:    ENDVERBATIM
	}
    }
    VERBATIM
        return 0;
    ENDVERBATIM
}	


FUNCTION fpoisrand(mean, noise) {

	: gets the interval from poisson distribution
	: mean = average values
	: noise = amount of randomness (must be between 0 and 1)

        if (mean <= 0.) {
                mean = .01 (ms) : I would worry if it were 0.
        }
        if (noise == 0) {
               fpoisrand = mean
        }else{
               fpoisrand = (1. - noise)*mean + noise*mean*exprand(1)
        }
}







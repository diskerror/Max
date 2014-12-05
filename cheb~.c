/**
*	Chebyshev filter implimentation
*	version 2004-09-02
*   
*   Copyright 2008 Reid A. Woodbury Jr.
*	
*	Part of this code was adapted from
*       "The Scientist and Engineer's Guide to Digital Signal Processing" 2nd edition
*       by Steven W. Smith
*       Chebyshev filter page 340, Table 20-4
*   
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*   
*      http://www.apache.org/licenses/LICENSE-2.0
*   
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#include "ext.h"  		// c74support/max-includes/ext.h
#include "z_dsp.h"		// c74support/max-includes/z_dsp.h
#include <math.h>
//#include <string.h>

//	from fp.h
#ifndef	pi
#define	pi			3.14159265358979323846
#endif

//#define	numCoeffs	23
//#define sampRate	44100.0

//	double	T		= 2.0 * tan(0.5);
#define	T			1.092604979687581
//	double	TT		= T * T;
#define	TT			1.193785641638099

struct t_chebCoef	// defines our object's internal variables for each instance in a patch
{
	t_object p_ob;		// object header - ALL objects MUST begin with this...
	double	omegah;
	double	sinhVXoKX;
	double	coshVXoKX;
	char	lowHIGH;
	short	poles;
	double	piPoles;
	double	piPoles2;
	double	*a;
	double	*b;
	double	*ta;
	double	*tb;
	void	**a_outlet;		//	*a_outlet[numCoeffs]
	void	**b_outlet;		//	*b_outlet[numCoeffs], [0] is never used
};
typedef struct t_chebCoef t_chebCoef;

void *chebCoef_class;		// global pointer to the object class - so max can reference the object 

// these are prototypes for the methods that are defined below
void *chebCoef_new	(long p, t_symbol *s);
void chebCoef_free	(t_chebCoef *x);
void chebCoef_assist(t_chebCoef *x, void *b, long m, long a, char *s);
void chebCoef_bang	(t_chebCoef *x);
void chebCoef_low	(t_chebCoef *x);
void chebCoef_high	(t_chebCoef *x);
void chebCoef_print	(t_chebCoef *x);
void chebCoef_cutoff(t_chebCoef *x, float c);
void chebCoef_cutoffFreq(t_chebCoef *x, long f);
void chebCoef_ripple(t_chebCoef *x, float r);
void chebCoef_calc	(t_chebCoef *x);
long strcmp		(char *a, char *b);

////////////////////////////////////////////////////////////////////////////////////////////////////
void main(void)
{
	setup((t_messlist **)&chebCoef_class, (method)chebCoef_new, (method)chebCoef_free, (short)sizeof(t_chebCoef), 0L, A_LONG, A_DEFSYM, 0);
	// setup() loads our external into Max's memory so it can be used in a patch
	
	addmess((method)chebCoef_assist, "assist", A_CANT, 0); // (optional) assistance method needs to be declared like this
	addbang((method)chebCoef_bang);			// the method it uses when it gets a bang in the left inlet
	addmess((method)chebCoef_low, "low", 0);
	addmess((method)chebCoef_high, "high", 0);
	addmess((method)chebCoef_print, "print", 0);
	addfloat((method)chebCoef_cutoff);		// the method for a float in the left inlet (inlet 0)
	addint((method)chebCoef_cutoffFreq);
	addftx((method)chebCoef_ripple, 1);		// the method for an int in the right inlet (inlet 1)
	
	post("cheb object loaded...",0);	// post any important info to the max window when our object is laoded
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void *chebCoef_new(long p, t_symbol *s)
{
	long i;
	
	t_chebCoef *x;				// local variable (pointer to a t_chebCoef data structure)
	x = (t_chebCoef *)newobject(chebCoef_class); // create a new instance of this object
	
	//	impose limits
	if ( p<0 || p>20 || p%2>0 )
	{
		post("WARNING: Number of poles must be even and between 2 and 20, inclusive.");
		post("   Setting value to 2.");
		p=2;
	}
	x->poles		= p;
	x->piPoles		= pi/p;
	x->piPoles2		= pi/(p*2.0);
	
	if ( strcmp(s->s_name, "low") && strcmp(s->s_name, "high") && s->s_name[0] )	//	zero means match
	{
		post("WARNING: Must specify either \"high\" or \"low\" pass.");
		post("   Setting value to \"low\".");
		x->lowHIGH = 0;
	}
	else
	{
		if (!strcmp(s->s_name, "low") || s->s_name[0] == 0 )
			x->lowHIGH = 0;
		else
			x->lowHIGH = 1;
	}
	
	x->a	= (double *)getbytes((p+3) * sizeof(double));
	x->b	= (double *)getbytes((p+3) * sizeof(double));
	x->ta	= (double *)getbytes((p+3) * sizeof(double));
	x->tb	= (double *)getbytes((p+3) * sizeof(double));
	x->a_outlet	= (void **)getbytes((p+1) * sizeof(void *));
	x->b_outlet	= (void **)getbytes((p+1) * sizeof(void *));
	
	//	create inlets and outlets
	floatin(x,1);				// create a second float inlet
	for (i=x->poles; i>=1; i--)	//	create outlets from right to left
	{
		x->b_outlet[i]	= floatout(x);
		x->a_outlet[i]	= floatout(x);
	}
	x->a_outlet[0]	= floatout(x);
	
	//	set default values
	x->omegah		= 0.0;		// set initial (default) left operand value in the instance's data structure
	x->sinhVXoKX	= 0.0;		// set initial (default) right operand value
	x->coshVXoKX	= 0.0;
	
	post(" new cheb object instance added to patch...",0); // post info to the max window when new instance is created
	
	return(x);					// return a reference to the object instance 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_free(t_chebCoef *x)
{
	freebytes(x->a, (x->poles+3) * sizeof(double));
	freebytes(x->b, (x->poles+3) * sizeof(double));
	freebytes(x->ta, (x->poles+3) * sizeof(double));
	freebytes(x->tb, (x->poles+3) * sizeof(double));
	freebytes(x->a_outlet, (x->poles+1) * sizeof(void *));
	freebytes(x->b_outlet, (x->poles+1) * sizeof(void *));
	post("chebCoef_free called",0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// arguments are always the same for assistance method
void chebCoef_assist(t_chebCoef *, void *, long m, long a, char *s)
{
	if (m == ASSIST_OUTLET)
	{
		if (a == 0)
			sprintf(s,"Coefficient output: a0.");
		else
		{
			a++;
			if (a%2 == 0)
				sprintf(s,"Coefficient output: a%d.", a/2);
			else
				sprintf(s,"Coefficient output: b%d.", a/2);
		}
	}
	else
	{
		switch (a)
		{	
			case 0:
			sprintf(s,"Cutoff: integer 0 to 22050 or float 0 to 0.5.");
			break;
			
			case 1:
			sprintf(s,"Percentage ripple. (zero to 29\%)");
			break;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_bang(t_chebCoef *x)		// x = reference to this instance of the object
{
	chebCoef_calc(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_low(t_chebCoef *x)
{
	x->lowHIGH = 0;
	chebCoef_calc(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_high(t_chebCoef *x)
{
	x->lowHIGH = 1;
	chebCoef_calc(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_print(t_chebCoef *x)
{
	long i;
	post("a[00] = % .10e", x->a[0]);
	for ( i=1; i <= x->poles; i++ )
		post("a[%02d] = % .10e   b[%02d] = % .10e", i, x->a[i], i, x->b[i]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_cutoffFreq(t_chebCoef *x, long f)
{
	chebCoef_cutoff(x, (double)f / sys_getsr());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_cutoff(t_chebCoef *x, float c)
{
	x->omegah	= c > 0.5 ? 0.5 : (c < 0.0 ? 0.0 : c);	//	use omegah at temporary storage
	x->omegah	*= pi;									//	now contains true omega/2 (omega-half)
	
	//post("W = %g, LH = %d", x->omegah*2, x->lowHIGH);
	
	chebCoef_calc(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_ripple(t_chebCoef *x, float r)
{
	double	epsInv, VX, KX;

	r	= r > 29.0 ? 29.0 : (r < 0.0 ? 0.0 : r );
	
	if ( r > 0.0 )
	{
		epsInv	= 1.0 / sqrt( pow(100.0/(100.0 - r), 2.0) - 1.0 );
		VX		= asinh(epsInv) / (double)x->poles;
		KX		= cosh( acosh(epsInv) / (double)x->poles );
		
		x->sinhVXoKX	= sinh(VX) / KX;
		x->coshVXoKX	= cosh(VX) / KX;
	}
	else
	{
		x->sinhVXoKX	= 0;
		x->coshVXoKX	= 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_calc(t_chebCoef *x)
{
	// holds the "a" & "b" coefficients upon program completion
	//double	a[numCoeffs], b[numCoeffs];
	// internal use for combining stages
	//double	ta[numCoeffs], tb[numCoeffs];
	
	double	K, KK;
	double	RP, IP, M, D, X0, X1, X2, Y1, Y2;
	double	A0, A1, A2, B1, B2, sa, sb, gain;
	long 	p, i;
	
	// INITIALIZE VARIABLES
	//for ( i=0; i < numCoeffs; i++ )
	for ( i=0; i < x->poles+3; i++ )
	{
		x->a[i] = 0.0;
		x->b[i] = 0.0;
	}

	x->a[2] = 1.0;
	x->b[2] = 1.0;
	
	// LP TO LP, or LP TO HP transform	(new calculation not needed when ripple changes)
	if ( x->lowHIGH )
		K = -cos(x->omegah + 0.5) / cos(x->omegah - 0.5);
	else
		K =  sin(0.5 - x->omegah) / sin(0.5 + x->omegah);
	
	KK = K * K;
		
	// LOOP FOR EACH POLE-PAIR
	for ( p=1; p <= x->poles/2; p++ )
	{
		// Calculate the pole location on the unit circle
		RP = -cos(x->piPoles2 + (p-1) * x->piPoles);
		IP =  sin(x->piPoles2 + (p-1) * x->piPoles);
		
		// Warp from a circle to an ellipse when ripple is greater than zero
		if ( x->sinhVXoKX != 0 && x->coshVXoKX != 0 )
		{
			RP *= x->sinhVXoKX;
			IP *= x->coshVXoKX;
		}
		
		//post("%d  RP = %g, IP = %g", p, RP, IP);
		
		// s-domain to z-domain conversion
		M	= RP*RP + IP*IP;
		D	= 4.0 - 4.0*RP*T + M*TT;
		X0	= TT/D;
		X1	= 2.0*X0;	//	X1	= 2.0*TT/D;
		X2	= X0;		//	X2	= TT/D;
		Y1	= (8.0 - 2.0*M*TT)/D;
		Y2	= (-4.0 - 4.0*RP*T - M*TT)/D;
		
		//post("%d  M = %g, D = %g, X0 = %g, X1 = %g, X2 = %g, Y1 = %g, Y2 = %g", p, M, D, X0, X1, X2, Y1, Y2);
		
		D = 1 + Y1*K - Y2*KK;
		
		A0	= (X0 - X1*K + X2*KK)/D;
		A1	= (-2*X0*K + X1 + X1*KK - 2*X2*K)/D;
		A2	= (X0*KK - X1*K + X2)/D;
		B1	= (2*K + Y1 + Y1*KK - 2*Y2*K)/D;
		B2	= (-KK - Y1*K + Y2)/D;
		
		if ( x->lowHIGH )
		{
			A1 = -A1;
			B1 = -B1;
		}
		
		//post("%d  A0 = %g, A1 = %g, A2 = %g, B1 = %g, B2 = %g", p, A0, A1, A2, B1, B2);
		
		// Add coefficients to the cascade
		//for ( i=0; i < numCoeffs; i++ )
		for ( i=0; i < x->poles+3; i++ )
		{
			x->ta[i] = x->a[i];
			x->tb[i] = x->b[i];
		}
		
		//for ( i=2; i < numCoeffs; i++ )
		for ( i=2; i < x->poles+3; i++ )
		{
			x->a[i] = A0*x->ta[i] + A1*x->ta[i-1] + A2*x->ta[i-2];
			x->b[i] = x->tb[i] - B1*x->tb[i-1] - B2*x->tb[i-2];
		}
	}
		
	// Finish combining coefficients
	x->b[2] = 0;
	//for ( i=0; i<numCoeffs-2; i++ )
	for ( i=0; i<x->poles+1; i++ )
	{
		x->a[i] = x->a[i+2];
		x->b[i] = -x->b[i+2];
	}
	
	// NORMALIZE THE GAIN
	sa = 0.0, sb = 0.0;
	//for ( i=0; i<numCoeffs-2; i++ )
	for ( i=0; i<x->poles+1; i++ )
	{
		if (x->lowHIGH)
		{
			sa += x->a[i] * pow(-1, i);
			sb += x->b[i] * pow(-1, i);
		}
		else
		{
			sa += x->a[i];
			sb += x->b[i];
		}
	}
	
	gain = sa / (1 - sb);
	
	//for ( i=0; i<numCoeffs-2; i++ )
	for ( i=0; i<x->poles+1; i++ )
		x->a[i] /= gain;
	
	//	send to outputs
	outlet_float(x->a_outlet[0], x->a[0]);
	for (i=1; i<=x->poles; i++)
	{
		outlet_float(x->a_outlet[i], x->a[i]);
		outlet_float(x->b_outlet[i], x->b[i]);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//	using this instead of string.h
//	keeps the compiled code smaller
long strcmp (char *a, char *b)
{
	long i=0;
	while (a[i] && b[i] && a[i] == b[i])
		i++;
	
	if (!a[i] && !b[i])
		return 0;
	else if (a[i] > b[i])
		return 1;
	else
		return -1;
}

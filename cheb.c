/**
*	Chebyshev coefficient generator
*	version 2004-08-11
*   
*   Copyright 2004 Reid A. Woodbury Jr.
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

#include "ext.h"  		// you must include this - it contains the external object's link to available Max functions
#include "z_dsp.h"
#include "ext_strings.h"

#include <math.h>

#ifndef	pi
//#define	pi		3.14159265358979323846d
#define	pi		3.1415926535897932384626433
#endif

//	double	T	= 2.0 * tan(0.5);
#define	T		1.0926049796875809683172064978862181305885
//	double	TT	= T * T;
#define	TT		1.1937856416380991930736854556016623973846

struct t_chebCoef	// defines our object's internal variables for each instance in a patch
{
	t_object p_ob;		// object header - ALL objects MUST begin with this...
	double	omegah;
	char	lowHIGH;
	char	poles;
	double	piPoles;
	double	piPoles2;
	float	ripple;
	double	sinhVXoKX;
	double	coshVXoKX;
	double	*a;
	double	*b;
	char	outOrder;		//	0 = abab, 1 = aabb
	void	*outlet;		//	list outlet
};
typedef struct t_chebCoef t_chebCoef;

void *chebCoef_class;		// global pointer to the object class - so max can reference the object 

// these are prototypes for the methods that are defined below
void *chebCoef_new	(t_symbol *s, long p, float r, t_symbol *o);
void chebCoef_free	(t_chebCoef *x);
void chebCoef_assist(t_chebCoef *x, void *b, long m, long a, char *s);
void chebCoef_bang	(t_chebCoef *x);
void chebCoef_low	(t_chebCoef *x);
void chebCoef_high	(t_chebCoef *x);
void chebCoef_aabab	(t_chebCoef *x);
void chebCoef_aaabb	(t_chebCoef *x);
void chebCoef_print	(t_chebCoef *x);
void chebCoef_cutoff(t_chebCoef *x, float c);
void chebCoef_cutoffInt(t_chebCoef *x, long l);
void chebCoef_ripple(t_chebCoef *x, float r);
void chebCoef_poles	(t_chebCoef *x, long p);

void chebCoef_rippleCalc	(t_chebCoef *x);
void chebCoef_calculate		(t_chebCoef *x);
void chebCoef_getPointers	(t_chebCoef *x);
void chebCoef_releasePtrs	(t_chebCoef *x);

////////////////////////////////////////////////////////////////////////////////////////////////////
void main(void)
{
	setup((t_messlist **)&chebCoef_class, (method)chebCoef_new, (method)chebCoef_free, (short)sizeof(t_chebCoef), 0L,
															A_DEFSYM, A_DEFLONG, A_DEFFLOAT, A_DEFSYM, 0);
	
	addmess((method)chebCoef_assist, "assist", A_CANT, 0);
	addbang((method)chebCoef_bang);
	addmess((method)chebCoef_low, "low", 0);
	addmess((method)chebCoef_high, "high", 0);
	addmess((method)chebCoef_print, "print", 0);
	addmess((method)chebCoef_aabab, "aabab", 0);
	addmess((method)chebCoef_aaabb, "aaabb", 0);
	addfloat((method)chebCoef_cutoff);		// the method for a float in the left inlet (inlet 0)
	addint((method)chebCoef_cutoffInt);
	addinx((method)chebCoef_poles, 1);
	addftx((method)chebCoef_ripple, 2);
	
	post("cheb object loaded...",0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void *chebCoef_new(t_symbol *s, long p, float r, t_symbol *o)
{
	t_chebCoef *x;				// local variable (pointer to a t_chebCoef data structure)
	x = (t_chebCoef *)newobject(chebCoef_class); // create a new instance of this object
	
	//	create list outlet
	x->outlet	= listout(x);
	
	//	create other inlets, HIGHEST TO LOWEST or RIGHT TO LEFT
	floatin(x,2);	//	ripple
	intin(x,1);		//	poles
	
    //	post message
	post("cheb [low|high] [#poles] [(float)%%ripple] [aabab|aaabb]");
	
	////////////////////	impose limits	///////////////////////////////
	//	high or low pass
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
	
	//	number of poles
	p = p==0 ? 2 : p;	//	default of zero is always changed to 2, then won't be considered an error
	if ( p<2 || p>20 || p%2>0 )
	{
		post("WARNING: Number of poles must be even and between 2 and 20, inclusive.");
		post("   Setting number of poles to 2.");
		p=2;
	}
	x->poles	= p;
	x->piPoles	= pi/p;
	x->piPoles2	= pi/(p*2.0);
	
	//	percentage ripple
	if ( r > 29.0 || r < 0.0 )
	{
		post("WARNING: Ripple value must be even and between 0 and 29%.");
		post("   Setting ripple percentage to zero (Butterworth).");
		r = 0.0;
	}
	x->ripple	= r;
	chebCoef_rippleCalc(x);
	
	//	output order
	if ( strcmp(o->s_name, "aabab") && strcmp(o->s_name, "aaabb") && o->s_name[0] )	//	zero means match
	{
		post("WARNING: Must specify an output order of either \"aabab\" or \"aaabb\".");
		post("   Setting value to \"aabab\".");
		x->outOrder = 0;
	}
	else
	{
		if (!strcmp(o->s_name, "abab") || !o->s_name[0] )
			x->outOrder = 0;
		else
			x->outOrder = 1;
	}
	
	//	start with result pointers == zero
	x->a = 0L; x->b = 0L;
	
	//	set up space for results; this depends on the number of poles
	chebCoef_getPointers(x);
	
	//	set default values
	x->omegah		= 0.03926990816987241d;	//	aprox 1100Hz at 44.1kHz
	
	return(x);					// return a reference to the object instance 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_free(t_chebCoef *x)
{
	chebCoef_releasePtrs(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_assist(t_chebCoef *, void *, long m, long a, char *s) // arguments are always the same for assistance method
{
	if (m == ASSIST_OUTLET)
	{
		sprintf(s,"Coefficient output.");
	}
	else
	{
		switch (a)
		{	
			case 0:
			sprintf(s,"Cutoff Frequency");
			break;
			
			case 1:
			sprintf(s,"Number of poles (2-20)");
			break;
			
			case 2:
			sprintf(s,"Percentage ripple (0-29)");
			break;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_bang(t_chebCoef *x)		// x = reference to this instance of the object
{
	chebCoef_calculate(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_low(t_chebCoef *x)
{
	x->lowHIGH = 0;
	chebCoef_calculate(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_high(t_chebCoef *x)
{
	x->lowHIGH = 1;
	chebCoef_calculate(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_aabab(t_chebCoef *x)
{
	x->outOrder = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_aaabb(t_chebCoef *x)
{
	x->outOrder = 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_print(t_chebCoef *x)
{
	long i;
	post("Fc = % .4f,  # poles = %d,  %% ripple = % .2f", x->omegah/pi, x->poles, x->ripple);
	post("a[00] = % .15e", x->a[0]);
	for ( i=1; i <= x->poles; i++ )
		post("a[%02d] = % .15e   b[%02d] = % .15e", i, x->a[i], i, x->b[i]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_cutoffInt(t_chebCoef *x, long l)
{
	chebCoef_cutoff(x, (double)l);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_cutoff(t_chebCoef *x, float c)
{
	//	change to fraction of the sample rate
	c	/=	sys_getsr();
	
	//	check the range and mulitply by pi (not 2¹); omegah contains true omega/2 (omega-half)
	x->omegah	= ( c > 0.5 ? 0.5 : (c < 0.0 ? 0.0 : c) ) * pi;
	
	//post("Fc = %g, LH = %d", x->omegah/pi, x->lowHIGH);
	
	chebCoef_calculate(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_ripple(t_chebCoef *x, float r)
{
	x->ripple	= r > 29.0 ? 29.0 : (r < 0.0 ? 0.0 : r );
	
	chebCoef_rippleCalc(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_poles(t_chebCoef *x, long p)
{
	p	= p > 20 ? 20 : (p < 0 ? 0 : p );
#pragma optimization_level 0
	p /= 2;
	p *= 2;
#pragma optimization_level reset
	
	chebCoef_releasePtrs	(x);	//	number of a and b values may have changed
	
	x->poles	= p;
	x->piPoles	= PI/p;
	x->piPoles2	= PI/(p*2.0);
	
	chebCoef_rippleCalc(x);
	chebCoef_getPointers	(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_rippleCalc(t_chebCoef *x)
{
	double	ESinv, VX, KX;
	
	if ( x->ripple > 0.0 )
	{
		ESinv	= 1.0 / sqrt( pow(100.0/(100.0 - (double)x->ripple), 2.0) - 1.0 );
		VX		= asinh(ESinv) / (double)x->poles;
		KX		= cosh( acosh(ESinv) / (double)x->poles );
		
		x->sinhVXoKX	= sinh(VX) / KX;
		x->coshVXoKX	= cosh(VX) / KX;
	}
	else
	{
		x->sinhVXoKX	= 0.0;
		x->coshVXoKX	= 0.0;
	}
	
	//post("ES = %g, VX = %g, KX = %g", 1/ESinv, VX, KX);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void chebCoef_calculate(t_chebCoef *x)
{
	// holds the "a" & "b" coefficients upon program completion
	//double	a[numCoeffs], b[numCoeffs];
	// internal use for combining stages
	//double	ta[numCoeffs], tb[numCoeffs];
	double	*ta	= (double *)getbytes((x->poles+3) * sizeof(double));
	double	*tb	= (double *)getbytes((x->poles+3) * sizeof(double));
	
	t_atom	*list;
	
	double	K, KK;
	double	RP, IP, M, D, X0, X1, X2, Y1, Y2;
	double	A0, A1, A2, B1, B2, sa, sb, gain;
	long 	p, i;
	
	// INITIALIZE VARIABLES
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
		// chebCoef_calculate the pole location on the unit circle
		RP = -cos(x->piPoles2 + (p-1) * x->piPoles);
		IP =  sin(x->piPoles2 + (p-1) * x->piPoles);
		
		// Warp from a circle to an ellipse when ripple is greater than zero
		if ( x->ripple > 0 )
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
		
		//post("%d  w = %g, M = %g, D = %g, X0 = %g, X1 = %g, X2 = %g, Y1 = %g, Y2 = %g", p, x->omegah*2, M, D, X0, X1, X2, Y1, Y2);
		
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
		for ( i=0; i < x->poles+3; i++ )
		{
			ta[i] = x->a[i];
			tb[i] = x->b[i];
		}
		
		for ( i=2; i < x->poles+3; i++ )
		{
			x->a[i] = A0*ta[i] + A1*ta[i-1] + A2*ta[i-2];
			x->b[i] = tb[i] - B1*tb[i-1] - B2*tb[i-2];
		}
	}
	freebytes(ta, (x->poles+3) * sizeof(double));
	freebytes(tb, (x->poles+3) * sizeof(double));
		
	// Finish combining coefficients
	x->b[2] = 0;
	for ( i=0; i<x->poles+1; i++ )
	{
		x->a[i] = x->a[i+2];
		x->b[i] = -x->b[i+2];
	}
	
	// NORMALIZE THE GAIN
	sa = 0.0, sb = 0.0;
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
	
	for ( i=0; i<x->poles+1; i++ )
		x->a[i] /= gain;
	
	
	//	send to outputs
	list	= (t_atom *)getbytes((x->poles*2+1) * sizeof(t_atom));
	SETFLOAT(list, x->a[0]);
	if (x->outOrder)
	{
		for (p=1; p<=x->poles; p++)
		{
			SETFLOAT(list+(p), x->a[p]);
			SETFLOAT(list+(x->poles+p), x->b[p]);
		}
	}
	else
	{
		for (p=1; p<=x->poles; p++)
		{
			SETFLOAT(list+(2*p-1), x->a[p]);
			SETFLOAT(list+(2*p), x->b[p]);
		}
	}
	
	outlet_list(x->outlet, 0L, x->poles*2+1, list);
	
	freebytes(list, (x->poles*2+1) * sizeof(t_atom));
}

///////////////////////////////////////////////
//	get and free memory functions
void chebCoef_getPointers	(t_chebCoef *x)
{
	if (!x->a && !x->b)	//	all must be zero
	{
		x->a	= (double *)getbytes((x->poles+3) * sizeof(double));
		x->b	= (double *)getbytes((x->poles+3) * sizeof(double));
	}
	else
		error("chebCoef_getPointers; one pointer was not zero.");
}

void chebCoef_releasePtrs	(t_chebCoef *x)
{
	if (x->a && x->b)	//	all must != zero
	{
		freebytes(x->a, (x->poles+3) * sizeof(double));
		freebytes(x->b, (x->poles+3) * sizeof(double));
		
		x->a = 0L; x->b = 0L;;
	}
	else
		error("chebCoef_releasePtrs; one pointer was already zero.");
}

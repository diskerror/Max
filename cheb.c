/**
*	Chebyshev coefficient generator
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

#include "ext.h"				// standard Max include, always required
#include "ext_obex.h"			// required for new style Max object
#include "z_dsp.h"				//	for sys_getsr(), t_double, t_float, t_vptr
#include "ext_strings.h"

#include <math.h>

#ifndef	pi
#define	pi		3.1415926535897932384626433
#endif

//	double	T	= 2.0 * tan(0.5);
#define	T		1.0926049796875809683172064978862181305885
//	double	TT	= T * T;
#define	TT		1.1937856416380991930736854556016623973846

#define MAX_CHEB_POLES	20

typedef struct _cheb
{
	t_object	p_ob;		// object header - ALL objects MUST begin with this...
	t_double	omegah;
	t_uint8		lowHIGH;
	t_uint8		poles;
	t_double	piPoles;
	t_double	piPoles2;
	t_double	ripple;
	t_double	sinhVXoKX;
	t_double	coshVXoKX;
	t_double	*a;
	t_double	*b;
	t_uint8		outOrder;	//	0 = abab, 1 = aabb
	t_vptr		outlet;		//	list outlet
} t_cheb;

void *cheb_class;

//// standard set
void *cheb_new(t_symbol *s, long argc, t_atom *argv);
void cheb_free(t_cheb *x);
void cheb_assist(t_cheb *x, void *b, long m, long a, char *s);

void cheb_bang(t_cheb *x);
void cheb_low(t_cheb *x);
void cheb_high(t_cheb *x);
void cheb_aabab(t_cheb *x);
void cheb_aaabb(t_cheb *x);
void cheb_print(t_cheb *x);
void cheb_cutoff(t_cheb *x, double c);
void cheb_cutoffInt(t_cheb *x, long l);
void cheb_poles(t_cheb *x, long p);
void cheb_ripple(t_cheb *x, double r);

void cheb_rippleCalc(t_cheb *x);
void cheb_calculate(t_cheb *x);
void cheb_getPointers(t_cheb *x);
void cheb_releasePtrs(t_cheb *x);

////////////////////////////////////////////////////////////////////////////////////////////////////
int C74_EXPORT main(void)
{	
	t_class *c;
	
	c = class_new("cheb", (method)cheb_new, (method)cheb_free, (long)sizeof(t_cheb), 0L, A_DEFSYM, A_DEFLONG, A_DEFFLOAT, A_DEFSYM, 0);
    class_addmethod(c, (method)cheb_assist, "assist",	A_CANT, 0);  
	
    class_addmethod(c, (method)cheb_bang, "bang", 0);
    class_addmethod(c, (method)cheb_low, "low", 0);
    class_addmethod(c, (method)cheb_high, "high", 0);
    class_addmethod(c, (method)cheb_print, "print", 0);
    class_addmethod(c, (method)cheb_aabab, "aabab", 0);
    class_addmethod(c, (method)cheb_aaabb, "aaabb", 0);
    class_addmethod(c, (method)cheb_cutoff, "float", A_FLOAT, 0);	// the method for a float in the left inlet (inlet 0)
    class_addmethod(c, (method)cheb_cutoffInt, "int", A_LONG, 0); 	// the method for a integer in the left inlet (inlet 0)
    class_addmethod(c, (method)cheb_poles, "in1", A_DEFLONG, 0);
    class_addmethod(c, (method)cheb_ripple, "ft2", A_DEFFLOAT, 0);  
	
	class_register(CLASS_BOX, c);
	cheb_class = c;

//	post("cheb object loaded...",0);
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//	cheb_new(t_symbol *s, long p, float r, t_symbol *o)
void *cheb_new(t_symbol *s, long argc, t_atom *argv)
{
	long p;
	t_float r;
	
	t_cheb *x = NULL;
    
	if ( (x = (t_cheb *)object_alloc(cheb_class)) )
    {
		//	create other inlets, HIGHEST TO LOWEST or RIGHT TO LEFT
		floatin(x, 2);	//	ripple
		intin(x, 1);	//	poles
	
		//	create list outlet
		x->outlet = listout(x);
	
		//	post message
		post("cheb [low|high] [#poles] [(float)%%ripple] [aabab|aaabb]");
	
		////////////////////	impose limits	///////////////////////////////
		//	high or low pass
		if( argc && !strcmp(atom_getsym(argv)->s_name, "high") )
			x->lowHIGH = 1;
		else
			x->lowHIGH = 0;
		
		//	number of poles
		x->poles = 2;
		if( argc > 1 ) 
		{
			p = atom_getlong(argv+1);
			
			//	floor() to the lower even number (ignore last bit)
			x->poles = ( p > MAX_CHEB_POLES ? MAX_CHEB_POLES : (p < 0 ? 0 : p ) ) & 0xFFFFFFFE;
		}
		x->piPoles	= pi/x->poles;
		x->piPoles2	= pi/(x->poles*2.0);
	
		//	start with result pointers == zero
		x->a = x->b = 0L;
	
		//	set up space for results; this depends on the number of poles
		cheb_getPointers(x);
	
		//	percentage ripple
		x->ripple = 0.0;
		if( argc > 2 ) 
		{
            r = atom_getfloat(argv+2);
			if ( r > 0.0 && r <= 29.0 )
			{
				x->ripple = r;
			}
		}
		cheb_rippleCalc(x);
	
		//	output order
		if( argc > 3 && !strcmp(atom_getsym(argv+3)->s_name, "aaabb") )
			x->outOrder = 1;
		else
			x->outOrder = 0;
	
		//	set default values
		x->omegah		= 0.03926990816987241;	//	aprox 1100Hz at 44.1kHz
	
		cheb_calculate(x);
	}
	return x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_free(t_cheb *x)
{
	cheb_releasePtrs(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_assist(t_cheb *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_OUTLET)
	{
		sprintf(s,"Coefficient output (list).");
	}
	else
	{
		switch (a)
		{	
			case 0:
			sprintf(s,"Cutoff Frequency (Hz)");
			break;
			
			case 1:
			sprintf(s,"Number of poles (2-20)");
			break;
			
			case 2:
			sprintf(s,"Percentage ripple (0-29%%)");
			break;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_bang(t_cheb *x)		// x = reference to this instance of the object
{
	long 	p;
	t_atom	*list;
	
	//	send to outputs
	list	= (t_atom *) getbytes((x->poles*2+1) * sizeof(t_atom));
	atom_setfloat(list, x->a[0]);
	if (x->outOrder)
	{
		for (p=1; p<=x->poles; p++)
		{
			atom_setfloat(list+(p), x->a[p]);
			atom_setfloat(list+(x->poles+p), x->b[p]);
		}
	}
	else
	{
		for (p=1; p<=x->poles; p++)
		{
			atom_setfloat(list+(2*p-1), x->a[p]);
			atom_setfloat(list+(2*p), x->b[p]);
		}
	}
	
	outlet_list(x->outlet, 0L, x->poles*2+1, list);
	
	freebytes(list, (x->poles*2+1) * sizeof(t_atom));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_low(t_cheb *x)
{
	x->lowHIGH = 0;
	cheb_calculate(x);
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_high(t_cheb *x)
{
	x->lowHIGH = 1;
	cheb_calculate(x);
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_aabab(t_cheb *x)
{
	x->outOrder = 0;
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_aaabb(t_cheb *x)
{
	x->outOrder = 1;
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_print(t_cheb *x)
{
	long i;
	post("Fc = % .4f,  # poles = %d,  %% ripple = %.2f", x->omegah/pi, x->poles, x->ripple);
	post("a[00] = % .15e", x->a[0]);
	for ( i=1; i <= x->poles; i++ )
		post("a[%02d] = % .15e   b[%02d] = % .15e", i, x->a[i], i, x->b[i]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_cutoffInt(t_cheb *x, long l)
{
	cheb_cutoff(x, (float)l);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_cutoff(t_cheb *x, double c)
{
	//	change to fraction of the sample rate
	c	/=	sys_getsr();
	
	//	check the range and mulitply by pi (not 2¹); omegah contains true omega/2 (omega-half)
	x->omegah	= ( c > 0.5 ? 0.5 : (c < 0.0 ? 0.0 : c) ) * pi;
	
// 	post("Fc = %g, LH = %d", x->omegah/pi, x->lowHIGH);
	
	cheb_calculate(x);
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_ripple(t_cheb *x, double r)
{
	x->ripple = (r>29.0) ? 29.0 : ((r<0.0) ? 0.0 : r);
	
	cheb_rippleCalc(x);
	cheb_calculate(x);
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_poles(t_cheb *x, long p)
{
	p = ( (p > MAX_CHEB_POLES) ? MAX_CHEB_POLES : ((p < 0) ? 0 : p ) ) & 0xFFFFFFFE;
	
	cheb_releasePtrs(x);	//	count of a and b values may change

	x->poles	= p;
	x->piPoles	= pi/p;
	x->piPoles2	= pi/(p*2.0);
	
	cheb_getPointers	(x);
	
	cheb_rippleCalc(x);
	cheb_calculate(x);
	cheb_bang(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_rippleCalc(t_cheb *x)
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
	
// 	post("ES = %g, VX = %g, KX = %g", 1/ESinv, VX, KX);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_calculate(t_cheb *x)
{
	// holds the "a" & "b" coefficients upon program completion
	//double	a[numCoeffs], b[numCoeffs];
	// internal use for combining stages
	//double	ta[numCoeffs], tb[numCoeffs];
	double	*ta	= (double *)getbytes((x->poles+3) * sizeof(double));
	double	*tb	= (double *)getbytes((x->poles+3) * sizeof(double));
	
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
		// cheb_calculate the pole location on the unit circle
		RP = -cos(x->piPoles2 + (p-1) * x->piPoles);
		IP =  sin(x->piPoles2 + (p-1) * x->piPoles);
		
		// Warp from a circle to an ellipse when ripple is greater than zero
		if ( x->ripple > 0 )
		{
			RP *= x->sinhVXoKX;
			IP *= x->coshVXoKX;
		}
		
// 		post("%d  RP = %g, IP = %g", p, RP, IP);
		
		// s-domain to z-domain conversion
		M	= RP*RP + IP*IP;
		D	= 4.0 - 4.0*RP*T + M*TT;
		X0	= TT/D;
		X1	= 2.0*X0;	//	X1	= 2.0*TT/D;
		X2	= X0;		//	X2	= TT/D;
		Y1	= (8.0 - 2.0*M*TT)/D;
		Y2	= (-4.0 - 4.0*RP*T - M*TT)/D;
		
// 		post("%d  w = %g, M = %g, D = %g, X0 = %g, X1 = %g, X2 = %g, Y1 = %g, Y2 = %g", p, x->omegah*2, M, D, X0, X1, X2, Y1, Y2);
		
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
		
// 		post("%d  A0 = %g, A1 = %g, A2 = %g, B1 = %g, B2 = %g", p, A0, A1, A2, B1, B2);
		
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
}

///////////////////////////////////////////////
//	get and free memory functions
void cheb_getPointers	(t_cheb *x)
{
	if (!x->a && !x->b)	//	all must be zero
	{
		x->a	= (double *)getbytes((x->poles+3) * sizeof(double));
		x->b	= (double *)getbytes((x->poles+3) * sizeof(double));
	}
	else
		error("cheb_getPointers; one pointer was not zero.");
}

void cheb_releasePtrs	(t_cheb *x)
{
	if (x->a && x->b)	//	all must != zero
	{
		freebytes(x->a, (x->poles+3) * sizeof(double));
		freebytes(x->b, (x->poles+3) * sizeof(double));
		
		x->a = x->b = 0L;
	}
	else
		error("cheb_releasePtrs; one pointer was already zero.");
}

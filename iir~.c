/**
*	This will do a IIR or recursive filter based on an input list of coefficients.
*	version 2004-09-11
*   
*   Copyright 2004 Reid A. Woodbury Jr.
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

#include "ext.h"
#include "ext_obex.h"
#include "ext_strings.h"
#include "z_dsp.h"
#include <math.h>

void *iir_class;

typedef struct
{
	t_pxobject l_obj;
	t_uint8 poles;		//	number of poles
	t_uint8 inOrder;	//	order of coefficients
	t_double *a, *b;	//	coefficients from input list
	t_double *x, *y;	//	delayed input and output values
} t_iir;

void *iir_new(t_symbol *o, long argc, t_atom *argv);
void iir_free(t_iir *x);
void iir_assist(t_iir *x, void *b, long m, long a, char *s);

void iir_aabab(t_iir *x);
void iir_aaabb(t_iir *x);
void iir_print(t_iir *x);
void iir_dsp(t_iir *x, t_signal **sp, short *count);
void iir_dsp64(t_iir *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
t_int *iir_perform(t_int *w);
void iir_perform64(t_iir *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void iir_clear(t_iir *x);
void iir_coeffs(t_iir *x, t_symbol *, short argc, t_atom *argv);
void iir_releasePtrs(t_iir *x);

int C74_EXPORT main(void)
{
	t_class *c;
	c = class_new("iir~", (method)iir_new, (method)iir_free, sizeof(t_iir), 0L, A_DEFSYM, 0);
	class_addmethod(c, (method)iir_assist, "assist", A_CANT, 0);

	class_addmethod(c, (method)iir_dsp, "dsp", A_CANT, 0);
	class_addmethod(c, (method)iir_dsp64, "dsp64", A_CANT, 0);
	class_addmethod(c, (method)iir_clear, "clear", 0);
	class_addmethod(c, (method)iir_aabab, "aabab", 0);
	class_addmethod(c, (method)iir_aaabb, "aaabb", 0);
	class_addmethod(c, (method)iir_print, "print", 0);
	class_addmethod(c, (method)iir_coeffs, "list", A_GIMME, 0);
	
	class_dspinit(c);
	class_register(CLASS_BOX, c);
	
	iir_class = c;
	
	return 0;
}

void *iir_new(t_symbol *o, long argc, t_atom *argv)
{
    t_iir *x = NULL;
    
    if( (x = (t_iir *)object_alloc(iir_class)) )
    {
		dsp_setup((t_pxobject *)x, 1);
	
		// one signal outlet
		outlet_new((t_object *)x, "signal");
	
		//	post message
		post("iir~ [aabab|aaabb]");
	
		//	output order
		if( argc && !strcmp(atom_getsym(argv)->s_name, "aaabb") )
			x->inOrder = 1;
		else
			x->inOrder = 0;

		//	start with all pointers == null
		x->a = x->b = x->x = x->y = NULL;
	}

    return (x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_free(t_iir *x)
{
	iir_releasePtrs(x);
	dsp_free((t_pxobject *)x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_assist(t_iir *x, void *b, long m, long a, char *s)
{
	if (m == 2)
		sprintf(s,"(signal) Output");
	else
	{
		sprintf(s,"(signal) Input, List Input");
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_aabab(t_iir *x)
{
	x->inOrder = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_aaabb(t_iir *x)
{
	x->inOrder = 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_print(t_iir *x)
{
	long i;
	if (x->a)
	{
		post("a[00] = % .15e", x->a[0]);
		for ( i=1; i <= x->poles; i++ )
			post("a[%02d] = % .15e   b[%02d] = % .15e", i, x->a[i], i, x->b[i]);
	}
	else
		post("a[00] = 0.0");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_dsp(t_iir *x, t_signal **sp, short *count)
{
	iir_clear(x);
	dsp_add(iir_perform, 6, sp[0]->s_vec, sp[3]->s_vec, x, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}


void iir_dsp64(t_iir *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	iir_clear(x);
	dsp_add64(dsp64, (t_object*)x, (t_perfroutine64)iir_perform64, 0, NULL);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
t_int *iir_perform(t_int *w)
{
	// assign from parameters
    t_float *in = (t_float *) w[1];
    t_float *out = (t_float *) w[2];
    t_iir *x = (t_iir *) w[3];
    long n = (long) w[4];
    register long p;
    
    if (x->l_obj.z_disabled)
		return (w+5);
	
	// DSP loops
	if (x->a && x->b && x->x && x->y)
	{
		while (n--)
		{
			//	delay old values one sample
			for (p=x->poles; p>0; p--)
			{
				x->x[p] = x->x[p-1];
				x->y[p] = x->y[p-1];
			}
			x->x[0] = (double)*in++;
			
			x->y[0] = x->x[0] * x->a[0];
			for (p=1; p<=x->poles; p++)
			{
				x->y[0] += x->x[p] * x->a[p];
				x->y[0] += x->y[p] * x->b[p];
			}
			
			*out++ = (float)x->y[0];
		}
	}
	else	//	if pointers are no good...
	{
		while (n--)
			*out++ = *in++;	//	...just copy input to output
	}
	
    return (w+5);
}


void iir_perform64(t_iir *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double *in = ins[0];
    t_double *out = outs[0];
    long n = sampleframes;
    register long p;
    
    if (x->l_obj.z_disabled)
		return;
	
	// DSP loops
	if (x->a && x->b && x->x && x->y)
	{
		while (n--)
		{
			//	delay old values one sample
			for (p=x->poles; p>0; p--)
			{
				x->x[p] = x->x[p-1];
				x->y[p] = x->y[p-1];
			}
			x->x[0] = *in++;
			
			x->y[0] = x->x[0] * x->a[0];
			for (p=1; p<=x->poles; p++)
			{
				x->y[0] += x->x[p] * x->a[p];
				x->y[0] += x->y[p] * x->b[p];
			}
			
			*out++ 	= x->y[0];
		}
	}
	else	//	if pointers are no good...
	{
		while (n--)
			*out++ = *in++;	//	...just copy input to output
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_clear(t_iir *x)
{
	long p;
	for ( p=1; p<=x->poles; p++ )
		x->y[p] = 0.0;	//	y[0] is not used
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_coeffs	(t_iir *x, t_symbol *s, short argc, t_atom *argv)
{
	long i, p;
	
	for(i=0; i<argc; i++)
	{
		if (argv[i].a_type != A_FLOAT)
		{
			post("WARNING: All list members must be of type float or double.");
			return;
		}
	}
	
	iir_releasePtrs(x);	//	does nothing if pointer are null
	x->poles = argc/2;	//	integer division
	
	//	now we need pointers for our new data
	x->a	= (double *)getbytes((x->poles+1) * sizeof(double));
	x->b	= (double *)getbytes((x->poles+1) * sizeof(double));
	x->x	= (double *)getbytes((x->poles+1) * sizeof(double));
	x->y	= (double *)getbytes((x->poles+1) * sizeof(double));
	
	//	clear the x and y space. Make them start as silence.
	for ( p=0; p<=x->poles; p++ )
		x->x[p] = x->y[p] = 0.0;

	x->a[0] = argv[0].a_w.w_float;	//	the first is always the same no matter the order
	if (x->inOrder)	//	aaabb
	{
		for (p=1; p<=x->poles; p++)
		{
			x->a[p] = (double)argv[p].a_w.w_float;
			x->b[p] = (double)argv[x->poles+p].a_w.w_float;
		}
	}
	else				//	aabab
	{
		for (p=1; p<=x->poles; p++)
		{
			x->a[p] = (double)argv[p*2-1].a_w.w_float;
			x->b[p] = (double)argv[p*2].a_w.w_float;
		}
	}
	
	if (!(argc % 2))	//	if the list is short...
		x->b[x->poles] = 0.0;	//	...then there's bogus data in the last b
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_releasePtrs(t_iir *x)
{
	if (x->a && x->b && x->x && x->y)	//	all must != zero
	{
		freebytes(x->a, (x->poles+1) * sizeof(double));
		freebytes(x->b, (x->poles+1) * sizeof(double));
		freebytes(x->x, (x->poles+1) * sizeof(double));
		freebytes(x->y, (x->poles+1) * sizeof(double));
		
		x->a = x->b = x->x = x->y = 0L;
	}
	//else
	//	error("iir_releasePtrs; one pointer was already zero.");
	
	//not an error
	//	assumes pointers were never allocated
}

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

#define IIR_MAX_POLES		20
#define IIR_COEF_MEM_SIZE	((IIR_MAX_POLES+1) * sizeof(double))

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
    long c;
    
    if( (x = (t_iir *)object_alloc(iir_class)) )
    {
		dsp_setup((t_pxobject *)x, 1);
	
		// one signal outlet
		outlet_new((t_object *)x, "signal");
	
		//	post message
		object_post((t_object *)x, "iir~ [aabab|aaabb]");
	
		//	output order
		if( argc && !strcmp(atom_getsym(argv)->s_name, "aaabb") )
			x->inOrder = 1;
		else
			x->inOrder = 0;
		
		//	now we need pointers for our new data
		x->a = (double *)getbytes(IIR_COEF_MEM_SIZE);	//	input coefs
		x->b = (double *)getbytes(IIR_COEF_MEM_SIZE);	//	output coefs
		x->x = (double *)getbytes(IIR_COEF_MEM_SIZE);	//	delayed input
		x->y = (double *)getbytes(IIR_COEF_MEM_SIZE);	//	delayed output
		
		for ( c=0; c<=IIR_MAX_POLES; c++ )
			x->x[c] = x->y[c] = 0.0;
		
		if (!x->a || !x->b || !x->x || !x->y)
			object_error((t_object *)x, "BAD INIT POINTER");
		
		x->poles = 0.0;
	}

    return (x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_free(t_iir *x)
{
	if(x->a) freebytes(x->a, IIR_COEF_MEM_SIZE);
	if(x->b) freebytes(x->b, IIR_COEF_MEM_SIZE);
	if(x->x) freebytes(x->x, IIR_COEF_MEM_SIZE);
	if(x->y) freebytes(x->y, IIR_COEF_MEM_SIZE);
	
	dsp_free((t_pxobject *)x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_assist(t_iir *x, void *b, long m, long a, char *s)
{
	if (m == 2)
		sprintf(s,"(signal) Output");
	else
		sprintf(s,"(signal) Input, List Input");
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
		object_post((t_object *)x, "a[00] = % .15e", x->a[0]);
		for ( i=1; i <= x->poles; i++ )
			object_post((t_object *)x, "a[%02d] = % .15e   b[%02d] = % .15e", i, x->a[i], i, x->b[i]);
	}
	else
		object_post((t_object *)x, "a[00] = 0.0");
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
    register t_double y0;	//	de-referenced "x->y[0]"
    t_iir *x = (t_iir *) w[3];
    long n = (long) w[4];
    register long p;
    
    if (x->l_obj.z_disabled)
		return (w+5);
	
	// DSP loops
	if (x->poles && x->a && x->b && x->x && x->y) {
		while (n--) {
			//	delay old values one sample
			for (p=x->poles; p>0; p--) {
				x->x[p] = x->x[p-1];
				x->y[p] = x->y[p-1];
			}
			x->x[0] = (double)*in++;
			
			y0 = x->x[0] * x->a[0];
			for (p=1; p<=x->poles; p++) {
				y0 += x->x[p] * x->a[p];
				y0 += x->y[p] * x->b[p];
			}
			
			*out++ = x->y[0] = (t_float)y0;
		}
	}
	else {	//	if pointers are no good...
		while (n--)
			*out++ = *in++;	//	...just copy input to output
	}
	
    return (w+5);
}


void iir_perform64(t_iir *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double *in = ins[0];
    t_double *out = outs[0];
    register t_double y0;	//	de-referenced "x->y[0]"
    long n = sampleframes;
    register long p;
    
    if (x->l_obj.z_disabled)
		return;
	
	// DSP loops
	if (x->poles && x->a && x->b && x->x && x->y) {
		while (n--) {
			//	delay old values one sample
			for (p=x->poles; p>0; p--) {
				x->x[p] = x->x[p-1];
				x->y[p] = x->y[p-1];
			}
			x->x[0] = *in++;
			
			y0 = x->x[0] * x->a[0];
			for (p=1; p<=x->poles; p++) {
				y0 += x->x[p] * x->a[p];
				y0 += x->y[p] * x->b[p];
			}
			
			*out++ = x->y[0] = y0;
		}
	}
	else {	//	if pointers are no good...
		while (n--)
			*out++ = *in++;	//	...just copy input to output
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_clear(t_iir *x)
{
	long c;
	for ( c=0; c<=IIR_MAX_POLES; c++ )
		x->y[c] = 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_coeffs	(t_iir *x, t_symbol *s, short argc, t_atom *argv)
{
	long i, p, poles;
	
	for(i=0; i<argc; i++) {
		if (argv[i].a_type != A_FLOAT) {
			object_post((t_object *)x, "WARNING: All list members must be of type float or double.");
			x->poles = 0.0;
			return;
		}
	}
	
	poles = argc/2;	//	integer division
	
	//	clear the x and y space. Make them start as silence.
// 	for ( p=0; p<=poles; p++ )
// 		x->x[p] = x->y[p] = 0.0;

	x->a[0] = argv[0].a_w.w_float;	//	the first is always the same no matter the order
	if (x->inOrder) {	//	aaabb
		for (p=1; p<=poles && p<=IIR_MAX_POLES; p++) {
			x->a[p] = (double)argv[p].a_w.w_float;
			x->b[p] = (double)argv[poles+p].a_w.w_float;
		}
	}
	else {				//	aabab
		for (p=1; p<=poles && p<=IIR_MAX_POLES; p++) {
			x->a[p] = (double)argv[p*2-1].a_w.w_float;
			x->b[p] = (double)argv[p*2].a_w.w_float;
		}
	}
	
	if( poles > IIR_MAX_POLES )
		x->poles = IIR_MAX_POLES;
	else {
		x->poles = poles;
		if (!(argc % 2))	//	if the list is short...
			x->b[poles] = 0.0;	//	...then there's bogus data in the last b
	}
}

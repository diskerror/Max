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

#include "ext_strings.h"

void *iir_class;

#include "ext.h"
#include "z_dsp.h"
#include <math.h>

typedef struct _iir
{
	t_pxobject l_obj;
	char	poles;		//	number of poles
	char	inOrder;	//	order of coefficients
	double	*a, *b;		//	coefficients from input list
	double	*x, *y;		//	delayed input and output values
} t_iir;

void *iir_new	(t_symbol *o);
void iir_dsp	(t_iir *iir, t_signal **sp, short *count);
t_int *iir_perform(t_int *w);
void iir_clear	(t_iir *iir);
void iir_assist	(t_iir *iir, void *b, long m, long a, char *s);
void iir_free	(t_iir *iir);
void iir_aabab	(t_iir *iir);
void iir_aaabb	(t_iir *iir);
void iir_coeffs	(t_iir *iir, t_symbol *, short argc, t_atom *argv);
void iir_print	(t_iir *iir);

void iir_releasePtrs	(t_iir *iir);

void main(void)
{
	setup((t_messlist **)&iir_class, (method)iir_new, (method)dsp_free, (short)sizeof(t_iir), 0L, A_DEFSYM, 0);
	addmess((method)iir_dsp, "dsp", A_CANT, 0);
	addmess((method)iir_assist, "assist", A_CANT, 0);
	addmess((method)iir_clear, "clear", 0);
	addmess((method)iir_aabab, "aabab", 0);
	addmess((method)iir_aaabb, "aaabb", 0);
	addmess((method)iir_print, "print", 0);
	addmess((method)iir_coeffs, "list", A_GIMME, 0);
	dsp_initclass();
}

void *iir_new(t_symbol *o)
{
	long i;
	
    t_iir *iir = (t_iir *)newobject(iir_class);
    dsp_setup((t_pxobject *)iir, 1);
	
    // one signal outlet
    outlet_new((t_object *)iir, "signal");
    
    //	post message
	post("iir~ [aabab|aaabb]");
	
	//	output order
	if ( strcmp(o->s_name, "aabab") && strcmp(o->s_name, "aaabb") && o->s_name[0] )	//	zero means match
	{
		post("WARNING: Must specify an input order of either \"aabab\" or \"aaabb\".");
		post("   Setting value to \"aabab\".");
		iir->inOrder = 0;
	}
	else
	{
		if (!strcmp(o->s_name, "aabab") || !o->s_name[0] )
			iir->inOrder = 0;
		else
			iir->inOrder = 1;
	}
	
	//	start with all pointers == null
	iir->a = iir->b = iir->x = iir->y = NULL;

    return (iir);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_aabab(t_iir *iir)
{
	iir->inOrder = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_aaabb(t_iir *iir)
{
	iir->inOrder = 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_print(t_iir *iir)
{
	long i;
	if (iir->a)
	{
		post("a[00] = % .15e", iir->a[0]);
		for ( i=1; i <= iir->poles; i++ )
			post("a[%02d] = % .15e   b[%02d] = % .15e", i, iir->a[i], i, iir->b[i]);
	}
	else
		post("a[00] = 0.0");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_free(t_iir *iir)
{
	iir_releasePtrs	(iir);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_assist(t_iir *iir, void *b, long m, long a, char *s)
{
	if (m == 2)
		sprintf(s,"(signal) Output");
	else
	{
		sprintf(s,"(signal) Input, List Input");
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
t_int *iir_perform(t_int *w)
{
	// assign from parameters
    t_float	*in		= (t_float *) w[1];
    t_iir	*iir	= (t_iir *) w[2];
    t_float *out	= (t_float *) w[3];
    long	n		= (long) w[4];
    register long p;
    
    if (!iir->l_obj.z_disabled)
    {
	    // DSP loops
	    if (iir->a && iir->b && iir->x && iir->y)
	    {
			while (n--)
		    {
		    	//	delay old values one sample
		    	for (p=iir->poles; p>0; p--)
		    	{
		    		iir->x[p] = iir->x[p-1];
		    		iir->y[p] = iir->y[p-1];
		    	}
		    	iir->x[0] = (double)*in++;
		    	
		    	iir->y[0] = iir->x[0] * iir->a[0];
		    	for (p=1; p<=iir->poles; p++)
		    	{
		    		iir->y[0] += iir->x[p] * iir->a[p];
		    		iir->y[0] += iir->y[p] * iir->b[p];
		    	}
		    	
		    	*out++	= (float)iir->y[0];
		    }
		}
		else	//	if pointers are no good...
		{
			while (n--)
				*out++ = *in++;	//	...just copy input to output
		}
	}
    return (w+5);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//	Tell host how to insert this plugin in DSP chain.
void iir_dsp(t_iir *iir, t_signal **sp, short *count)
{
	iir_clear(iir);

	dsp_add(iir_perform, 4, sp[0]->s_vec, iir, sp[1]->s_vec, sp[0]->s_n);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//	Clear output. Needed if filter blows up.
void iir_clear(t_iir *iir)
{
	long p;
	for ( p=1; p<=iir->poles; p++ )
		iir->y[p] = 0.0;	//	y[0] is not used
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_coeffs	(t_iir *iir, t_symbol *, short argc, t_atom *argv)
{
	long i, p;
	
	for(i=0; i<argc; i++)
	{
		if (argv[i].a_type != A_FLOAT)
		{
			post("WARNING: All list members must be of type float.");
			return;
		}
	}
	
	iir_releasePtrs(iir);	//	does nothing if pointer are null
	iir->poles = argc/2;	//	integer division
	
	//	now we need pointers for our new data
	iir->a	= (double *)getbytes((iir->poles+1) * sizeof(double));
	iir->b	= (double *)getbytes((iir->poles+1) * sizeof(double));
	iir->x	= (double *)getbytes((iir->poles+1) * sizeof(double));
	iir->y	= (double *)getbytes((iir->poles+1) * sizeof(double));
	
	//	clear the x and y space. Make them start as silence.
	for ( p=0; p<=iir->poles; p++ )
		iir->x[p] = iir->y[p] = 0.0;

	iir->a[0] = argv[0].a_w.w_float;	//	the first is always the same no matter the order
	if (iir->inOrder)	//	aaabb
	{
		for (p=1; p<=iir->poles; p++)
		{
			iir->a[p] = (double)argv[p].a_w.w_float;
			iir->b[p] = (double)argv[iir->poles+p].a_w.w_float;
		}
	}
	else				//	aabab
	{
		for (p=1; p<=iir->poles; p++)
		{
			iir->a[p] = (double)argv[p*2-1].a_w.w_float;
			iir->b[p] = (double)argv[p*2].a_w.w_float;
		}
	}
	
	if (!(argc % 2))	//	if the list is short...
		iir->b[iir->poles] = 0.0;	//	...then there's bogus data in the last b
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void iir_releasePtrs	(t_iir *iir)
{
	if (iir->a && iir->b && iir->x && iir->y)	//	all must != zero
	{
		freebytes(iir->a, (iir->poles+1) * sizeof(double));
		freebytes(iir->b, (iir->poles+1) * sizeof(double));
		freebytes(iir->x, (iir->poles+1) * sizeof(double));
		freebytes(iir->y, (iir->poles+1) * sizeof(double));
		
		iir->a = 0L; iir->b = 0L; iir->x = 0L; iir->y = 0L;
	}
	//else
	//	error("iir_releasePtrs; one pointer was already zero.");
	
	//not an error
	//	assumes pointers were never allocated
}

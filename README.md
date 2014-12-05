MAX Plug-ins
================

cheb~.c
-------
An implementation of a Chebyshev recursive filter for MAX/MSP.

This was my first attempt at writing a plug-in for the MAX programming environment in 2004 for MAX/MSP (version 3? 4?). It was build with Metrowerks Codewarrior and is an implementation of the algorithm presented by Stephen W. Smith in his book “The Scientist and Engineer's Guide to Digital Signal Processing” 2nd edition.

It has a particular “sound” that's interesting and useful. It needs additional code to prevent the “zipper” effect or to smooth the transitions between settings.

iir~.c
------
This will do a IIR or recursive filter based on an input list of coefficients. This was also created for the 2004 version of MAX/MSP.

cheb.c
------
Code extracted from “cheb~.c” and only generates a list of coefficients that can then be sent to an “iir~” object. It's also interesting to watch how the coefficients change with cuttoff frequency and ripple settings.

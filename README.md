MAX Plug-ins
================

cheb
-------
An implementation of a Chebyshev recursive filter for MAX 6.1/7. It only generates a list of coefficients that can then be sent to an “iir~” object. The output can also be sent to a Max multi-slider to watch how the coefficients change with cutoff frequency, poles, and ripple settings.

This is an implementation of the algorithm presented by Stephen W. Smith in his book “The Scientist and Engineer's Guide to Digital Signal Processing” 2nd edition.

iir~
------
This will do a IIR or recursive filter based on an input list of coefficients. This is in the process of being updated for Max 6.1/7. There are some instabilities being worked out.

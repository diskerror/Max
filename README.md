Cycling74 Max Externals (v6.1/v7)
================

cheb
-------
An implementation of a Chebyshev recursive filter. It only generates a list of coefficients that can then be sent to an “iir~” object. The output can also be sent to a Max multi-slider to watch how the coefficients change with cutoff frequency, poles, and ripple settings. It can output a coefficient list in one of two orders, either "aabab…" or "aaa…bb…".

This is an implementation of the algorithm presented by Stephen W. Smith in his book “The Scientist and Engineer's Guide to Digital Signal Processing” 2nd edition.

iir~
------
This will do a IIR or recursive convolution based on an input list of float or double precision coefficients. Coefficients can be in one of two orders, either "aabab…" or "aaa…bb…". Handles both 32- and 64-bit MSP streams.

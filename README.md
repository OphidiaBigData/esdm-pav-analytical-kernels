# ESDM-PAV Analytical Kernels

This module provides a set of analytical kernels supporting in-flight analytics operations on top of the ESDM read streaming interface. 

### Requirements

In order to compile this module the following libraries are required:

- [ESDM](https://github.com/ESiWACE/esdm)

### Installation from sources

If you are building from git, you also need automake, autoconf and libtool. To install the libraries run:

```
$ git clone https://github.com/OphidiaBigData/esdm-pav-analytical-kernels.git
$ cd esdm-pav-analytical-kernels
$ ./bootstrap 
$ ./configure --prefix=PREFIX -with-esdm-path=PATH-TO-ESDM-LIBS
$ make 
$ make install
```

### List of supported functions

- Statitical operations: *maximum, minimum, average, sum, standard deviation, variance*
- Arithmetical operations: *scalar sum, scalar multiplication, absolute value, square root, square, ceil, floor, round, power, exponential, logarithmic, reciprocal value, negation*
- Trigonometrical operations: *sine, cosine, tangent, arcsine, arccosine, arctangent, hyperbolic sine, hyperbolic cosine, hyperbolic tangent*

### Acknowledgement

This software has been developed in the context of the *[ESiWACE2](http://www.esiwace.eu)* project: the *Centre of Excellence in Simulation of Weather and Climate in Europe phase 2*. ESiWACE2 has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No. 823988.

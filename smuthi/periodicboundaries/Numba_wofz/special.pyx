cimport scipy.special.cython_special
cimport numpy as np

cdef api wofz(double in1_real,double in1_imag,double *out_real,double *out_imag):
  cdef double complex z
  z.real=in1_real
  z.imag=in1_imag
  
  cdef double complex out=scipy.special.cython_special.wofz(z)
  
  out_real[0]=out.real
  out_imag[0]=out.imag
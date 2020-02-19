from __future__ import absolute_import

cpdef hyp2f1(double a, double b, double c, z):
  if isinstance(z, float):
    return abk_hyp2f1_real(a,b,c,z)
  if isinstance(z, complex):
    return abk_hyp2f1_complex(a,b,c,z)
  return -1.0

cdef double complex abk_hyp2f1_complex(double a, double b, double c, double complex z):
  return z*z

cdef double abk_hyp2f1_real(double a, double b, double c, double z):
  return z + 1.0
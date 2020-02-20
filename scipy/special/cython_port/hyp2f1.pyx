cimport libc.math

include "constants.pyx"

cpdef hyp2f1(double a, double b, double c, z):
  if isinstance(z, float):
    return abk_hyp2f1_real(a,b,c,z)
  if isinstance(z, complex):
    return abk_hyp2f1_complex(a,b,c,z)
  return -1.0

cdef complex abk_hyp2f1_complex(double a, double b, double c, double complex z):
  return z*z

cdef double abk_hyp2f1_real(double a, double b, double c, double z):
  #  Handle some of the easy cases when one of arguments is zero.
  if z == 0.0: return 1.0
  if (a == 0.0 or b == 0.0) and (c != 0.0): return 1.0

  # Check to see if one of a or b is a negative integer.
  cdef double nearest_int_a, nearest_int_b
  nearest_int_a = round(a)
  nearest_int_b = round(b)

  cdef bint is_negative_int_a, is_negative_int_b
  is_negative_int_a = False
  is_negative_int_b = False
  
  if a <= 0.0 and libc.math.fabs(a - nearest_int_a) < EPS:
    is_negative_int_a = True

  if b <= 0.0 and libc.math.fabs(b - nearest_int_b) < EPS:
    is_negative_int_b = True

  return nearest_int_a
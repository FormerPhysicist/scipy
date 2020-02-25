from libc.math cimport fabs, pow
import numpy as np

include "constants.pyx"

cpdef hyp2f1(double a, double b, double c, z):
  if isinstance(z, float):
    return _hyp2f1_real(a,b,c,z)
  if isinstance(z, complex):
    return _hyp2f1_complex(a,b,c,z)
  return -1.0

cdef complex _hyp2f1_complex(double a, double b, double c, double complex z):
  return z*z

cdef double _hyp2f1_real(double a, double b, double c, double z):
  #  Handle some of the easy cases when one of arguments is zero.
  if z == 0.0: return 1.0
  if (a == 0.0 or b == 0.0) and (c != 0.0): return 1.0

  # Compute some intermediate values
  cdef double d, s, int_d
  s = 1.0 - z
  d = c - a - b
  int_d = round(d)

  # Check to see if one of a or b is a negative integer.
  cdef double int_a, int_b
  int_a = round(a)
  int_b = round(b)

  cdef bint is_neg_int_a, is_neg_int_b
  is_neg_int_a = a <= 0.0 and fabs(a - int_a) < EPS
  is_neg_int_b = b <= 0.0 and fabs(b - int_b) < EPS

  # Check to see if d is a negative integer.
  cdef bint is_neg_int_d
  is_neg_int_d = d <= -1 and fabs(d - int_d) < EPS

  # Abramowitz and Stegun 15.3.3
  if is_neg_int_d and s > 0.0 and not (is_neg_int_a or is_neg_int_b):
    return pow(s, d) * _hyp2f1_real(c - a, c - b, c, z)

  # series diverges at z = 1.0 when d < 0.0
  if d < 0.0 and z == 1.0 and not (is_neg_int_a or is_neg_int_b):
    return np.inf

  # Handle the special cases when a = c or b = c
  cdef bint ac_equal, bc_equal
  if fabs(z) < 1.0 or z == -1.0:
    ac_equal = fabs(a - c) < EPS
    bc_equal = fabs(b - c) < EPS

    # Abramowitz and Stegun 15.4.2
    if ac_equal and is_neg_int_a:
      return _hyp2f1_neg_c_equal(b, a, z)

    if bc_equal and is_neg_int_b:
      return _hyp2f1_neg_c_equal(a, b, z)

    # Abramowitz and Stegun 15.1.8
    if ac_equal:
      return pow(s, -b)
    
    if bc_equal:
      return pow(s, -a)

  # Check to see if c is a negative integer. This means hyp2f1 diverges
  # unless the series terminates before c can blow everything up.
  cdef double int_c
  cdef bint is_neg_int_c, terminates_early
  
  int_c = round(c)
  is_neg_int_c = c <= 0.0 and fabs(c - int_c) < EPS
  terminates_early = (is_neg_int_a and int_a > int_c) or (is_neg_int_b and int_b > int_c)
  if is_neg_int_c and not terminates_early:
    return np.inf

  # If a or b is a negative integer then hyp2f1 is a polynomial.
  if is_neg_int_a or (is_neg_int_a and is_neg_int_b and a > b):
    return _hyp2f1_poly(-a, b, c, z)
  elif is_neg_int_b:
    return _hyp2f1_poly(-b, a, c, z)

  return 1.0

cdef double _hyp2f1_neg_c_equal(double a, double b, double z):
  pass

cdef double _hyp2f1_poly(double m, double b, double c, double z):
  pass
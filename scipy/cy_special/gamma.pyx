from libc.math cimport fabs, pow, exp, floor, sin
import numpy as np
from .polyeval cimport polyeval

include "constants.pyx"

cdef double P[7]
cdef double Q[8]
cdef double STIR[5]

P[:] = [
  	1.60119522476751861407E-4,
		1.19135147006586384913E-3,
		1.04213797561761569935E-2,
		4.76367800457137231464E-2,
		2.07448227648435975150E-1,
		4.94214826801497100753E-1,
		9.99999999999999996796E-1
]

Q[:] = [
    -2.31581873324120129819E-5,
		5.39605580493303397842E-4,
		-4.45641913851797240494E-3,
		1.18139785222060435552E-2,
		3.58236398605498653373E-2,
		-2.34591795718243348568E-1,
		7.14304917030273074085E-2,
		1.00000000000000000320E0
]

STIR[:] = [
    7.87311395793093628397E-4,
		-2.29549961613378126380E-4,
		-2.68132617805781232825E-3,
		3.47222221605458667310E-3,
		8.33333333333482257126E-2,
]

cdef double stirf(double x):
  cdef double y, w, v
  
  if x >= MAXGAM:
    return np.inf

  w = 1.0/x
  w = 1.0 + w * polyeval(w, STIR, 4)
  y = exp(x)

  # Avoid overflow in pow()
  if x > MAXSTIR:
    v = pow(x, 0.5 * x - 0.25)
    y = v * (v / y)
  else:
    y = pow(x, x - 0.5) / y

  return SQTPI * y * w

cdef double gamma(double x):
  cdef double p, q, z, i
  cdef int sgngam
  i = floor(x)
  sgngam = 1

  if np.isinf(x):
    return x
  
  # Check to see if x is a zero or a negative integer.
  if x == 0.0 or (x == i and i < 0):
    return np.inf

  # Use stirlings formula for large x.
  q = fabs(x)
  if x > 33.0:
    return stirf(x)

  if x < -33.0:
    p = round(q) 

    # If the integer part is ODD then the sign is negative.
    if floor(q) % 2 == 0:
      sgngam = -1
    z = q - p
    z = q * sin(PI * z)
    if z == 0.0:
      return sgngam * np.inf
    z = fabs(z)
    return sgngam * PI / (z * stirf(q))
  
  # Use gamma recursion relation to decrease absolute magnitude of argument.
  z = 1.0
  while x >= 3.0:
    x -= 1.0
    z *= x

  while x < -1.e-9:
    z /= x
    x += 1.0

  # If x has a very small magnitude handle that first.
  if fabs(x) < 1.e-9:
    return z / ((1.0 + 0.5772156649015329 * x) * x)

  # If here x has to be between 1e-9 and 3. Make sure x is between 2 and 3
  while x < 2.0:
    z /= x
    x += 1.0
  
  # Now evaluate gamma for the updated x value.
  if x == 2.0:
    return z

  x -= 2.0
  p = polyeval(x, P, 6)
  q = polyeval(x, Q, 7)
  return z * p / q
  
cimport libc.math

cdef double round(double x):
  cdef double i, f

  # The largest integer still smaller than x.
  i = libc.math.floor(x)

  # The fractional part.
  f = x - i

  # Round up.
  if f > 0.5:
    return i + 1.0

  if f == 0.5:
    f = i - 2.0 * libc.math.floor(0.5 * i)
    if f == 1.0:
      return i + 1.0

  # Round down.
  return i
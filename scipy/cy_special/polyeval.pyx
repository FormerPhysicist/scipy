cdef double polyeval(double x, double* coef, int N):
  cdef double value
  value = coef[0]
  for i in range(1, N):
    value = value * x + coef[i]
  return value
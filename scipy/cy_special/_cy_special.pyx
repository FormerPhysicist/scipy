from __future__ import absolute_import

from .round cimport roundy

cpdef double abk(double y):
  return y + 1.0

cpdef double round_wrapper(double y):
  return roundy(y)

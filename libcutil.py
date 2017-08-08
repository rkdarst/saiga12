from ctypes import *

_libraries = {}
_libraries['/home/richard/research/saiga12/_cutil.so'] = CDLL('/home/richard/research/saiga12/_cutil.so')


calcFs_lmp = _libraries['/home/richard/research/saiga12/_cutil.so'].calcFs_lmp
calcFs_lmp.restype = c_double
calcFs_lmp.argtypes = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_double), c_int, POINTER(c_double), c_int, c_int, POINTER(c_double), POINTER(c_double), c_int]
__all__ = ['calcFs_lmp']

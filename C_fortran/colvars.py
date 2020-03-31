import os, sys
import ctypes

cv = ctypes.cdll.LoadLibrary('./libcolvars_C.so')

print(cv.allocate_Colvars('Allocated from Python through ctypes'))

print(cv.get_colvarscript_result())

print(cv.run_colvarscript_command('version'))

print(cv.get_colvarscript_result())

import os, sys
import ctypes

# Convert python string to byte object for passing to C functions expecting char *
def b(s):
  return bytes(s, 'utf-8')

def run_cmd(cmd):
  words = [b(w) for w in cmd.split()]
  objc = len(words) + 1
  objv = (ctypes.c_char_p * objc) (b('cv'), *words)
  return cv.run_colvarscript_command(objc, objv)

def run_cmds(*cmds):
  words = [b(w) for w in cmds]
  objc = len(words) + 1
  objv = (ctypes.c_char_p * objc) (b('cv'), *words)
  return cv.run_colvarscript_command(objc, objv)

def cv_result():
  return cv.get_colvarscript_result().decode('utf-8')

def run_echo(cmd):
  print(f'Calling cv {cmd}')
  run_cmd(cmd)
  print(cv_result())


MAX_ARRAY_SIZE = 10
cmd_type = ctypes.c_char_p * MAX_ARRAY_SIZE


cv = ctypes.cdll.LoadLibrary('./libcolvars_C.so')

cv.allocate_Colvars(b('Allocated from Python through ctypes'))

cv.get_colvarscript_result.restype = ctypes.c_char_p

if run_cmd('version'):
  print('Error running cv version')

print(f'Detected cv version {cv_result()}')

run_cmds('config', 'units electron')
run_echo('help')
run_echo('getconfig')
run_echo('delete')
run_echo('reset')

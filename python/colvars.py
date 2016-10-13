import os, sys
import ctypes

class colvarscript:
    # load all symbols from the NAMD executable
    _interp = ctypes.cdll.LoadLibrary('')
    def run(self, cmd):
        args = ['cv'] + cmd.split()
        self._interp.run_colvarscript_command.argtypes = \
            [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        self._interp.get_colvarscript_result.restype = ctypes.c_char_p
        if (self._interp.run_colvarscript_command(len(args), \
             (ctypes.c_char_p * len(args))(*args)) == 0):
            return self._interp.get_colvarscript_result()
        else:
            return 'Colvarscript error.'

if __name__ == '__main__':
    # Example use
    cv = colvarscript()
    cv.run('update')
    for b in cv.run('list biases').split():
        energy = float(cv.run("bias %s energy" % b))
        sys.stdout.write("Energy[%s] = %f\n" % (b, energy))


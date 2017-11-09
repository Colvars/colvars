from __future__ import print_function
import ctypes


class cvscript(object):
    '''
    Interface class between Python and the Colvars module.
    Example use:
        cv = cvscript()
        cv.run('update')
        for b in cv.run('list biases').split():
            energy = float(cv.run("bias %s energy" % b))
            print("Energy[%s] = %f" % (b, energy))
    '''

    _interp = ctypes.cdll.LoadLibrary('')

    _interp.cvscript_n_commands.restype = ctypes.c_int
    _interp.cvscript_n_commands.argtypes = []
    _interp.cvscript_command_names.restype = ctypes.POINTER(ctypes.c_char_p)
    _interp.cvscript_command_names.argtypes = [ctypes.c_int]
    _interp.cvscript_help.restype = ctypes.c_char_p
    _interp.cvscript_help.argtypes = [ctypes.c_int]

    n_commands = _interp.cvscript_n_commands()
    commands_arr = ctypes.cast(_interp.cvscript_command_names(),
                               ctypes.POINTER(c_char_p))

    for c in range(n_commands):

        cname = commands_arr[c]
        chelp = _interp.cvscript_help(cname).splitlines()
        cargs = ""
        arghelp = ""

        for line in chelp[1:]:
            if (len(cargs)): 
                cargs += ", "
            words = line.split(' ')
            cargs += words[0]
            if ('optional' in words):
                cargs += "=None"
            arghelpelems = line.split(' - ')
            arghelp += """
    %s
        %s
""" % (arghelpelems[0], arghelpelems[1])

        exec("""
def %s(%s):
    '''
    %s

    Parameters
    ----------
    %s
    '''
""" % (cname, cargs, chelp, arghelp))


if __name__ == '__main__':
    # Example use: print energy of each bias
    cv = cvscript()
    cv.run('update')
    for b in cv.run('list biases').split():
        energy = float(cv.run("bias %s energy" % b))
        print("Energy[%s] = %f\n" % (b, energy))


# -*- tcl -*-

python {
import numpy as np

def calc_cvsq(comps, xi):
    # comp[0] is the first (and only) component
    x = comps[0]
    return float(x*x)
def calc_cvsq_gradient(comps, xi, grads):
    x = comps[0]
    grads[0] = 2.0*x
    return grads
}



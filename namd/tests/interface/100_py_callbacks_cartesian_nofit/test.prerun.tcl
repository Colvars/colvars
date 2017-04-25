# -*- tcl -*-

python {
from __future__ import print_function
import numpy as np

def calc_gyration(comps, xi):
    xyz = comps[0] # comp[0] is the first (and only) component
    n = xyz.size/3
    return np.sqrt(np.sum(np.square(xyz))/n)
def calc_gyration_gradient(comps, xi, grads):
    xyz = comps[0]
    n = xyz.size/3
    np.multiply(xyz, 1.0/(xi*n), out=grads[0])
    return


def calc_distancevec(comps, xi):
    xyz1 = comps[0]
    n1 = (xyz1.size/3)
    xyz1 = xyz1.reshape((n1, 3))

    xyz2 = comps[1]
    n2 = (xyz2.size/3)
    xyz2 = xyz2.reshape((n2, 3))

    return xyz2.mean(axis=0) - xyz1.mean(axis=0)

def calc_distance(comps, xi):
    distv = calc_distancevec(comps, xi)
    print("v = ", distv)    
    dist2 = np.sum(np.square(distv))
    return np.sqrt(dist2)

def calc_distance_gradient(comps, xi, grads):
    distv = calc_distancevec(comps, xi)
    print("xi = ", np.array_str(xi, precision=16))
    print("v = ", np.array_str(distv, precision=16))
    print("u = ", np.array_str(distv/xi, precision=16))
    xyz1 = comps[0]
    n1 = (xyz1.size/3)
    xyz1.reshape((n1*3))

    xyz2 = comps[1]
    n2 = (xyz2.size/3)
    xyz2.reshape((n2*3))

    np.multiply(-1.0/float(n1)/xi, 
                np.tensordot(np.ones((n1)), distv, axes=0).reshape((n1*3)), 
                out=grads[0])
    np.multiply(1.0/float(n2)/xi,
                np.tensordot(np.ones((n2)), distv, axes=0).reshape((n2*3)), 
                out=grads[1])

    print(np.array_str(grads[0], precision=16))
    print(np.array_str(grads[1], precision=16))
    return
    
}



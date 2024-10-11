#!/usr/bin/env python3
import torch
import re
import numpy as np

torch.set_printoptions(precision=10)


# function to extract grad
def set_grad(var):
    def hook(grad):
        var.grad = grad
    return hook

mass1 = torch.tensor([14.0070, 12.0110, 12.0110, 15.999], dtype=torch.float64)

mass2 = torch.tensor([14.0070, 12.0110, 12.0110, 15.999], dtype=torch.float64)

ref_fitting_group = torch.tensor(
    [[ 5.69469   , -1.80418   , -5.22382   ],
       [ 6.37169   , -1.08118   , -4.10382   ],
       [ 6.61169   ,  0.365824  , -4.58082   ],
       [ 5.62269   , -1.09518   , -2.77982   ],
       [ 4.39969   , -0.986176  , -2.72882   ],
       [ 6.36569   , -1.21018   , -1.66282   ],
       [ 5.83469   , -1.37018   , -0.318824  ],
       [ 7.00969   , -1.45018   ,  0.678176  ],
       [ 4.84269   , -0.296176  ,  0.126176  ],
       [ 3.77369   , -0.596176  ,  0.648176  ],
       [ 5.15669   ,  0.991824  , -0.113824  ],
       [ 4.27869   ,  2.10682   ,  0.190176  ],
       [ 5.03669   ,  3.42282   , -0.0748235 ],
       [ 2.96669   ,  2.09282   , -0.596824  ],
       [ 1.89969   ,  2.39182   , -0.0678235 ],
       [ 3.01469   ,  1.69682   , -1.88382   ],
       [ 1.84869   ,  1.57082   , -2.73482   ],
       [ 2.30269   ,  1.40782   , -4.19882   ],
       [ 0.962686  ,  0.399824  , -2.31982   ],
       [-0.255314  ,  0.522824  , -2.22882   ],
       [ 1.57569   , -0.758176  , -2.00882   ],
       [ 0.900686  , -1.92218   , -1.46882   ],
       [ 1.91969   , -3.07218   , -1.33582   ],
       [ 0.223686  , -1.65518   , -0.123824  ],
       [-0.930314  , -2.01818   ,  0.0891765 ],
       [ 0.918686  , -0.964176  ,  0.800176  ],
       [ 0.372686  , -0.514176  ,  2.06618   ],
       [ 1.49069   ,  0.164824  ,  2.88418   ],
       [-0.814314  ,  0.435824  ,  1.92118   ],
       [-1.84631   ,  0.257824  ,  2.56118   ],
       [-0.707314  ,  1.45482   ,  1.04518   ],
       [-1.78431   ,  2.37982   ,  0.745176  ],
       [-1.24331   ,  3.48882   , -0.178824  ],
       [-3.00931   ,  1.72082   ,  0.109176  ],
       [-4.14731   ,  2.02882   ,  0.455176  ],
       [-2.80331   ,  0.769824  , -0.821824  ],
       [-3.86731   , -0.00917647, -1.42282   ],
       [-3.28731   , -0.820176  , -2.59882   ],
       [-4.57231   , -0.946176  , -0.440824  ],
       [-5.78931   , -1.10018   , -0.475824  ],
       [-3.81431   , -1.58718   ,  0.469176  ],
       [-4.35131   , -2.44718   ,  1.50418   ],
       [-3.19331   , -3.25818   ,  2.12018   ],
       [-5.10031   , -1.70118   ,  2.61018   ],
       [-6.15131   , -2.13818   ,  3.07218   ],
       [-6.33131   ,  1.19882   ,  3.62918   ],
       [-7.08531   ,  1.74482   ,  4.44018   ],
       [-6.45431   ,  1.44382   ,  2.30618   ],
       [-4.56431   , -0.557176  ,  3.07518   ],
       [-5.19131   ,  0.240824  ,  4.11318   ],
       [-4.10131   ,  1.05882   ,  4.83518   ]], dtype=torch.float64)

# fitting_group = atom_positions
# main_group = atom_positions[[0, 1, 2, 3]]


def get_optimal_rotation(fitting_atoms: torch.tensor, reference_atoms: torch.tensor):
    fitting_atoms_centered = fitting_atoms - fitting_atoms.mean(dim=0)
    reference_atoms_centered = reference_atoms - reference_atoms.mean(dim=0)
    mat_R = torch.linalg.matmul(fitting_atoms_centered.T, reference_atoms_centered)
    F00 = mat_R[0][0] + mat_R[1][1] + mat_R[2][2]
    F01 = mat_R[1][2] - mat_R[2][1]
    F02 = mat_R[2][0] - mat_R[0][2]
    F03 = mat_R[0][1] - mat_R[1][0]
    F10 = F01
    F11 = mat_R[0][0] - mat_R[1][1] - mat_R[2][2]
    F12 = mat_R[0][1] + mat_R[1][0]
    F13 = mat_R[0][2] + mat_R[2][0]
    F20 = F02
    F21 = F12
    F22 = -mat_R[0][0] + mat_R[1][1] - mat_R[2][2]
    F23 = mat_R[1][2] + mat_R[2][1]
    F30 = F03
    F31 = F13
    F32 = F23
    F33 = -mat_R[0][0] - mat_R[1][1] + mat_R[2][2]
    row0 = torch.stack((F00, F01, F02, F03))
    row1 = torch.stack((F10, F11, F12, F13))
    row2 = torch.stack((F20, F21, F22, F23))
    row3 = torch.stack((F30, F31, F32, F33))
    F = torch.stack((row0, row1, row2, row3))
    w, v = torch.linalg.eigh(F)
    q = v[:, -1]
    q.register_hook(set_grad(q))
    R00 = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
    R01 = 2.0 * (q[1] * q[2] - q[0] * q[3])
    R02 = 2.0 * (q[1] * q[3] + q[0] * q[2])
    R10 = 2.0 * (q[1] * q[2] + q[0] * q[3])
    R11 = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3]
    R12 = 2.0 * (q[2] * q[3] - q[0] * q[1])
    R20 = 2.0 * (q[1] * q[3] - q[0] * q[2])
    R21 = 2.0 * (q[2] * q[3] + q[0] * q[1])
    R22 = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]
    row0 = torch.stack((R00, R01, R02))
    row1 = torch.stack((R10, R11, R12))
    row2 = torch.stack((R20, R21, R22))
    R = torch.stack((row0, row1, row2))
    return w, v, R, q


def center_of_mass(pos, mass):
    return torch.sum(pos * mass[:, None], dim=0) / torch.sum(mass)


def read_data_from_traj(filename):
    group1_positions = []
    group2_positions = []
    fitting_group1_positions = []
    fitting_group2_positions = []
    forces = []
    with open(filename, 'r') as f_input:
        for line in f_input:
            s = line.strip()
            if s.startswith('#'):
                # skip comment lines
                continue
            if len(s) > 0:
                # skip empty lines
                # match the data in parentheses
                result = re.findall(r'\(([^)]+)\)', s)
                # get the force applied on the orientation
                forces.append(torch.tensor(list(map(float, result[1].split(','))), dtype=torch.float64))
                # get the group1 position
                tmp = np.array(list(map(float, result[2].split(',')))).reshape(-1, 3)
                group1_positions.append(torch.tensor(tmp, dtype=torch.float64, requires_grad=True))
                # get the group2 position
                tmp = np.array(list(map(float, result[3].split(',')))).reshape(-1, 3)
                group2_positions.append(torch.tensor(tmp, dtype=torch.float64, requires_grad=True))
                # get the fitting group1 position
                tmp = np.array(list(map(float, result[4].split(',')))).reshape(-1, 3)
                fitting_group1_positions.append(torch.tensor(tmp, dtype=torch.float64, requires_grad=True))
                fitting_group2_positions.append(torch.tensor(tmp, dtype=torch.float64, requires_grad=True))
    return forces, group1_positions, group2_positions, fitting_group1_positions, fitting_group2_positions


# Use Colvars definition: the distance between two unit vectors is angle between them
def dist2_distance_dir(v1, v2):
    v1_normalized = v1 / torch.linalg.norm(v1)
    v2_normalized = v2 / torch.linalg.norm(v2)
    dist = torch.acos(torch.sum(v1_normalized * v2_normalized))
    return dist * dist


def scan_force_from_log(filename):
    main_group_forces = []
    fitting_group_forces = []
    with open(filename, 'r') as f_input:
        new_main_group_started = False
        main_group_force = []
        new_fitting_group_started = False
        fitting_group_force = []
        for line in f_input:
            if new_main_group_started:
                result = re.findall(r'\(([^)]+)\)', line)
                if result:
                    main_group_force.append(list(map(float, result[0].split(','))))
                else:
                    new_main_group_started = False
                    main_group_forces.append(np.array(main_group_force))
                    main_group_force = []
            if new_fitting_group_started:
                result = re.findall(r'\(([^)]+)\)', line)
                if result:
                    fitting_group_force.append(list(map(float, result[0].split(','))))
                else:
                    new_fitting_group_started = False
                    fitting_group_forces.append(np.array(fitting_group_force))
                    fitting_group_force = []
            if line.startswith('colvars:   Applying force on the fitting group of main group:'):
                new_fitting_group_started = True
                continue
            if line.startswith('colvars:   Applying force on main group'):
                new_main_group_started = True
                continue
    return main_group_forces, fitting_group_forces


def test(group1, group2, fitting_group_1, fitting_group_2):
    # Fit group 1
    w1, v1, R1, q1 = get_optimal_rotation(fitting_group_1, ref_fitting_group)
    rpg_cog1 = fitting_group_1.mean(dim=0)
    # group1_c = group1.mean(dim=0)
    group1_centered = group1 - rpg_cog1
    group1_rotated = (R1 @ group1_centered.T).T
    # Fit group 2
    w2, v2, R2, q2 = get_optimal_rotation(fitting_group_2, ref_fitting_group)
    rpg_cog2 = fitting_group_2.mean(dim=0)
    # group2_c = group2.mean(dim=0)
    group2_centered = group2 - rpg_cog2
    group2_rotated = (R2 @ group2_centered.T).T
    # distanceDir
    c2 = center_of_mass(group2_rotated, mass2)
    c1 = center_of_mass(group1_rotated, mass1)
    d = c2 - c1
    d_dir = d / torch.linalg.norm(d)
    ref = torch.tensor([1.0, 0.1, 0.2], dtype=torch.float64)
    energy = 0.5 * 10.0 * dist2_distance_dir(d_dir, ref) / 0.25
    energy.backward(retain_graph=True)
    return -group1.grad.cpu().detach().numpy(), -group2.grad.cpu().detach().numpy(), \
        -fitting_group_1.grad.cpu().detach().numpy(), -fitting_group_2.grad.cpu().detach().numpy()


if __name__ == '__main__':
    print('Compare AutoDiff/test.colvars.traj and AutoDiff/test.colvars.out')
    # read the forces from the debug Colvars log file
    main_group_forces_colvars, fitting_group_forces_colvars = scan_force_from_log('AutoDiff/test.colvars.out')
    group1_forces_colvars = main_group_forces_colvars[::2]
    group2_forces_colvars = main_group_forces_colvars[1::2]
    fitting_group1_forces_colvars = fitting_group_forces_colvars[::2]
    fitting_group2_forces_colvars = fitting_group_forces_colvars[1::2]
    # read the positions from .colvars.traj
    _, group1_positions, group2_positions, fitting_group1_positions, fitting_group2_positions = read_data_from_traj('AutoDiff/test.colvars.traj')
    for i in range(len(group1_positions)):
        # calculate forces from pytorch
        gf1_pytorch, gf2_pytorch, fgf1_pytorch, fgf2_pytorch = test(
            group1_positions[i], group2_positions[i], fitting_group1_positions[i], fitting_group2_positions[i])
        # error with respect to colvars
        error_group1 = np.mean(np.abs(group1_forces_colvars[i] - gf1_pytorch) / np.abs(gf1_pytorch))
        error_group2 = np.mean(np.abs(group2_forces_colvars[i] - gf2_pytorch) / np.abs(gf2_pytorch))
        error_fitting_group1 = np.mean(np.abs(fitting_group1_forces_colvars[i] - fgf1_pytorch) / np.abs(fgf1_pytorch))
        error_fitting_group2 = np.mean(np.abs(fitting_group2_forces_colvars[i] - fgf2_pytorch) / np.abs(fgf2_pytorch))
        print(f'Error of forces: group1 = {error_group1:12.5e}, group2 = {error_group2:12.5e}, fitting_group1 = {error_fitting_group1:12.5e}, fitting_group2 = {error_fitting_group2:12.5e}')

    # Do the same thing with the restart file
    print('Compare AutoDiff/test.restart.colvars.traj and AutoDiff/test.restart.colvars.out')
    # read the forces from the debug Colvars log file
    main_group_forces_colvars, fitting_group_forces_colvars = scan_force_from_log('AutoDiff/test.restart.colvars.out')
    group1_forces_colvars = main_group_forces_colvars[::2]
    group2_forces_colvars = main_group_forces_colvars[1::2]
    fitting_group1_forces_colvars = fitting_group_forces_colvars[::2]
    fitting_group2_forces_colvars = fitting_group_forces_colvars[1::2]
    # read the positions from .colvars.traj
    _, group1_positions, group2_positions, fitting_group1_positions, fitting_group2_positions = read_data_from_traj('AutoDiff/test.restart.colvars.traj')
    for i in range(len(group1_positions)):
        # calculate forces from pytorch
        gf1_pytorch, gf2_pytorch, fgf1_pytorch, fgf2_pytorch = test(
            group1_positions[i], group2_positions[i], fitting_group1_positions[i], fitting_group2_positions[i])
        # error with respect to colvars
        error_group1 = np.mean(np.abs(group1_forces_colvars[i] - gf1_pytorch) / np.abs(gf1_pytorch))
        error_group2 = np.mean(np.abs(group2_forces_colvars[i] - gf2_pytorch) / np.abs(gf2_pytorch))
        error_fitting_group1 = np.mean(np.abs(fitting_group1_forces_colvars[i] - fgf1_pytorch) / np.abs(fgf1_pytorch))
        error_fitting_group2 = np.mean(np.abs(fitting_group2_forces_colvars[i] - fgf2_pytorch) / np.abs(fgf2_pytorch))
        print(f'Error of forces: group1 = {error_group1:12.5e}, group2 = {error_group2:12.5e}, fitting_group1 = {error_fitting_group1:12.5e}, fitting_group2 = {error_fitting_group2:12.5e}')

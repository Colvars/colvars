#!/usr/bin/env python3
import torch
import re
import numpy as np

torch.set_printoptions(precision=10)

ref_main_group = torch.tensor(
    [[6.3061, -1.592, -3.585],
     [5.4031, -0.784, 0.045] ,
     [3.8061, 2.804, -0.129] ,
     [1.2441, 1.799, -2.726] ,
     [0.5861, -1.803, -1.516],
     [-0.2169, -0.926, 2.133],
     [-2.1679, 2.13, 1.258]  ,
     [-4.4689, 0.434, -1.263],
     [-4.8649, -2.25, 1.46]  ,
     [-5.6269, 0.188, 4.323]], dtype=torch.float64)

ref_fitting_group = torch.tensor(
    [[5.69469, -1.80418, -5.22382]    ,
     [6.37169, -1.08118, -4.10382]    ,
     [6.61169, 0.365824, -4.58082]    ,
     [5.62269, -1.09518, -2.77982]    ,
     [4.39969, -0.986176, -2.72882]   ,
     [6.36569, -1.21018, -1.66282]    ,
     [5.83469, -1.37018, -0.318824]   ,
     [7.00969, -1.45018, 0.678176]    ,
     [4.84269, -0.296176, 0.126176]   ,
     [3.77369, -0.596176, 0.648176]   ,
     [5.15669, 0.991824, -0.113824]   ,
     [4.27869, 2.10682, 0.190176]     ,
     [5.03669, 3.42282, -0.0748235]   ,
     [2.96669, 2.09282, -0.596824]    ,
     [1.89969, 2.39182, -0.0678235]   ,
     [3.01469, 1.69682, -1.88382]     ,
     [1.84869, 1.57082, -2.73482]     ,
     [2.30269, 1.40782, -4.19882]     ,
     [0.962686, 0.399824, -2.31982]   ,
     [-0.255314, 0.522824, -2.22882]  ,
     [1.57569, -0.758176, -2.00882]   ,
     [0.900686, -1.92218, -1.46882]   ,
     [1.91969, -3.07218, -1.33582]    ,
     [0.223686, -1.65518, -0.123824]  ,
     [-0.930314, -2.01818, 0.0891765] ,
     [0.918686, -0.964176, 0.800176]  ,
     [0.372686, -0.514176, 2.06618]   ,
     [1.49069, 0.164824, 2.88418]     ,
     [-0.814314, 0.435824, 1.92118]   ,
     [-1.84631, 0.257824, 2.56118]    ,
     [-0.707314, 1.45482, 1.04518]    ,
     [-1.78431, 2.37982, 0.745176]    ,
     [-1.24331, 3.48882, -0.178824]   ,
     [-3.00931, 1.72082, 0.109176]    ,
     [-4.14731, 2.02882, 0.455176]    ,
     [-2.80331, 0.769824, -0.821824]  ,
     [-3.86731, -0.00917647, -1.42282],
     [-3.28731, -0.820176, -2.59882]  ,
     [-4.57231, -0.946176, -0.440824] ,
     [-5.78931, -1.10018, -0.475824]  ,
     [-3.81431, -1.58718, 0.469176]   ,
     [-4.35131, -2.44718, 1.50418]    ,
     [-3.19331, -3.25818, 2.12018]    ,
     [-5.10031, -1.70118, 2.61018]    ,
     [-6.15131, -2.13818, 3.07218]    ,
     [-6.33131, 1.19882, 3.62918]     ,
     [-7.08531, 1.74482, 4.44018]     ,
     [-6.45431, 1.44382, 2.30618]     ,
     [-4.56431, -0.557176, 3.07518]   ,
     [-5.19131, 0.240824, 4.11318]    ,
     [-4.10131, 1.05882, 4.83518]], dtype=torch.float64)

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
    return w, v, R


def test(force, main_group, fitting_group):
    # Using all atoms as the fitting group
    w, v, R = get_optimal_rotation(fitting_group, ref_fitting_group)
    main_group_centered = main_group - main_group.mean(dim=0)
    ref_main_group_centered = ref_main_group - ref_main_group.mean(dim=0)
    # Rotate the main group
    main_group_r = (R @ main_group_centered.T).T
    # main_group_r.register_hook(set_grad(main_group_r))
    w2, v2, R2 = get_optimal_rotation(ref_main_group_centered, main_group_r)
    q = v2[:, -1]
    ref_quat = torch.tensor([1.0, 0.0, 0.0, 0.0], dtype=torch.float64)
    q_inner_ref_quat = torch.sum(ref_quat * q)
    if (q_inner_ref_quat > 0):
        final_q = q
    else:
        final_q = -1.0 * q
    # print('Orientation:\n', final_q)
    # Force on q
    # force = torch.tensor([3.40697e-06, -8.24548e-05, -5.71708e-05, 5.95895e-05], dtype=torch.float64)
    # Apply force
    # torch.sum(force * final_q).backward()
    (force[0] * final_q[0]).backward(retain_graph=True)
    (force[1] * final_q[1]).backward(retain_graph=True)
    (force[2] * final_q[2]).backward(retain_graph=True)
    (force[3] * final_q[3]).backward(retain_graph=True)
    # print('main group force:\n', main_group.grad)
    # print('main group (rotated) force:\n', main_group_r.grad)
    # print('fitting group force:\n', fitting_group.grad)
    return main_group.grad.cpu().detach().numpy(), fitting_group.grad.cpu().detach().numpy()


def read_data_from_traj(filename):
    main_group_positions = []
    fitting_group_positions = []
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
                # get the main group position
                tmp = np.array(list(map(float, result[2].split(',')))).reshape(-1, 3)
                main_group_positions.append(torch.tensor(tmp, dtype=torch.float64, requires_grad=True))
                # get the fitting group position
                tmp = np.array(list(map(float, result[3].split(',')))).reshape(-1, 3)
                fitting_group_positions.append(torch.tensor(tmp, dtype=torch.float64, requires_grad=True))
    return forces, main_group_positions, fitting_group_positions


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


if __name__ == '__main__':
    print('Compare AutoDiff/test.colvars.traj and AutoDiff/test.colvars.out')
    # read positions from Colvars traj
    forces, main_group_positions, fitting_group_positions = read_data_from_traj('AutoDiff/test.colvars.traj')
    main_group_forces_pytorch = []
    fitting_group_forces_pytorch = []
    # calculate the forces by pytorch
    for f, pm, pf in zip(forces, main_group_positions, fitting_group_positions):
        f1, f2 = test(f, pm, pf)
        main_group_forces_pytorch.append(f1)
        fitting_group_forces_pytorch.append(f2)
    # read the forces from the debug Colvars log file
    main_group_forces_colvars, fitting_group_forces_colvars = scan_force_from_log('AutoDiff/test.colvars.out')
    # compare the force (calculate the relative error)
    for i in range(len(main_group_forces_colvars)):
        error_main = np.mean(np.abs(main_group_forces_colvars[i] - main_group_forces_pytorch[i]) / np.abs(main_group_forces_colvars[i]))
        error_fitting = np.mean(np.abs(fitting_group_forces_colvars[i] - fitting_group_forces_pytorch[i]) / np.abs(fitting_group_forces_colvars[i]))
        print(f'Frame {i}: error of the main group forces {error_main:12.5e} , fitting group forces {error_fitting:12.5e}')

    # Do the same thing with the restart file
    print('Compare AutoDiff/test.restart.colvars.traj and AutoDiff/test.restart.colvars.out')
    forces, main_group_positions, fitting_group_positions = read_data_from_traj('AutoDiff/test.restart.colvars.traj')
    main_group_forces_pytorch = []
    fitting_group_forces_pytorch = []
    # calculate the forces by pytorch
    for f, pm, pf in zip(forces, main_group_positions, fitting_group_positions):
        f1, f2 = test(f, pm, pf)
        main_group_forces_pytorch.append(f1)
        fitting_group_forces_pytorch.append(f2)
    # read the forces from the debug Colvars log file
    main_group_forces_colvars, fitting_group_forces_colvars = scan_force_from_log('AutoDiff/test.restart.colvars.out')
    # compare the force (calculate the relative error)
    for i in range(len(main_group_forces_colvars)):
        error_main = np.mean(np.abs(main_group_forces_colvars[i] - main_group_forces_pytorch[i]) / np.abs(main_group_forces_colvars[i]))
        error_fitting = np.mean(np.abs(fitting_group_forces_colvars[i] - fitting_group_forces_pytorch[i]) / np.abs(fitting_group_forces_colvars[i]))
        print(f'Frame {i}: error of the main group forces {error_main:12.5e} , fitting group forces {error_fitting:12.5e}')

from numpy import cos, sin, eye, ones
from math import pi
import numpy as np
from classes.data.data_structs import RotorData
from utils import skew
def  fun_blades2rotor(blade_aero, blade_struc, n_pitch, a_pitch, n_rotor, a_rotor, n_tilt, a_tilt, x_r):
    # translating, pitching and rotating of blade to bring in rotor position
    # rotation matrix around pitch axis
    assert n_pitch.ndim == 2, "n_pitch must be 1x3"
    assert n_rotor.ndim == 2, "n_rotor must be 1x3"
    assert n_tilt.ndim == 2, "n_tilt must be 1x3"
    Rp = cos(a_pitch*pi/180)*eye(3) + \
                sin(a_pitch*pi/180)*skew(n_pitch)+ \
                (1-cos(a_pitch*pi/180))*n_pitch.T@n_pitch
    # rotation matrix around rotor axis
    Rr = cos(a_rotor*pi/180)*eye(3) + \
                sin(a_rotor*pi/180)*skew(n_rotor)+ \
                (1-cos(a_rotor*pi/180))*n_rotor.T @ n_rotor;
    # rotation matrix around tilt
    Rt = cos(a_tilt*pi/180)*eye(3) +  \
               sin(a_tilt*pi/180)*skew(n_tilt)+ \
                (1-cos(a_tilt*pi/180))*n_tilt.T @ n_tilt;

    # print(f'Rr = {Rr}')
    # print(f'Rp = {Rp}')
    # print(f'Rt = {Rt}')
    # calculate new position vector for aero grid in global cos
    nnodes = (blade_aero.M+1)*(blade_aero.N+1);
    temp = ones((nnodes,1))
    X_R = np.hstack([ temp * x_r[0],
                        temp * x_r[1],
                        temp * x_r[2]])
    blade_ro_cos = RotorData()
    blade_ro_cos.aero.X_C  = (Rt @ Rr @ (Rp @ blade_aero.X_C.T + X_R.T)).T

    nnodes = (blade_aero.M+1)*(2*blade_aero.N+1);
    temp = ones((nnodes,1))
    X_R = np.hstack([ temp * x_r[0],
                        temp * x_r[1],
                        temp * x_r[2]])
    blade_ro_cos.aero.X_W  = (Rt @ Rr @ (Rp @ blade_aero.X_W.T + X_R.T)).T

    nnodes = (blade_struc.M+1);
    temp = ones((nnodes,1))
    X_R = np.hstack([ temp * x_r[0],
                        temp * x_r[1],
                        temp * x_r[2]])
    
    blade_ro_cos.struc.X_RE = (Rt @ Rr @                                 (Rp @ blade_struc.arr_coordinates[:,:3].T                                      + X_R.T)).T
    blade_ro_cos.struc.D1 = (Rt @ Rr @ Rp @                                 blade_struc.arr_coordinates[:, 3:6].T).T
    blade_ro_cos.struc.D2 = (Rt @ Rr @ Rp @                                 blade_struc.arr_coordinates[:, 6:9].T).T
    blade_ro_cos.struc.D3 = (Rt @ Rr @ Rp @                                 blade_struc.arr_coordinates[:,9:12].T).T
    return blade_ro_cos

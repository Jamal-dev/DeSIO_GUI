from __future__ import annotations
from dataclasses import dataclass
from dataclasses import field
from typing import List
import os
import numpy as np
import yaml
from pathlib import Path
import json
with open(Path(os.getcwd()) / Path(f'simu_settings.json'),'r') as f:
    simu_settings = json.load(f)
    simu_settings['d_vinf'] = np.array(simu_settings['d_vinf'])

@dataclass
class HubNacceleSettings:
    """Class for keeping track hub and ncelle settings."""
    type: str = 'nacelle'
    center_mass: np.array = np.array([0,0,0])
    D1: np.array = np.array([0,0,0])
    D2: np.array = np.array([0,0,0])
    D3: np.array = np.array([0,0,0])
    mass_matrix: np.array = np.eye(3)


@dataclass
class SimulationSettings:
    """Class for keeping track of simulation settings."""
    @dataclass
    class PointMass:
        mass:float = 0.0
        node:int = 0
    @dataclass
    class Constraint:
        type:str = 'fixed'
        nodes:int = 0
        phi1:np.array = np.empty((3,))
        phi2:np.array = np.empty((3,))
        dir:np.array = np.empty((3,))
    @dataclass
    class MaterialBeams:
        cmat:np.array = np.eye(3)
        mmat:np.array = np.eye(3)
        diss:np.array = np.eye(3)
        strname:str = 'steel'
    pointmass:PointMass = PointMass()
    constraint:Constraint = Constraint()
    materialbeam:MaterialBeams = MaterialBeams()
    rb:list[HubNacceleSettings] = field(default_factory=lambda: [])
    pointmass12:list[PointMass] = field(default_factory=lambda: [])
    constraints:list[Constraint] = field(default_factory=lambda: [])
    matbeam:list[MaterialBeams] = field(default_factory=lambda: [])
    currDir:str = simu_settings['currDir']
    strfilename:str = simu_settings['strfilename']
    time:float = simu_settings["time"]
    deltat:float = simu_settings["deltat"]
    cutoff:float = simu_settings["cutoff"]
    density:float = simu_settings["density"]
    i_vinf:float = simu_settings["i_vinf"]
    d_vinf:np.array = np.array(simu_settings["d_vinf"])
    tol:float = simu_settings["tol"]
    niter:int = simu_settings["niter"]
    grav:float = simu_settings["grav"]
    grav_vec:np.array = np.array(simu_settings["grav_vector"])
    loads:np.array = np.array([])
    simType:str = 'dynamic'
    sort:str = 'file'
    wind_field_file:str = ''
    grid_center:np.array = np.array([])    
    flag_matbeam:bool = True
    jobname:str = ''


@dataclass
class FSISettings:
    """Class for keeping track of FSI settings."""
    radius_rbf:float = 0.0
    lf_type:str = 'constant'
    lf_duration:float = 10.0
    flag_linearization:bool = False
    currDir:str = simu_settings['currDir']
    strfilename:str = simu_settings['strfilename']
    data:List = field(default_factory=list)
    data_input_node_fsi_radius:List = field(default_factory=list)

class SimulationSetting_speical(SimulationSettings):
    @dataclass
    class Boundary:
        sort:str = ''
        intensity:np.array = np.zeros(3)
        time:float = 0.0
        constraintID:int = 0
    def __init__(self, *args, **kwargs):
        super(SimulationSetting_speical,self).__init__(*args, **kwargs)
        self.type = ''
        self.settings = np.array([])
        self.remove_aero_attributes()
    def remove_aero_attributes(self):
        att_remove_list = ["cutoff", "density", 
                            "i_vinf", "d_vinf",
                            "sort", "wind_field_file",
                            "grid_center"]
        for att in att_remove_list:
            if hasattr(self, att):
                delattr(self, att)
    def add_dynamic_attributes(self):
        self.boundary12 = self.Boundary()
        self.prevjobname = ''

    def copy(self, other):
        
        for k, v in other.__dict__.items():
            if hasattr(self, k) and hasattr(other, k):
                setattr(self, k, v)
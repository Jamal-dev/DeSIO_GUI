from dataclasses import dataclass
from dataclasses import field
import numpy as np
from copy import deepcopy
from typing import List, Any


@dataclass
class MatProp:
    n_gl: np.array = np.array([0,0,0])
    inz : np.array = np.array([0,0,0])
    Cmat : np.array = np.array([0,0,0])
    mmat : np.array = np.array([0,0,0])
    diss : np.array = np.array([0,0,0])

@dataclass
class MeshInfo:
    """Class for keeping track of mesh information."""
    M:int = 0
    N:int = 0
    connectivity: np.array = np.array([0,0,0])
    def init_aero_params(self):
        self.X = np.array([0,0,0])
        self.node_fsi_radius = np.array([])
    def init_struc_params(self):
        self.matBeam = deepcopy(MatProp())
        self.D1 = np.array([1,0,0])
        self.D2 = np.array([0,1,0])
        self.D3 = np.array([0,0,1])
        self.X_RE = np.array([0,0,0])

@dataclass
class Grid:
    """Class for keeping track of grid settings."""
    struc: MeshInfo = deepcopy(MeshInfo())
    aero : MeshInfo = deepcopy(MeshInfo())
@dataclass
class Monopile:
    """Class for keeping track of monopile mesh."""
    grid: Grid = deepcopy(Grid())

@dataclass
class MeshStruct:
    mx:int = 0
    nn:int = 0
    connectivity: np.array = np.array([[0,0,0]])
    imat:int = 0
    nodes: np.array = np.array([0,0,0])
    strname: str = ''

@dataclass
class Tower:
    """Class for keeping track of monopile mesh."""
    grid: Grid = Grid()


@dataclass
class RotorData:
    @dataclass
    class RotorDataAero:
        X_C:np.array = np.array([])
        X_W:np.array = np.array([])
    @dataclass
    class RotorDataStruc:
        X_RE:np.array = np.array([])
        D1:np.array = np.array([])
        D2:np.array = np.array([])
        D3:np.array = np.array([])
        matBeam:MatProp = MatProp()
    aero:RotorDataAero = RotorDataAero(  X_C = np.array([])
                         , X_W  = np.array([]))
    struc:RotorDataStruc = RotorDataStruc(X_RE = np.array([]), 
                            D1 = np.array([]), 
                            D2 = np.array([]), 
                            D3 = np.array([]),  
                            matBeam = MatProp())

@dataclass
class WakeObj:
    inf:np.array = np.array([])
    nsurf:int = 0
    nsegments:int = 0
    nproperty:int = 0
@dataclass
class Wake:
    
    wakes:List[WakeObj] = field(default_factory=lambda: [WakeObj()])
    property:np.array = np.empty((1,3)) 

class Blade:
    """Class for blade data"""
    def __init__(self) -> None:
        struc = deepcopy(MeshInfo(M=0,
                         N=0,
                         connectivity=np.array([])))
        self.grid: Grid = deepcopy(Grid(struc=struc, aero=struc))
    def __repr__(self) -> str:
        return f"Blade: {self.grid}"

@dataclass
class BladeAero:
    X:np.array = np.array([])
    connectivity:np.array = np.array([])
    M:int = 0
    N:int = 0
@dataclass
class Blade2:

    @dataclass
    class BladeGrid:
        aero:BladeAero = BladeAero(X = np.array([]),
                                  connectivity = np.array([]),
                                  M = 0,
                                  N = 0)
    temp = BladeAero(X = np.array([]),
                                  connectivity = np.array([]),
                                  M = 0,
                                  N = 0)
    grid:BladeGrid = BladeGrid(aero = temp)
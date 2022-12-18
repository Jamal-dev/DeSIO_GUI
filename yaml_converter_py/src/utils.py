import yaml
from yaml.loader import SafeLoader
from scipy.interpolate import interp1d
import numpy as np
from typing import List, Any

def yaml_loader(filepath: str) -> dict:
    with open(filepath, 'r') as file_descriptor:
        data = yaml.load(file_descriptor, Loader=SafeLoader)
    return data

def if_assign(field, subkey, xq, sep_values=None):
    """
        This function checks if the subkey is available in the field.keys
        if it does the it does linear interpolation between grid and values
        Inputs:
            field: dict
            subkey: str
            xq: querry points
            sep_values: If this is defined it will treated as y for the interpolation
    """
    if sep_values is None:
        if subkey in field.keys():
            x = field[subkey]["grid"]
            y = field[subkey]["values"]
            fun = interp1d(x,y)
            return fun(xq)
    else:
        if subkey in field.keys():
            x = field[subkey]["grid"]
            y = sep_values
            fun = interp1d(x,y)
            return fun(xq)
def interp1D(array, xq):
    """
        This function does linear interpolation between grid and values
        Inputs:
            array: array
            xq: querry points
    """
    x = array[:,0]
    y = array[:,1]
    fun = interp1d(x,y)
    return fun(xq)


def skew(a):
    if a.ndim == 2:
        a = a.flatten()
    return np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])

class Helpers:
    def __init__(self):
        pass
    @staticmethod
    def hstack( temp, x_r):
        if x_r.ndim ==2:
            x_r = x_r.flatten()
        return np.hstack([ temp * x_r[0],                             temp * x_r[1],                                 temp * x_r[2]])
    @staticmethod
    def try_append(obj:List, idx:int, value:Any):
        try: 
            obj[idx] = value
        except IndexError:
            obj.append(value)
        return obj

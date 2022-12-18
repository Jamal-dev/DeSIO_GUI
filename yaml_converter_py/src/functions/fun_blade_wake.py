import logging
from copy import deepcopy, copy
from classes.data.data_structs import Wake
import numpy as np
from utils import skew
def   fun_blade_wake(wake, wakesurf, obj, tnrows, nrows, wake_diffusion):
    logging.info('inside fun_blade_wake')
    # define seperation edge
    if not wake:
        w = 0
        wake = Wake()
        wake.wakes.append(deepcopy(wake.wakes[0]))
        wake.property = np.append(wake.property, [[-1,-1,-1]], axis=0)
        logging.debug(f'wake was empty, creating new wake')
        logging.debug(f'wake index: w = {w}')
        logging.debug(f'type(wake.wakes) = {type(wake.wakes)}')
        logging.debug(f'len(wake.wakes) = {len(wake.wakes)}')
    else:
        w = len(wake.wakes)
        logging.debug(f'wake was not empty, appending new wake')
        logging.debug(f'wake index: w = {w}')
        wake.wakes.append(deepcopy(wake.wakes[0]))
        wake.property = np.append(wake.property, [[-1,-1,-1]], axis=0)
        wake.wakes.append(deepcopy(wake.wakes[0]))
        wake.property = np.append(wake.property, [[-1,-1,-1]], axis=0)
        logging.debug(f'wake property.shape: {wake.property.shape}')
    
    
    
    a = 0;
    # logging.debug(f'(obj.mx+1) = {obj.mx+1}')
    # logging.debug(f'(obj.mx+1)*(obj.my+1)-1 = {(obj.mx+1)*(obj.my+1)-1}')
    # logging.debug(f'obj.connectivity.shape = {obj.connectivity.shape}')
    # logging.info('For loop started')
    iter_ = np.arange( (obj.mx+1)-1,(obj.mx+1)*(obj.my+1)-1, (obj.mx+1))
    wake.wakes[w].inf = np.zeros((len(iter_), 3))
    connectivity = copy(obj.connectivity)
    connectivity = connectivity.flatten('F')
    for j in iter_:
        
        # logging.debug(f'j = {j}, j+(obj.mx+1) = {j+(obj.mx+1)}')
        temp = np.argwhere(obj.connectivity == j)
        row1 = temp[:,0]
        # logging.debug(f'obj.connectivity[row1,:] = {obj.connectivity[row1,:]}')
        temp = np.argwhere(obj.connectivity[row1,:] == j+(obj.mx+1))
        row2 = temp[:,0]
        
        if not row2: 
            row2 = 0
        # logging.debug(f'j= {j}, row1 = {row1}')
        # logging.debug(f'j= {j}, row2 = {row2}')
        # logging.debug(f'j= {j}, row1[row2] = {row1[row2]}')
        wake.wakes[w].inf[a,:] = np.array([j,j+obj.mx+1, row1[row2]]) 
        a = a + 1;
    
    # logging.info('For loop ended')
    wake.wakes[w].nsurf     = wakesurf;
    wake.wakes[w].nsegments = a;
    wake.wakes[w].nproperty = w +1;

    # define wake properties
    wake.property[w,:3] = [tnrows,nrows, wake_diffusion];
    # wake at blade tip
    a = 0;
    iter_ = np.arange((obj.mx+1)*(obj.my+1)-1,(obj.mx+1)*obj.my+1+1,-1)
    # logging.debug(f'wake.wakes[w+1].inf = {wake.wakes[w+1].inf}')
    # logging.debug(f'type(wake.wakes[w+1].inf) = {type(wake.wakes[w+1].inf)}')
    wake.wakes[w+1].inf = np.zeros((len(iter_), 3))
    # logging.debug(f'obj.connectivity[:,-1] = {obj.connectivity[:,-1]}')
    for j in iter_:
        # logging.debug(f'j = {j}')
        # logging.debug(f'j = {j}, j+(obj.mx+1) = {j+(obj.mx+1)}')
        temp = np.argwhere(obj.connectivity == j)
        row1 = temp[:,0]
        # logging.debug(f'obj.connectivity[row1,:] = {obj.connectivity[row1,:]}')
        temp = np.argwhere(obj.connectivity[row1,:] == j-1)
        row2 = temp[:,0]
        if not row2: 
            row2 = 0
        # logging.debug(f'row1 = {row1}')
        # logging.debug(f"[j,j-1, row1[row2]] = {[j,j-1, row1[row2]]}")
        wake.wakes[w+1].inf[a,:] = np.array([j,j-1, row1[row2]])
        a = a + 1;
    wake.wakes[w+1].nsurf     = wakesurf;
    wake.wakes[w+1].nsegments = a;
    wake.wakes[w+1].nproperty = w + 1;
    wake.property[w+1,:3] = [tnrows,nrows, wake_diffusion];
    logging.info('fun_blade_wake done')
    return wake
import numpy as np
from utils import if_assign, interp1D
from scipy.interpolate import interp1d
import logging
class UVLM:
    def __init__(self, strName,model_uvlm,airfoils):
        # span-wise and chord-wise discretization
        self.M = model_uvlm["M_aero"]
        self.N = model_uvlm["N_aero"]
        
        # natural coordinates of span-and chord-wise discretization
        self.xhi_x  = np.arange(0,1 + 1/self.M,1/self.M) # span-wise
        self.xhi_y_c  = np.arange(0,1 + 1/self.N,1/self.N) # chord-wise

        self.connectivity_c = []
        self.connectivity_w = []
        
        
    

        nbr             = 0;                    # number of blade root cross-section
        arr_inz_airfoil = [];                   # index array for locating airfoils in span-wise direction
        
        # interpolating coordinates of reference axis according to span-wise discretization
        self.arr_xre_x = if_assign(model_uvlm["reference_axis"], "x", self.xhi_x)
        self.arr_xre_y = if_assign(model_uvlm["reference_axis"], "y", self.xhi_x)
        self.arr_xre_z = if_assign(model_uvlm["reference_axis"], "z", self.xhi_x)
        arr_val_c = []
        arr_val_w = []
        
        self.arr_node_fsi_radius_w = []
        self.arr_node_fsi_radius = []
        # if-statement according to type of surface cross-section
        if "blade" in strName:
            self.airfoil = model_uvlm["airfoil_position"]
            nbr              = model_uvlm["blade_root_position"]

            # interpolating chord length in span-wise direction
            self.arr_c = if_assign(model_uvlm, "chord", self.xhi_x)

            # interpolating twist angle according to span-wise discretization
            self.arr_twist = np.zeros((self.M,1))
            
            self.arr_twist = if_assign(model_uvlm, "twist", self.xhi_x)
            
            # interpolating airfoil in span-wise direction. This is important
            # to identify the blade root and blade's lifting surfaces
            airfoil_grid    = self.airfoil["grid"]
            # TODO:
#             airfoil_values  = [1:length(uvlm_obj.airfoil.grid)]
#             print(self.airfoil)
#             print('*'*80)
#             print(len(airfoils))
#             print(airfoils[0].keys())
            airfoil_values  = range(len(self.airfoil["grid"]))
            fun = interp1d(airfoil_grid,airfoil_values)
            arr_inz_airfoil = np.fix(fun(self.xhi_x))

            arr_pitch_ax = np.ones((self.M+1,1))*1.0;
            self.arr_pitch_ax = if_assign(model_uvlm, "pitch_axis",self.xhi_x)
            
            # natural coordinates of chord-wise discretization for whole
            # surface of cross-section
            xhi_y_w = self.xhi_y_c;
        elif "pipe" in strName:
            logging.debug("UVLM class:\nInside pipe:")
            # interpolating chord length in span-wise direction
            temp = np.asarray(model_uvlm["outer_diameter"]["values"])                 - np.asarray(model_uvlm["thickness"]["values"])
            logging.debug('model_uvlm["outer_diameter"]["values"]: ' + str(model_uvlm["outer_diameter"]["values"]) )
            logging.debug('model_uvlm["thickness"]["values"]: ' + str(model_uvlm["thickness"]["values"]) )
            self.arr_c          = if_assign(model_uvlm,
                                            "outer_diameter",
                                            self.xhi_x, # these are query points
                                            temp # this is y for interpolation
                                           );
            self.arr_twist      = np.zeros((self.M+1,1)) # zero twist
            self.arr_pitch_ax   = np.ones((self.M+1,1))*0.5; # location of pitch axis
            self.airfoil = {}
            self.airfoil["labels"] = ['circular','circular']  # artificial "airfoil" sections for pipe
            self.airfoil["grid"]   = [0.0, 1.0];               # artificial "airfoil" sections grid for pipe

            # mesh discretization around circular cross-section
            temp = np.linspace(0,np.pi, self.N+1)
            xhi_y_w = 0.5*(1-np.cos(temp))
        # loop over airfoils in the model to calculate coordinates of
        # cross-section surfaces
        for i in range(len(self.airfoil["labels"])):
            airfoil_name = self.airfoil["labels"][i]

            # searching current airfoil from airfoil-list
            for j in range(len(airfoils)):
                if len(airfoils) == 1:
                    if airfoil_name in airfoils[j]["name"]:
                        airfoilj = airfoils;
                        break

                else:
                    if airfoil_name in airfoils[j]["name"]  :
                        airfoilj = airfoils[j]
                        break




            # extracting coordinates for upper and lower airfoil 
            airfoil_coord = np.vstack([airfoilj["coordinates"]["x"],airfoilj["coordinates"]["y"]]).T
            logging.debug(f"airfoil_coord.shape: {airfoil_coord.shape}")
            inz0 = np.nonzero(airfoil_coord[:,0]==0) 
            inz0 = inz0[0][0] + 1
            logging.debug(f"inz0: {inz0}")
            airfoil_coord_u = airfoil_coord[:inz0,:]; 
            inzsort = np.argsort(airfoil_coord_u[:,0])
            airfoil_coord_u = airfoil_coord_u[inzsort,:]
            airfoil_coord_l = airfoil_coord[inz0-1:,:]; 
            inzsort = np.argsort(airfoil_coord_l[:,0]) 
            airfoil_coord_l = airfoil_coord_l[inzsort,:]

            # interpolating values in airfoil thickness direction according to
            # discretization in chord-wise direction to calculate camber
            # surface coordinates
#             print("airfoil_coord_u: ",airfoil_coord_u)
            xhi_airf_u = interp1D(airfoil_coord_u,self.xhi_y_c);
            xhi_airf_l = interp1D(airfoil_coord_l,self.xhi_y_c);
            # coordinates of camber surface
#             print(xhi_airf_u.shape, xhi_airf_l.shape)
            xhi_airf_c = (xhi_airf_u + xhi_airf_l)/2;

            # interpolating values in airfoil thickness direction according to
            # discretization in chord-wise direction to calculate whole
            # cross-section surface coordinates
            xhi_airf_u = interp1D(airfoil_coord_u,xhi_y_w);
            xhi_airf_l = interp1D(airfoil_coord_l,xhi_y_w);
            # coordinates of whole cross-section surface
            xhi_airf_w = np.append(np.append(xhi_airf_u[:-1],(xhi_airf_u[-1] + xhi_airf_l[-1])/2),xhi_airf_l[-2::-1])
#             print("xhi_airf_w.shape: ",xhi_airf_w.shape)
            
            # if-statement to detect, if blade root s. In case blade root
            # , switching from whole surface to camber surface
            if i == nbr+1:
                if nbr != 0:
                    xhi_airf_w = np.append(xhi_airf_c[:-1],xhi_airf_c[::-1])
#                     print("inside: xhi_airf_w.shape: ",xhi_airf_w.shape)



#             print("xhi_airf_c.shape :", xhi_airf_c.shape)
            arr_val_c.append( xhi_airf_c)      # camber surface coordinates 
            arr_val_w.append( xhi_airf_w) # whole surface coordinates
        arr_val_c = np.asarray(arr_val_c).T
        arr_val_w = np.asarray(arr_val_w).T
        

        # interpolating camber surface coordinates in span-wise direction 
        arr_xhi_airf_c = np.zeros((self.xhi_x.size, arr_val_c.shape[0]))
        for i in range(arr_val_c.shape[0]):
            x = np.asarray(self.airfoil["grid"])
            y = arr_val_c[i,:]
            # print("line 165: (x.shape, y.shape): ",x.shape, y.shape)
            fun = interp1d(x,y)
            arr_val = fun(self.xhi_x)
            arr_xhi_airf_c[:,i] = arr_val


        # interpolating whole surface coordinates in span-wise direction
        arr_xhi_airf_w =  np.zeros((self.xhi_x.size, arr_val_w.shape[0]))
        for i in range(arr_val_w.shape[0]):
            x = self.airfoil["grid"]
            y = arr_val_w[i,:]
            fun = interp1d(x,y)
            arr_val = fun(self.xhi_x)
            
            arr_xhi_airf_w[:,i] = arr_val


        if nbr!=0:
            inz  = np.nonzero(arr_inz_airfoil<=nbr)
            inz = inz[0]
            temp_coo = arr_xhi_airf_c[inz[-1]+1,:]
            arr_xhi_airf_w[inz[-1]+1,:] = np.append(temp_coo,temp_coo[-2::-1])
        self.arr_xhi_airf_c  = arr_xhi_airf_c;
        self.arr_xhi_airf_w  = arr_xhi_airf_w;
        self.xhi_y_w         = xhi_y_w;
        self.arr_xhi_x       = self.xhi_x
        self.nbr             = nbr;
        self.arr_inz_airfoil = arr_inz_airfoil;

        # elements in x (span) and y (chord) direction
        mx = len(self.arr_xhi_x) - 1
        my = self.N
        max_dim = self.arr_c.shape[0]
        dim_w = max_dim*(2*my+1)
        dim_c = max_dim*(my+1)
        self.X_W = np.zeros((dim_w,3))
        self.X_C = np.zeros((dim_c,3))
        self.X_0 = np.zeros((dim_c,3))
        self.X_0_W = np.zeros((dim_w,3))
    def __repr__(self) -> str:
        output = {"M": self.M,
                  "N": self.N,
                  "arr_xre_x": self.arr_xre_x,
                  "arr_xre_y": self.arr_xre_y,
                  "arr_xre_z": self.arr_xre_z,
                  "airfoil": self.airfoil,
                  "arr_c": self.arr_c,
                  "arr_twist": self.arr_twist,
                  "arr_pitch_ax": self.arr_pitch_ax,
                  "arr_xhi_airf_c": self.arr_xhi_airf_c,
                  "arr_xhi_airf_w": self.arr_xhi_airf_w,
                  "xhi_y_w": self.xhi_y_w,
                  "arr_xhi_x": self.arr_xhi_x,
                  "nbr": self.nbr,
                  "arr_inz_airfoil": self.arr_inz_airfoil,
                  "X_W": self.X_W,
                  "X_C": self.X_C,
                  "X_0": self.X_0,
                  "X_0_W": self.X_0_W,
                  "arr_node_fsi_radius_w": self.arr_node_fsi_radius_w,
                  "arr_node_fsi_radius": self.arr_node_fsi_radius,
                  "connectivity_c": self.connectivity_c,
                  "connectivity_w": self.connectivity_w,}

        return ', \n'.join('{}: {}'.format(k, v) for k, v in output.items())

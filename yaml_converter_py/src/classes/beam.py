from scipy.interpolate import interp1d
from math import pi, cos, sin
import numpy as np

class Beam:
    def __init__(self,strName,model_beam, materials):
        
        self.strName = strName
        self.materials = materials
        # span wise discritization
        self.model_beam = model_beam
        self.M = model_beam["M_struc"]
        # initializing stiffness and mass matrix
        self.arr_stiff_matrix = np.zeros((self.M,21))
        self.arr_mass_matrix = np.zeros((self.M,6))
        self.O = np.zeros((self.M,1));
        
        # nodal natural coordinates in span-wise direction
        self.xhi_x  = np.arange(0,1 + 1/self.M,1/self.M)
        # print(f"xhi_x.shape : {self.xhi_x.shape}")
        
        # element natural coordinates in span-wise direction
        self.xhi_x_elem = (self.xhi_x[0:-1]+self.xhi_x[1:])/2
        # print(f"xhi_x_elem.shape : {self.xhi_x_elem.shape}")
        self.interpolate_ref_values()
        
        self.connectivity = np.zeros((self.M,3))
        self.arr_coordinates = [] 
    def __repr__(self) -> str:
        output_str = self.strName + "("
        output = dict()
        if "blade" in self.strName:
            output = {"M": self.M, 
                      "arr_stiff_matrix" : self.arr_stiff_matrix,
                      "arr_mass_matrix" : self.arr_mass_matrix, 
                    "arr_xre_x" : self.arr_xre_x, 
                    "arr_xre_y" : self.arr_xre_y,
                    "arr_xre_z" : self.arr_xre_z,
                    "arr_twist" : self.arr_twist,
                    "dissipation" : self.dissipation,
                    "arr_coordinates" : self.arr_coordinates,
                    "connectivity" : self.connectivity}
        elif "pipe" in self.strName:
            output = {"M": self.M, 
                        "arr_stiff_matrix" : self.arr_stiff_matrix,
                        "arr_mass_matrix" : self.arr_mass_matrix, 
                        "arr_xre_x" : self.arr_xre_x, 
                        "arr_xre_y" : self.arr_xre_y,
                        "arr_xre_z" : self.arr_xre_z,
                        "arr_outer_diameter" : self.arr_outer_diameter,
                        "arr_thickness" : self.arr_thickness,
                        "arr_EA" : self.arr_EA,
                        "arr_GA1" : self.arr_GA1,
                        "arr_GA2" : self.arr_GA2,
                        "arr_EI1" : self.arr_EI1,
                        "arr_EI2" : self.arr_EI2,
                        "arr_GI3" : self.arr_GI3,
                        "arr_ES1" : self.arr_ES1,
                        "arr_ES2" : self.arr_ES2,
                        "arr_GS2" : self.arr_GS2,
                        "arr_GS1" : self.arr_GS1,
                        "arr_EI12" : self.arr_EI12,
                        "arr_rhoA" : self.arr_rhoA,
                        "arr_rhoI1" : self.arr_rhoI1,
                        "arr_rhoI2" : self.arr_rhoI2,
                        "arr_rhoI12" : self.arr_rhoI12,
                        "arr_rhoS1" : self.arr_rhoS1,
                        "arr_rhoS2" : self.arr_rhoS2,
                        "arr_twist" : self.arr_twist,
                        "dissipation" : self.dissipation,
                        "arr_coordinates" : self.arr_coordinates,
                        "connectivity" : self.connectivity}
        output_str = output_str +                     ', \n'.join(f'{key}:{value}' for key, value in output.items())                     + ")"
                
        return output_str
    def ifassign(self,field, subkey, xq):
        """
            This function checks if the subkey is available in the field.keys
            if it does the it does linear interpolation between grid and values
            Inputs:
                field: dict
                subkey: str
                xq: querry points
        """
        if subkey in field.keys():
            x = field[subkey]["grid"]
            y = field[subkey]["values"]
            fun = interp1d(x,y)
            return fun(xq)
    def interpolate_cross_sectional(self):
        # if-statement according to type of surface cross-section
        if "blade" in self.strName:
            # interpolating along span-wise direction
            if 'six_x_six' in self.model_beam["elastic_properties_mb"]:
                x_el_phy = self.model_beam["elastic_properties_mb"]["six_x_six"]["stiff_matrix"]["grid"]
                
                self.model_beam["elastic_properties_mb"]["six_x_six"]["stiff_matrix"]["values"] = np.asarray(self.model_beam["elastic_properties_mb"]["six_x_six"]["stiff_matrix"]["values"])
                for i in range(21):
                    
                    
                    x_el_ref = self.model_beam["elastic_properties_mb"]["six_x_six"]["stiff_matrix"]["values"][:,i]
                    
                    el_model = interp1d(x_el_phy, x_el_ref)
                    self.arr_stiff_matrix[:,i] = el_model(self.xhi_x_elem)
                x_in_phy = self.model_beam["elastic_properties_mb"]["six_x_six"]["inertia_matrix"]["grid"]
                self.model_beam["elastic_properties_mb"]["six_x_six"]["inertia_matrix"]["values"] = np.asarray(self.model_beam["elastic_properties_mb"]["six_x_six"]["inertia_matrix"]["values"])
                for i in range(6):
                    x_in_ref = self.model_beam["elastic_properties_mb"]["six_x_six"]["inertia_matrix"]["values"][:,i]
                    in_model = interp1d(x_in_phy, x_in_ref)
                    self.arr_mass_matrix[:,i] = in_model(self.xhi_x_elem)
                self.dissipation = np.asarray([self.model_beam["dissipation"]["alpha_s"], self.model_beam["dissipation"]["alpha_v"]])
            else:
                self.arr_EA  = np.zeros((self.M,1)); 
                self.arr_GA1 = np.zeros((self.M,1)); 
                self.arr_GA2 = np.zeros((self.M,1));
                self.arr_EI1 = np.zeros((self.M,1)); 
                self.arr_EI2 = np.zeros((self.M,1)); 
                self.arr_GI3 = np.zeros((self.M,1));
                self.arr_ES1 = np.zeros((self.M,1)); 
                self.arr_ES2 = np.zeros((self.M,1)); 
                self.arr_EI12 = np.zeros((self.M,1));
                self.arr_GS1 = np.zeros((self.M,1)); 
                self.arr_GS2 = np.zeros((self.M,1));

                self.arr_rhoA  = np.zeros((self.M,1)); self.arr_rhoI1 = np.zeros((self.M,1)); self.arr_rhoI2  = np.zeros((self.M,1));
                self.arr_rhoS1 = np.zeros((self.M,1)); self.arr_rhoS2 = np.zeros((self.M,1)); self.arr_rhoI12 = np.zeros((self.M,1));

                if 'EA' in self.model_beam["elastic_properties_mb"].keys():
                    x = self.model_beam["elastic_properties_mb"]["EA"]["grid"]
                    y = self.model_beam["elastic_properties_mb"]["EA"]["values"]
                    fun = interp1d(x,y)
                    self.arr_EA   = fun(self.xhi_x_elem);


                
                self.arr_GA1  = self.ifassign(self.model_beam["elastic_properties_mb"],'GA1',self.xhi_x_elem)
                self.arr_GA2  = self.ifassign(self.model_beam["elastic_properties_mb"],'GA2',self.xhi_x_elem)
                self.arr_EI1  = self.ifassign(self.model_beam["elastic_properties_mb"],'EI1',self.xhi_x_elem)
                self.arr_EI2  = self.ifassign(self.model_beam["elastic_properties_mb"],'EI2',self.xhi_x_elem)
                self.arr_GI3  = self.ifassign(self.model_beam["elastic_properties_mb"],'GI3',self.xhi_x_elem)
                self.arr_ES1  = self.ifassign(self.model_beam["elastic_properties_mb"],'ES1',self.xhi_x_elem)
                self.arr_ES2  = self.ifassign(self.model_beam["elastic_properties_mb"],'ES2',self.xhi_x_elem)
                self.arr_GS1  = self.ifassign(self.model_beam["elastic_properties_mb"],'GS1',self.xhi_x_elem)
                self.arr_GS2  = self.ifassign(self.model_beam["elastic_properties_mb"],'GS2',self.xhi_x_elem)
                self.arr_EI12  = self.ifassign(self.model_beam["elastic_properties_mb"],'EI12',self.xhi_x_elem)
                self.arr_rhoA  = self.ifassign(self.model_beam["elastic_properties_mb"],'rhoA',self.xhi_x_elem)
                self.arr_rhoI1  = self.ifassign(self.model_beam["elastic_properties_mb"],'rhoI1',self.xhi_x_elem)
                self.arr_rhoI2  = self.ifassign(self.model_beam["elastic_properties_mb"],'rhoI2',self.xhi_x_elem)
                self.arr_rhoS1  = self.ifassign(self.model_beam["elastic_properties_mb"],'rhoS1',self.xhi_x_elem)
                self.arr_rhoS2  = self.ifassign(self.model_beam["elastic_properties_mb"],'rhoS2',self.xhi_x_elem)
                self.arr_rhoI12  = self.ifassign(self.model_beam["elastic_properties_mb"],'rhoI12',self.xhi_x_elem)

                
                # GA1 GA2 EA EI1 EI2 GI3 0 0 0 GS2 -GS1 0 0 0 0 0 ES1 -EI12 -ES2 0 0
                self.arr_stiff_matrix[:,:] = [self.arr_GA1, self.arr_GA2, self.arr_EA, self.arr_EI1, self.arr_EI2, self.arr_GI3,                                                 self.O, self.O, self.O, self.arr_GS2, -self.arr_GS1, self.O, self.O, self.O, self.O, self.O, self.arr_ES1, -self.arr_EI12,                                              -self.arr_ES2, self.O, self.O];

                # rhoA, rhoI2, rhoI1, rhoI12, rhoS1, rhoS2
                self.arr_mass_matrix[:,:]   = [self.arr_rhoA, self.arr_rhoI2, self.arr_rhoI1, self.arr_rhoI12, self.arr_rhoS1, self.arr_rhoS2];
                self.dissipation         = np.asarry([self.model_beam["dissipation"]["alpha_s"], self.model_beam["dissipation"]["alpha_v"]])
        
        if 'pipe' in self.strName:
            # interpolating cross-section properties along span-wise direction
            E = 0.0; G = 0.0; rho = 0.0;

            # self.materials is the list of dictionaries of cross-section properties
            if 'material' in self.model_beam.keys():
                strmaterial = self.model_beam["material"];
                for i in range(len( self.materials )):
                    if strmaterial in self.materials[i]["name"]:
                        if len(self.materials[i]["name"]) == len(strmaterial):
                            E   = self.materials[i]["E"]
                            nu  = self.materials[i]["nu"]
                            rho = self.materials[i]["rho"]
                            G   = E/(2.0*(1.0+nu))




            if 'elastic_properties_mb' in self.model_beam.keys():
                E   = self.model_beam["elastic_properties_mb"]["E"];
                G   = E/(2.0*(1+self.model_beam["elastic_properties_mb"]["nu"]));
                rho = self.model_beam["elastic_properties_mb"]["rho"];


            k1  = 1.0;
            k2  = 1.0;
            if 'shear_factor' in self.model_beam.keys():
                k1 = self.model_beam["shear_factor"]["k1"];
                k2 = self.model_beam["shear_factor"]["k2"];


            self.arr_outer_diameter = self.ifassign(self.model_beam,"outer_diameter",self.xhi_x_elem)
            
            self.arr_thickness      = self.ifassign(self.model_beam,"thickness",self.xhi_x_elem)
            self.arr_EA             = E*pi/4*(self.arr_outer_diameter**2 - (self.arr_outer_diameter-2*self.arr_thickness)**2);
            self.arr_GA1            = k1*G*pi/4*(self.arr_outer_diameter**2 - (self.arr_outer_diameter-2*self.arr_thickness)**2);
            self.arr_GA2            = k2*G*pi/4*(self.arr_outer_diameter**2 - (self.arr_outer_diameter-2*self.arr_thickness)**2);
            self.arr_EI1            = E*pi/64*(self.arr_outer_diameter**4 - (self.arr_outer_diameter-2*self.arr_thickness)**4);
            self.arr_EI2            = E*pi/64*(self.arr_outer_diameter**4 - (self.arr_outer_diameter-2*self.arr_thickness)**4);
            self.arr_GI3            = G*pi/32*(self.arr_outer_diameter**4 - (self.arr_outer_diameter-2*self.arr_thickness)**4);
            self.arr_ES1            = self.arr_EA*0;
            self.arr_ES2            = self.arr_EA*0;
            self.arr_GS2            = self.arr_EA*0;
            self.arr_GS1            = self.arr_EA*0;
            self.arr_EI12           = self.arr_EA*0;

            self.arr_rhoA           = rho*pi/4*(self.arr_outer_diameter**2 - (self.arr_outer_diameter-2*self.arr_thickness)**2);
            self.arr_rhoI1          = rho*pi/64*(self.arr_outer_diameter**4 - (self.arr_outer_diameter-2*self.arr_thickness)**4);
            self.arr_rhoI2          = rho*pi/64*(self.arr_outer_diameter**4 - (self.arr_outer_diameter-2*self.arr_thickness)**4);
            self.arr_rhoI12         = self.arr_rhoA*0;
            self.arr_rhoS1          = self.arr_rhoA*0;
            self.arr_rhoS2          = self.arr_rhoA*0;

            # print(type(self.arr_GA1), self.arr_GA1.shape)
            # print(type(self.arr_GA2), self.arr_GA2.shape)
            # print(type(self.arr_EA), self.arr_EA.shape)
            # print(type(self.arr_EI1), self.arr_EI1.shape)
            # print(type(self.arr_EI2), self.arr_EI2.shape)
            # print(type(self.arr_GI3), self.arr_GI3.shape)
            # print(type(self.arr_ES1), self.arr_ES1.shape)
            # print(type(self.arr_ES2), self.arr_ES2.shape)
            # print(type(self.arr_GS1), self.arr_GS1.shape)
            # print(type(self.arr_GS2), self.arr_GS2.shape)
            # print(type(self.arr_EI12), self.arr_EI12.shape)
            # print(type(self.O), self.O.shape)

            O = self.O[:,0]
            # GA1 GA2 EA EI1 EI2 GI3 0 0 0 GS2 -GS1 0 0 0 0 0 ES1 -EI12 -ES2 0 0
            self.arr_stiff_matrix[:,:] = np.vstack([self.arr_GA1, self.arr_GA2, self.arr_EA, self.arr_EI1, self.arr_EI2, self.arr_GI3,                                             O, O, O,                                                 self.arr_GS2, -self.arr_GS1,                                                     O, O, O, O, O,                                                         self.arr_ES1, -self.arr_EI12, -self.arr_ES2,                                                              O, O]).T

            self.arr_mass_matrix[:,:] = np.vstack([self.arr_rhoA, self.arr_rhoI2,                                                     self.arr_rhoI1, self.arr_rhoI12,                                                         self.arr_rhoS1, self.arr_rhoS2]).T
            self.dissipation       = np.asarray([ self.model_beam["dissipation"]["alpha_s"], self.model_beam["dissipation"]["alpha_v"]])                  
        
    def interpolate_ref_values(self):
        x_phy = self.model_beam["reference_axis"]["x"]["grid"]
        y_phy = self.model_beam["reference_axis"]["y"]["grid"]
        z_phy = self.model_beam["reference_axis"]["z"]["grid"]
        x_ref = self.model_beam["reference_axis"]["x"]["values"]
        y_ref = self.model_beam["reference_axis"]["y"]["values"]
        z_ref = self.model_beam["reference_axis"]["z"]["values"]
        fx = interp1d(x_phy, x_ref)
        fy = interp1d(y_phy, y_ref)
        fz = interp1d(z_phy, z_ref)
        self.arr_xre_x = fx(self.xhi_x)
        self.arr_xre_y = fy(self.xhi_x)
        self.arr_xre_z = fz(self.xhi_x)
        
        # interpolating twist angle according to span-wise discretization
        if 'twist' in self.model_beam.keys():
            phi_phy = self.model_beam["twist"]["grid"]
            phi_ref = self.model_beam["twist"]["values"]
            phi_interp = interp1d(phi_phy, phi_ref)
            self.arr_twist = phi_interp(self.xhi_x);
        self.interpolate_cross_sectional()

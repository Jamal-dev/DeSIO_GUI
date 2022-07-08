import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
import numpy as np

def sind(deg):
    return np.sin(np.deg2rad(deg))
def cosd(deg):
    return np.cos(np.deg2rad(deg))
def Rz(deg):
    """
        Rotation matrix around z-axis
    """
    return np.array([[cosd(deg), -sind(deg), 0],
                     [sind(deg), cosd(deg), 0],
                     [0, 0, 1]])

class RigidBody:
    def __init__(self,num_rigid_bodies:int,num_cross_prop:int):
        """
            num_rigid_bodies: number of rigid bodies
            num_cross_prop: number of cross-properties
        """
        self.num_rigid_bodies = num_rigid_bodies
        self.num_cross_prop = num_cross_prop
        self.masses = []
        self.positions = {}
        self.intertia_tensor_cm = {}
        self.directors_x = {}
        self.directors_y = {}
        self.directors_z = {}
        self.Rz = Rz(180)
        # assume the body's center of mass is at the origin
        self.x_bar = 0.0
        self.y_bar = 0.0
        self.z_bar = 0.0 # this will lead to 0 values of S1, S2, and S3 static first order moments
        # relating the property id annd number of rigid body
        self.prop_id_rigid_bodies_map = {}
        self.file_name =path_main / Path("io/rigidbodyinput.txt")
    
    def set_prop_id_rigid_bodies_map(self,prop_id_rigid_bodies_map:dict=None):
        """
            prop_id_rigid_bodies_map: dictionary of property id and number of rigid body
        """
        if prop_id_rigid_bodies_map is None:
            # we should set then one to one correspondence
            assert self.num_cross_prop == self.num_rigid_bodies, "If the map of rigid body to corresectioal properties is not given. The, the Number of cross-properties and number of rigid bodies should be the same"
            self.prop_id_rigid_bodies_map = {i:i for i in range(self.num_cross_prop)}
        else:
            self.prop_id_rigid_bodies_map = prop_id_rigid_bodies_map
    def set_directors_x(self,directors_x:dict):
        """
            directors_x: dictionary of directors x
        """
        self.directors_x = directors_x
    def set_directors_y(self,directors_y:dict):
        """
            directors_y: dictionary of directors y
        """
        self.directors_y = directors_y
    def set_directors_z(self,directors_z:dict):
        """
            directors_z: dictionary of directors z
        """
        self.directors_z = directors_z
    
    def set_masses(self,masses:list):
        """
            masses: list of masses
        """
        self.masses = masses
    def set_positions(self,positions:dict):
        """
            positions: dictionary of positions
        """
        self.positions = positions
    def set_inertias(self,intertia_tensor_cm:dict):
        """
            intertia_tensor: dictionary of inertia tensors 3x3
            This is defined at the center of mass
        """
        self.intertia_tensor_cm = intertia_tensor_cm
    def transfer_interia(self,interia, mass, position):
        """
            interia: ndarray of inertia tensor
            mass: mass
            position: ndarray of position
            Returns: ndarray of inertia tensor
                    at the r=[r1,r2,r3] away from origin
        """
        I11cm = interia[0,0]
        I22cm = interia[1,1]
        I33cm = interia[2,2]
        I12cm = interia[0,1]
        I13cm = interia[0,2]
        I23cm = interia[1,2]

        I11 = I11cm + mass*(position[1]**2 + position[2]**2)
        I22 = I22cm + mass*(position[0]**2 + position[2]**2)
        I33 = I33cm + mass*(position[0]**2 + position[1]**2)
        I12 = I12cm - mass*position[0]*position[1]
        I13 = I13cm - mass*position[0]*position[2]
        I23 = I23cm - mass*position[1]*position[2]

        return np.array([[I11, I12, I13],
                            [I12, I22, I23],
                            [I13, I23, I33]])
    def transfer_polar(self, interia):
        """
            interia: ndarray of inertia tensor 3x3
            Returns: ndarray of polar inertia tensors 3x3
        """
        I11 = interia[0,0]
        I22 = interia[1,1]
        I33 = interia[2,2]
        I12 = interia[0,1]
        I13 = interia[0,2]
        I23 = interia[1,2]

        I11p = 0.5 * (I22 + I33 - I11)
        I22p = 0.5 * (I11 + I33 - I22)
        I33p = 0.5 * (I11 + I22 - I33)
        I12p = - I12
        I13p = - I13
        I23p = - I23
        

        return np.array([[I11p, I12p, I13p],
                            [I12p, I22p, I23p],
                            [I13p, I23p, I33p]])


    def rotate_interia(self,interia):
        """
            interia: ndarray of inertia tensor
            Returns: ndarray of inertia tensor
                    rotated by Rz(180)
        """
        return (self.Rz.T @ interia) @ self.Rz
    def cal_intertial_tensor_transformed(self):
        """
            intertia_tensor_cm: ndarray
            This is defined at the center of mass
            Returns: nndarray of inertia tensors
                    at the r=[r1,r2,r3] away from origin
        """
        interia_tensor_transfrmed = {}
        polar_interia_tensor_transfrmed = {}
        polar_inerterial_rotated = {}
        for prop_id, interia in self.intertia_tensor_cm.items():
            mass = self.masses[prop_id]
            position = self.positions[prop_id]
            interia_tensor_transfrmed[prop_id] = self.transfer_interia(interia, mass, position)
            polar_interia_tensor_transfrmed[prop_id] = self.transfer_polar(interia)
            polar_inerterial_rotated[prop_id] = self.rotate_interia(polar_interia_tensor_transfrmed[prop_id])
        self.inertia_tensor_transformed = interia_tensor_transfrmed
        self.polar_interial_rotated = polar_inerterial_rotated
    def cal_cross_sectional_properties(self):
        """
            Returns: ndarray of cross-sectional properties
        """
        S1 = 0.0
        S2 = 0.0
        S3 = 0.0
        S23 = 0.0
        S13 = 0.0
        cross_prop = {}
        for prop_id, interia in self.polar_interial_rotated.items():
            arg7 = self.masses[prop_id]
            arg8 = interia[0,0]
            arg9 = interia[1,1]
            arg10 = interia[2,2]
            arg11 = S23
            arg12 = S13
            arg13 = S3
            arg14 = S2
            arg15 = S1
            arg16 = interia[0,1]
            cross_prop[prop_id] = np.asarray([arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16])
        self.cross_proprities = cross_prop
    def cal_position__directors(self):
        """
            This will create a dictioary where position, directors in x direction
            , directors in y direction, ad directors in z direction are stored. 
        """
        position_directors = {}
        for prop_id, arg2 in self.positions.items():
            arg3 = self.directors_x[prop_id]
            arg4 = self.directors_y[prop_id]
            arg5 = self.directors_z[prop_id]

            position_directors[prop_id] = np.hstack([arg2,arg3,arg4,arg5])
        self.position_directors = position_directors
    @staticmethod
    def string(array,if_int= False):
        
        def checkThreshold(X):
            threshold = 1e-11
            X[abs(X)<threshold] = 0.0
            return X
        def convert(X):
            data = ""
            # d0 is added because of the precison
            # for element in X:
                # element = X[i,:]
            element = X
            if if_int:
                data += ' '.join(map(str, element))+ "\n"
            else:    
                data += 'd0 '.join(map(str, element))+ "d0\n"
            return data
        array = checkThreshold(array)
        return convert(array)
    
    def write_data(self):
        """
            Writes the data to a file
        """
        data = ""
        data += "!! rigid body input written by automatic Wind Energy Converter generator\n" 
        data += "!!\n" 
        data += "!! number of rigid bodies (1), number of body properties (2).\n"
        data += f"{self.num_rigid_bodies}\t {self.num_cross_prop}\n"
        for prop_id in range(self.num_rigid_bodies):
            body_num = prop_id + 1
            cross_num = self.prop_id_rigid_bodies_map[prop_id] +1
            data += "!!\n"
            data += "!!\n"
            data += f"!! body {body_num} - : phi (1, 2, 3), d1 (4, 5, 6), d2 (7, 8, 9), d3 (10, 11, 12).\n"
            data += RigidBody.string(self.position_directors[prop_id])
            data += "!!\n"
            data += "!!\n"
            data += f"!! body {body_num} - : property (1)\n"
            data += f"{cross_num}\n"
        for prop_id in range(self.num_rigid_bodies):
            body_num = prop_id + 1
            cross_num = self.prop_id_rigid_bodies_map[prop_id] +1
            data += "!!\n"
            data += "!!\n"
            data += f"!! property body {cross_num} - : mass (1), I_11 (2), I_22 (3), I_33 (4), S_23 (5), S_13 (6), S_3 (7), S_2 (8), S_1 (9), I_12 (10)\n"
            data += RigidBody.string( self.cross_proprities[prop_id])
        with open(self.file_name, 'w') as f:
            f.write(data)
            print(f"Data written to {self.file_name}") 
    def check_properties(self):
        # TODO: Change the print to errors
        if not self.prop_id_rigid_bodies_map:
            # if the map is not created, create it
            print("prop_id_rigid_bodies_map is not created yet")
            return False
        if not self.positions:
            print("positions is not created yet")
            return False
        if not self.directors_x:
            print("directors_x is not created yet")
            return False
        if not self.directors_y:
            print("directors_y is not created yet")
            return False
        if not self.directors_z:
            print("directors_z is not created yet")
            return False
        if not self.masses:
            print("masses is not created yet")
            return False
        if not self.intertia_tensor_cm:
            print("intertia_tensor_cm is not created yet")
            return False
        return True

    def main(self):
        if self.check_properties():
            self.cal_intertial_tensor_transformed()
            self.cal_position__directors()
            self.cal_cross_sectional_properties()
            self.write_data()
        else:
            print("Properties are not created yet")
            return False
        
        

        
if __name__ == "__main__":
    # test
    num_rigid_bodies = 2
    num_cross_prop = 2
    rigid_body = RigidBody(num_rigid_bodies=num_rigid_bodies, num_cross_prop=num_cross_prop)
    
    ############ Data from GUI ##############
    # Enter position of each rigid body
    positions = {0:np.array([4.68800,	0.0	,148.822]),1:np.array( [-10.93417,	0.0,	150.07473])}
    # Enter directors of each rigid body in x direction
    directors_x = {0:np.array([1.0,	0.0,	0.0]),1:np.array([0.707107, 0.707107, 0.0])}
    # Enter directors of each rigid body in y direction
    directors_y = {0:np.array([0.0,	1.0,	0.0]),1:np.array([-0.707107, 0.707107, 0.0])}
    # Enter directors of each rigid body in z direction
    directors_z = {0:np.array([0.0,	0.0,	1.0]),1:np.array([0.0, 0.0, 1.0])}
    # Enter masses of each rigid body
    masses = [820888.0, 190000.0]
    # Enter second moment of interia of each rigid body
    I11_cm = 9437478.00823881
    I22_cm = 9808155.98764257
    I33_cm = 9780449.1838345
    I12_cm = 446585.459454508
    I13_cm = 1960561.61974822
    I23_cm = -232499.00841237
    II_cm = np.asarray([[I11_cm, I12_cm, I13_cm],[I12_cm, I22_cm, I23_cm],[I13_cm, I23_cm, I33_cm]])

    inertias = {0:II_cm,1:II_cm}
    # map of each rigid body to its cross sectional properties # This will change when we have more cross sectional properties than the rigid bodies
    prop_id_rigid_bodies_map = {0:0,1:1}

    ##############################################

    ############## After getting the data these functions needed to be called ##############
    rigid_body.set_positions(positions)
    rigid_body.set_directors_x(directors_x)
    rigid_body.set_directors_y(directors_y)
    rigid_body.set_directors_z(directors_z)
    rigid_body.set_masses(masses)
    rigid_body.set_inertias(inertias)
    rigid_body.set_prop_id_rigid_bodies_map(prop_id_rigid_bodies_map)
    ######################################################################################
    rigid_body.main()

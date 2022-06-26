import numpy as np
class Constraints:
    def __init__(self,num_constraint,
               num_rigid_supports=1,
               num_revolute_joints =1, 
               num_rigid_connection=4):
        
        self.num_constraints = num_constraint
        self.num_rigid_supports = num_rigid_supports
        self.num_revolute_joints = num_revolute_joints
        self.num_rigid_connections = num_rigid_connection
        self.jointTypes = {0:"rigidsupport",1:"rigidconnection",2:"revolutejoint"}
        self.initialize_lists()
        self.internal_nodes_set_by_external = False
    def initialize_lists(self):
        self.from_rse = []
        self.to_rse = []
        # list for revolute joints
        self.from_rje = []
        self.to_rje = []
        self.phi1_rje = []
        self.phi2_rje = []
        self.dir_rje = []
        # lists for rigid connection elements
        self.from_rce = []
        self.to_rce = []
    def setInternalNodes(self,nodes):
        self.internal_nodes_set_by_external= True
        self.internal_nodes = nodes
    def addRigidSupportElement(self, start_node, end_node):
        # it needs information from starting node and the ending node
        if not (isinstance(start_node,int) and isinstance(end_node,int)):
            raise ValueError("Starting node of the elment and ending node should be int")
        self.from_rse.append(start_node)
        self.to_rse.append(end_node)
    
    def addRigidConnectionElement(self, start_node, end_node):
        # it needs information from starting node and the ending node
        if not (((isinstance(start_node,int) and isinstance(end_node,int))) or (isinstance(start_node,np.int64) and isinstance(end_node,np.int64))):
            raise ValueError("Starting node of the elment and ending node should be int")
        self.from_rce.append(start_node)
        self.to_rce.append(end_node)
    def addRevoluteElement(self, start_node, end_node,phi1:list,phi2:list,dir_:list):
        # it needs information from starting node and the ending node
        # phi1, phi2, and dir_ are lists of length 3
        if not (isinstance(start_node,int) and isinstance(end_node,int)):
            raise ValueError("Starting node of the elment and ending node should be int")
        if not (isinstance(phi1,list) and isinstance(phi2,list) and isinstance(dir_,list)):
            raise ValueError("phi1, phi2, and dir should be list of 3 elements")
        if not (len(phi1)==3 and len(phi2)==3 and len(dir_)==3):
            raise ValueError("Length of phi1, phi2, and dir should be 3")
        self.from_rje.append(start_node)
        self.to_rje.append(end_node)
        self.phi1_rje.append(phi1)
        self.phi2_rje.append(phi2)
        self.dir_rje.append(dir_)
        
    def call_header(self):
        return r"!!constraint 12 input for the NREL 15 MW Reference Wind Turbine tower convergence analysis" \
        + "\n!!\n" \
        + r"!!number of constraints for nodes with 12 coordinates" + '\n' \
        + f"{self.num_constraints}" \
        + "\n!!\n!!\n" \
        + r"!!constraints for nodes with 12 coordinates: sort (1), nodes (2, 3), phi1 (4, 5, 6), phi2 (7, 8, 9), dir (10, 11, 12)" + "\n"
        
    def genrate_internal(self):
        self.body_struc = {}
        if self.internal_nodes_set_by_external:
            for n in self.internal_nodes:
                self.body_struc[n] = f"internal\t{n}\t0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\n" 
        else:
            for c in range(self.num_constraints):
                self.body_struc[c+1] = f"internal\t{c+1}\t0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\t0.0d0\n" 
    def default_0Vaues(self,num_elements):
        phi1 = [[0.0,0.0,0.0] for _ in range(num_elements)]
        phi2 = [[0.0,0.0,0.0] for _ in range(num_elements)]
        dir_ = [[0.0,0.0,0.0] for _ in range(num_elements)]
        return phi1, phi2, dir_
    def generate_rigid_support_elements(self):
        # rigid support elements
        phi1, phi2, dir_ = self.default_0Vaues(self.num_rigid_supports)
        self.rigid_support_elements = {}
        for i,(f,t) in enumerate(zip(self.from_rse,self.to_rse)):
            self.rigid_support_elements[f] = {"from":f,"to":t,"phi1":f"{phi1[i][0]}d0\t{phi1[i][1]}d0\t{phi1[i][2]}d0","phi2":f"{phi1[i][0]}d0\t{phi1[i][1]}d0\t{phi1[i][2]}d0","dir":f"{dir_[i][0]}d0\t{dir_[i][1]}d0\t{dir_[i][2]}d0"}
    def generate_revolute_joint_elements(self):
        # revolt joint elements
        self.revolt_joint_elements = {}
        for i,(f,t) in enumerate(zip(self.from_rje,self.to_rje)):
            self.revolt_joint_elements[f] = {"from":f,"to":t,"phi1":f"{self.phi1_rje[i][0]}d0\t{self.phi1_rje[i][1]}d0\t{self.phi1_rje[i][2]}d0","phi2":f"{self.phi1_rje[i][0]}d0\t{self.phi1_rje[i][1]}d0\t{self.phi1_rje[i][2]}d0","dir":f"{self.dir_rje[i][0]}d0\t{self.dir_rje[i][1]}d0\t{self.dir_rje[i][2]}d0"}
    def generate_rigid_connection_elements(self):
        # rigid connection elements
        phi1, phi2, dir_ = self.default_0Vaues(self.num_rigid_connections)
        self.rigid_connection_elements = {}
        for i,(f,t) in enumerate(zip(self.from_rce,self.to_rce)):
            self.rigid_connection_elements[i] = {"from":f,"to":t,"phi1":f"{phi1[i][0]}d0\t{phi1[i][1]}d0\t{phi1[i][2]}d0","phi2":f"{phi1[i][0]}d0\t{phi1[i][1]}d0\t{phi1[i][2]}d0","dir":f"{dir_[i][0]}d0\t{dir_[i][1]}d0\t{dir_[i][2]}d0"}
    def dic2line(self,v):
        line = ""
        for e in v.values():
            line += str(e) + "\t"
        line += "\n"
        return line
    def checkSizes(self):
        if len(self.from_rse) != len(self.to_rse) and len(self.from_rse) != self.num_rigid_supports:
            raise ValueError("Defined number of rigid support elements is not equal to added rigid support elements")
        if len(self.from_rje) != len(self.to_rje) and len(self.from_rje) != self.num_revolute_joints:
            raise ValueError("Defined number of revolt joint elements is not equal to added revolt joint elements")
        if len(self.from_rce) != len(self.to_rce) and len(self.from_rce) != self.num_rigid_connections:
            raise ValueError("Defined number of rigid connection elements is not equal to added rigid connection elements")
    def setup(self):
        # create first all internal nodes
        self.genrate_internal()
        # check if the correct number of elements have been added
        self.checkSizes()
        # genrate elements
        self.generate_rigid_support_elements()
        self.generate_revolute_joint_elements()
        self.generate_rigid_connection_elements()
        # replacing where rigid support elements are defined
        for node,v in self.rigid_support_elements.items():
            self.body_struc[node] = self.jointTypes[0] + "\t" + self.dic2line(v)

        # replacing where revolt joint elements are defined
        for node,v in self.revolt_joint_elements.items():
            self.body_struc[node] = self.jointTypes[2] + "\t" + self.dic2line(v)

        self.header = self.call_header()
        # appending header with rigid connection elements
        for _,v in self.rigid_connection_elements.items():
            self.header += self.jointTypes[1] + "\t" + self.dic2line(v)
    def write(self,filename:str):
        self.setup()
        body = ""
        for _,v in self.body_struc.items():
            body += v
        self.file_content = self.header + body
        with open(filename,'w') as f:
            f.write(self.file_content)
        print(f"Total constraints: {self.num_constraints}")
        print(f"Total rigid support constrains: {self.num_rigid_supports}")
        print(f"Total rigid connections: {self.num_rigid_connections}")
        print(f"Total revolt joints: {self.num_revolute_joints}")
        print(f"Total internal constraints: {self.num_constraints-self.num_revolute_joints-self.num_rigid_supports}")
        print(f"{filename} is successfully written!")





def test():
    cn = Constraints(num_constraint=253,
                num_rigid_supports=1,
                num_revolute_joints =1, 
                num_rigid_connection=4) 
    # adding rigid support
    cn.addRigidSupportElement(3,0)
    # adding revolute element
    cn.addRevoluteElement(1,2,phi1 = [-15.622170,0.000000,1.252730]
                            ,phi2 = [0.0,0.0,0.0]
                            ,dir_ = [-0.994521895368273,    0.000000000000000,    0.104528463267653])
    # adding rigid connection
    cn.addRigidConnectionElement(52,1)
    cn.addRigidConnectionElement(2,53)
    cn.addRigidConnectionElement(2,103)
    cn.addRigidConnectionElement(2,153)
    # generate file
    cn.write('test.txt')

if __name__ == "__main__":
    test()
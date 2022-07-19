from collections import namedtuple
import sys
from pathlib import Path
import os

from sqlalchemy import true
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")/ Path("..")
sys.path.append(str(path_main))
from beam.segment import Segment
from constraints.Constraint import Constraints
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# import networkx as nx

def sind(deg):
    return np.sin(np.deg2rad(deg))
def cosd(deg):
    return np.cos(np.deg2rad(deg))
def tand(deg):
    return np.tan(np.deg2rad(deg))
def atand(x):
    return np.rad2deg(np.arctan(x))
class interpCrossSection():
    """
        This class is for the interpolation between segments which are defiend for a single beam
        inputs:
                segments_list: list of Segment objects
                L: length of the beam
                num_elements: number of elements in the beam
                beam_name: name of the beam, such as SB, S1, S2, etc.
    """
    def __init__(self, segments_list, L,num_elements, beam_name = None, global_id = 1):
        # L is the length of the beam
        self.num_segments = len(segments_list)
        self.segments_list = segments_list
        self.L = L
        self.beam_name = beam_name
        self.num_elements = num_elements
        self.f = None
        self.global_id = global_id
        self.getSegmentData()
        self.get_interp_funcs()
        self.get_values_mid_elements()
        

    def getSegmentData(self):
        """
            It's for the extraction of the data from the Segment objects
        """
        # Returns a list of lists of the form:
        # [[segment_id, length_ratio, diameter_start, diameter_end, thickness_start, thickness_end, density, e, g, alpha_s, alpha_v, scf_start, scf_end],...]
        
        length_ratios = []
        diameter_starts = []
        diameter_ends = []
        thickness_starts = []
        thickness_ends = []
        densities = []
        e_values = []
        g_values = []
        alpha_s_values = []
        alpha_v_values = []
        scf_starts = []
        scf_ends = []

        for segment in self.segments_list:
            length_ratios.append(segment.length_ratio)
            diameter_starts.append(segment.diameter_start)
            diameter_ends.append(segment.diameter_end)
            thickness_starts.append(segment.thickness_start)
            thickness_ends.append(segment.thickness_end)
            densities.append(segment.density)
            e_values.append(segment.e)
            g_values.append(segment.g)
            alpha_s_values.append(segment.alpha_s)
            alpha_v_values.append(segment.alpha_v)
            scf_starts.append(segment.scf_start)
            scf_ends.append(segment.scf_end)
        self.length_ratios = np.asarray(length_ratios, dtype = np.float32)
        self.L = float(self.L)
        self.length_each_segment = self.L * self.length_ratios
        x_start = [0.0]
        x_end = [self.length_each_segment[0]]
        for i in range(self.num_segments-1):
            x_start.append(x_end[i])
            x_end.append(x_end[i] + self.length_each_segment[i+1])
        self.x_start = np.asarray(x_start, dtype = np.float32)
        self.x_end = np.asarray(x_end, dtype = np.float32)
        x = []
        for x1,x2 in zip(self.x_start, self.x_end):
            x.append(x1)
            x.append(x2)
        self.x = np.asarray(x, dtype = np.float32)

        self.diameter_starts = np.asarray(diameter_starts, dtype = np.float32)
        self.diameter_ends = np.asarray(diameter_ends, dtype = np.float32)
        # outer radius
        self.r_out_starts = self.diameter_starts/2.0
        self.r_out_ends = self.diameter_ends/2.0
        



        self.thickness_starts = np.asarray(thickness_starts, dtype = np.float32)
        self.thickness_ends = np.asarray(thickness_ends, dtype = np.float32)
        # internal radius for the tube cross section
        self.r_int_start = self.r_out_starts - self.thickness_starts
        self.r_int_end = self.r_out_ends - self.thickness_ends
        

        self.densities = np.asarray(densities, dtype = np.float32)
        self.e_values = np.asarray(e_values, dtype = np.float32)
        self.g_values = np.asarray(g_values,    dtype = np.float32)
        self.alpha_s_values = np.asarray(alpha_s_values, dtype = np.float32)
        self.alpha_v_values = np.asarray(alpha_v_values, dtype = np.float32)
        self.scf_starts = np.asarray(scf_starts, dtype = np.float32)
        self.scf_ends = np.asarray(scf_ends, dtype = np.float32)

        
    def get_interp_funcs(self):
        
        self.outer_r_func = self.func2(self.r_out_starts, self.r_out_ends)
        self.inner_r_func = self.func2(self.r_int_start, self.r_int_end)
        self.density_func = self.func(self.densities)
        self.e_func = self.func(self.e_values)
        self.g_func = self.func(self.g_values)
        self.alpha_s_func = self.func(self.alpha_s_values)
        self.alpha_v_func = self.func(self.alpha_v_values)        
        # self.alpha_func = self.func2(self.alpha_s_values,self.alpha_v_values)
        self.scf_func = self.func2(self.scf_starts, self.scf_ends)
    def eval_interp_funcs(self, x):
        return [self.outer_r_func(x), self.inner_r_func(x), self.density_func(x), self.e_func(x), self.g_func(x), self.alpha_s_func(x), self.alpha_v_func(x), self.scf_func(x)]
    
    @staticmethod
    def deriveCrossSectionProperty(elem_prop):
        """
            This is for the calculation of the cross sectional properties
        """
        
        rho = elem_prop[2]	# Mass Density
        E = elem_prop[3]	# Young's Modulus of Elasticity
        G = elem_prop[4]	# Shear Modulus
        alpha_s = elem_prop[-3] if (elem_prop[-3] != None) else 0	# Stress induced dissipation
        alpha_v = elem_prop[-2] if (elem_prop[-2] != None) else 0	# Velocity induced dissipation

        outer_r = elem_prop[0]	# Outer radius
        inner_r =  elem_prop[1]	# Inner radius: outer radius - thickness (th.: start = end)

        # Formulae for the cylinder coordinates of a ring
        area = (np.pi)*(outer_r**2 - inner_r**2)
        # Ixx second moment of interia
        i_1 = (np.pi/4)*(outer_r**4 - inner_r**4)
        # Iyy second moment of interia
        i_2 = i_1
        # centeriod of the ring s_1, and s_2
        s_1 = 0
        s_2 = 0
        # Ixy second moment of interia
        i_12 = 0

        arg11 = E * area
        arg12 = G * area # area1 and area2 are equal
        arg13 = G * area 
        arg14 = E * i_1
        arg15 = E * i_2
        arg16 = G * (i_1 + i_2) # Ip is the polar moment of the inertia
        arg17 = E * s_1
        arg18 = E * s_2
        arg19 = G * s_1
        arg20 = G * s_2
        arg21 = E * i_12
        arg22 = rho * area
        arg23 = rho * i_1
        arg24 = rho * i_2
        arg25 = rho * s_1
        arg26 = rho * s_2
        arg27 = rho * i_12
        arg28 = alpha_s 
        arg29 = alpha_v
        # Property shown as output
        updated_property = [
            [arg11 ,arg12 ,arg13 ,arg14 ,arg15 ,arg16 ,arg17 ,arg18 ,arg19 ,arg20 ,arg21],
            [arg22 ,arg23 ,arg24 ,arg25 ,arg26 ,arg27],
            [arg28, arg29]
        ]

        # Value correction: Convert all negative zeroes to positive zeroes
        for line_ind in range(len(updated_property)):
            for val_ind in range(len(updated_property[line_ind])):
                if float(updated_property[line_ind][val_ind]) == -0.0:
                    updated_property[line_ind][val_ind] = 0.0

        # Return the updated cross section property
        return updated_property
    def __repr__(self) -> str:
        if self.cross_section_properties is None:
            return "No cross section properties loaded"
        else:
            return str(self.cross_section_properties)
    def get_values_mid_elements(self):
        """
            This function returns the values of the cross section properties at the mid-location of the elements
            cross_prop: dict
                it has the following keys
                'num_elements', 'beam_name', 'num_nodes', 1, 2, 'local_node_nums'
                1, 2 is the key for the cross sectional properties for the elements 1 and 2
                'local_node_nums' is the key for the local node numbers of the elements. Morover, this is appended by the global 
                element id. 
        """
        pnts = np.linspace(0, self.L, self.num_elements+1)
        cross_prop = {"num_elements": self.num_elements, "beam_name": self.beam_name, "num_nodes":self.num_elements + 1}
        
        for e, (x1, x2) in enumerate(zip(pnts[:-1], pnts[1:])):
            if self.L == 0:
                continue
            # e is local element number
            x = (x1+x2)/2.0
            
            values = self.eval_interp_funcs(x)
            cross_prop[e+1] = self.deriveCrossSectionProperty(values)
        # giving global element ids
        pnts_num = range(1,cross_prop["num_nodes"]+1)
        local_node_nums = []
        for pt1 , pt2 in zip(pnts_num[:-1],pnts_num[1:]):
            local_node_nums.append([pt1,pt2, self.global_id])
            self.global_id += 1
        cross_prop["local_node_nums"] = np.asarray(local_node_nums, dtype = np.int32)

        self.cross_section_properties = cross_prop
            

    def func(self, var):
        gen = []
        for v1,v2 in zip(var, var):
            gen.append(v1)
            gen.append(v2)
        values = np.asarray(gen)
        return interp1d(self.x, values)
    def func2(self, var1, var2):
        gen = []
        for v1,v2 in zip(var1, var2):
            gen.append(v1)
            gen.append(v2)
        values = np.asarray(gen)
        return interp1d(self.x, values)
    
class Jacket:
    def __init__(self, num_legs, num_bays,L, losg, r_base,r_top, l_tp=0.0, file_name = 'beaminput.txt', axes=None):
        self.num_legs = num_legs
        self.num_bays = num_bays
        # total height of the jacket
        self.L = L
        # height from ground to the first bay
        self.losg = losg
        # height from top bay to the total height
        self.l_tp = l_tp
        self.r_base = r_base
        self.r_top = r_top

        self.oblique_angle_plane = atand(self.L/(self.r_base-self.r_top))
        self.angle_pos = 360/self.num_legs
        self.r = lambda h: self.r_base + (self.r_top-self.r_base)/(self.L - 0 ) * (h - 0)
        # total length of jacket stand
        self.len_jacket_stand = np.sqrt((self.r_base-self.r_top)**2 + self.L**2)
        
        self.beams_info = None
        self.segment_data = None
        self.bays_height_data = None # list of heights of the bays
        self.lengths_bays = {}
        if axes is None:
            self.ax = plt.figure().add_subplot(projection='3d')
        else:
            self.ax = axes
        self.file_name = path_main/Path(f"io/{file_name}")
        

    def calc_abs_height(self):
        """
            This function calculates the absolute height of each stand. 
            By absolute height we mean the height of the stand from the ground.
            position_stands is a dictionary. Its first key is the name of the beam such
            as "SB", "S1", "S2", etc. 
            Then it has another dictionary which is correspoinding to each point in the leg
            For instance if the number of elements for SB is defined to be 4 ,and there are three
            legs of the jacket structure, it will pop the following keys:
            >>> position_stands["SB"].keys()  
            dict_keys(['num_pnts', (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (4, 2)])
            
            num_pnts is correspoing to the number of points defined for this beam
            (0, 0) is the first point of the first leg
            (1, 0) is the second point of the first leg
            In the tuple the first element is correspoind to the number of points and the second is for the 
            leg number. For instance here if you want to extract the points for 
            the second point of the third leg
            >>> position_stands["SB"][(1, 2)]
            {'x': -1.4375000000000013, 'y': -2.4898230358802604, 'z': 0.25}
            
            By this key you will get the x, y, z coordinates of the point.
        """
        position_stands = {}
        abs_height_radius = {}
        h = 0.0
        for s in range(self.heights_stands.shape[0]):
            current_height = self.heights_stands[s]
            current_beam_name = self.stand_class_names[s]
            num_elements = self.beams_info[current_beam_name]['no_elements']
            num_pnts = num_elements + 1
            # vector of heights
            local_heights = np.linspace(0,current_height, num_pnts)
            abs_height = local_heights + h
            h += current_height
            # abs_radius is vector of the same length as abs_height
            abs_radius = self.r(abs_height)
            abs_height_radius[current_beam_name] = {"h": abs_height, "r":abs_radius,"len":num_pnts}
            # it containts the positions of all points in the beam
            position_stands[current_beam_name] = self.convert_hr_position( abs_height_radius[current_beam_name] , self.num_legs)
        self.abs_height_radius = abs_height_radius
        # for each class we have the position of the points
        self.stands_pnts = position_stands
        
    def doubleNum_elements(self, dict):
        """
            This function doubles the number of elements beams which are defined at the bays location.
            The idea is that instead of defining 4 beams for one X structure at the bay. We define two
            beams but with doubled number of elements.
        """
        new_dict = {}
        for beam_name, values in dict.items():
            if "C" in beam_name:
                dict[f"{beam_name}"] = values
                num_elements = values["no_elements"]
                num_segments = values["no_segments"]
                new_dict[f"{beam_name}"] = {'no_segments': num_segments, 'no_elements': num_elements*2}
            else:
                new_dict[beam_name] = values
        
        return new_dict    
    
    def setBeamInfo(self,beam_info, beam_class_names, num_beam_classes):
        """
            This is the double dictionary data. First dictionary contains the names of the beams. The beam names
            can be such as SB which is for the bottom stand. The sequence of names is SB, S1, S2, ..., SN, ST, 
            C1L, C1U, C2L, C2U, ... CNU.
            Corresponding to each name we have the second dictionary which contains the information for the number 
            of segments and number of elements.
            For 2 number of comartments, the data would be:
            {'SB': {'no_segments': 2, 'no_elements': 2}, 
            'S1': {'no_segments': 2, 'no_elements': 2}, 
            'S2': {'no_segments': 2, 'no_elements': 2}, 
            'ST': {'no_segments': 2, 'no_elements': 2}, 
            'C1L': {'no_segments': 2, 'no_elements': 2}, 
            'C1U': {'no_segments': 2, 'no_elements': 2}, 
            'C2L': {'no_segments': 2, 'no_elements': 2}, 
            'C2U': {'no_segments': 2, 'no_elements': 2}}

            SB, S1, S2 is defined for the individual leg. For the other legs the data would be assumed to be same.
            C1L stands for the lower compartment 1. By compartment we mean the bay.    
        """
        self.beams_info = self.doubleNum_elements(beam_info)
        self.num_beam_classes = num_beam_classes
        # names of all beams classes
        self.beam_class_names = beam_class_names
        # bays names
        self.bays_names  = [beam_name for beam_name in self.beam_class_names if "C" in beam_name]
        # names of stands
        self.stand_class_names = [stand_name for stand_name in self.beam_class_names if "S" in stand_name]
        # stands names which are appearing with bays
        self.bay_stands_names = [stand_name for stand_name in self.stand_class_names if not "T" in stand_name and not "B"  in stand_name]
    def setBeamSegment(self,segment_data):
        """
            segment_data is a dictionary which containts the information of segment data. The segment data is 
            defined by the Segment object which contains properties such as E, diameter_start, diameter_end,....
            The segment data for 2 bays and 2 no segments is defined as follows:
            {'SB': [Segment ID: 0, Segment ID: 1],
             'S1': [Segment ID: 2, Segment ID: 3],
             'S2': [Segment ID: 4, Segment ID: 5],
             'ST': [Segment ID: 6, Segment ID: 7],
             'C1L': [Segment ID: 8, Segment ID: 9], 'C1U': [Segment ID: 10, Segment ID: 11],
             'C2L': [Segment ID: 12, Segment ID: 13], 'C2U': [Segment ID: 14, Segment ID: 15]}
            Corresponding to each beam name we have the list of segments
        """
        self.segment_data = segment_data
    def setBaysHeight(self,bays_height_data):
        """
            bays_height_data is a list which contains the height of each bay
        """
        self.bays_height_data = bays_height_data
       

    def calc_heights(self):
        # this will be called in the initialization
        """
            This function calculates the heights of each stand.
        """
        # height data will come from the height dialog
        # self.compartments_height_data = self.dlgCompartment.load() will output the list of heights for each compartment
        heights_stands = np.zeros( (self.num_bays +2,))
        heights_stands[0] = self.losg
        # height_comp_stand = (self.L - self.losg - self.l_tp)/(self.num_bays)
        h = self.losg
        for i in range(self.num_bays):
            heights_stands[i+1] =  self.bays_height_data[i]
            h += heights_stands[i+1]
        heights_stands[-1] = self.L - h
        self.l_tp = heights_stands[-1]
        # TODO: Implement this in try block so that program does not crash if the user does not enter the correct data
        if heights_stands[-1] < 0:
            raise ValueError("The sum of heights of the jackets is greater than the heights of all the jacket")
        self.heights_stands = heights_stands
    def calc_lengths(self):
        """
            This function calculates the lengths of each stand. Height can be considered
            as the vertical distance from the ground, and length is the slanted dsistance.
            Each height and the length is of the member of the stand.
        """
        lenghts_stands = [h/self.L*self.len_jacket_stand for h in self.heights_stands]    
        self.lenghts_stands = np.asarray(lenghts_stands)
        
        for leg in range(self.num_legs):
            for i,l in enumerate(self.lenghts_stands):
                self.lengths_bays[(self.stand_class_names[i], leg)] = l
            
    def convert_hr_position(self, _dic:dict, number_legs:int):
        """
            This function converts the absolute height and radius to the position of the stand
            to the position of the stand in the leg.
            returns:
                    pnts: a dictionary which contains the position of the stand in the leg.
                    For a point 1 and the leg 1, we will have another dictionary which will containt the points
                    x, y, z
                     {'(0,0)': {'x': 0.0, 'y': 0.0, 'z': 0.0}}
        """
        heights = _dic["h"]
        radius = _dic["r"]
        num_pnts = _dic["len"]
        angle_pos = 360/number_legs
        pnts = {"num_pnts":num_pnts}
        for k in range(num_pnts):
            h = heights[k]
            r = radius[k]
            for l in range(number_legs):
                x = r*cosd(l*angle_pos)
                y = r*sind(l*angle_pos)
                z = h
                # we should give names as per the leg and also as per the stand
                pnts[(k,l)] = {'x':x, 'y':y, 'z':z}
        return pnts

    def line3D(self,lower_pnt,mag,dir):
        """
            This function returns the points of a line in 3D space.
            lower_pnt: the lower point of the line
            mag: the magnitude of the line
            dir: the direction of the line
        """
        pnts = np.zeros((len(mag),3))
        for pnt_num in range(len(mag)):
            pnts[pnt_num,:] = lower_pnt + mag[pnt_num] * dir
        return pnts
    def  calc_bay_pnts(self,position_stands,number_legs):
        """
            This function calculates the points of the bay.
            position_stands: a dictionary which contains the position of the stands in the leg.
            bay_pnts: a dictionary which contains the points of the bay.
                For the leg combinition leg1 and leg2 we have the bay and the X type of the structure.
                The first key for bays_pnts is leg number, and this leg number is corresponding to the combinition 
                of the legs.
                The second key is the name of the bay which C1L or C1U which will give all points for the C1U 
        """
        wrapN = lambda i: ( 1 + (i-1) % number_legs)
        #  bays will be between leg1, leg2 ; leg2, leg1; leg3, leg1
        # bays_pnts = cell(number_legs,1);
        
        bays_pnts = {}
        for leg in range(number_legs):
            bays_pnts[leg] = {} #cell(number_bays,2);
            l1 = leg
            l2 = leg +1
            if l2 == number_legs:
                l2 = 0
            bay_idx = 1
            for bay in self.bay_stands_names:
                temp = position_stands[bay]
                lower_pnt1 = temp[(0,l1)]["x"] * np.array([1,0,0]) + temp[(0,l1)]["y"] * np.array([0,1,0]) + temp[(0,l1)]["z"] * np.array([0,0,1])
                lower_pnt2 = temp[(0,l2)]["x"] * np.array([1,0,0]) + temp[(0,l2)]["y"] * np.array([0,1,0]) + temp[(0,l2)]["z"] * np.array([0,0,1])
                idx = temp["num_pnts"] - 1
                upper_pnt1 = temp[(idx,l1)]["x"] * np.array([1,0,0]) + temp[(idx,l1)]["y"] * np.array([0,1,0]) + temp[(idx,l1)]["z"] * np.array([0,0,1])
                upper_pnt2 = temp[(idx,l2)]["x"] * np.array([1,0,0]) + temp[(idx,l2)]["y"] * np.array([0,1,0]) + temp[(idx,l2)]["z"] * np.array([0,0,1])

                
                comp_1_name_1 = f"C{bay_idx}L"
                comp_2_name_1 = f"C{bay_idx}U"
                # for finding the direction
                comp_1 = upper_pnt2 - lower_pnt1
                comp_2 = upper_pnt1 - lower_pnt2
                len_1 = np.linalg.norm(comp_1,2)
                len_2 = np.linalg.norm(comp_2,2)

                self.lengths_bays[(comp_1_name_1,leg)] = len_1
                self.lengths_bays[(comp_2_name_1,leg)] = len_2
                
                num_elements_comp_1 = self.beams_info[comp_1_name_1]["no_elements"]
                num_elements_comp_2 = self.beams_info[comp_2_name_1]["no_elements"]

                ratio_1 = np.linspace(0,len_1,num_elements_comp_1+1);
                ratio_2 = np.linspace(0,len_2,num_elements_comp_2+1);
                
                dir_1 = comp_1 /len_1
                dir_2 = comp_2 /len_2

                pnts_1 = self.line3D(lower_pnt1,ratio_1, dir_1)
                pnts_2 = self.line3D(lower_pnt2,ratio_2, dir_2)
                bays_pnts[leg][comp_1_name_1] = pnts_1
                bays_pnts[leg][comp_2_name_1] = pnts_2
                bay_idx = bay_idx + 1
        
        return bays_pnts
    def calc_cross_sectional_properties(self):
        cross_sectional_properties = {}
        global_id = 1
        total_num_cross_sectional_properties = 0
        for (beam_name,leg), L in self.lengths_bays.items():
            num_elements = self.beams_info[beam_name]["no_elements"]
            total_num_cross_sectional_properties += num_elements
            seg_list = self.segment_data[beam_name]
            temp = interpCrossSection(segments_list=seg_list, L=L ,num_elements=num_elements, beam_name =beam_name, global_id=global_id)
            global_id = temp.global_id
            cross_sectional_properties[(beam_name,leg)] = temp.cross_section_properties
        total_num_beams = len(cross_sectional_properties.keys())

        self.total_num_beams = total_num_beams
        self.total_num_cross_sectional_properties = total_num_cross_sectional_properties
        
        return cross_sectional_properties

    @staticmethod
    def string(array,if_int= False):
        def checkThreshold(X):
            threshold = 1e-11
            X[abs(X)<threshold] = 0.0
            return X
        def convert(X):
            data = ""
            # d0 is added because of the precison
            for i in range(X.shape[0]):
                elements = X[i,:]
                if if_int:
                    data += ' '.join(map(str, elements))+ "\n"
                else:    
                    data += 'd0 '.join(map(str, elements))+ "d0\n"
            return data
        array = checkThreshold(array)
        return convert(array)
    def listVal2String(self,list_val):
        threshold = 1e-11
        for i,elements in enumerate(list_val):
            for j,c in enumerate(elements):
                if abs(c) < threshold:
                    list_val[i][j] = 0.0
        data = ""
        # d0 is added because of the precison
        for i,elements in enumerate(list_val):
            data += 'd0 '.join(map(str, elements))+ "d0\n"
        return data
    def constraints_interactions(self,origin_seq, destiniation_seq,if_only_mid_connections=True):
        threshold = 1e-11
        if if_only_mid_connections:
            mid_pnt1 = self.beams_global_nodes_info[("C1L",0)]["pnts"][int(self.beams_info["C1L"]["no_elements"]/2)]
            mid_pnt2 = self.beams_global_nodes_info[("C1U",0)]["pnts"][int(self.beams_info["C1U"]["no_elements"]/2)]
            threshold = np.linalg.norm(mid_pnt1-mid_pnt2,2) +.01
        rigid_connections = []
        rigid_support_constraints = []
        explored_set = set()
        for beam_name_origin in origin_seq:
            for leg_origin in range(self.num_legs):
                # print(beam_name_origin, ", leg=",leg_origin)
                name_origin = f"{beam_name_origin} ,leg={leg_origin}"
                info_origin = self.beams_global_nodes_info[(beam_name_origin,leg_origin)]
                pts_origin = info_origin["pnts"]
                # print(f"total_ne = ",self.beams_info[beam_name_origin]["no_elements"])
                mid_index_origin = int(self.beams_info[beam_name_origin]["no_elements"]/2)
                nodes_origin = info_origin["global_nodes_id"]
                middle_node_origin = nodes_origin[mid_index_origin]
                if beam_name_origin == "SB":
                    rigid_support_constraints.append(nodes_origin[0])

                starting_pnt_origin = pts_origin[0,:]
                ending_pnt_origin = pts_origin[-1,:]
                mid_pnt_origin = (starting_pnt_origin + ending_pnt_origin)/2
                # print(f"mid_pnt_origin = {mid_pnt_origin}, middle_pnt = {pts_origin[mid_index_origin,:]}")
                for beam_name_dest in destiniation_seq:
                    for leg_dest in range(self.num_legs):
                        name_dest = f"{beam_name_dest} ,leg={leg_dest}"
                        info_dest = self.beams_global_nodes_info[(beam_name_dest,leg_dest)]
                        if name_dest == name_origin:
                            continue
                        if {name_origin,name_dest} in explored_set:
                            continue
                        pts_dest = info_dest["pnts"]
                        nodes_dest = info_dest["global_nodes_id"]
                        mid_index_dest = int(self.beams_info[beam_name_dest]["no_elements"]/2)
                        middle_node_dest = nodes_dest[mid_index_dest]
                        starting_pnt_dest = pts_dest[0,:]
                        ending_pnt_dest = pts_dest[-1,:]
                        mid_pnt_dest = (starting_pnt_dest + ending_pnt_dest)/2
                        flag = False
                        
                        if if_only_mid_connections:
                            # print(f"{name_origin} and {name_dest}: dist = {np.linalg.norm(mid_pnt_origin - mid_pnt_dest,2)}")
                            if np.linalg.norm(mid_pnt_origin - mid_pnt_dest,2) <= threshold :
                                # print(f"{name_origin} and {name_dest} are connected by their mid points")
                                rigid_connections.append([middle_node_origin, middle_node_dest])
                                flag = True

                        else:
                            if np.linalg.norm(ending_pnt_origin - ending_pnt_dest,1) < threshold :
                                print(f"{name_origin} and {name_dest} are connected by their ending points")
                                rigid_connections.append([nodes_origin[-1], nodes_dest[-1]])
                                flag = True
                            
                            if np.linalg.norm(ending_pnt_origin - starting_pnt_dest,1) < threshold :
                                print(f"{name_origin} and {name_dest} are connected. ")
                                rigid_connections.append([nodes_origin[-1], nodes_dest[0]])
                                flag = True
                        if flag:    
                            explored_set.add(frozenset([name_dest, name_origin]))
        return rigid_connections, rigid_support_constraints
    def find_constraints(self):
        threshold = 1e-11
        rigid_connections = []
        rigid_support_constraints = []
        explored_set = set()
        total_num_nodes = self.beams_global_nodes_info[(f"C{self.num_bays}U", self.num_legs-1)]["global_nodes_id"][-1]
        print(f"total number of nodes = {total_num_nodes}")
        seq_beam_origin = self.stand_class_names
        # for (beam_name_origin,leg_origin), info_origin in self.beams_global_nodes_info.items():
        for beam_name_origin in seq_beam_origin:
            for leg_origin in range(self.num_legs):
                # print(beam_name_origin, ", leg=",leg_origin)
                name_origin = f"{beam_name_origin}_{leg_origin}"
                info_origin = self.beams_global_nodes_info[(beam_name_origin,leg_origin)]
                pts_origin = info_origin["pnts"]
                nodes_origin = info_origin["global_nodes_id"]
                if beam_name_origin == "SB":
                    rigid_support_constraints.append(nodes_origin[0])

                starting_pnt_origin = pts_origin[0,:]
                ending_pnt_origin = pts_origin[-1,:]
                for (beam_name_dest,leg_dest), info_dest in self.beams_global_nodes_info.items():
                    name_dest = f"{beam_name_dest}_{leg_dest}"
                    if name_dest == name_origin:
                        continue
                    if {name_origin,name_dest} in explored_set:
                        continue
                    pts_dest = info_dest["pnts"]
                    nodes_dest = info_dest["global_nodes_id"]
                    starting_pnt_dest = pts_dest[0,:]
                    ending_pnt_dest = pts_dest[-1,:]
                    flag = False
                    def nearest(n1,n2):
                        pn1 = self.beams_global_nodes_pos[n1]
                        pn2 = self.beams_global_nodes_pos[n2]
                        dist1 = sum(pn1*pn1)
                        dist2 = sum(pn2*pn2)
                        if dist1 < dist2:
                            return [n1,n2]
                        elif dist1 > dist2:
                            return [n2,n1]
                        else:
                            return [n1,n2]
                    # if np.linalg.norm(starting_pnt_origin - starting_pnt_dest,1) < threshold : 
                    #     # print(f"{name} and {name2} are connected by their starting points")
                    #     rigid_connections.append(nearest(nodes_origin[0], nodes_dest[0]))
                    #     flag = True
                    if np.linalg.norm(ending_pnt_origin - ending_pnt_dest,1) < threshold :
                        # print(f"{name_origin} and {name_dest} are connected by their ending points")
                        rigid_connections.append([nodes_origin[-1], nodes_dest[-1]])
                        flag = True
                    # if np.linalg.norm(starting_pnt_origin - ending_pnt_dest,1) < threshold :
                    #     # print(f"{name} and {name2} are connected. The starting point of {name} is connected to the ending point of {name2}")
                    #     rigid_connections.append(nearest(nodes_origin[0], nodes_dest[-1]))
                    #     flag = True
                    if np.linalg.norm(ending_pnt_origin - starting_pnt_dest,1) < threshold :
                        # print(f"{name_origin} and {name_dest} are connected. ")
                        rigid_connections.append([nodes_origin[-1], nodes_dest[0]])
                        flag = True
                    if flag:
                        
                        explored_set.add(frozenset([name_dest, name_origin]))
                        # explored1.append(name2)
                    
        origin_seq = [name for name in self.bays_names if "L" in name]
        dest_seq = [name for name in self.bays_names if "U" in name]
        rigid_connections_2,_ = self.constraints_interactions(origin_seq, dest_seq)
        rigid_connections.extend(rigid_connections_2)
        rigid_connections = np.asarray(rigid_connections, dtype=int)
        
        # rigid_connections_2 = np.asarray(rigid_connections_2, dtype=int)
        # print(f"rigid_connections_2 = {rigid_connections_2.shape[0]}")
        # Graph Visualization
        # G = nx.Graph()
        # for r in rigid_connections:
        #     G.add_edge(r[0], r[1])
        # plt.figure()
        # nx.draw(G, with_labels=True)
        total_num_rigid_connections = rigid_connections.shape[0]
        total_num_rigid_supports = len(rigid_support_constraints)
        nodes = set(range(1,total_num_nodes+1))
        non_internal_points_rc = set()
        # internal nodes for the rigid connection
        internal_nodes_rc= set(rigid_connections[:,0])
        non_internal_points_rc = set(rigid_connections[:,1])
        flag = True
        for ii in internal_nodes_rc:
            if ii in non_internal_points_rc:
                flag = False
                print("P internal nodes: ",ii)
        if flag:
            print("sucess! in generating origin nodes and destination nodes for rigid connections")
        else:
            print("error! in generating origin nodes and destination nodes for rigid connections")
        # print(f"len(internal_nodes_rc) = {len(internal_nodes_rc)}, len(non_internal_points_rc) = {len(non_internal_points_rc)}")
        internal_nodes_remaining = internal_nodes_rc - non_internal_points_rc
        non_internal_nodes_remaining = non_internal_points_rc - internal_nodes_rc
        
        # print(f"len(non_internal_points_rc): {len(non_internal_points_rc)}, len(internal_nodes_remaining): {len(internal_nodes_remaining)}")
        # print(f"len(non_internal_nodes_remaining): {len(non_internal_nodes_remaining)}")
        
        
        # print("len(explored): ",len(explored_set))
        # print("total rigid connections:", len(rigid_connections))
        # print("total rigid support constraints:", len(rigid_support_constraints))
        internal_nodes = nodes - non_internal_nodes_remaining
        internal_nodes = internal_nodes - set(rigid_support_constraints)
        total_num_internal_nodes = len(internal_nodes)
        total_num_constraints = total_num_internal_nodes + total_num_rigid_connections + total_num_rigid_supports
        
        # print(f"total number of constraints = {total_num_constraints}", f"total internal constraints: {total_num_internal_nodes}", f"total num of rigid connections: {total_num_rigid_connections}", f"total num of rigid supports: {total_num_rigid_supports}")
        
        
        cn = Constraints(num_constraint= total_num_constraints,
               num_rigid_supports=total_num_rigid_supports,
               num_revolute_joints =0, 
               num_rigid_connection=total_num_rigid_connections)
        cn.setInternalNodes(internal_nodes)
        for n in rigid_support_constraints:
            cn.addRigidSupportElement(n,0)
        for n1,n2 in rigid_connections:
            cn.addRigidConnectionElement(n1,n2)
        cn.write(Path('io/constraint12input.txt'))
        
    
    def write_data(self):
        data = "!! beam input for NREL 15 MW reference wind turbine (based on HAWC2 input files)\n" +"!!\n" +"!! number of beams (1), number of cross-section properties (2)\n"
        data += f"{self.total_num_beams} {self.total_num_cross_sectional_properties}\n" 
        # connection information. This will store the global ids of the nodes
        beams_global_nodes_info = {}
        # global nodes position
        beams_global_nodes_pos = {}
        current_node = int(1)
        for (beam_name,leg), cross_sectional_properties in self.cross_sectional_properties.items():
            num_nodes = cross_sectional_properties["num_nodes"]
            beams_global_nodes_info[(beam_name,leg)] = {"global_nodes_id":list(range(current_node,current_node+num_nodes))}
            current_node += num_nodes
            num_elements = cross_sectional_properties["num_elements"]
            local_node_nums = cross_sectional_properties["local_node_nums"]
            data += f'!!\n'
            data += f'!!\n'
            data += f"!! jacket beam ({beam_name}) , leg = {leg}: number of nodes (1), number of elements (2)\n"
            data += f"{num_nodes} {num_elements}\n"
            data += f'!!\n'
            data += f'!!\n'
            data += f'!! jacket beam  ({beam_name}) , leg = {leg}: nodal phi (1, 2, 3), nodal d1 (4, 5, 6), nodal d2 (7, 8, 9), nodal d3 (10, 11, 12)\n'
            pts = self.directors_data[(beam_name,leg)]
            beams_global_nodes_info[(beam_name,leg)]["pnts"] = pts[:,:3]
            for node, pnt in zip(beams_global_nodes_info[(beam_name,leg)]["global_nodes_id"],beams_global_nodes_info[(beam_name,leg)]["pnts"]):
                beams_global_nodes_pos[node] = pnt
            data += Jacket.string(pts)
            data += f'!!\n'
            data += f'!!\n'
            data += f'!! jacket beam  ({beam_name}) , leg = {leg}: connectivities (1, 2), cross-section property (3)\n'
            
            data += Jacket.string(local_node_nums, if_int=True)
        self.beams_global_nodes_info = beams_global_nodes_info
        self.beams_global_nodes_pos = beams_global_nodes_pos
        self.find_constraints()
        
        
        for (beam_name,leg), cross_sectional_properties in self.cross_sectional_properties.items():
            local_node_nums = cross_sectional_properties["local_node_nums"]
            num_elements = cross_sectional_properties["num_elements"]
            for e in range(1, num_elements+1):    
                data += f'!!\n'
                data += f'!!\n'
                data += f'!! jacket beam  ({beam_name}) , leg = {leg}: cross-section properties {local_node_nums[e-1,2]}\n'
                data += self.listVal2String(cross_sectional_properties[e])
                # data += f'!!\n'
                # data += f'!!\n'
        with open(self.file_name, 'w') as f:
            f.write(data)
            print(f"Data written to {self.file_name}")
        return data



    def cal_dirs(self, beam_pnts, if_plot=False):
        """
            This calculates the directors for each point in the beam
            result: will have (num_points, 12) shape
            where first 3 will be the position of the point
            and the other 9 will be the director of the point
        """
        num_pnts = beam_pnts.shape[0]
        pnts_num_list = range(num_pnts)
        result = np.zeros((num_pnts,12))
        for i, (pt1, pt2) in enumerate(zip(pnts_num_list[:-1], pnts_num_list[1:])):
            dir = beam_pnts[pt2,:] - beam_pnts[pt1,:]
            dir = dir / np.linalg.norm(dir,2)
            x = dir
            y = np.asarray([-x[1], x[0], 0.0], dtype=np.float32)
            z = np.cross(x, y)
            result[i,:] = np.concatenate((beam_pnts[pt1,:],x,y,z))
        result[-1,3:] = result[-2,3:]
        result[-1,:3] = beam_pnts[-1,:]
        if if_plot:
            plt = self.ax
            plt.quiver(result[:,0], result[:,1], result[:,2], result[:,3], result[:,4], result[:,5], color="r")
            plt.quiver(result[:,0], result[:,1], result[:,2], result[:,6], result[:,7], result[:,8], color="g")
            plt.quiver(result[:,0], result[:,1], result[:,2], result[:,9], result[:,10], result[:,11], color="b")
        return result
    def calc_directors(self, if_plot=False):
        """
            This calculates the directors for each point in the beam
            result: will have (num_points, 12) shape
            where first 3 will be the position of the point
            and the other 9 will be the director of the point
        """
        result = {}
        for beam_name, attr in self.stands_pnts.items():
            num_pnts = attr["num_pnts"]
            for l in range(self.num_legs):
                x = []
                y = []
                z = []
                for i in range(num_pnts):
                    x.append(attr[(i,l)]["x"])
                    y.append(attr[(i,l)]["y"])
                    z.append(attr[(i,l)]["z"])
                    
                x = np.asarray(x).reshape((-1,1))
                y = np.asarray(y).reshape((-1,1))
                z = np.asarray(z).reshape((-1,1))
                pnts = np.concatenate((x,y,z), axis=1)
                
                dirs = self.cal_dirs(pnts, if_plot)
                result[(beam_name,l)] = dirs
                
        for leg, beam in self.bay_pnts.items():
            for beam_name, pnts in beam.items(): 
                dirs = self.cal_dirs(pnts, if_plot)
                result[(beam_name,leg)] = dirs
        
        return result
        
    def plotJacket(self):
        linewidth = 2.5
        for beam_name, attr in self.stands_pnts.items():
            num_pnts = attr["num_pnts"]
            for leg in range(self.num_legs):
                x = []
                y = []
                z = []
                for i in range(num_pnts):
                    x.append(attr[(i,leg)]["x"])
                    y.append(attr[(i,leg)]["y"])
                    z.append(attr[(i,leg)]["z"])
                self.ax.plot(x,y,z,color="k",linewidth=linewidth)
                self.ax.text(np.mean(x),np.mean(y),np.mean(z),f"{beam_name} {leg}",color="k",fontsize=10)
                
        for leg, beam in self.bay_pnts.items():
            for beam_name, pnts in beam.items():
                
                x = pnts[:,0]
                y = pnts[:,1]
                z = pnts[:,2]

                eps = 0.1
                if "L" in beam_name:
                    self.ax.text(np.mean(x)+eps,np.mean(y)+eps,np.mean(z)+eps,f"{beam_name} {leg}",color="r",fontsize=10)
                    self.ax.plot(x,y,z,color="r",linewidth=linewidth)
                else:
                    self.ax.text(np.mean(x),np.mean(y),np.mean(z),f"{beam_name} {leg}",color="g",fontsize=10)
                    self.ax.plot(x,y,z,color="g",linewidth=linewidth)
        

    def calc_pnts(self):
        """
            This function calculates the points of the bay.
             
        """
        self.calc_heights()
        self.calc_abs_height()
        self.calc_lengths()
        self.bay_pnts = self.calc_bay_pnts(self.stands_pnts,self.num_legs)
        self.cross_sectional_properties = self.calc_cross_sectional_properties()
        self.plotJacket()
        self.directors_data = self.calc_directors(if_plot=False)
        self.write_data()
        # TODO: remove this plt.show
        # plt.show()
        plt.draw()
        return true
        


if __name__ == "__main__":
    jack = Jacket(num_legs = 3, 
                num_bays = 2,
                L = 4, 
                losg = 1, 
                l_tp = 1, 
                r_base = 3,
                r_top = 2)
    jack.setBaysHeight( [1.0,1.0])
    beam_info = {'SB': {'no_segments': 2, 'no_elements': 2}, 
            'S1': {'no_segments': 2, 'no_elements': 2}, 
            'S2': {'no_segments': 2, 'no_elements': 2}, 
            'ST': {'no_segments': 2, 'no_elements': 2}, 
            'C1L': {'no_segments': 2, 'no_elements': 2}, 
            'C1U': {'no_segments': 2, 'no_elements': 2}, 
            'C2L': {'no_segments': 2, 'no_elements': 2}, 
            'C2U': {'no_segments': 2, 'no_elements': 2}}
    segment_data = {'SB': [Segment( segment_id= 0), Segment( segment_id = 1)],
            'S1': [Segment( segment_id= 2), Segment( segment_id = 3)],
            'S2': [Segment( segment_id= 4), Segment( segment_id = 5)],
            'ST': [Segment( segment_id= 6), Segment( segment_id = 7)],
            'C1L': [Segment( segment_id= 8), Segment( segment_id = 9)],
            'C1U': [Segment( segment_id= 10), Segment( segment_id = 11)],
            'C2L': [Segment( segment_id= 12), Segment( segment_id = 13)],
            'C2U': [Segment( segment_id= 14), Segment( segment_id = 15)]}
    beam_class_name = ['SB','S1','S2','ST','C1L','C1U','C2L','C2U']
    num_beam_classes = 8
    jack.setBeamInfo(beam_info,beam_class_name,num_beam_classes)
    jack.setBeamSegment(segment_data)
    jack.calc_pnts()

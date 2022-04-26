import sys
from pathlib import Path
import os
from turtle import color
import numpy as np
import matplotlib.pyplot as plt

file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir) / Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from Utils.geometry import Geometry
from structure.tower.tower import Tower
from beam.beam import Beam
from beam.segment import Segment
from Utils.segments_userInterface import segments_ui, main_beamInfo

from PyQt5.QtWidgets import QDialog, QWidget
from segment_table import *
import traceback


class SegmentsByBeams(QDialog):
    def __init__(self, beam_info, beam_heights):
        super().__init__()
        # self.segments is a dictionary: segments[class_name][segment_number-1] where segment_number starts from 1
        # elements of segments is a Segment object
        self.segments = {}
        self.beams_data = beam_info
        self.beam_heights = beam_heights
        self.beams_classes_names = self.get_beam_classes_names()
        self.get_class_data()
        # classes of segment_ui which will keep the information of buttons
        self.segment_btns = {}
        self.vertical_pose = 0
        self.segments_valid_outputs = {}
        self.segments_valid_id = []
        self.all_segments_valid = False
    
    def get_beam_classes_names(self):
        return list(self.beams_data.keys())
    def get_class_data(self):
        self.class_number_elements = []
        self.class_number_segments = []
        
        for _,class_setting in self.beams_data.items():
            self.class_number_elements.append(class_setting['number_of_elements'])
            self.class_number_segments.append(class_setting['number_of_segments'])
    def compartmentData2ID(self):
        # it contains the a unique id int which is corresponding to the classname and segment number
        # classname is the key and segment number is the value
        self.cd2ID = {}
        # ID to classname and segment number
        self.ID2cd = {}
        # corresponding classname we have the global number of segment and the local number of the segments both are lists
        self.cd2IDlist = {}
        i = 0
        for cls_name in self.beams_classes_names:
            global_segment_id_list = []
            segments_local_list = []
            for seg_id in range(self.class_number_segments[i]):
                self.cd2ID[cls_name,seg_id+1] = i*self.class_number_segments[i]+seg_id
                self.ID2cd[i*self.class_number_segments[i]+seg_id] = cls_name,seg_id+1
                global_segment_id_list.append(i*self.class_number_segments[i]+seg_id)
                segments_local_list.append(seg_id + 1)
                # print(f"{cls_name} {seg_id+1} {i*self.class_number_segments[i]+seg_id}")
            self.cd2IDlist[cls_name] = global_segment_id_list,segments_local_list
            i+=1
        

            
    # def __getitem__(self, key):
    #     return self.segments[key]

    # def __len__(self):
    #     return len(self.segments)

    # def __iter__(self):
    #     return iter(self.segments)

    # def __repr__(self):
    #     return 'SegmentsByBeams(segments={}, beams={})'.format(
    #         self.segments, self.beams)

    # def __str__(self):
    #     return 'SegmentsByBeams(segments={}, beams={})'.format(
    #         self.segments, self.beams)

    # def __eq__(self, other):
    #     return (self.segments == other.segments and
    #             self.beams == other.beams)

    # def __ne__(self, other):
    #     return not self == other

    
    
    
    # def checkAllSegments(self):
    #     for id,seg in enumerate(self.segments):
    #         if not self.checkDialogInput(seg):
    #             util.showErrorCustomMsg(f"Segment {id+1} field values were not correct","Please correct them as per suggestions. Output file is not created yet")
    #             return False
    #     return True
    
    
      

    def dispBtnBox(self,parent,no_of_segments,class_name):
        # display botton box
        # taking values from segment table
        self.segments[class_name] = []
        self.cur_class_name = class_name
        self.segment_btns[class_name] = segments_ui(parent, no_of_segments,btnGrpIds=self.cd2IDlist[class_name])
        self.segment_btns[class_name].setupUI_layoutStyle(pos_y = self.vertical_pose)
        # self.segment_btns[class_name].ui_prop()
        self.vertical_pose += 80+10 #self.segment_btns[class_name].height()
        self.segment_btns[class_name].label.setText(f"Enter the segment data for {class_name}")
        print(f"Enter the segment data for {class_name}")
        

        for id in range(no_of_segments):
            self.segments[class_name].append(Segment(segment_id=id))
            self.cur_segment_id = id
            print(self.segment_btns[class_name])
            
            # self.segment_btns[class_name].btns[id].clicked.connect(self.dispDialogSegmentTable)
            
            # self.segment_btns[class_name].btns[id].clicked.connect(self.dispCheck)
            self.segment_btns[class_name].btn_grp.buttonClicked.connect(self.dispCheck)
            self.segment_btns[class_name].btns[id].show()
            self.segment_btns[class_name].btns[id].setEnabled(True)
    
    def dispCheck(self):
        print("hello")
        print(f"Checking {self.cur_class_name} {self.cur_segment_id}")
        

    def setup_ui(self):
        self.compartmentData2ID()
        if self.segments_valid_id:
            self.segments, status_accept, self.segments_valid_id = main_beamInfo(self.beams_data,self.cd2IDlist,self.ID2cd,self.segments,self.segments_valid_id)    
        else:
            self.segments, status_accept, self.segments_valid_id = main_beamInfo(self.beams_data,self.cd2IDlist,self.ID2cd)
            
        
        if status_accept:
            b = list(np.nonzero(self.segments_valid_id))
            b = b[0]
            c = set(range(len(self.segments_valid_id))); 
            b = np.asarray(list(c.difference(set(b))))
            
            if b.size>0:
                for k in b:
                    cl_name,seg_num = self.ID2cd[k]
                    util.showErrorCustomMsg(f"Segment {seg_num} for {cl_name} is not valid","Please correct it as per suggestions. Output file is not created yet")
                    break
                self.setup_ui()
            else:
                self.all_segments_valid = True
                print('All segments data for each is valid') 
                return self.segments
            

    def load(self,parent=None):
        # DialogBeamData = QtWidgets.QDialog(parent=parent)
        return self.setup_ui()
        # DialogBeamData.show()
        # DialogBeamData.exec_()
        


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    beam_info = {"SB":{'number_of_segments':2,'number_of_elements':2}}
    beam_info["S1"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["S2"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["ST"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL1"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU1"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL2"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU2"] = {"number_of_segments":2,"number_of_elements":2}


    # beam_info["CL3"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU3"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL4"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU4"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL5"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU5"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL6"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU6"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL7"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU7"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CL8"] = {"number_of_segments":2,"number_of_elements":2}
    # beam_info["CU8"] = {"number_of_segments":2,"number_of_elements":2}

    

    heights= [3,4]
    sbb = SegmentsByBeams(beam_info,heights)
    sbb.load()
    # sbb.show()
    # sbb.close()
    sys.exit(sbb.exec_())
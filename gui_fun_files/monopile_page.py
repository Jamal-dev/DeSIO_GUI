import sys
from pathlib import Path
import numpy as np
import os

file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir) / Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from structure.monopile.monopile import Monopile
# from structure.jacket.jacket import Jacket
from structure.jacket.jacket_contained_3 import Jacket
from beam.beam import Beam
from beam.segment import Segment
from beam.stand import Stand
from Utils.segments_userInterface import segments_ui
from Utils.geometry import Geometry
from Utils.compartment import Compartment
from Utils.compartments_userInterface import DialogCompartment as dlgCompartment
from Utils.beamsInfo_userInterface import DialogBeamInfo as dlgBeamsInfo

from Utils.segmentsByBeams_userinterface import SegmentsByBeams as dlgSbB


import traceback

from PyQt5.QtWidgets import QDialog
from segment_table import *
from scipy.io import savemat
import matplotlib.pyplot as plt


class MonopilePage:
    def __init__(self, parent=None):
        self.ui = parent

        self.ui.comboBox_Modal.currentIndexChanged.connect(self.load_default_picture)
        # lineStructureMono_StandLength this came from object name inside Qt designer
        self.ui.lineStructureMono_StandLength.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureMono_NoOfSegments.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureMono_NoOfElements.textChanged.connect(self.disableGenbtn)
        self.mpl = self.ui.widStructureMono_mpl
        self.ax = self.mpl.canvas.axes

        self.ui.lineStructureJ_StandLength.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureJ_NoOfBays.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureJ_LOSG.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureJ_Rhead.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureJ_Rfoot.textChanged.connect(self.disableGenbtn)
        
        
        self.beam_types = Beam.getBeamTypes() 	# 0: Stand | 1: Compartment

        # these btns will be used in the monopile
        self.segment_btns = None

    def getValuesMonoPile(self):
        self.mono_page_fields = {}
        self.mono_page_fields["stand_length"] = self.ui.lineStructureMono_StandLength.text()
        self.mono_page_fields["no_segments"] = self.ui.lineStructureMono_NoOfSegments.text()
        self.mono_page_fields["no_elements"] = self.ui.lineStructureMono_NoOfElements.text()
        

    def getValuesJacket(self):
        self.jacket_page_fields = {}
        # Grab J3 text values
        self.jacket_page_fields["stand_length"] = self.ui.lineStructureJ_StandLength.text()
        self.jacket_page_fields["num_bays"] = self.ui.lineStructureJ_NoOfBays.text()
        self.jacket_page_fields["LOSG"] = self.ui.lineStructureJ_LOSG.text()
        self.jacket_page_fields["Rhead"] = self.ui.lineStructureJ_Rhead.text()
        self.jacket_page_fields["Rfoot"] = self.ui.lineStructureJ_Rfoot.text()
        self.jacket_page_fields["num_legs"] = self.ui.lineStructureNumLegs.text()
    
    
    def disableGenbtn(self):
        self.ui.btnStructureTowerGGenrateFile_2.setEnabled(False)

    
    def checkJacketPageinputs(self):
        # stand length J3
        v = util.errorMsg_greaterthan0_float(self.jacket_page_fields["stand_length"], "Stnad length")
        if not v:
            return False
        # no of compartments J3
        v = util.errorMsg_greaterThan0_int(self.jacket_page_fields["num_bays"], "Number of bays")
        if not v:
            return False
        
        # gaps from below
        v = util.errorMsg_greaterOrequal0_float(self.jacket_page_fields["LOSG"], "LOSG")
        if not v:
            return False
        
        # distance above
        v = util.errorMsg_greaterOrequal0_float(self.jacket_page_fields["Rhead"],"Head radius")
        if not v:
            return False
        
        # # distance below
        v = util.errorMsg_greaterOrequal0_float(self.jacket_page_fields["Rfoot"],"Foot radius")
        if not v:
            return False
        # number of legs
        v = util.errorMsg_greaterThan0_int(self.jacket_page_fields["num_legs"], "Number of legs")
        if not v:
            return False
        # convert J3 fields
        self.jacket_page_fields = util.convertEmpty2zero(self.jacket_page_fields)
        self.jacket_page_fields["stand_length"] = float(self.jacket_page_fields["stand_length"])
        self.jacket_page_fields["num_bays"] = int(self.jacket_page_fields["num_bays"])
        self.jacket_page_fields["LOSG"] = float(self.jacket_page_fields["LOSG"])
        self.jacket_page_fields["Rhead"] = float(self.jacket_page_fields["Rhead"])
        self.jacket_page_fields["Rfoot"] = float(self.jacket_page_fields["Rfoot"])
        self.jacket_page_fields["num_legs"] = int(self.jacket_page_fields["num_legs"])
        self.num_compartments = self.jacket_page_fields["num_bays"]
        return True

    def checkMonoPageInputs(self):
        # stand length
        v = util.errorMsg_greaterThan0_int(self.mono_page_fields["stand_length"], "Stnad length")
        if not v:
            return False
        # no of segments
        v = util.errorMsg_greaterThan0_int(self.mono_page_fields["no_segments"], "Number of segments")
        if not v:
            return False

        # no of elements
        v = util.errorMsg_greaterThan0_int(self.mono_page_fields["no_elements"], "Number of elements")
        if not v:
            return False

        
        

        # convert tower fields
        self.mono_page_fields = util.convertEmpty2zero(self.mono_page_fields)
        self.mono_page_fields["stand_length"] = int(self.mono_page_fields["stand_length"])
        self.mono_page_fields["no_segments"] = int(self.mono_page_fields["no_segments"])
        self.mono_page_fields["no_elements"] = int(self.mono_page_fields["no_elements"])
        self.num_compartments = self.mono_page_fields["no_segments"]

        return True

    def checkDialogInput(self, input_fields):
        if isinstance(input_fields, dict):
            values = input_fields
        if isinstance(input_fields, Segment):
            values = {}
            values["id"] = input_fields.segment_id
            values["length_ratio"] = input_fields.length_ratio
            values["Dstar"] = input_fields.diameter_start
            values["Dend"] = input_fields.diameter_end
            values["Tstar"] = input_fields.thickness_start
            values["Tend"] = input_fields.thickness_end
            values["rho"] = input_fields.density
            values["E"] = input_fields.e
            values["G"] = input_fields.g
            values["alpha_s"] = input_fields.alpha_s
            values["alpha_v"] = input_fields.alpha_v
            values["scfStart"] = input_fields.scf_start
            values["scfEnd"] = input_fields.scf_end

        # length_ratio
        v = util.errorMsg_greaterOrequal0_float(values["length_ratio"], "Length ratio")
        if not v:
            return False
        if not values["length_ratio"]:
            values["length_ratio"] = 0.0
        else:
            values["length_ratio"] = float(values["length_ratio"])
        if values["length_ratio"] > 1:
            util.showErrorMsg("Length ratio", "It should not be greater than 1. The sum of all length ratios must be 1")
            return False
        # E
        v = util.errorMsg_greaterthan0_float(values["E"], "E")
        if not v:
            return False

        # G
        v = util.errorMsg_greaterthan0_float(values["G"], "G")
        if not v:
            return False
        # calculating poision ratio E = 2G( 1 + nu)
        nu = float(values["E"]) / (2 * float(values["G"])) - 1
        if nu <= -1 or nu >= 0.5:
            util.showErrorCustomMsg("E and G values are not in correct range",
                                    f"Poison ratio is {nu} which is not between -1 and 0.5")
            return False

        # Dstar
        v = util.errorMsg_greaterOrequal0_float(values["Dstar"], "D_start")
        if not v:
            return False

        # Dend
        v = util.errorMsg_greaterOrequal0_float(values["Dend"], "D_end")
        if not v:
            return False

        # alpha_s
        v = util.errorMsg_greaterOrequal0_float(values["alpha_s"], "alpha_s")
        if not values["alpha_s"]:
            values["alpha_s"] = 0.0
        if not v:
            return False
        elif float(values["alpha_s"]) < 0 or float(values["alpha_s"]) > 1:
            util.showErrorMsg("alpha_s", "It should be greater or equal to 0 and less or equal to 1")
            return False

        # alpha_v

        v = util.errorMsg_greaterOrequal0_float(values["alpha_v"], "alpha_v")
        if not values["alpha_v"]:
            values["alpha_v"] = 0.0
        if not v:
            return False
        elif float(values["alpha_v"]) < 0 or float(values["alpha_v"]) > 1:
            util.showErrorMsg("alpha_v", "It should be greater or equal to 0 and less or equal to 1")
            return False

        # Tstar
        v = util.errorMsg_greaterOrequal0_float(values["Tstar"], "t_start")
        if not v:
            return False

        # Tend
        v = util.errorMsg_greaterOrequal0_float(values["Tend"], "t_end")
        if not v:
            return False

        # scfStart
        v = util.errorMsg_greaterOrequal0_float(values["scfStart"], "SCF_start")
        if not v:
            return False

        # scfEnd
        v = util.errorMsg_greaterOrequal0_float(values["scfEnd"], "SCF_end")
        if not v:
            return False

        # rho
        v = util.errorMsg_greaterOrequal0_float(values["rho"], "Density")
        if not v:
            return False

        return True

    def checkAllSegments(self):
        for id, seg in enumerate(self.segments):
            if not self.checkDialogInput(seg):
                util.showErrorCustomMsg(f"Segment {id + 1} field values were not correct",
                                        "Please correct them as per suggestions. Output file is not created yet")
                return False
        return True

    def dispDialogSegmentTable(self, id=None):
        # current segment

        if id is None or type(id) == bool:
            id = -1 * (self.segment_btns.btn_grp.checkedId() + 2)

        # id = self.cur_segment_id
        cur_segment = self.segments[id]
        Dialog = QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        # Set the labels txts
        ui.lblTitle.setText("Enter data for segment " + str(id + 1))
        ui.lineLengthRatio.setText(str(cur_segment.length_ratio))
        ui.lineDStar.setText(str(cur_segment.diameter_start))
        ui.lineDEnd.setText(str(cur_segment.diameter_end))
        ui.lineTStar.setText(str(cur_segment.thickness_start))
        ui.lineTend.setText(str(cur_segment.thickness_end))
        ui.lineRho.setText(str(cur_segment.density))
        ui.lineE.setText(str(cur_segment.e))
        ui.lineG.setText(str(cur_segment.g))
        ui.lineAlphaS.setText(str(cur_segment.alpha_s))
        ui.lineAlphaV.setText(str(cur_segment.alpha_v))
        ui.lineScfStart.setText(str(cur_segment.scf_start))
        ui.lineScfEnd.setText(str(cur_segment.scf_end))

        Dialog.show()
        resp = Dialog.exec_()

        values = {}
        values["id"] = int(id)
        self.segments_valid_output = False
        if resp == QDialog.Accepted:
            values["length_ratio"] = ui.lineLengthRatio.text()
            values["Dstar"] = ui.lineDStar.text()
            values["Dend"] = ui.lineDEnd.text()
            values["Tstar"] = ui.lineTStar.text()
            values["Tend"] = ui.lineTend.text()
            values["rho"] = ui.lineRho.text()
            values["E"] = ui.lineE.text()
            values["G"] = ui.lineG.text()
            values["alpha_s"] = ui.lineAlphaS.text()
            values["alpha_v"] = ui.lineAlphaV.text()
            values["scfStart"] = ui.lineScfStart.text()
            values["scfEnd"] = ui.lineScfEnd.text()

            cur_segment.segment_id = values["id"]
            cur_segment.length_ratio = values["length_ratio"]
            cur_segment.diameter_start = values["Dstar"]
            cur_segment.diameter_end = values["Dend"]
            cur_segment.thickness_start = values["Tstar"]
            cur_segment.thickness_end = values["Tend"]
            cur_segment.density = values["rho"]
            cur_segment.e = values["E"]
            cur_segment.g = values["G"]
            cur_segment.alpha_s = values["alpha_s"]
            cur_segment.alpha_v = values["alpha_v"]
            cur_segment.scf_start = values["scfStart"]
            cur_segment.scf_end = values["scfEnd"]

            if self.checkDialogInput(values):
                # everything is good now
                # values = util.convertEmpty2zero(values)
                # values = util.convertDic2float(values)
                cur_segment.convertFields2numeric()
                Dialog.close()
                self.segments_valid_output = True

            else:
                # give use another chance for correct inputs
                self.dispDialogSegmentTable(id)

        else:
            print("Cancel is pressed")

    def selectPage(self, text):
        cur_txt = text
        if self.segment_btns:
            # if it exists hide it
            self.segment_btns.hide()


        if cur_txt == 'Please select Input':
            self.ui.stackedWidget.setCurrentWidget(self.ui.Main_page)
            self.select = False # reset the select flag
        elif cur_txt == 'Monopile':
            self.ui.stackedWidget.setCurrentWidget(self.ui.Mono_page)
            self.select = 1 #Monopile
            self.num_stands = 1
        elif 'Jacket' in cur_txt:
            self.ui.stackedWidget.setCurrentWidget(self.ui.Jacket_page)
            self.select = 2 #Jacket structure
        

    def get_compartment_data(self):
        # get the compartment data
        # this data contains id and height of each compartment
        if self.select == 2:
            self.num_compartments = int(self.jacket_page_fields["num_bays"])

        
        self.dlgCompartment = dlgCompartment(self.num_compartments)
        self.compartments_height_data = self.dlgCompartment.load()
        # print("self.compartments_height_data: ",self.compartments_height_data)
    def get_beam_ns_ne(self):
        # get the beam's number of segments and number of elements
        # this data contains number of segments and number of elements for each beam
        self.dlgBeamNsNe = dlgBeamsInfo(self.num_compartments)
        self.beam_info = self.dlgBeamNsNe.load()
        
        
        self.num_beam_classes = self.dlgBeamNsNe.number_beamClasses
        self.beam_class_names = self.dlgBeamNsNe.beam_class_names
        
        if self.beam_info is None:
            # user cancelled
            return False
        else:
            # user entered data
            return True
    def get_beam_segments_data(self):
        # get the beam segments data
        # this data contains the segments data for each beam
        self.dlgBeamSegments = dlgSbB(self.beam_info, self.compartments_height_data)
        self.beam_segments_data = self.dlgBeamSegments.load()
        if self.beam_segments_data is None:
            # user cancelled
            return False
        else:
            # user entered data
            return True
    
    
    def load_jacket(self):
       
        if self.select == 2:
            self.getValuesJacket()
            if not self.checkJacketPageinputs():
                return
            
            num_bays = self.jacket_page_fields["num_bays"]
            L = self.jacket_page_fields["stand_length"]
            r_head = self.jacket_page_fields["Rhead"]
            r_foot = self.jacket_page_fields["Rfoot"]
            losg = self.jacket_page_fields["LOSG"]
            num_legs = self.jacket_page_fields["num_legs"]
            
        
        # Define Jacket object with initial parameters
        # TODO
        self.mpl = self.ui.widStructureJ3_mpl
        self.ax = self.mpl.canvas.axes
        self.ui.J_pic.hide()
        self.ui.widStructureJ3_mpl.show()
        self.mpl.canvas.figure.delaxes(self.mpl.canvas.axes)
        self.mpl.canvas.axes = self.mpl.canvas.figure.add_subplot(111, projection='3d')
        self.ax = self.mpl.canvas.axes
        self.ax.clear()
        
        jacket = Jacket(num_legs = num_legs, 
                num_bays = num_bays,
                L = L, 
                losg = losg, 
                l_tp = 0, 
                r_base = r_foot,
                r_top = r_head,
                axes = self.ax)
        
        # Get the height data from the compartments sheet
        self.get_compartment_data()
        jacket.setBaysHeight(self.compartments_height_data)
        # Get the number of segments and number of elements for each beam class
        status = self.get_beam_ns_ne()
        if not status:
            # user did not enter any data
            return None
        
        
        # Get the segments data for each beam class
        status = self.get_beam_segments_data()
        if not status:
            # user did not enter any data
            return None

        # setting data to the jacket class
        # print(self.beam_info)
        jacket.setBeamInfo(self.beam_info, self.beam_class_names, self.num_beam_classes)
        jacket.setBeamSegment(self.beam_segments_data)
        
        
        return jacket
        
        
        
        
    
    def generate_jacket(self):
        # Load the jacket object
        self.jacket = self.load_jacket()
        if self.jacket is None:
            return
        status = self.jacket.calc_pnts()
        if status:
            self.mpl.canvas.draw()
            util.showInfoMsg(tit="Jacket Input Files Generated", message="Jacket Input Files have been successfully generated.!")
            return True
        else:
            return False
        
    
    def generate_monopile_input_files(self):
        self.monopile_stand.setBeams(self.beams)
        self.monopile_stand.beams[0].setSegments(self.segments)
        # Generate log file containing entered segments input data
        # if self.ui.comboBox_Modal.activated[str].connect(self.display) == "Monopile":  # Monopile

        # Create Monopile object
        self.monopile = Monopile(self.monopile_stand)

        # Raise exception if length ratios not valid
        lengthRatiosValidity = self.monopile.checkLengthRatiosValidity()
        if not lengthRatiosValidity:
            util.showErrorCustomMsg(
                "SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.",
                "Output file is not generated")
            # raise Exception(f'SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.')
            return

        # Generate Beam Input Data (Beam, Segments, Nodes, Elements, etc.) of Monopile
        beamInputGenerated = self.monopile.generateBeamInputData()
        if beamInputGenerated:
            # Write Beam Input File
            self.monopile.writeBeamInput()
            # Write Log File
            self.monopile.writeLogFile()
            # Info Message box stating that input files have successfully been generated
            util.showInfoMsg(tit="Input Files Generated",
                                message="The Monopile input files (Beam Input and Log) have successfully been generated.")
        else:
            pass

    def main_bts_monopile(self):
        self.getValuesMonoPile()
        if not self.checkMonoPageInputs():
            return
        if not self.checkAllSegments():
            return
        self.generate_monopile_input_files()
        self.mpl.canvas.figure.delaxes(self.mpl.canvas.axes)
        self.mpl.canvas.axes = self.mpl.canvas.figure.add_subplot(111)
        self.ax = self.mpl.canvas.axes
        self.visualize_MonopileData()
    def main_bts_jacket(self):
        self.selectPage(self.ui.comboBox_Modal.currentText())
        if self.select == 2:
            # Jacket3 case
            self.getValuesJacket()
            
            if not self.checkJacketPageinputs():
                return
            L = float(self.jacket_page_fields["stand_length"])

            title = f'3D-Simulation (Jacket {self.jacket_page_fields["num_legs"]})'
        else:
            return
        
        self.generate_jacket()
        # TODO: adding graphing capability
        # if self.jacket:
        #     if self.jacket.all_coordinates:
        #         save_path = path_main / Path('test_codes/jacket_input.mat')
        #         savemat(save_path, {"jacket_corrd":self.jacket.all_coordinates, "jacket_line_end_points":self.jacket.all_line_end_points})
        #         self.update_graph(self.jacket.all_coordinates, self.jacket.all_line_end_points, title, scalling_parameter=L)
        #     else:
        #         return

    def visualize_MonopileData(self):
        n_pnts_seg = len(self.segments) + 1
        y_values = np.linspace(start=0, stop=self.mono_page_fields["stand_length"], num=n_pnts_seg)
        seg_y_val = []
        for st, en in zip(y_values, y_values[1:]):
            seg_y_val.append([st, en])
        x = []
        y = []
        for seg, (yst, yen) in zip(self.segments, seg_y_val):
            x.append(seg.diameter_start / 2)
            x.append(seg.diameter_end / 2)
            y.append(yst)
            y.append(yen)
        x = np.asarray(x)
        y = np.asarray(y)
        self.ax.clear()
        self.ax.plot(x, y, color='blue', linewidth=4)
        self.ax.plot(-x, y, color='blue', linewidth=4)
        self.ax.scatter(x, y, s=80, marker='o', c="g", alpha=0.5)
        self.ax.scatter(-x, y, s=80, marker='o', c="g", alpha=0.5)
        self.ax.set_title('Monopile beam Front View')
        self.mpl.canvas.draw()

    def load_default_picture(self):
        self.selectPage(self.ui.comboBox_Modal.currentText())
        if self.select == 1:
            # Monopile
            self.mpl = self.ui.widStructureMono_mpl
            self.ax = self.mpl.canvas.axes
            # image = plt.imread(str(path_main / Path('desio/tower_beam.png')))
            self.ui.Mono_pic.show()  
            self.ui.widStructureMono_mpl.hide()   
    
        # self.mpl.canvas.draw()
        elif self.select == 2:
            # J3 structure
            self.mpl = self.ui.widStructureJ3_mpl
            self.ax = self.mpl.canvas.axes
            # image = plt.imread(str(path_main / Path('desio/J3.png')))
            self.ui.J_pic.show()
            self.ui.widStructureJ3_mpl.hide()
        else:

            return
        # self.mpl.canvas.axes.imshow(image, interpolation = 'nearest', aspect='auto')
    
    def plotDataInAxes(self, scatter_coordinates, line_end_points, windowTitle='3D-Simulation', orthoBaseLength=0.5,
                       scatterEnabled=True, axesOn=True, orthonormalBaseOn=True, scalling_parameter = None):
        # Display scatter points
        if scatterEnabled:
            x_scatter, y_scatter, z_scatter = list(), list(), list()
            for coordinate in scatter_coordinates:
                x_scatter.append(coordinate[0])
                y_scatter.append(coordinate[1])
                z_scatter.append(coordinate[2])
            scalling = list(scalling_parameter * 10 * np.ones(shape=np.shape(x_scatter)))
            # for x,y,z in zip(x_scatter, y_scatter, z_scatter):
            #     print(x,y,z)
            self.ax.scatter(x_scatter, y_scatter, z_scatter, c='b', marker='o', s=scalling)

        # Display Lines connecting edge points
        print("shape(line_end_points): ", np.shape(line_end_points))
        idx_iterator = 0
        for coord_pair in line_end_points:
            
            start_point = coord_pair[0]
            end_point = coord_pair[1]
            if idx_iterator == 0:
                print("shape(coord_pair): ", np.shape(coord_pair))
                print("start_point: ", np.shape(start_point))
                print("end_point: ", np.shape(end_point))
                idx_iterator += 1

            # Display 3D-line connecting start and end points
            self.ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
                         [start_point[2], end_point[2]], 'blue')
            # self.ax.plot3D([start_point[0], end_point[0]], [start_point[1], end_point[1]],
            #     [start_point[2], end_point[2]], 'blue')

        highest_val = -2
        # Display the orthonormal bases
        if orthonormalBaseOn:
            for coordinate in scatter_coordinates:
                # Slice the coordinates into their constituents
                point = coordinate[:3]
                point = list(map(float, point))
                x_base = coordinate[3:6]
                y_base = coordinate[6:9]
                z_base = coordinate[9:]

                # Finds the orthonormal base points from the given point, direction vector and base length
                x_base_point = Geometry.findPointOnVector(point, x_base, orthoBaseLength)
                y_base_point = Geometry.findPointOnVector(point, y_base, orthoBaseLength)
                z_base_point = Geometry.findPointOnVector(point, z_base, orthoBaseLength)

                # Display the orthonormal bases as 3D-lines
                # self.ax.plot3D([point[0], x_base_point[0]], [point[1], x_base_point[1]],
                #     [point[2], x_base_point[2]], 'red')
                # self.ax.plot3D([point[0], y_base_point[0]], [point[1], y_base_point[1]],
                #     [point[2], y_base_point[2]], 'red')
                # self.ax.plot3D([point[0], z_base_point[0]], [point[1], z_base_point[1]],
                #     [point[2], z_base_point[2]], 'red')
                head_ratio = 0.2
                self.ax.quiver(point[0], point[1], point[2], x_base_point[0] - point[0], x_base_point[1] - point[1],
                               x_base_point[2] - point[2], color='red', arrow_length_ratio=head_ratio)
                self.ax.quiver(point[0], point[1], point[2], y_base_point[0] - point[0], y_base_point[1] - point[1],
                               y_base_point[2] - point[2], color='green', arrow_length_ratio=head_ratio)
                self.ax.quiver(point[0], point[1], point[2], z_base_point[0] - point[0], z_base_point[1] - point[1],
                               z_base_point[2] - point[2], color='blue', arrow_length_ratio=head_ratio)
                if z_base_point[2] > highest_val:
                    highest_val = z_base_point[2]
                # self.ax.plot([point[0], x_base_point[0]], [point[1], x_base_point[1]],
                #     [point[2], x_base_point[2]], 'red')
                # self.ax.plot([point[0], y_base_point[0]], [point[1], y_base_point[1]],
                #     [point[2], y_base_point[2]], 'green')
                # self.ax.plot([point[0], z_base_point[0]], [point[1], z_base_point[1]],
                #     [point[2], z_base_point[2]], 'blue')

        # self.ax.set_xlim3d(-0.3, 0.55)
        # self.ax.set_ylim3d(-0.3, 0.55)
        # self.ax.set_zlim3d(-0.01, highest_val)
        # Axes can be disabled
        if not axesOn:
            self.ax.set_axis_off()

    def update_graph(self, scatter_coordinates, line_end_points, windowTitle='3D-Simulation', orthoBaseLength=0.5,
                     scatterEnabled=True, axesOn=True, orthonormalBaseOn=True, scalling_parameter = None):
        if self.select == 1:
            # Monopile
            self.mpl = self.ui.widStructureMono_mpl
            self.ax = self.mpl.canvas.axes
            self.ui.Mono_pic.hide()  
            self.ui.widStructureMono_mpl.show()
            self.ui.widStructureMono_mpl.showMaximized()
            
        elif self.select == 2:
            # J3 structure
            self.mpl = self.ui.widStructureJ3_mpl
            self.ax = self.mpl.canvas.axes
            self.ui.J_pic.hide()
            self.ui.widStructureJ3_mpl.show()
        else:
            util.showErrorCustomMsg('Error', 'Please select a structure')
        
        
        # Saving the input parameters for later use
        self.scatter_coordinates = scatter_coordinates
        self.line_end_points = line_end_points
        self.windowTitle = windowTitle
        self.orthoBaseLength = orthoBaseLength
        self.scatterEnabled = scatterEnabled
        self.axesOn = axesOn
        self.orthonormalBaseOn = orthonormalBaseOn

        # Plot data into the figure
        # self.mpl.set3d()
        self.mpl.canvas.figure.delaxes(self.mpl.canvas.axes)
        self.mpl.canvas.axes = self.mpl.canvas.figure.add_subplot(111, projection='3d')
        self.ax = self.mpl.canvas.axes
        self.ax.clear()
        self.plotDataInAxes(scatter_coordinates, line_end_points, windowTitle, orthoBaseLength,
                            scatterEnabled, axesOn, orthonormalBaseOn, scalling_parameter = scalling_parameter)
        self.mpl.canvas.draw()

        # Geometry.plotDataInAxes(self.ax, scatter_coordinates, line_end_points, windowTitle, orthoBaseLength,
        # scatterEnabled, axesOn, orthonormalBaseOn)

    def load_monopile(self):
        self.getValuesMonoPile()
        
        
        self.selectPage(self.ui.comboBox_Modal.currentText())
        if self.select == 1:
            if not self.checkMonoPageInputs():
                return
        elif self.select == 2:
            if not self.checkJacketPageinputs():
                return
        
        else:
            return
        self.beam_classes = Beam.getBeamClasses(self.num_compartments)
        self.stand_beam_classes = Beam.getStandBeamClasses(self.num_compartments)
        self.comp_beam_classes = Beam.getCompBeamClasses(self.num_compartments)
        # Initialization
        self.monopile_stand = Stand(1, self.mono_page_fields["stand_length"])
        self.monopile_beam = Beam(1, "Stand", self.mono_page_fields["no_segments"], self.mono_page_fields["no_elements"])

        # Append monopile beam to the list of beams (single item list)
        self.beams = list()
        self.beams.append(self.monopile_beam)

        # taking values from segment table
        self.segments = []
        self.segment_btns = segments_ui(self.ui.tabSupportStructure, self.mono_page_fields["no_segments"])
        self.segment_btns.ui_prop()

        # Initializ segment Table base on no_segment
        for id in range(self.mono_page_fields["no_segments"]):
            self.segments.append(Segment(segment_id=id))
            self.cur_segment_id = id
            self.segment_btns.btns[id].clicked.connect(self.dispDialogSegmentTable)
        self.ui.btnStructureTowerGGenrateFile_2.setEnabled(True)

        for id in range(self.mono_page_fields["no_segments"]):
            print("Data entering for segment =", id)
            self.dispDialogSegmentTable(id)
            if not self.segments_valid_output:
                util.showWarningMsg("Beam file is not generated!", 'Generated File')
                return
        self.generate_monopile_input_files()
        self.update_graph(self.monopile.coordinates, self.monopile.line_end_points, '3D-Simulation (Monopile beam)', scalling_parameter=self.mono_page_fields["stand_length"])



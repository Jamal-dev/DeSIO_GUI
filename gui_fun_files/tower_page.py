
import sys
from pathlib import Path
import os

file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir) / Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from structure.tower.tower import Tower
from beam.beam import Beam
from beam.segment import Segment
from Utils.segments_userInterface import segments_ui

from PyQt5.QtWidgets import QDialog
from segment_table import *


class TowerPage:
    def __init__(self,parent=None) :
        # ,stand_length,no_segments=None,no_elements=None,distance_above=None,distance_below=None,
        # self.tower_page_fields = {}
        # self.tower_page_fields["stand_length"] = stand_length
        # self.tower_page_fields["no_segments"] = no_segments
        # self.tower_page_fields["no_elements"] = no_elements
        # self.tower_page_fields["distance_above"] = distance_above
        # self.tower_page_fields["distance_below"] = distance_below
        self.ui = parent
        self.ui.lineStructureTower_NoOfSegments.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureTower_NoOfElements.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureTower_StandLength.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureTower_DistanceAbove.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureTower_DistanceBelow.textChanged.connect(self.disableGenbtn)
    def getValuesFromParent(self):
        self.tower_page_fields = {}
        self.tower_page_fields["stand_length"] = self.ui.lineStructureTower_StandLength.text()
        self.tower_page_fields["no_segments"] = self.ui.lineStructureTower_NoOfSegments.text()
        self.tower_page_fields["no_elements"] = self.ui.lineStructureTower_NoOfElements.text()
        self.tower_page_fields["distance_above"] = self.ui.lineStructureTower_DistanceAbove.text()
        self.tower_page_fields["distance_below"] = self.ui.lineStructureTower_DistanceBelow.text()

    def disableGenbtn(self):
        self.ui.btnStructureTowerGGenrateFile.setEnabled(False)
    
    def checkTowerPageInputs(self):
        # stand length
        v = util.errorMsg_greaterThan0_int(self.tower_page_fields["stand_length"],"Stnad length")
        if not v:
            return False
        # no of segments
        v = util.errorMsg_greaterThan0_int(self.tower_page_fields["no_segments"],"Number of segments")
        if not v:
            return False
        
        # no of elements
        v = util.errorMsg_greaterThan0_int(self.tower_page_fields["no_elements"],"Number of elements")
        if not v:
            return False

        # distance above
        v = util.errorMsg_greaterOrequal0_float(self.tower_page_fields["distance_above"],"Distance above")
        if not v:
            return False
        # distance below
        v = util.errorMsg_greaterOrequal0_float(self.tower_page_fields["distance_below"],"Distance below")
        if not v:
            return False
        # convert tower fields
        self.tower_page_fields = util.convertEmpty2zero(self.tower_page_fields)
        self.tower_page_fields["stand_length"] = int(self.tower_page_fields["stand_length"])
        self.tower_page_fields["no_segments"] = int(self.tower_page_fields["no_segments"])
        self.tower_page_fields["no_elements"] = int(self.tower_page_fields["no_elements"])
        self.tower_page_fields["distance_above"] = float(self.tower_page_fields["distance_above"])
        self.tower_page_fields["distance_below"] = float(self.tower_page_fields["distance_below"])
        
        return True
    
    def checkDialogInput(self,input_fields):
        if isinstance(input_fields,dict):
            values = input_fields
        if isinstance(input_fields,Segment):
            values = {}
            values["id"] = input_fields.segment_id 
            values["length_ratio"] = input_fields.length_ratio
            values["Dstar"] = input_fields.diameter_start
            values["Dend"]  = input_fields.diameter_end
            values["Tstar"]  = input_fields.thickness_start
            values["Tend"] = input_fields.thickness_end
            values["rho"] = input_fields.density
            values["E"] = input_fields.e
            values["G"]  = input_fields.g
            values["alpha_s"] = input_fields.alpha_s
            values["alpha_v"]  = input_fields.alpha_v
            values["scfStart"]  = input_fields.scf_start
            values["scfEnd"] = input_fields.scf_end
        
        # length_ratio
        v = util.errorMsg_greaterOrequal0_float(values["length_ratio"],"Length ratio")
        if not v:
            return False
        if not values["length_ratio"]:
            values["length_ratio"] = 0.0
        else:
            values["length_ratio"] = float(values["length_ratio"])
        if values["length_ratio"]>1:
            util.showErrorMsg("Length ratio","It should not be greater than 1. The sum of all length ratios must be 1")
            return False
        # E
        v = util.errorMsg_greaterthan0_float(values["E"],"E")
        if not v:
            return False
        
        # G
        v = util.errorMsg_greaterthan0_float(values["G"],"G")
        if not v:
            return False
        # calculating poision ratio E = 2G( 1 + nu) 
        nu = float(values["E"])/(2*float(values["G"])) - 1
        if nu<=-1 or nu>=0.5:
            util.showErrorCustomMsg("E and G values are not in correct range",f"Poison ratio is {nu} which is not between -1 and 0.5")
            return False

        # Dstar
        v = util.errorMsg_greaterOrequal0_float(values["Dstar"],"D_start")
        if not v:
            return False
        
        # Dend
        v = util.errorMsg_greaterOrequal0_float(values["Dend"],"D_end")
        if not v:
            return False

        # alpha_s
        v = util.errorMsg_greaterOrequal0_float(values["alpha_s"],"alpha_s")
        if not values["alpha_s"]:
            values["alpha_s"] = 0.0
        if not v:
            return False
        elif float(values["alpha_s"])<0 or float(values["alpha_s"])>1:
            util.showErrorMsg("alpha_s","It should be greater or equal to 0 and less or equal to 1")
            return False
        
        # alpha_v
        
        v = util.errorMsg_greaterOrequal0_float(values["alpha_v"],"alpha_v")
        if not values["alpha_v"]:
            values["alpha_v"] = 0.0
        if not v:
            return False
        elif float(values["alpha_v"])<0 or float(values["alpha_v"])>1:
            util.showErrorMsg("alpha_v","It should be greater or equal to 0 and less or equal to 1")
            return False

        # Tstar
        v = util.errorMsg_greaterOrequal0_float(values["Tstar"],"t_start")
        if not v:
            return False
        
        # Tend
        v = util.errorMsg_greaterOrequal0_float(values["Tend"],"t_end")
        if not v:
            return False
        
        # scfStart
        v = util.errorMsg_greaterOrequal0_float(values["scfStart"],"SCF_start")
        if not v:
            return False
        
        # scfEnd
        v = util.errorMsg_greaterOrequal0_float(values["scfEnd"],"SCF_end")
        if not v:
            return False

        # rho
        v = util.errorMsg_greaterOrequal0_float(values["rho"],"Density")
        if not v:
            return False
        
        return True
    
    def checkAllSegments(self):
        for id,seg in enumerate(self.segments):
            if not self.checkDialogInput(seg):
                util.showErrorCustomMsg(f"Segment {id+1} field values were not correct","Please correct them as per suggestions. Output file is not created yet")
                return False
        return True
    
    def dispDialogSegmentTable(self,id=None):
        # current segment
        
        if id is None or type(id)==bool:
            id = -1*(self.segment_btns.btn_grp.checkedId() + 2)
            
        # id = self.cur_segment_id
        cur_segment = self.segments[id]
        Dialog = QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        # Set the labels txts
        ui.lblTitle.setText(  "Enter data for segment " + str(id+1))
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
            values["Dstar"]= ui.lineDStar.text()
            values["Dend"]= ui.lineDEnd.text()
            values["Tstar"]= ui.lineTStar.text()
            values["Tend"]= ui.lineTend.text()
            values["rho"]= ui.lineRho.text()
            values["E"] = ui.lineE.text()
            values["G"] = ui.lineG.text()
            values["alpha_s"] = ui.lineAlphaS.text()
            values["alpha_v"] = ui.lineAlphaV.text()       
            values["scfStart"]= ui.lineScfStart.text()
            values["scfEnd"]= ui.lineScfEnd.text()
            
            cur_segment.segment_id=values["id"] 
            cur_segment.length_ratio=values["length_ratio"]
            cur_segment.diameter_start=values["Dstar"]
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
            
            if  self.checkDialogInput(values):
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
    
    def writeBeamFile(self):
        # Setting intermediate parameters
        self.tower.setBeams(self.beams)
        self.tower.beams[0].setSegments(self.segments)

        # Raise exception if length ratios not valid
        lengthRatiosValidity = self.tower.checkLengthRatiosValidity()
        if not lengthRatiosValidity:
            util.showErrorCustomMsg("SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.","Output file is not generated")
            # raise Exception(f'SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.')
            return
        # Generate Beam Input Data (Beam, Segments, Nodes, Elements, etc.) of Tower
        beamInputGenerated = self.tower.generateBeamInputData()
        if beamInputGenerated:
            # Write Beam Input File
            self.tower.writeBeamInput()
            # Write Log File
            self.tower.writeLogFile()
            # Info Message box stating that input files have successfully been generated
            util.showInfoMsg(tit="Input Files Generated", message="The Tower input files have successfully been generated!")
            # # Display Monopile graph in image frame
            # self.update_graph(tower.coordinates, tower.line_end_points, '3D-Simulation (Tower)')
        else:
            # Throw exception
            #raise Exception(f'Error: Beam Input Data (Tower) not generated.')	# DEBUG_TEST
            pass # DEBUG_TEST
    
    def main_bts(self):
        self.getValuesFromParent()
        if not self.checkTowerPageInputs():
            return
        if not self.checkAllSegments():
            return
        self.writeBeamFile()
        
        
    def main(self):
        self.getValuesFromParent()
        if not self.checkTowerPageInputs():
            return
        # Initialization
        self.tower = Tower(self.tower_page_fields["stand_length"])
        self.tower_beam = Beam(1, "Stand", self.tower_page_fields["no_segments"], self.tower_page_fields["no_elements"])
        
        # Append tower beam to the list of beams (single item list)
        self.beams = list()
        self.beams.append(self.tower_beam)
        # taking values from segment table
        self.segments = []
        self.segment_btns = segments_ui(self.ui.tabTower, self.tower_page_fields["no_segments"])
        self.segment_btns.ui_prop()
        

        for id in range(self.tower_page_fields["no_segments"]):
            self.segments.append(Segment(segment_id=id))
            self.cur_segment_id = id
            self.segment_btns.btns[id].clicked.connect(self.dispDialogSegmentTable)
        self.ui.btnStructureTowerGGenrateFile.setEnabled(True)
        for id in range(self.tower_page_fields["no_segments"]):
            print("Data entering for segment =",id)
            self.dispDialogSegmentTable(id)
            if  not self.segments_valid_output :
                util.showWarningMsg("Beam file is not generated!",'Generated File')
                return
        self.writeBeamFile()
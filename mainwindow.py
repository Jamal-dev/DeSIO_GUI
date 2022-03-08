
import sys
from pathlib import Path
import os

from sklearn import utils
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from structure.tower.tower import Tower
from beam.beam import Beam
from beam.segment import Segment


# PyQt libs
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QDialog, QMainWindow, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
# import PyQt5.QtCore as QtCore
# from PySide6 import QtCore
from subprocess import call
import os
from form import *
from segment_table import *

def create_py_gui():
    cur_dir = os.getcwd()
    file_name = 'form'
    cmd_input = ["pyuic5", file_name +'.ui' ,
            "-o", file_name +'.py']
    call(cmd_input,timeout=None)

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent=parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.btnStructureTowerSegmentsTable.clicked.connect(self.btnTowerSegmentTable)
        self.ui.btnStructureMonoSegmentsTable.clicked.connect(self.dispDialogSegmentTable)
        self.ui.btnStructureJ3SegmentsTable.clicked.connect(self.dispDialogSegmentTable)

        # self.ui.btnClickMe.clicked.connect(self.dispmsg)
        # QtCore.QObject.connect(self.ui.btnClickMe,QtCore.SIGNAL('clicked()'),self.dispmsg)
    
    def btnTowerSegmentTable(self):
        values = {}
        values["stand_length"] = self.ui.lineStructureTower_StandLength.text()
        values["no_segments"] = self.ui.lineStructureTower_NoOfSegments.text()
        values["no_elements"] = self.ui.lineStructureTower_NoOfElements.text()
        values["distance_above"] = self.ui.lineStructureTower_DistanceAbove.text()
        values["distance_below"] = self.ui.lineStructureTower_DistanceBelow.text()
        def checkInput(values):
            # stand length
            v = util.errorMsg_greaterThan0_int(values["stand_length"],"Stnad length")
            if not v:
                return False
            # no of segments
            v = util.errorMsg_greaterThan0_int(values["no_segments"],"Number of segments")
            if not v:
                return False
            
            # no of elements
            v = util.errorMsg_greaterThan0_int(values["no_elements"],"Number of elements")
            if not v:
                return False

            # distance above
            v = util.errorMsg_greaterOrequal0_float(values["distance_above"],"Distance above")
            if not v:
                return False
            # distance below
            v = util.errorMsg_greaterOrequal0_float(values["distance_below"],"Distance below")
            if not v:
                return False
            
            return True
        if  checkInput(values):
            # everything is good now
            values = util.convertEmpty2zero(values)
            values["stand_length"] = int(values["stand_length"])
            values["no_segments"] = int(values["no_segments"])
            values["no_elements"] = int(values["no_elements"])
            values["distance_above"] = float(values["distance_above"])
            values["distance_below"] = float(values["distance_below"])
        else:
            return
        
        # Initialization
        tower = Tower(values["stand_length"])
        tower_beam = Beam(1, "Stand", values["no_segments"], values["no_elements"])
        segments = list()
        # Append tower beam to the list of beams (single item list)
        beams = list()
        beams.append(tower_beam)
        # taking values from segment table
        segments = []
        for id in range(values["no_segments"]):
            print("Data entering for segment =",id)
            val, status = self.dispDialogSegmentTable(id)
            if status :
                segments.append(val)
            else:
                util.showWarningMsg("Beam file is not generated!",'Generated File')
                return
        # Setting intermediate parameters
        tower.setBeams(beams)
        tower.beams[0].setSegments(segments)

        # Raise exception if length ratios not valid
        lengthRatiosValidity = tower.checkLengthRatiosValidity()
        if not lengthRatiosValidity:
            raise Exception(f'SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.')
        # Generate Beam Input Data (Beam, Segments, Nodes, Elements, etc.) of Tower
        beamInputGenerated = tower.generateBeamInputData()
        if beamInputGenerated:
            # Write Beam Input File
            tower.writeBeamInput()
            # Write Log File
            tower.writeLogFile()
            # Info Message box stating that input files have successfully been generated
            util.showInfoMsg(tit="Input Files Generated", message="The Tower input files have successfully been generated!")
            # # Display Monopile graph in image frame
            # self.update_graph(tower.coordinates, tower.line_end_points, '3D-Simulation (Tower)')
        else:
            # Throw exception
            #raise Exception(f'Error: Beam Input Data (Tower) not generated.')	# DEBUG_TEST
            pass # DEBUG_TEST

            



    def dispDialogSegmentTable(self,id):
        Dialog = QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        ui.lblTitle.setText(  "Enter data for segment " + str(id+1))
        Dialog.show()
        resp = Dialog.exec_()
        def checkInput(values):
            # length_ratio
            v = util.errorMsg_greaterOrequal0_float(values["length_ratio"],"Length ratio")
            if not v:
                return False
            # E
            v = util.errorMsg_greaterthan0_float(values["E"],"E")
            if not v:
                return False
            
            # G
            v = util.errorMsg_greaterthan0_float(values["G"],"G")
            if not v:
                return False

            # Dstar
            v = util.errorMsg_greaterOrequal0_float(values["Dstar"],"D_start")
            if not v:
                return False
            
            # Dend
            v = util.errorMsg_greaterOrequal0_float(values["Dend"],"D_end")
            if not v:
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
        values = {}
        valid_output = False
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
            values["alpah_v"] = ui.lineAlphaV.text()       
            values["scfStart"]= ui.lineScfStart.text()
            values["scfEnd"]= ui.lineScfEnd.text()
            
            
            if  checkInput(values):
                # everything is good now
                values = util.convertEmpty2zero(values)
                values = util.convertDic2float(values)
                values["id"] = int(id)
                Dialog.close()
                valid_output = True
                seg = Segment( segment_id=values["id"], 
                                length_ratio=values["length_ratio"], 
                                diameter_start=values["Dstar"], 
                                diameter_end = values["Dend"], 
                                thickness_start = values["Tstar"], 
                                thickness_end = values["Tend"], 
                                density = values["rho"], 
                                e = values["E"], 
                                g = values["G"], 
                                alpha_s = values["alpha_s"], 
                                alpha_v = values["alpah_v"], 
                                scf_start = values["scfStart"], 
                                scf_end = values["scfEnd"])
                return seg, valid_output
            else:
                # give use another chance for correct inputs
                values, valid_output = self.dispDialogSegmentTable(id)
                return values, valid_output
            
        else:
            print("Cancel is pressed")
            return values, valid_output
            # self.dispDialogSegmentTable(id)

    def dispmsg(self):
        self.ui.lblInfo.setText("Welocme " + self.ui.lineUserName.text() + " to QT")


if __name__ == "__main__":
#    create_py_gui()
    app = QApplication(sys.argv)
#    myapp = MyForm()
    myapp = MainWindow()
    myapp.show()
    sys.exit(app.exec())


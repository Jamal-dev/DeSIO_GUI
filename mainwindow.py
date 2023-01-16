
import sys
from pathlib import Path
import os
from matplotlib.pyplot import step
import matplotlib.pyplot as plt
import numpy as np

from sklearn import utils
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from structure.tower.tower import Tower
from beam.beam import Beam
from beam.segment import Segment
from Utils.segments_userInterface import segments_ui
from gui_fun_files.tower_page import TowerPage
from gui_fun_files.monopile_page import MonopilePage
from Utils.utilities import BrowseLineEdit


# PyQt libs
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QDialog, QMainWindow, QMessageBox, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
# import PyQt5.QtCore as QtCore
# from PySide6 import QtCore
from subprocess import call
import subprocess
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
        self.setWindowIcon(QtGui.QIcon(':newPrefix/desio/J3.png'))
        self.ui.comboBox_Modal.activated[str].connect(self.display)
        self.ui.stackedWidget.setCurrentWidget(self.ui.Main_page)

        self.towerPage = TowerPage(parent=self.ui)
        self.monopilePage = MonopilePage(parent=self.ui)
        self.ui.btnStructureTowerSegmentsTable.clicked.connect(self.towerPage.main)
        self.ui.btnStructureTowerGGenrateFile.clicked.connect(self.towerPage.main_bts)
        self.ui.btnStructureMonoSegmentsTable.clicked.connect(self.monopilePage.load_monopile)
        self.ui.btnStructureJSegmentsTable.clicked.connect(self.monopilePage.main_bts_jacket)
        self.ui.btnStructureTowerGGenrateFile_2.clicked.connect(self.monopilePage.main_bts_monopile)
        self.ui.actionImport.triggered.connect(self.mainmenu_import_browse)
        self.ui.btnSimuSettingGenerate.clicked.connect(self.btnGenerate)
        
        self.ui.widStructureMono_mpl.setWindowModality(QtCore.Qt.ApplicationModal)
        self.ui.widStructureMono_mpl.setWindowTitle('Monopile')

        
        # self.mpl = self.ui.widStructureTower_mpl
        
        # self.mpl.canvas.axes.clear()
   
        # image = plt.imread(str(path_main / Path('desio/tower_beam.png')))
        # self.mpl.canvas.axes.imshow(image)
    
        # self.mpl.canvas.draw()

         
        

        # self.ui.btnStructureMonoSegmentsTable.clicked.connect(self.dispDialogSegmentTable)
        # self.ui.btnStructureJ3SegmentsTable.clicked.connect(self.dispDialogSegmentTable)

        # self.ui.btnClickMe.clicked.connect(self.dispmsg)
        # QtCore.QObject.connect(self.ui.btnClickMe,QtCore.SIGNAL('clicked()'),self.dispmsg)

        self.ui.lineDesioPath = BrowseLineEdit()
        self.ui.lineDesioPath.show()
        self.ui.lineDesioPath.setPlaceholderText("Select the desio path")
        # self.ui.lineDesioPath.setReadOnly(True)
        self.ui.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.ui.lineDesioPath)
        self.ui.lineDesioPath.setObjectName("lineDesioPath")

        self.ui.lineMteePath = BrowseLineEdit()
        self.ui.lineMteePath.show()
        self.ui.lineMteePath.setPlaceholderText("Select the mtee path")
        # self.ui.lineMteePath.setReadOnly(True)
        self.ui.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.ui.lineMteePath)
        self.ui.lineMteePath.setObjectName("lineMteePath")
    
    
    def mainmenu_import_browse(self):
        fname = QFileDialog.getOpenFileName(self, 'Open yaml file', str(path_main),"Image files (*.yaml *.yml)")
        if fname[0]:
            self.yaml_file_path = fname[0]
    def btnGenerate(self):
        '''
           This function is connected with Simulation setting generate button
           This will check the inputs for this page and then subsequently generate the output files
        '''
        # check path of yaml file
        if not hasattr(self, 'yaml_file_path'):
            QMessageBox.about(self, "Error", "Please select the yaml file from File->Import")
            return
        # check if the yaml file is valid
        if not os.path.isfile(self.yaml_file_path):
            QMessageBox.about(self, "Error", "Please select the valid yaml file")
            return
        # check if jobname is not empty
        if not util.errorMsg_empty(self.ui.lineJobName.text(), "Job name"):
            return
        jobname = self.ui.lineJobName.text()
        # check if scaling of blade is float
        if not util.errorMsg_greaterthan0_float(self.ui.lineBladeScaling.text(), "Scaling of blade"):
            return
        blade_scaling = float(self.ui.lineBladeScaling.text())
        # check if scaling of tower is float
        if not util.errorMsg_greaterthan0_float(self.ui.lineTowerScaling.text(), "Scaling of tower"):
            return
        tower_scaling = float(self.ui.lineTowerScaling.text())
        # check if vector sum of yaw axis is 1 and its values are between 0 and 1
        if not util.errorMsg_0to1_vectorSum1(self.ui.lineYawAxis_1.text(),
                                                self.ui.lineYawAxis_2.text(),
                                                self.ui.lineYawAxis_3.text(),
                                                "Yaw axis"):
            return
        yaw_axis = [float(self.ui.lineYawAxis_1.text()),
                    float(self.ui.lineYawAxis_2.text()),
                    float(self.ui.lineYawAxis_3.text())]
        # check if vector sum of tilt axis is 1 and its values are between 0 and 1
        if not util.errorMsg_0to1_vectorSum1(self.ui.lineTiltAxis_1.text(),
                                                self.ui.lineTiltAxis_2.text(),
                                                self.ui.lineTiltAxis_3.text(),
                                                "Tilt axis"):
            return
        tilt_axis = [float(self.ui.lineTiltAxis_1.text()),
                    float(self.ui.lineTiltAxis_2.text()),
                    float(self.ui.lineTiltAxis_3.text())]
        # check if vector MSL poisition is float
        if not util.errorMsg_vectorFloat(self.ui.lineMSLpos_1.text(),
                                    self.ui.lineMSLpos_2.text(),
                                    self.ui.lineMSLpos_3.text(),
                                    "MSL position"):
            return
        msl_pos = [float(self.ui.lineMSLpos_1.text()),
                    float(self.ui.lineMSLpos_2.text()),
                    float(self.ui.lineMSLpos_3.text())]
        # all inputs are valid
        # TODO: add paths for batch files
        path_desio = self.ui.lineDesioPath.text()
        path_mtee = self.ui.lineMteePath.text()
        # read check boxes
        if self.ui.cbBlade.isChecked():
            flag_blade = '--flag_blade'
        else:
            flag_blade = '--no-flag_blade'
        if self.ui.cbBladeAero.isChecked():
            flag_blade_aero = '--flag_blade_aero'
        else:
            flag_blade_aero = '--no-flag_blade_aero'
        if self.ui.cbBladerStructure.isChecked():
            flag_blade_structure = '--flag_blade_struc'
        else:
            flag_blade_structure = '--no-flag_blade_struc'
        if self.ui.cbTower.isChecked():
            flag_tower = '--flag_tower'
        else:
            flag_tower = '--no-flag_tower'
        if self.ui.cbTowerAero.isChecked():
            flag_tower_aero = '--flag_tower_aero'
        else:
            flag_tower_aero = '--no-flag_tower_aero'
        if self.ui.cbTowerStructure.isChecked():
            flag_tower_structure = '--flag_tower_struc'
        else:
            flag_tower_structure = '--no-flag_tower_struc'

        if self.ui.cbFoundationrStructure.isChecked():
            flag_foundation_structure = '--flag_foundation_struc'
        else:
            flag_foundation_structure = '--no-flag_foundation_struc'
        if self.ui.cbFoundationAero.isChecked():
            flag_foundation_aero = '--flag_foundation_aero'
        else:
            flag_foundation_aero = '--no-flag_foundation_aero'
        if self.ui.cbNacelle.isChecked():
            flag_nacelle = '--flag_nacelle'
        else:
            flag_nacelle = '--no-flag_nacelle'
        if self.ui.cbHub.isChecked():
            flag_hub = '--flag_hub'
        else:
            flag_hub = '--no-flag_hub'

        # run python file
        main_file = path_main/Path('yaml_converter_py/read_yaml.py')
        main_file = str(main_file)
        subprocess.run(['python', main_file,
                        '--desio_file_path', path_desio,
                        '--mtee_file_path', path_mtee,
                        '--jobname', jobname,
                        '--yaml_file_path', self.yaml_file_path,
                        '--scaling_blade', str(blade_scaling),
                        '--scaling_tower', str(tower_scaling),
                        '--n_yaw', str(yaw_axis[0]), str(yaw_axis[1]), str(yaw_axis[2]),
                        '--n_tilt', str(tilt_axis[0]), str(tilt_axis[1]), str(tilt_axis[2]),
                        '--cos_msl', str(msl_pos[0]), str(msl_pos[1]), str(msl_pos[2]),
                        flag_blade,
                        flag_blade_aero,
                        flag_blade_structure,
                        flag_tower,
                        flag_tower_aero,
                        flag_tower_structure,
                        flag_foundation_structure,
                        flag_foundation_aero,
                        flag_nacelle,
                        flag_hub,])


    
    def dispmsg(self):
        self.ui.lblInfo.setText("Welocme " + self.ui.lineUserName.text() + " to QT")

    def display(self, text):
        cur_txt = text
        if cur_txt == 'Please select Input':
            self.ui.stackedWidget.setCurrentWidget(self.ui.Main_page)
        elif cur_txt == 'Monopile':
            self.ui.stackedWidget.setCurrentWidget(self.ui.Mono_page)
        elif cur_txt == 'Jacket 3-Stand':
            self.ui.stackedWidget.setCurrentWidget(self.ui.J3_page)
        
if __name__ == "__main__":
#    create_py_gui()
    app = QApplication(sys.argv)
#    myapp = MyForm()
    myapp = MainWindow()
    myapp.show()
    sys.exit(app.exec())



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
from Utils.segments_userInterface import segments_ui
from gui_fun_files.tower_page import TowerPage
from gui_fun_files.monopile_page import MonopilePage


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
        self.ui.comboBox_Modal.activated[str].connect(self.display)
        self.ui.stackedWidget.setCurrentWidget(self.ui.Main_page)

        self.towerPage = TowerPage(parent=self.ui)
        self.monopilePage = MonopilePage(parent=self.ui)
        self.ui.btnStructureTowerSegmentsTable.clicked.connect(self.towerPage.main)
        self.ui.btnStructureTowerGGenrateFile.clicked.connect(self.towerPage.main_bts)
        self.ui.btnStructureMonoSegmentsTable.clicked.connect(self.monopilePage.load_monopile)
        self.ui.btnStructureTowerGGenrateFile_2.clicked.connect(self.monopilePage.main_bts)

         
        

        # self.ui.btnStructureMonoSegmentsTable.clicked.connect(self.dispDialogSegmentTable)
        # self.ui.btnStructureJ3SegmentsTable.clicked.connect(self.dispDialogSegmentTable)

        # self.ui.btnClickMe.clicked.connect(self.dispmsg)
        # QtCore.QObject.connect(self.ui.btnClickMe,QtCore.SIGNAL('clicked()'),self.dispmsg)
    
    
   
    
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
        else:
            self.ui.stackedWidget.setCurrentWidget(self.ui.J4_page)

if __name__ == "__main__":
#    create_py_gui()
    app = QApplication(sys.argv)
#    myapp = MyForm()
    myapp = MainWindow()
    myapp.show()
    sys.exit(app.exec())


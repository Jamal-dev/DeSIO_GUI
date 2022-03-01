
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QDialog, QMainWindow
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
        self.ui.btnStructureTowerSegmentsTable.clicked.connect(self.dispDialogSegmentTable)
        # self.ui.btnClickMe.clicked.connect(self.dispmsg)
        # QtCore.QObject.connect(self.ui.btnClickMe,QtCore.SIGNAL('clicked()'),self.dispmsg)
    def dispDialogSegmentTable(self):
        Dialog = QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        Dialog.show()
        Dialog.exec_()

    def dispmsg(self):
        self.ui.lblInfo.setText("Welocme " + self.ui.lineUserName.text() + " to QT")


if __name__ == "__main__":
#    create_py_gui()
    app = QApplication(sys.argv)
#    myapp = MyForm()
    myapp = MainWindow()
    myapp.show()
    sys.exit(app.exec())


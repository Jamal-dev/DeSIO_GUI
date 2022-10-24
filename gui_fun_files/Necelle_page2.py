
import sys
from pathlib import Path
import os
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox, QDialog
from rigid_body import RigidBody
from Necelle_RB_3DInputPara import Ui_Dialog1

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(385, 369)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label_RB = QtWidgets.QLabel(self.centralwidget)
        self.label_RB.setGeometry(QtCore.QRect(30, 40, 161, 16))
        self.label_RB.setObjectName("label_RB")
        self.label_CP = QtWidgets.QLabel(self.centralwidget)
        self.label_CP.setGeometry(QtCore.QRect(30, 70, 191, 31))
        self.label_CP.setObjectName("label_CP")
        self.input_RB = QtWidgets.QLineEdit(self.centralwidget)
        self.input_RB.setGeometry(QtCore.QRect(230, 40, 81, 21))
        self.input_RB.setObjectName("lineEdit")
        self.input_CP = QtWidgets.QLineEdit(self.centralwidget)
        self.input_CP.setGeometry(QtCore.QRect(230, 80, 81, 21))
        self.input_CP.setObjectName("lineEdit_2")
        self.bt1_genFile = QtWidgets.QPushButton(self.centralwidget)
        self.bt1_genFile.setGeometry(QtCore.QRect(10, 290, 231, 51))
        self.bt1_genFile.setObjectName("bt1_genFile")
        
        self.bt2_genTable = QtWidgets.QPushButton(self.centralwidget)
        self.bt2_genTable.setGeometry(QtCore.QRect(100, 150, 171, 41))
        self.bt2_genTable.setDefault(True)
        self.bt2_genTable.setObjectName("bt2_GenTable")
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.ui_2 = Ui_Dialog1(ob_RB)
        self.Dialog1 = QtWidgets.QDialog()
        # Run the .setupUi() method to show the GUI
        self.ui_2.setupUi(self.Dialog1)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.bt2_genTable.clicked.connect(self.TableDlg)

        self.bt1_genFile.clicked.connect(ob_RB.write_data)


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "NecellePage"))
        self.label_RB.setText(_translate("MainWindow", "<html><head/><body><p>Number of Rigid Bodies</p></body></html>"))
        self.label_CP.setText(_translate("MainWindow", "<html><head/><body><p>Number of Cross-Properties</p></body></html>"))
        self.bt1_genFile.setText(_translate("MainWindow", "Generate File"))
        self.bt2_genTable.setText(_translate("MainWindow", "Generate Table"))

    def TableDlg(self):
        # Creating an instance of the GUI
        self.Dialog1.show()



if __name__ == "__main__":
    import sys
    ob_RB = RigidBody(2,2)
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

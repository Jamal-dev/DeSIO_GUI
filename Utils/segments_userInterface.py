import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))

from PyQt5 import QtCore, QtGui, QtWidgets
# from mainwindow import MainWindow
class segments_ui(QtWidgets.QWidget): # (QtWidgets.QMainWindow)
    def __init__(self,tabTower,number_segments) :
        super(segments_ui, self).__init__()
        self.tabTower = tabTower
        self.number_segments = number_segments
        self.btns = []
        self.btn_clicked_ID = 0
    def ui_prop(self):
        self.scrollArea = QtWidgets.QScrollArea(self.tabTower)
        self.scrollArea.setGeometry(QtCore.QRect(0, 290, 161, 57))
        self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 159, 41))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.scrollAreaWidgetContents)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        
        self.btn_grp = QtWidgets.QButtonGroup()
        self.btn_grp.setExclusive(True)
        
        for i in range(self.number_segments):
            self.btns.append( QtWidgets.QPushButton(self.scrollAreaWidgetContents) )
            self.btns[i].setStyleSheet("border-radius : 20;\n"
                                                    "border : 2px solid black")
            self.btns[i].setObjectName(f"btnStructureTowerSegment{i}")
            self.btns[i].setText(str(i+1))
            self.horizontalLayout.addWidget(self.btns[i])
            self.btns[i].setCheckable(True)
            self.btns[i].show()
            self.btn_grp.addButton(self.btns[i])
        
        self.btn_grp.buttonClicked.connect(self.on_click)
        
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.scrollArea.show()

        self.label = QtWidgets.QLabel(self.tabTower)
        self.label.setGeometry(QtCore.QRect(0, 270, 171, 17))
        self.label.setObjectName("lblStructureTowerSegment")
        self.label.setText("Edit the segments here:")
        self.label.show()
    def on_click(self, btn):
        self.btn_clicked_ID = self.btn_grp.checkedId()

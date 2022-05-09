import sys
from pathlib import Path
import os

from matplotlib.widgets import Widget
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))

from Utils.utilities import Utilities as util
from PyQt5 import QtCore, QtGui, QtWidgets
import PyQt5
from beam.segment import Segment
# from mainwindow import MainWindow
class segments_ui(QtWidgets.QWidget): # (QtWidgets.QMainWindow)
    def __init__(self,parent,number_segments):
        super(segments_ui, self).__init__()
        self.parent = parent
        self.number_segments = number_segments
        self.btns = []
        self.btn_clicked_ID = 0
<<<<<<< HEAD
       
        
    def __repr__(self) -> str:
        return f"Button_Container : it has {len(self.btns)} buttons"
=======

>>>>>>> 97fa1d3f0f6c147af990f054c86924e8b0c29830
    def ui_prop(self):
        self.scrollArea = QtWidgets.QScrollArea(self.parent)
        self.scrollArea.setGeometry(QtCore.QRect(0, 300, 161, 57))
        self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        # self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 159, 41))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.verticalLayout = QtWidgets.QHBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout.setObjectName("horizontalLayout_2")
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
        
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.scrollArea.show()

        self.label = QtWidgets.QLabel(self.parent)
        # self.label.setGeometry(QtCore.QRect(0, 280, 171, 17))
        self.label.setGeometry(QtCore.QRect(0, 280, 350, 25))
        self.label.setObjectName("lblStructureTowerSegment")
        self.label.setText("Edit the segments here:")
        
        self.label.show()

   
     
    def on_click(self, btn):
        self.btn_clicked_ID = self.btn_grp.checkedId()
        

class segments_ui_beamData(QtWidgets.QWidget): # (QtWidgets.QMainWindow)
    def __init__(self,parent,number_segments, btnGrpIds=None,classNamesByBtnIds=None, segments_valid_id = None):
        super(segments_ui_beamData, self).__init__()
        self.parent = parent
        self.number_segments = number_segments
        self.btns = []
        self.btn_clicked_ID = 0
        self.segments_id_local_list = None
        self.segments_id_global_list = None
        self.segments_valid_id = None
        if btnGrpIds is not None:
            self.segments_id_global_list,self.segments_id_local_list = btnGrpIds
        if classNamesByBtnIds is not None:
            self.classNamesByBtnIds = classNamesByBtnIds
        if segments_valid_id is not None:
            self.segments_valid_id = segments_valid_id


    def setupUI_instance(self,widget):
        pos_x,pos_y,width,height=0,0,350,80
        widget.resize(width, height)
        widget.setMinimumSize(QtCore.QSize(width, height))
        self.verticalLayout = QtWidgets.QVBoxLayout(widget)
        
        self.scrollArea = QtWidgets.QScrollArea(widget)
        self.scrollArea.setWidgetResizable(True)

        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 203, 221))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        
        # horizontalLayout
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.scrollAreaWidgetContents)

        self.btn_grp = QtWidgets.QButtonGroup(self.scrollAreaWidgetContents)
        self.btn_grp.setObjectName(f"btn_grp_{pos_x}_{pos_y}")
        self.btn_grp.setExclusive(True)

        # loop on all the segments
        for i in range(self.number_segments):
            self.btns.append( QtWidgets.QPushButton(self.scrollAreaWidgetContents) )
            self.btns[i].setStyleSheet("QPushButton"
                                                    "{"
                                                    "background-color : lightblue;"
                                                    "}"
                                                    "QPushButton::pressed"
                                                    "{"
                                                    "background-color : red;"
                                                    "}")
            self.btns[i].setObjectName(f"btnStructureTowerSegment{i}_{pos_x}_{pos_y}")
            if self.segments_id_local_list is not None:
                self.btns[i].setText(str(self.segments_id_local_list[i]))
            else:
                self.btns[i].setText(str(i+1))
            self.btns[i].setMinimumSize(QtCore.QSize(30, 20))
            self.horizontalLayout.addWidget(self.btns[i])
            self.btns[i].setCheckable(True)
            self.btns[i].show()
            if self.segments_id_global_list is not None:
                self.btn_grp.addButton(self.btns[i],self.segments_id_global_list[i])
            else:
                self.btn_grp.addButton(self.btns[i],i)
        
        self.btn_grp.buttonClicked.connect(self.on_click)
        self.label = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.label.setObjectName(f"lblStructureTowerSegment_{pos_x}_{pos_y}")
        self.label.setText("Edit the segments here:")
        self.label.show()
        
        self.verticalLayout.addWidget(self.label)
        
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        
        self.scrollArea.show()
        QtCore.QMetaObject.connectSlotsByName(widget)
        
    
    
    def setupUI(self,pos_x=0,pos_y=0,width=350,height=80):
        # parent = self.tabTower
        self.scrollArea = QtWidgets.QScrollArea(self.parent)
        self.scrollArea.setGeometry(QtCore.QRect(pos_x, pos_y, width, height))
        # self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget(self.parent)
        
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(pos_x, pos_y, width, height))
        
        self.scrollAreaWidgetContents.setObjectName(f"scrollAreaWidgetContents_{pos_x}_{pos_y}")
        # self.parent.setObjectName(f"scrollAreaWidgetContents_{pos_x}_{pos_y}")
        
        self.verticalLayout = QtWidgets.QVBoxLayout(self.scrollAreaWidgetContents)
        # self.verticalLayout = QtWidgets.QVBoxLayout(self.parent)
        self.verticalLayout.setObjectName(f"verticalLayout_{pos_x}_{pos_y}")
        
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName(f"horizontalLayout_{pos_x}_{pos_y}")
        
        # self.btn_grp = QtWidgets.QButtonGroup(self.scrollAreaWidgetContents)
        self.btn_grp = QtWidgets.QButtonGroup(self.scrollAreaWidgetContents)
        self.btn_grp.setObjectName(f"btn_grp_{pos_x}_{pos_y}")
        self.btn_grp.setExclusive(True)
        
        for i in range(self.number_segments):
            self.btns.append( QtWidgets.QPushButton(self.btn_grp.parent()) )
            # self.btns.append( QtWidgets.QPushButton(self.parent) )
            self.btns[i].setStyleSheet("QPushButton"
                                                    "{"
                                                    "background-color : lightblue;"
                                                    "}"
                                                    "QPushButton::pressed"
                                                    "{"
                                                    "background-color : red;"
                                                    "}")
            self.btns[i].setObjectName(f"btn{self.segments_id_global_list[i]}")
            
            # self.btns[i].clicked.connect(self.on_click)
            if self.segments_id_local_list is not None:
                self.btns[i].setText(str(self.segments_id_local_list[i]))
            else:
                self.btns[i].setText(str(i+1))
            self.btns[i].setMinimumSize(QtCore.QSize(30, 20))
            self.horizontalLayout.addWidget(self.btns[i])
            self.btns[i].setCheckable(True)
            self.btns[i].show()
            if self.segments_id_global_list is not None:
                self.btn_grp.addButton(self.btns[i],self.segments_id_global_list[i])
            else:
                self.btn_grp.addButton(self.btns[i],i)
        
        self.btn_grp.buttonClicked.connect(self.on_click)
        
        self.label = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        # self.label.setGeometry(QtCore.QRect(0, 280, 171, 17))
        # self.label.setGeometry(QtCore.QRect(pos_x, pos_y-height, width, height))
        self.label.setObjectName(f"lblStructureTowerSegment_{pos_x}_{pos_y}")
        self.label.setText("Edit the segments here:")
        self.label.setStyleSheet("font-weight: bold")
        self.label.show()
        
        self.verticalLayout.addWidget(self.label)
        
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        # self.scrollArea.setWidget(self.parent)
        self.scrollArea.show()
        
        return self.scrollAreaWidgetContents
    
        
     
    def on_click(self):
        # print("Clicked: ")
        self.btn_clicked_ID = self.btn_grp.checkedId()
        class_name, segment_number = self.classNamesByBtnIds[self.btn_clicked_ID]
        current_segment = segments[class_name][segment_number-1]
        status = util.dispDialogSegmentTable(cur_segment = current_segment,
                                    title=f"Enter the data for {class_name}, segment {segment_number}",
                                    id   = self.btn_clicked_ID)
        if status:
            print(f"{class_name}, segment {segment_number} updated")
            if self.segments_valid_id is not None:
                self.segments_valid_id[self.btn_clicked_ID-1] = True
        else:
            print(f"{class_name}, segment {segment_number} not updated")
        # print("Clicked: ",self.btn_clicked_ID)


def genrate_multiple_uis_by_previous_data(parent,beams_data,cl2ID,ID2cl,segments_dic,segments_valid_id):
    global segments
    segments = segments_dic
    y_pos = 0
    qwidegts = [QtWidgets.QWidget(parent) for _ in range(len(beams_data))]
    uis = []
    # TODO: change the width and height of qwiidget to the parent
    width = 350
    height = 80
    for nc, (class_name,class_setting) in enumerate(beams_data.items()):
        
        qwidegts[nc].setGeometry(QtCore.QRect(0, y_pos, width, height))
        qwidegts[nc].setMinimumHeight(height)
        qwidegts[nc].setObjectName(f"qw_{class_name}")
        ne = class_setting['number_of_elements']
        ns = class_setting['number_of_segments']
        # segments_id_global, segments_id_local = cl2ID[class_name]
        # segments[class_name] = [Segment(segment_id=id) for id in segments_id_global]
        uis.append(segments_ui_beamData(parent=qwidegts[nc], number_segments=ns,
                                btnGrpIds=cl2ID[class_name],
                                classNamesByBtnIds=ID2cl,
                                segments_valid_id=segments_valid_id))
        # uis[nc].setupUI_layoutStyle(pos_y=y_pos)
        uis[nc].setupUI(pos_y=y_pos)
        uis[nc].label.setText(f"Edit the segments for {class_name}:")
        segments_valid_id = uis[nc].segments_valid_id
        # y_pos += 80
    return qwidegts,uis, segments_valid_id

def generate_multiple_uis(beams_data:dict, cl2ID:dict, ID2cl:dict,parent=None):
    y_pos = 0
    uis = []
    global segments
    segments = {}
    segments_valid_id = [False for _ in range(len(ID2cl))]
    qwidegts = [QtWidgets.QWidget(parent) for _ in range(len(beams_data))]
    # TODO: change the width and height of qwiidget to the parent
    width = 350
    height = 80
    for nc, (class_name,class_setting) in enumerate(beams_data.items()):
        
        qwidegts[nc].setGeometry(QtCore.QRect(0, y_pos, width, height))
        qwidegts[nc].setMinimumHeight(height)
        qwidegts[nc].setObjectName(f"qw_{class_name}")
        ne = class_setting['number_of_elements']
        ns = class_setting['number_of_segments']
        segments_id_global, segments_id_local = cl2ID[class_name]
        segments[class_name] = [Segment(segment_id=id) for id in segments_id_global]
        uis.append(segments_ui_beamData(parent=qwidegts[nc], number_segments=ns,
                                btnGrpIds=cl2ID[class_name],
                                classNamesByBtnIds=ID2cl,
                                segments_valid_id=segments_valid_id))
        
        uis[nc].setupUI(pos_y=y_pos)
        uis[nc].label.setText(f"Edit the segments for {class_name}:")
        
        segments_valid_id = uis[nc].segments_valid_id
        # y_pos += 80
    return qwidegts,uis, segments_valid_id


class Communicate(QtCore.QObject):
    accepted = QtCore.pyqtSignal()
    rejected = QtCore.pyqtSignal()    

class BeamInfoDataWindow(QtWidgets.QDialog) :
    def __init__(self,beams_data:dict, cl2ID:dict, ID2cl:dict, segments=None, segments_valid_id=None, parent=None):
        super().__init__()
        width = 350
        self.setGeometry(300, 300, width, 400)
        self.setWindowTitle("Enter Segments data for each beam")

        self.segments = segments

        self.widget = QtWidgets.QWidget(parent=parent)
        
        self.scrollArea = QtWidgets.QScrollArea(parent=parent)
        self.layout = QtWidgets.QVBoxLayout(self.widget)
                
        
        self.scrollArea.setWidgetResizable(True) 
        self.widget.resize(width, 95*len(beams_data))
        
        
        

        if segments is None and segments_valid_id is None:
            self.qwidgets,self.uis ,self.segments_valid_id = generate_multiple_uis(beams_data=beams_data,
                                                                cl2ID=cl2ID,ID2cl=ID2cl, 
                                                                parent=self.widget)
        else:
            self.qwidgets,self.uis, self.segments_valid_id = genrate_multiple_uis_by_previous_data(beams_data=beams_data,
                                                                cl2ID=cl2ID,ID2cl=ID2cl, 
                                                                segments_dic=segments,
                                                                segments_valid_id=segments_valid_id,
                                                                parent=self.widget)
                
        # Ok and Cancel button
        self.buttonBox = QtWidgets.QDialogButtonBox(self.widget)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBoxOK")
        self.buttonBox.accepted.connect(self.accept) # type: ignore
        self.buttonBox.rejected.connect(self.reject) # type: ignore

        self.keyPressEvent = self.keyPressEvent
        self.closeEvent = self.closeEvent
        
        self.buttonBox.setGeometry(QtCore.QRect(0, self.widget.height()-50, self.widget.width(), 50))

        # self.setCentralWidget(self.scroll)
        
        
        for i in range(len(self.qwidgets)):
            self.layout.addWidget(self.qwidgets[i])         
        self.layout.addWidget(self.buttonBox)
        

        self.scrollArea.setWidget(self.widget)
        self.layout2 = QtWidgets.QVBoxLayout(self)
        self.layout2.addWidget(self.scrollArea)
        self.setLayout(self.layout2)

        self.c = Communicate()
        self.c.accepted.connect(self.acceptedTrue)
        self.c.rejected.connect(self.rejectTrue)
        self.accepted = False
        self.rejected = False
    def acceptedTrue(self):
        self.accepted = True
    def closeEvent(self,e):
        e.accept()

    def rejectTrue(self):
        self.rejected = True    
    
    def keyPressEvent(self, e):
        # print("event", e)
        if e.key()  == QtCore.Qt.Key_Return :
            self.accept()
        elif e.key() == QtCore.Qt.Key_Enter :
            self.accept()
            
        elif e.key() == QtCore.Qt.Key_Escape :
            self.reject()
            
    def accept(self):
        self.c.accepted.emit()
        self.close()
        # print("Accepted")
        # segments is the global variable
        self.segments = segments
        
    def reject(self):
        self.c.rejected.emit()
        self.close()
        self.accepted = False
        # print("Rejected")
        self.segments = segments
        
      
        


def test_beamInfo():
    beam_info = {"SB":{'number_of_segments':2,'number_of_elements':2}}
    beam_info["S1"] = {"number_of_segments":2,"number_of_elements":2}
    beam_info["S2"] = {"number_of_segments":2,"number_of_elements":2}
    
    cl2ID ={}
    cl2ID["SB"] = [0,1], [1,2]
    cl2ID["S1"] = [2,3], [1,2]
    cl2ID["S2"] = [4,5], [1,2]

    ID2cl = {}
    ID2cl[0] = "SB", 1
    ID2cl[1] = "SB", 2
    ID2cl[2] = "S1", 1
    ID2cl[3] = "S1", 2
    ID2cl[4] = "S2", 1
    ID2cl[5] = "S2", 2

    return beam_info, cl2ID, ID2cl

def main_beamInfo(beam_data, cl2ID, ID2cl,segments = None, segments_valid_id = None,parent=None):
    # app = QtWidgets.QApplication(sys.argv)
    # dlg = QtWidgets.QDialog(parent=parent)
    if segments is None and segments_valid_id is None:
        beamInfoWindow = BeamInfoDataWindow(beam_data, cl2ID,ID2cl,parent=parent)
    else:
        beamInfoWindow = BeamInfoDataWindow(beam_data, cl2ID,ID2cl,segments,segments_valid_id,parent=parent)
    beamInfoWindow.show()
    if  beamInfoWindow.exec_() == QtWidgets.QDialog.Accepted:
        print("yes")
    else:
        print("no")
    # dlg.show()
    # run = dlg.exec()
    aceepted_true = False
    if beamInfoWindow.accepted:
        print("Accepted the signal" )
        aceepted_true = True
    if beamInfoWindow.rejected:
        print("Rejected the signal")
    # sys.exit(run)
    segments_valid_id = beamInfoWindow.segments_valid_id
    # segments is a dictionary: segments[class_name][segment_number-1] where segment_number starts from 1
    segments = beamInfoWindow.segments
    del beamInfoWindow
    # del run
    # del dlg
    return segments, aceepted_true, segments_valid_id
    

if __name__ == "__main__":
    # import sys
    # app = QtWidgets.QApplication(sys.argv)
    # MainWindow = QtWidgets.QMainWindow()
    
    # # ui = segments_ui(MainWindow,16)
    # # ui2 = segments_ui(MainWindow,16)
    
    # # ui.setupUI_layoutStyle(pos_y=0)
    # # ui2.setupUI_layoutStyle(pos_y=100)
    # main(number_classes=9, parent=MainWindow)

    # MainWindow.show()
    # sys.exit(app.exec_())

    app = QtWidgets.QApplication(sys.argv)
    beam_data, cl2ID, ID2cl = test_beamInfo()
    main_beamInfo(beam_data, cl2ID, ID2cl)
    sys.exit(app.exec_())
    
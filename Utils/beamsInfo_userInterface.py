import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from PyQt5 import QtCore, QtGui, QtWidgets
# from mainwindow import MainWindow
class beamsinfo_ui(QtWidgets.QWidget): # (QtWidgets.QMainWindow)
    def __init__(self,number_compartments) :
        super(beamsinfo_ui, self).__init__()
        self.number_compartments = number_compartments
        self.number_beamClasses = 3 * number_compartments + 2
        self.no_segments_beam = []
        self.no_elements_beam = []
        self.lbls = []
        self.setBeamNames()
    def setBeamNames(self):
        self.beamNames = []
        self.beamNameDescriptions = {}
        self.id2beamName = {}
        self.beamName2id = {}
        flag = True
        for i in range(self.number_beamClasses):
            if i==0:
                self.beamNames.append("SB")
                self.beamNameDescriptions[self.beamNames[i]] = "Stand Bottom"
                self.id2beamName[i] = self.beamNames[i]
                self.beamName2id[self.beamNames[i]] = i
                continue
            if i>0 and i<=self.number_compartments:
                self.beamNames.append("S"+str(i))
                self.beamNameDescriptions[self.beamNames[i]] = "Stand "+str(i)
                self.id2beamName[i] = self.beamNames[i]
                self.beamName2id[self.beamNames[i]] = i
                continue
            if i==self.number_compartments+1:
                self.beamNames.append("ST")
                self.beamNameDescriptions[self.beamNames[i]] = "Stand Top"
                self.id2beamName[i] = self.beamNames[i]
                self.beamName2id[self.beamNames[i]] = i
                continue
            if i>self.number_compartments+1 and i<=self.number_beamClasses:
                if flag:
                    counter = 1
                    counterL = 1
                    counterU = 1
                    flag = False
                if not counter%2 == 0:
                    # odd case
                    self.beamNames.append(f"C{counterL}L")
                    self.beamNameDescriptions[self.beamNames[i]] = "Compartment "+str(counterL)+" (Lower)"
                    counter += 1
                    counterL += 1
                else:
                    self.beamNames.append(f"C{counterU}U")
                    self.beamNameDescriptions[self.beamNames[i]] = "Compartment "+str(counterU)+" (Upper)"
                    counter += 1
                    counterU += 1
                self.id2beamName[i] = self.beamNames[i]
                self.beamName2id[self.beamNames[i]] = i
                continue
            
        
    def ui_prop(self,DialogBeamInfo):
        DialogBeamInfo.setObjectName("DialogBeamInfo")
        DialogBeamInfo.resize(450, 300)
        
        self.gridLayout = QtWidgets.QGridLayout(DialogBeamInfo)
        self.gridLayout.setObjectName("gridLayout")
        self.scrollArea = QtWidgets.QScrollArea(DialogBeamInfo)
        # self.scrollArea.setGeometry(QtCore.QRect(0, 300, 161, 57))
        # self.scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 203, 221))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        
        self.gl_hor = QtWidgets.QGridLayout(self.scrollAreaWidgetContents)
        self.gl_hor.setContentsMargins(0, 0, 0, 0)
        self.gl_hor.setObjectName("gl_hor")
        self.gridLayoutlblsAndTxt = QtWidgets.QGridLayout()
        self.gridLayoutlblsAndTxt.setObjectName("gridLayoutlblsAndTxt")
        
        

        fontHeader = QtGui.QFont()
        fontHeader.setPointSize(11)
        fontHeader.setBold(True)

        self.lblID = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblID.setAlignment(QtCore.Qt.AlignCenter)
        self.lblID.setObjectName("lblID")
        self.lblID.setText("ID")

        self.lblBeamClass = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblBeamClass.setAlignment(QtCore.Qt.AlignCenter)
        self.lblBeamClass.setObjectName("lblBeamClass")
        self.lblBeamClass.setText("Beam Class")

        self.lblDescription = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblDescription.setAlignment(QtCore.Qt.AlignCenter)
        self.lblDescription.setObjectName("lblDescription")
        self.lblDescription.setText("Description")

        self.lblNoSegments = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblNoSegments.setAlignment(QtCore.Qt.AlignCenter)
        self.lblNoSegments.setObjectName("lblNoSegments")
        self.lblNoSegments.setText("No of Segments")

        self.lblNoElements = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblNoElements.setAlignment(QtCore.Qt.AlignCenter)
        self.lblNoElements.setObjectName("lblNoElements")
        self.lblNoElements.setText("No of Elements")

        self.lblID.setFont(fontHeader)
        self.lblID.setAlignment(QtCore.Qt.AlignCenter)
        self.lblBeamClass.setFont(fontHeader)
        self.lblBeamClass.setAlignment(QtCore.Qt.AlignCenter)
        self.lblDescription.setFont(fontHeader)
        self.lblDescription.setAlignment(QtCore.Qt.AlignCenter)
        self.lblNoSegments.setFont(fontHeader)
        self.lblNoSegments.setAlignment(QtCore.Qt.AlignCenter)
        self.lblNoElements.setFont(fontHeader)
        self.lblNoElements.setAlignment(QtCore.Qt.AlignCenter)

        self.gridLayoutlblsAndTxt.addWidget(self.lblID, 0, 0, 1, 1)
        self.gridLayoutlblsAndTxt.addWidget(self.lblBeamClass, 0, 1, 1, 1)
        self.gridLayoutlblsAndTxt.addWidget(self.lblDescription, 0, 2, 1, 1)
        self.gridLayoutlblsAndTxt.addWidget(self.lblNoSegments, 0, 3, 1, 1)
        self.gridLayoutlblsAndTxt.addWidget(self.lblNoElements, 0, 4, 1, 1)
        
        
        
        
        for i in range(self.number_beamClasses):
            self.no_segments_beam.append( QtWidgets.QLineEdit(self.scrollAreaWidgetContents) )
            self.no_elements_beam.append( QtWidgets.QLineEdit(self.scrollAreaWidgetContents) )
            self.lbls.append( QtWidgets.QLabel(self.scrollAreaWidgetContents) )
            lblBeamClass = QtWidgets.QLabel(self.scrollAreaWidgetContents)
            lblBeamDescription = QtWidgets.QLabel(self.scrollAreaWidgetContents)
            self.no_segments_beam[i].setObjectName(f"no_segments_beam{i}")
            self.no_elements_beam[i].setObjectName(f"no_elements_beam{i}")
            self.lbls[i].setObjectName(f"lblID{i}")
            self.lbls[i].setText(str(i+1))
            lblBeamClass.setText(self.beamNames[i])
            lblBeamDescription.setText(self.beamNameDescriptions[self.beamNames[i]])
            self.gridLayoutlblsAndTxt.addWidget(self.lbls[i], i+1, 0, 1, 1)
            self.gridLayoutlblsAndTxt.addWidget(lblBeamClass, i+1, 1, 1, 1)
            self.gridLayoutlblsAndTxt.addWidget(lblBeamDescription, i+1, 2, 1, 1)
            self.gridLayoutlblsAndTxt.addWidget(self.no_segments_beam[i], i+1, 3, 1, 1)
            self.gridLayoutlblsAndTxt.addWidget(self.no_elements_beam[i], i+1, 4, 1, 1)
            
            
            self.no_segments_beam[i].show()
            self.no_elements_beam[i].show()
            self.lbls[i].show()
            lblBeamClass.show()
            lblBeamDescription.show()
        
        self.gl_hor.addLayout(self.gridLayoutlblsAndTxt, 0, 0, 1, 1)
        
        
        
        
        
        
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.scrollArea.show()

        self.gridLayout.addWidget(self.scrollArea, 1, 4, 1, 4)
        self.buttonBox = QtWidgets.QDialogButtonBox(DialogBeamInfo)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 3, 7, 1, 1)
        self.lblTitle = QtWidgets.QLabel(DialogBeamInfo)
        
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)

        
        
        self.lblTitle.setFont(font)
        self.lblTitle.setAlignment(QtCore.Qt.AlignCenter)
        self.lblTitle.setObjectName("lblTitle")
        self.gridLayout.addWidget(self.lblTitle, 0, 6, 1, 2)

        self.retranslateUi(DialogBeamInfo)
        self.buttonBox.accepted.connect(DialogBeamInfo.accept) # type: ignore
        self.buttonBox.rejected.connect(DialogBeamInfo.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(DialogBeamInfo)

    def checkInputs(self):
        """
        Checks the inputs of the dialog.
        """
        
        self.getBeamsInfo()
        for i in range(self.number_beamClasses):
            
            if self.no_segments_beam[i].text() == "":
                self.no_segments_beam[i].setText("0")
            if self.no_elements_beam[i].text() == "":
                self.no_elements_beam[i].setText("0")
            
            sg = util.errorMsg_greaterThan0_int(self.no_segments_beam[i].text(),f"Number of segments for {self.beamNames[i]}")
            ne = util.errorMsg_greaterThan0_int(self.no_elements_beam[i].text(),f"Number of elements for {self.beamNames[i]}")
            if not (sg and ne):
                return False
        return True
    def getBeamsInfo(self):
        """
        Get the inputs of the dialog.
        """
        self.beam_info = {}
        self.beam_no_segments_data = []
        self.beam_no_elements_data = []
        for i in range(self.number_beamClasses):
            self.beam_info[self.beamNames[i]] = {}
            if self.no_segments_beam[i].text() == "":
                self.no_segments_beam[i].setText("0")
            if self.no_elements_beam[i].text() == "":
                self.no_elements_beam[i].setText("0")
            self.beam_info[self.beamNames[i]]["no_segments"] = int(self.no_segments_beam[i].text())
            self.beam_info[self.beamNames[i]]["no_elements"] = int(self.no_elements_beam[i].text())
            self.beam_no_segments_data.append(int(self.no_segments_beam[i].text()))
            self.beam_no_elements_data.append(int(self.no_elements_beam[i].text()))
        return self.beam_info, self.beam_no_segments_data, self.beam_no_elements_data
    def setBeamData(self,beam_no_segments_data,beam_no_elements_data):

        """
        Set the inputs of the dialog.
        """
        self.getBeamsInfo()
        for i in range(self.number_beamClasses):
            self.no_segments_beam[i].setText(str(beam_no_segments_data[i]))
            self.no_elements_beam[i].setText(str(beam_no_elements_data[i]))
        
            
        
    
    
    def retranslateUi(self, DialogCompartmentTable):
        _translate = QtCore.QCoreApplication.translate
        DialogCompartmentTable.setWindowTitle(_translate("DialogBeamInfo", "Beams Information"))
        self.lblTitle.setText(_translate("DialogBeamInfo", "Beams Data"))

class DialogBeamInfo:
    """
    Dialog to set the number of segments and elements for each beam.
    """
    def __init__(self, num_compartments):
        self.num_compartments = num_compartments
        self.beam_info = {}
        self.beam_no_segments_data = []
        self.beam_no_elements_data = []
        self.id2beamName = {}
        self.beamName2id = {}
        self.number_beamClasses = 0
    def setupUi(self):
        self.DialogBeams = QtWidgets.QDialog()
        self.ui = beamsinfo_ui(self.num_compartments)
        self.ui.ui_prop(self.DialogBeams)

        self.id2beamName = self.ui.id2beamName
        self.beamName2id = self.ui.beamName2id
        self.number_beamClasses = self.ui.number_beamClasses

        self.ui.no_segments_beam[0].setFocus()
        if self.beam_no_elements_data or self.beam_no_segments_data:
            self.setBeamInfo()
        self.resp = self.DialogBeams.exec_()
        if self.resp == QtWidgets.QDialog.Accepted:
            return self.update_cells()
        
        
    def update_cells(self):
        self.beam_info,self.beam_no_segments_data,self.beam_no_elements_data = self.ui.getBeamsInfo()
        if self.ui.checkInputs():
            return True
        else:
            self.DialogBeams.close()
            print("Error in the inputs")
            return False          
    
    def load(self):
        status = self.setupUi()
        if isinstance(status,bool) and status:
            print("Beam info loaded")
        elif isinstance(status,bool) and not status:
            self.load()
        else:
            print("cancel is pressed")
            

    def getBeamInfo(self):
        return self.beam_info
    
    def setBeamInfo(self):
        self.ui.setBeamData(self.beam_no_segments_data,self.beam_no_elements_data)
        
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    # DialogBeamsInfo = QtWidgets.QDialog()
    # number_compartments = 3
    # ui = beamsinfo_ui(number_compartments)
    # ui.ui_prop(DialogBeamsInfo)
    # ui.no_segments_beam[0].setFocus()
    # DialogBeamsInfo.show()
    # resp = DialogBeamsInfo.exec_()
    # if resp == QtWidgets.QDialog.Accepted:
    #     if ui.checkInputs():
    #         print("Accepted")
    #     else:
    #         print("Rejected")
    # sys.exit(app.exec_())
    dlg = DialogBeamInfo(3)
    dlg.load()
    dlg.close()
    sys.exit(app.exec_())
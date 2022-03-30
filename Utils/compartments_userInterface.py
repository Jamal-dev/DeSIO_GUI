import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))

from Utils.utilities import Utilities as util

from PyQt5 import QtCore, QtGui, QtWidgets
# from mainwindow import MainWindow
class compartments_ui(QtWidgets.QWidget): # (QtWidgets.QMainWindow)
    def __init__(self,number_compartments) :
        super(compartments_ui, self).__init__()
        self.number_compartments = number_compartments
        self.txts = []
        self.lbls = []
        
    def ui_prop(self,DialogCompartmentTable):
        DialogCompartmentTable.setObjectName("DialogCompartmentTable")
        DialogCompartmentTable.resize(223, 300)
        
        self.gridLayout = QtWidgets.QGridLayout(DialogCompartmentTable)
        self.gridLayout.setObjectName("gridLayout")
        self.scrollArea = QtWidgets.QScrollArea(DialogCompartmentTable)
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
        self.gridLayoutlblsAndTxt.setObjectName("glLbl")
        
        

        fontHeader = QtGui.QFont()
        fontHeader.setPointSize(11)
        fontHeader.setBold(True)

        self.lblID = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblID.setAlignment(QtCore.Qt.AlignCenter)
        self.lblID.setObjectName("lblID")
        self.lblID.setText("ID")

        self.lblHeight = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.lblHeight.setAlignment(QtCore.Qt.AlignCenter)
        self.lblHeight.setObjectName("lblHeight")
        self.lblHeight.setText("Height")

        self.lblID.setFont(fontHeader)
        self.lblID.setAlignment(QtCore.Qt.AlignCenter)
        self.lblHeight.setFont(fontHeader)
        self.lblHeight.setAlignment(QtCore.Qt.AlignCenter)

        self.gridLayoutlblsAndTxt.addWidget(self.lblID, 0, 0, 1, 1)
        self.gridLayoutlblsAndTxt.addWidget(self.lblHeight, 0, 1, 1, 1)
        
        
        
        
        for i in range(self.number_compartments):
            self.txts.append( QtWidgets.QLineEdit(self.scrollAreaWidgetContents) )
            self.lbls.append( QtWidgets.QLabel(self.scrollAreaWidgetContents) )
            self.txts[i].setObjectName(f"txtHeight{i}")
            self.lbls[i].setObjectName(f"lblID{i}")
            self.lbls[i].setText(str(i+1))
            self.gridLayoutlblsAndTxt.addWidget(self.lbls[i], i+1, 0, 1, 1)
            self.gridLayoutlblsAndTxt.addWidget(self.txts[i], i+1, 1, 1, 1)
            
            
            self.txts[i].show()
            self.lbls[i].show()
        
        self.gl_hor.addLayout(self.gridLayoutlblsAndTxt, 0, 0, 1, 1)
        
        
        
        
        
        
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.scrollArea.show()

        self.gridLayout.addWidget(self.scrollArea, 1, 4, 1, 4)
        self.buttonBox = QtWidgets.QDialogButtonBox(DialogCompartmentTable)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 3, 7, 1, 1)
        self.lblTitle = QtWidgets.QLabel(DialogCompartmentTable)
        
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)

        
        
        self.lblTitle.setFont(font)
        self.lblTitle.setAlignment(QtCore.Qt.AlignCenter)
        self.lblTitle.setObjectName("lblTitle")
        self.gridLayout.addWidget(self.lblTitle, 0, 6, 1, 2)

        self.retranslateUi(DialogCompartmentTable)
        self.buttonBox.accepted.connect(DialogCompartmentTable.accept) # type: ignore
        self.buttonBox.rejected.connect(DialogCompartmentTable.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(DialogCompartmentTable)

    def retranslateUi(self, DialogCompartmentTable):
        _translate = QtCore.QCoreApplication.translate
        DialogCompartmentTable.setWindowTitle(_translate("DialogCompartmentTable", "Compartment Table"))
        self.lblTitle.setText(_translate("DialogCompartmentTable", "Compartments Table"))
    
    def get_comp_values(self):
        self.comp_values = []
        for i in range(self.number_compartments):
            if self.txts[i].text() == "":
                self.txts[i].setText("0")
            self.comp_values.append(self.txts[i].text())
        return self.comp_values
        
    def check_values(self):
        self.get_comp_values()
        for i in range(self.number_compartments):
            if self.txts[i].text() == "":
                self.txts[i].setText("0")
            elif not self.txts[i].text().isnumeric():
                util.showErrorCustomMsg(f"Compartment {i+1} height must be a number")
                return False
            v = util.errorMsg_greaterthan0_float(self.txts[i].text(),f"Compartment {i+1}")
            if not v:
                return False
        return True
    
class DialogCompartment:
    def __init__(self,num_compartments):
        self.num_compartments = num_compartments
        self.height_values = []
        
    def setupUi(self):
        self.DialogCompartment = QtWidgets.QDialog()
        self.dialog = compartments_ui(self.num_compartments)
        self.dialog.ui_prop(self.DialogCompartment)

        self.dialog.txts[0].setFocus()
        
        if self.height_values:
            
            for i in range(self.num_compartments):
                self.dialog.txts[i].setText(str(self.height_values[i]))
        self.resp = self.DialogCompartment.exec_()
        if self.resp == QtWidgets.QDialog.Accepted:
            if self.dialog.check_values():
                self.height_values = self.dialog.get_comp_values()
                self.height_values = list(map(float,self.height_values))
                self.DialogCompartment.close()
                return True
            else:
                self.height_values = self.dialog.get_comp_values()
                self.height_values = list(map(float,self.height_values))
                self.DialogCompartment.close()
                print("Error in height values")
                return False
    def load(self):
        status = self.setupUi()
        
        if isinstance(status,bool) and status:
            print("Heights are loadded")
            return self.height_values
        elif isinstance(status,bool) and not status:
            self.load()
        else:
            print("cancel is pressed")
        


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    # DialogCompartment = QtWidgets.QDialog()
    # number_compartments = 11
    # ui = compartments_ui(number_compartments)
    # ui.ui_prop(DialogCompartment)
    # ui.txts[0].setFocus()
    # DialogCompartment.show()
    dlg = DialogCompartment(11)
    dlg.load()
    dlg.close()
    sys.exit(app.exec_())
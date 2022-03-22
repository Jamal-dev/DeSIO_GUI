from PyQt5.QtWidgets import*
from matplotlib import projections
from mpl_toolkits.mplot3d import Axes3D

from  matplotlib.backends.backend_qt5agg  import  FigureCanvas

from  matplotlib.figure  import  Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar

    
class mplwidget(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvasQTAgg(Figure())
        
        toolbar = NavigationToolbar(self.canvas,self)
        self.vertical_layout = QVBoxLayout()
        self.vertical_layout.addWidget(toolbar)
        self.vertical_layout.addWidget(self.canvas)
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(self.vertical_layout)
    def set3d(self):
        self.canvas.figure.delaxes(self.canvas.axes)
        self.canvas.axes = self.canvas.figure.add_subplot(111,projection='3d')
        # self.setLayout(self.vertical_layout)
    def set2d(self):
        self.canvas.figure.delaxes(self.canvas.axes)
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        # self.setLayout(self.vertical_layout)
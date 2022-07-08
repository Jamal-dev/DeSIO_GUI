import sys
from pathlib import Path
import os
from turtle import color
import numpy as np
import matplotlib.pyplot as plt

file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir) / Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from Utils.geometry import Geometry
from structure.tower.tower import Tower
from beam.beam import Beam
from beam.segment import Segment
from Utils.segments_userInterface import segments_ui

from PyQt5.QtWidgets import QDialog
from segment_table import *
from tower_page import *


class NecellePage:

    #def checkNecellePageInputs

import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities

class Node():

	def __init__(self, local_num, global_num, coordinates):
		self.local_num = local_num
		self.global_num = global_num
		self.coordinates = coordinates

	def getNodeAsStr(self):
		formatted_coords = Utilities.getFormattedCoordList(self.coordinates)
		return f"[Local Num: {self.local_num}, Global Num: {self.global_num},\nCoordinates: {formatted_coords}]"

	@staticmethod
	def getNodesAsStr(nodes):
		nodesStr = f"["
		if nodes != None:
			i = 0
			for node in nodes:
				nodesStr += f"\n{node.getNodeAsStr()}"
				if (i < len(nodes)-1):
					nodesStr += f","
					i += 1
				else:
					nodesStr += f"\n"
					break
		nodesStr += f"]"
		return nodesStr

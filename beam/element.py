import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from beam.node import Node
from Utils.crossSectionProperty import CrossSectionProperty

class Element():

	def __init__(self, local_num, global_num):
		self.local_num = local_num
		self.global_num = global_num
		self.start_node = None
		self.end_node = None
		self.cross_section_property = None

	def getElementAsStr(self):
		start_node_str = self.start_node.getNodeAsStr()
		end_node_str = self.end_node.getNodeAsStr()
		csp_str = self.cross_section_property.getPropertyAsStr()
		return f"[Local Num: {self.local_num}, Global Num: {self.global_num},\nStart Node: {start_node_str},\nEnd Node: {end_node_str},\nCross Section Property: {csp_str}\n]"
	
	def setStartAndEndNodes(self, start_node, end_node):
		self.start_node = start_node
		self.end_node = end_node
	
	def setCrossSectionProperty(self, cross_section_property):
		self.cross_section_property = cross_section_property

	@staticmethod
	def getElementsAsStr(elements):
		elemsStr = f"["
		if elements != None:
			i = 0
			for elem in elements:
				elemsStr += f"\n{elem.getElementAsStr()}" 
				if (i < len(elements)-1):
					elemsStr += f","
					i += 1
				else:
					elemsStr += f"\n"
					break
		elemsStr += f"]"
		return elemsStr

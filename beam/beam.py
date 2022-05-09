import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from beam.segment import Segment
from beam.node import Node
from beam.element import Element
from Utils.geometry import Geometry

class Beam():

	def __init__(self, beam_id, beam_type, n_segments, n_elements):
		self.beam_id = beam_id
		self.beam_type = beam_type
		self.n_segments = n_segments
		self.n_elements = n_elements
		self.segments = None
		self.nodes = None
		self.elements = None
		self.start_point = None
		self.end_point = None
		self.beam_length = None
		self.beam_class = None

	def getBeamAsStr(self):
		segmentsStr = "[]"
		nodesStr = "[]"
		elementsStr = "[]"
		if self.segments != None:
			segmentsStr = Segment.getSegmentsAsStr(self.segments)
		if self.nodes != None:
			nodesStr = Node.getNodesAsStr(self.nodes)
		if self.elements != None:
			elementsStr = Element.getElementsAsStr(self.elements)
		return f"[Beam ID: {int(float(self.beam_id))}, Beam Type: {self.beam_type}, Beam Class: {self.getBeamClass()}, Beam Length: {float(self.getBeamLength()):.3f}, No. of segments: {self.n_segments}, No. of elements: {self.n_elements},\nSegments: {segmentsStr},\nNodes: {nodesStr},\nElements: {elementsStr}\n]"

	@staticmethod
	def getBeamsAsStr(beams):
		beamsStr = f"["
		if beams != None:
			i = 0
			for beam in beams:
				beamsStr += f"\n{beam.getBeamAsStr()}"
				if (i < len(beams)-1):
					beamsStr += f","
					i += 1
				else:
					beamsStr += f"\n"
					break
		beamsStr += f"]"
		return beamsStr

	def setSegmentbyDic(self,dics):
		segments =[]
		for dic in dics:
			val =[]
			for _,v in dic.items():
				val.append(v)
			segments.append(tuple(val))
		if len(segments) != self.n_segments:
			self.segments = None
			return False
		# segments are ordered from top to bottom w.r.t their positions
		self.segments = segments
		return True

	def setSegments(self, segments):
		"""
			:param segments: list of Segment objects
		"""
		if len(segments) != self.n_segments:
			self.segments = None
			return False
		# segments are ordered from top to bottom w.r.t their positions
		self.segments = segments
		return True

	def getBeamLength(self):
		if self.start_point == None or self.end_point == None:
			return -1.0
		self.beam_length = Geometry.findPointsDistance(self.start_point, self.end_point)
		return self.beam_length

	def getBeamClass(self):
		if self.beam_class == None:
			return ""
		return self.beam_class

	def setStartAndEndPoints(self, start_point, end_point):
		self.start_point = start_point
		self.end_point = end_point

	def setNodes(self, nodes):
		self.nodes = nodes

	def setElements(self, elements):
		self.elements = elements

	def setBeamClass(self, beam_class):
		self.beam_class = beam_class

	def generateSegmentsStartAndEndPoints(self):
		if self.segments != None:
			curr_point = self.start_point
			for segment in self.segments:
				
				
				segment_length = float(segment.length_ratio) * float(self.getBeamLength())
				# print(segment,", segment.length_ratio =",segment.length_ratio,", self.getBeamLength() = ",self.getBeamLength())
				
				dir_vector = Geometry.findDirectionVector(curr_point, self.end_point)
				next_point = Geometry.findPointOnVector(curr_point, dir_vector, segment_length)
				segment.setStartAndEndPoints(curr_point, next_point)
				
				# print("stand and end point = ",segment.start_point, segment.end_point)
				curr_point = next_point

	@staticmethod
	def getBeamTypes():
		return ["Stand", "Compartment"]

	# Returns the beam classes and their respective descriptions
	@staticmethod
	def getBeamClasses(n_comps):
		beam_classes = list()

		# Append stand and compartment classes together into a single list
		stand_beam_classes = Beam.getStandBeamClasses(n_comps)
		comp_beam_classes = Beam.getCompBeamClasses(n_comps)

		beam_classes = stand_beam_classes + comp_beam_classes
		return beam_classes
	
	# Returns the stand beam classes and their respective descriptions
	@staticmethod
	def getStandBeamClasses(n_comps):
		"""
			It gives the beam class name and its description for stands
			This particularly is genrated for only stand
			>>> n_comps = 4
			>>> beam_classes = Beam.getStandBeamClasses(n_comps)
			>>> [("SB","Stand (Bottom)"), ("S1","Stand 1"), ("S2","Stand 2"), ("S3","Stand 3"), ("S4","Stand 4"), ("ST","Stand (Top)")]
		"""
		stand_beam_classes = list()

		# Append stand classes
		stand_beam_classes.append((f"SB", f"Stand (Bottom)"))
		for i in range(n_comps):
			curr_class_stand = f"S{i+1}"
			curr_descr_stand = f"Stand {i+1}"
			stand_beam_classes.append((curr_class_stand, curr_descr_stand))
		stand_beam_classes.append((f"ST", f"Stand (Top)"))

		return stand_beam_classes
	
	# Returns the compartment beam classes and their respective descriptions
	@staticmethod
	def getCompBeamClasses(n_comps):
		"""
			It gives the beam class name and its description for compartments
			This particularly is genrated for only compartments
			>>> n_comps = 2
			>>> beam_classes = Beam.getCompBeamClasses(n_comps)
			>>> [("C1L","Comp 1 (Lower)"),("C1U","Comp 1 (Upper)"), ("C2L","Comp 2 (Lower)"),("C2U","Comp 2 (Upper)")]
		"""
		
		comp_beam_classes = list()

		# Append compartment classes
		for i in range(n_comps):
			curr_class_lower_comp = f"C{i+1}L"
			curr_descr_lower_comp = f"Comp {i+1} (Lower)"

			curr_class_upper_comp = f"C{i+1}U"
			curr_descr_upper_comp = f"Comp {i+1} (Upper)"

			comp_beam_classes.append((curr_class_lower_comp, curr_descr_lower_comp))
			comp_beam_classes.append((curr_class_upper_comp, curr_descr_upper_comp))

		return comp_beam_classes

	@staticmethod
	def getInfoTableHeaders():
		return ["Beam Class", "Description", "No. of Segments", "No. of Elements"]

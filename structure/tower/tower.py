import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")/ Path("..")
sys.path.append(str(path_main))
from beam.beam import Beam
from Utils.geometry import Geometry
from beam.node import Node
from beam.element import Element
from Utils.crossSectionProperty import CrossSectionProperty
from Utils.utilities import Utilities
import traceback

class Tower():

	def __init__(self, height, tower_id=1):
		"""
		init is the Default Python Constructor. Height and Tower_id are method arguments passed.
		'self' refers to the current object
		
		height                  : Refers to height of tower
		tower_id                : Refers to unique identification of tower
		beams                   : Refers to 
		start_point             : Refers to initial coordinate of tower
		end_point               : Refers to final coordinate of tower
		n_stands                : Refers to legs of tower
		coordinates             : Refers to a specific point in the 3D coordinate system
		line_end_points         : 
		cross_section_properties: Refers to properties of bay in tower
		"""
		self.tower_id = tower_id
		self.height = height
		self.beams = None
		self.start_point = None
		self.end_point = None

		self.n_stands = 1
		self.coordinates = None
		self.line_end_points = None
		self.cross_section_properties = None
	
	def setBeams(self, beams):
		self.beams = beams
	
	def setStartAndEndPoints(self, start_point, end_point):
		self.start_point = start_point
		self.end_point = end_point
	
	def getTowerAsStr(self):
		"""
		beamStr is initialized. It will take the value of 'beam' from class beam if it has some value
		other than 'None'.

		Return: float and list data types
		"""
		beamStr = "[]"
		if self.beams != None:
			beamStr = Beam.getBeamsAsStr(self.beams)
		return f"[Tower ID: {int(float(self.tower_id))}, Height: {float(self.height):.3f},\nBeams: {beamStr}\n]"
	
	def generateBeamInputData(self):
		try:
			# Initialization
			beam = self.beams[0]

			n_elements = beam.n_elements 
			# print(f"Tower File: n_elements={n_elements}")
			n_nodes = n_elements + 1
			node_global_num = 0
			elem_global_num = 0
			nodes = list()
			elements = list()

			# Set beam class
			beam.setBeamClass("S")
			# Tower length/height (same value for respective monopile beam too)
			height = self.height
			# Origin: Meeting point of Tower (top) and Support Structure (bottom)
			origin = [0, 0, 0]
			# Fixing Start point (bottom) in 3d-coordinate system with respect to origin
			start_point = origin
			# Fixing End point in positive z-axis. Value given by tower height
			end_point = [start_point[0], start_point[1], height]
			# Setting start and end points of tower
			self.setStartAndEndPoints(start_point, end_point)

			# Set Start and End Points of Beam
			beam.setStartAndEndPoints(start_point, end_point)
			# Generate Start And End Points for Monopile Beam Segments
			beam.generateSegmentsStartAndEndPoints()

			# Get coordinates from start and end points (and no. of elements)
			coordinates = Geometry.getCoordinates(start_point, end_point, n_elements)
			line_end_points = [[start_point, end_point]]

			# Store the lists of coordinates and line end points for later graph display
			self.coordinates = coordinates
			self.line_end_points = line_end_points

			# Generate Cross Section Properties
			self.generateCrossSectionProperties()

			# Generate nodes
			for node_local_num in range(0, n_nodes):
				curr_coordinates = coordinates[node_local_num]
				node_global_num += 1

				node = Node(node_local_num + 1, node_global_num, curr_coordinates)
				nodes.append(node)

			# Generate elements
			for elem_local_num in range(0, n_elements):
				elem_global_num += 1

				# Create Element
				elem = Element(elem_local_num + 1, elem_global_num)
				# Set Start And End Nodes for each Element
				start_node = nodes[elem_local_num]
				end_node = nodes[elem_local_num+1]
				elem.setStartAndEndNodes(start_node, end_node)
				# Set Cross-Section Property to each Element
				elem.setCrossSectionProperty(self.cross_section_properties[elem_local_num])
				# Append Element to list of elements
				elements.append(elem)

			# Set nodes and elements to monopile beam
			beam.setNodes(nodes)
			beam.setElements(elements)

			# Return True if process successful (i.e. if no exceptions occur)
			return True
		except Exception:
			print(str(traceback.print_exc()))	# DEBUG_TEST
			return False

	# Creates cross section properties for the Tower beam
	def generateCrossSectionProperties(self):
		cross_section_properties = list()	# initializing the properties list
		prop_id = 1	# the property id which is to be incremented

		# Creating cross section properties for the Tower beam class
		beam = self.beams[0]
		beam_length = beam.getBeamLength()
		n_elems = beam.n_elements
		segments = beam.segments
		# print("Tower.py:",beam_length,n_elems,len(segments))

		# Get list of cross section properties
		cross_section_properties = CrossSectionProperty.getCrossSectionProperties(prop_id, 
					beam_length, n_elems, segments)

		# Store the list containing the cross section properties
		self.cross_section_properties = cross_section_properties
	
	# Gets list of all cross section properties sorted by prop id 
	def getCrossSectionProperties(self):
		cross_section_properties = list()
		for property in self.cross_section_properties:
			cross_section_properties.append(property)
		# Sort all cross section properties by prop id
		cross_section_properties.sort(key=lambda x: x.prop_id)
		# Return sorted list of cross section properties
		return cross_section_properties

	def checkLengthRatiosValidity(self):
		# Check if the length ratios sum upto 1 for the segments of the tower beam
		length_ratio_sum = 0
		beam = self.beams[0]
		for segment in beam.segments:
			length_ratio_sum += float(segment.length_ratio)
		length_ratio_valid = (length_ratio_sum > 0.999) and (length_ratio_sum <= 1)
		if not length_ratio_valid:
			return False
		# Return True if all length ratios valid
		return True
	
	def writeBeamInput(self):
		# Data to be written
		data = f""

		data += f"!! tower beam input for NREL 15 MW reference wind turbine (based on HAWC2 input files)\n!! \n"
		data += f"!! number of beams (1), number of cross-section properties (2)\n"

		# Calculate total number of cross section properties
		num_csps = len(self.cross_section_properties)

		data += f"1 {num_csps}\n!! \n!! \n"

		# Gets corresponding list of coordinates
		beam = self.beams[0]
		nodes = beam.nodes
		coordinates = list()
		for node in nodes:
			coordinates.append(node.coordinates)

		data += f"!! beam tower: number of nodes (1), number of elements (2)\n"
		data += f"{len(coordinates)} {len(coordinates)-1}\n!! \n!! \n"

		# Fill node data
		data += f"!! beam tower: nodal phi (1, 2, 3), nodal d1 (4, 5, 6), nodal d2 (7, 8, 9), nodal d3 (10, 11, 12)\n"
		# Iterates through the coordinates for the given beam and gets them as lines
		for coordinate in coordinates:
			data += Utilities.getCoordsAsLine(coordinate)
		data += f"!! \n!! \n"

		# Fill elements data
		data += f"!! beam tower: connectivities (1, 2), cross-section property (3)\n"
		for elem in beam.elements:
			data += f"{elem.start_node.local_num} {elem.end_node.local_num} {elem.cross_section_property.prop_id}\n"
		data += "!! \n!! \n"

		# Iterates through all cross section properties to extract their data
		cross_section_properties = self.getCrossSectionProperties()
		for i in range(len(cross_section_properties)):
			csp = cross_section_properties[i]
			data += f"!!cross-section property #{csp.prop_id} (tower element {csp.prop_id})\n"
			data += csp.getPropertiesAsLine()
			if i < len(cross_section_properties)-1:	# add footer, if not last element
				data += "!! \n!! \n"

		# Write Tower Beam Input File
		Utilities.writeBeamInput(data, filename="tower_beaminput.txt")
	
	def writeLogFile(self):
		# Data to be written
		data = f"==================\nTOWER LOG:\n"
		data += f"\n==================\nGENERAL INFORMATION:\n\nNo. of stands: {self.n_stands}\n"
		data += f"\n==================\nTOWER:\n\n{self.getTowerAsStr()}\n"
		data += f"\n=================="

		# Write Log File
		Utilities.writeLogFile(data, filename_start="tower-log-data")

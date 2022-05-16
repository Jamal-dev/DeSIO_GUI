import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")/ Path("..")
sys.path.append(str(path_main))
from beam.stand import Stand
from Utils.compartment import Compartment
from beam.beam import Beam
from beam.segment import Segment
from beam.node import Node
from beam.element import Element
from Utils.crossSectionProperty import CrossSectionProperty
from Utils.geometry import Geometry
from Utils.utilities import Utilities
import traceback

class Jacket():

	def __init__(self, n_stands, n_compartments, stand_length, stands_dist_above, stands_dist_below, gap_from_below=0):
		self.n_stands = n_stands	# can be 3 or 4 stands
		self.num_bays = n_compartments	# max 5
		self.stand_length = stand_length
		self.stands_Rhead = stands_dist_above
		self.stands_Rfoot = stands_dist_below
		self.LOSG = gap_from_below
		# TODO: why is this 4?
		self.n_beams_per_comp_side = 4
		self.stands = None
		self.compartments = None
		self.stand_beam_n_elems = None
		self.stand_beam_segments_settings = None
		self.stand_beam_lengths = list()
		self.upper_comp_beam_n_elems = None
		self.upper_comp_beam_segments_settings = None
		self.upper_comp_beam_lengths = list()
		self.lower_comp_beam_n_elems = None
		self.lower_comp_beam_segments_settings = None
		self.lower_comp_beam_lengths = list()
		self.all_coordinates = None
		self.all_line_end_points = None
		self.cross_section_properties = None
		self.node_global_num = 0
		self.elem_global_num = 0

	def setStands(self, stands):
		if len(stands) != self.n_stands:
			self.stands = None
			return False
		self.stands = stands
		return True

	def setCompartments(self, compartments):
		if len(compartments) != self.num_bays:
			self.compartments = None
			return False
		self.compartments = compartments
		return True

	def setCompBeamsData(self, stand_beam_n_elems, stand_beam_segments_settings, upper_comp_beam_n_elems, upper_comp_beam_segments_settings, lower_comp_beam_n_elems, lower_comp_beam_segments_settings):
		"""
			stand_beam_n_elems: list of ints, number of elements in each Stand
			self.stand_beam_segments_settings: list of list of Segments for each stand
			For example: SB: [Segment1, Segment2, Segment3], ST: [Segment1, Segment2, Segment3], S1: [Segment1, Segment2, Segment3]
			[ [Segment1, Segment2, Segment3], [Segment1, Segment2, Segment3], [Segment1, Segment2, Segment3] ]
			Here Segment1,... is the Segment object
		"""
		
		self.stand_beam_n_elems = stand_beam_n_elems
		self.stand_beam_segments_settings = stand_beam_segments_settings
		self.upper_comp_beam_n_elems = upper_comp_beam_n_elems
		self.upper_comp_beam_segments_settings = upper_comp_beam_segments_settings
		self.lower_comp_beam_n_elems = lower_comp_beam_n_elems
		self.lower_comp_beam_segments_settings = lower_comp_beam_segments_settings

	def generateBeamInputData(self):
		try:
			# Stores the beam and segment ids globally
			global curr_beam_id
			self.curr_beam_id = 1
			global curr_segment_id
			self.curr_segment_id = 1

			# Stores the lists of beam classes globally
			global beam_classes
			# beam_classes is the list of tuple which contains the class name and the class description
			self.beam_classes = Beam.getBeamClasses(self.num_bays)
			global stand_beam_classes
			# stand_beam_classes is the list of tuple which contains the class name and the class description for all stands SB, S1, S2, S3, S4,..., ST
			self.stand_beam_classes = Beam.getStandBeamClasses(self.num_bays)
			global comp_beam_classes
			# comp_beam_classes is the list of tuple which contains the class name and the class description for all compartments C1L, C1U, ..
			self.comp_beam_classes = Beam.getCompBeamClasses(self.num_bays)

			# print("*"*50)
			# print(self.beam_classes, self.stand_beam_classes, self.comp_beam_classes)
			# print("*"*50)

			# Generate Start And End Points for all Stands
			self.generateStandEdgePoints()
			### output is checked until here
			
			# Generate stand beam border points for all compartments
			self.generateStandCompInteractions()

			# Generate Beams and Segments for Stands
			self.generateStandBeams()

			# Generate Beams and Segments for Compartments
			self.generateCompBeams()
			
			# Generate Cross Section Properties
			self.generateCrossSectionProperties()
			
			# Generate Nodes and Elements for all Beams (for both Stands and Compartments)
			self.generateNodesAndElements()

			# Return True if process successful (i.e. if no exceptions occur)
			return True
		except Exception:
			print(str(traceback.print_exc()))	# DEBUG_TEST
			return False

	def generateStandEdgePoints(self):
		# Initial position and height settings
		origin = [0, 0, 0]
		jacket_height = self.findJacketHeight()
		
		# Find above and below radii for the circles encompassing the three stands 
		radius_above = Geometry.findRadius(self.n_stands, self.stands_Rhead)
		radius_below = Geometry.findRadius(self.n_stands, self.stands_Rfoot)
		
		# Return an invalid output if radii or height values negative
		if (radius_above < 0) or (radius_below < 0) or (jacket_height < 0): 
			return -1.0

		# TODO: change the origin to bottom so then above coordinates become positive
		# Determine the centers of the upper and lower circles
		center_above = origin
		center_below = [origin[0], origin[1], origin[2]-jacket_height]

		# Find equidistant points on circle
		beam_points_above = Geometry.findAdjacentCirclePoints(self.n_stands, center_above, radius_above)
		beam_points_below = Geometry.findAdjacentCirclePoints(self.n_stands, center_below, radius_below)

		# Iterate through all overall stands
		for stand_ind in range(self.n_stands):
			curr_stand = self.stands[stand_ind]
			curr_stand_start_point = beam_points_below[stand_ind]
			curr_stand_end_point = beam_points_above[stand_ind]
			# Set the start and end points for each stand
			curr_stand.setStartAndEndPoints(curr_stand_start_point, curr_stand_end_point)
			#print(f"Stand {stand_ind+1}: Start Point: {curr_stand.start_point:.2f}, End Point: {curr_stand.end_point:.2f}")

	def generateStandCompInteractions(self):
		# Find the jacket height
		jacket_height = self.findJacketHeight()

		# Iterate through all stands
		for stand_ind in range(self.n_stands):
			# Initialize list of stand-compartment interactions
			stand_comp_interactions = list()

			# Initial Settings
			stand = self.stands[stand_ind]
			start_point = stand.start_point
			end_point = stand.end_point
			dist = Geometry.findPointsDistance(start_point, end_point)
			jacket_dist_vector = Geometry.findDirectionVector(end_point, start_point)

			# Calculating gap from below points w.r.t origin (parallel to jacket height vector)
			curr_stand_comp_lower_pos = None
			curr_stand_comp_upper_pos = -jacket_height + self.LOSG

			# Iterate through all compartments and generate stand-compartment interactions
			for comp_ind in range(self.num_bays):
				# Setting the current lower and upper positions on the vertical jacket height
				curr_comp = self.compartments[comp_ind]
				curr_stand_comp_lower_pos = curr_stand_comp_upper_pos
				curr_stand_comp_upper_pos += curr_comp.height
				
				curr_lower_dist_ratio = curr_stand_comp_lower_pos/(-jacket_height)
				curr_upper_dist_ratio = curr_stand_comp_upper_pos/(-jacket_height)

				# Computing the lower and upper points on the stand with respect to the vertical length ratio
				stand_comp_lower_point = Geometry.findPointOnVectorFromRatio(end_point, jacket_dist_vector, curr_lower_dist_ratio)
				stand_comp_upper_point = Geometry.findPointOnVectorFromRatio(end_point, jacket_dist_vector, curr_upper_dist_ratio)

				# Computing the segment radii by which the lower and upper points are to be moved
				curr_lower_segments = self.lower_comp_beam_segments_settings[comp_ind]
				curr_upper_segments = self.upper_comp_beam_segments_settings[comp_ind]

				curr_lower_radius = (curr_lower_segments[0].diameter_start)/2.0
				curr_upper_radius = (curr_upper_segments[len(curr_upper_segments)-1].diameter_end)/2.0

				# Calculating the stand distance vectors and ratios corresponding to the segment radii
				lower_stand_dist_vector = Geometry.findDirectionVector(stand_comp_lower_point, end_point)
				upper_stand_dist_vector = Geometry.findDirectionVector(stand_comp_upper_point, start_point)

				lower_stand_dist_ratio = curr_lower_radius/Geometry.findPointsDistance(stand_comp_lower_point, end_point)
				upper_stand_dist_ratio = curr_upper_radius/Geometry.findPointsDistance(stand_comp_upper_point, start_point)

				# Finding the new points on the stand after adding the segment radii to the current positions
				stand_comp_lower_point = Geometry.findPointOnVectorFromRatio(stand_comp_lower_point, lower_stand_dist_vector, lower_stand_dist_ratio)
				stand_comp_upper_point = Geometry.findPointOnVectorFromRatio(stand_comp_upper_point, upper_stand_dist_vector, upper_stand_dist_ratio)

				# Updating the current upper reference points on the vector perpendicular to the jacket height
				new_upper_stand_dist_ratio = curr_upper_radius/Geometry.findPointsDistance(stand_comp_upper_point, end_point)
				new_upper_stand_dist_vector = Geometry.findDirectionVector(stand_comp_upper_point, end_point)
				updated_stand_comp_upper_point = Geometry.findPointOnVectorFromRatio(stand_comp_upper_point, new_upper_stand_dist_vector, new_upper_stand_dist_ratio)
				updated_upper_dist_ratio = Geometry.findPointsDistance(end_point, updated_stand_comp_upper_point)/dist
				curr_stand_comp_upper_pos = -jacket_height*updated_upper_dist_ratio

				# Add the found points to the list of stand-compartment interactions
				stand_comp_interactions.append({"lower": stand_comp_lower_point, "upper": stand_comp_upper_point})

			# Sets the Stand-Compartment Interactions for each stand
			stand.setStandCompInteractions(stand_comp_interactions)
	
	def generateStandBeams(self):
		# Generate Beams and Segments for Stands
		for stand_ind in range(self.n_stands):
			# Initialize list of all beams for current stand
			curr_stand_beams = list()

			# Current stand and corresponding stand-compartment interactions
			curr_stand = self.stands[stand_ind]
			curr_interactions = curr_stand.stand_comp_interactions

			# Generate Stand Beam (bottom-most portion)
			curr_beam = self.generateStandBeam(0)
			# Generate Start and End points for current stand beam and its corresponding segments
			curr_beam.setStartAndEndPoints(curr_stand.start_point, curr_interactions[0]["lower"])
			# Appends and store the beam length of current stand beam
			self.stand_beam_lengths.append(curr_beam.getBeamLength())
			curr_beam.generateSegmentsStartAndEndPoints()
			# Append current stand beam to list
			curr_stand_beams.append(curr_beam)

			# Iterate through all compartments to generate the stand beams corresponding to their interactions
			for comp_ind in range(self.num_bays):
				# Generate Stand Beam for current compartment (middle portion)
				curr_beam = self.generateStandBeam(comp_ind+1)
				# Generate Start and End points for current stand beam and its corresponding segments
				curr_beam.setStartAndEndPoints(curr_interactions[comp_ind]["lower"], curr_interactions[comp_ind]["upper"])
				# Appends and store the beam length of current stand beam
				self.stand_beam_lengths.append(curr_beam.getBeamLength())
				curr_beam.generateSegmentsStartAndEndPoints()
				# Append current stand beam to list
				curr_stand_beams.append(curr_beam)
			
			# Generate Stand Beam (top-most portion)
			curr_beam = self.generateStandBeam(len(self.stand_beam_classes)-1)
			# Generate Start and End points for current stand beam and its corresponding segments
			curr_beam.setStartAndEndPoints(curr_interactions[len(curr_interactions)-1]["upper"], curr_stand.end_point)
			# Appends and store the beam length of current stand beam
			self.stand_beam_lengths.append(curr_beam.getBeamLength())
			curr_beam.generateSegmentsStartAndEndPoints()
			# Append current stand beam to list
			curr_stand_beams.append(curr_beam)
			
			# Set list of generated Beams to Stand
			curr_stand.setBeams(curr_stand_beams)

	def generateStandBeam(self, beam_class_ind):
		# Get segments and elements settings corresponding to each stand portion
		n_segments = len(self.stand_beam_segments_settings[beam_class_ind])
		n_elems = self.stand_beam_n_elems[beam_class_ind]
		segments = self.stand_beam_segments_settings[beam_class_ind]

		# Create Beam (Stand portion) and set its beam class
		beam = Beam(self.curr_beam_id, "Stand", n_segments, n_elems)
		beam.setBeamClass(self.stand_beam_classes[beam_class_ind][0])

		# Set Segments to Beam (Stand portion)
		new_segments = list()
		for segment_ind in range(len(segments)):
			# Update the segment ids of the current segments
			new_segments.append(Segment.getSegmentFromSettings(self.curr_segment_id, segments[segment_ind]))
			# Increment Segment ID
			self.curr_segment_id += 1
		# Add updated Segments to Beam
		beam.setSegments(new_segments)

		# Increment Beam ID
		self.curr_beam_id += 1

		# Return the generated beam (Stand portion)
		return beam
	
	def generateCompBeams(self):
		# Generate Beams and Segments for Compartments
		for comp_ind in range(self.num_bays):
			# Initial segments and elements settings
			curr_upper_n_segments = len(self.upper_comp_beam_segments_settings[comp_ind])
			curr_upper_n_elems = self.upper_comp_beam_n_elems[comp_ind]
			curr_upper_segments = self.upper_comp_beam_segments_settings[comp_ind]
			curr_lower_n_segments = len(self.lower_comp_beam_segments_settings[comp_ind])
			curr_lower_n_elems = self.lower_comp_beam_n_elems[comp_ind]
			curr_lower_segments = self.lower_comp_beam_segments_settings[comp_ind]

			# Initalize list of current compartment beams
			curr_comp_beams = list()

			# Iterate through all compartment sides
			for prev_stand_ind in range(self.n_stands):
				# prev_stand_ind is the stand index of the previous stand which is the long leg of the tower
				next_stand_ind = None
				if prev_stand_ind < (self.n_stands-1):
					next_stand_ind = prev_stand_ind + 1
				else:
					next_stand_ind = 0

				# Get previous and following stands to extract Compartment Borders for current side
				prev_stand = self.stands[prev_stand_ind]
				next_stand = self.stands[next_stand_ind]

				# Get stand beam border points for current compartment
				prev_stand_lower_point = prev_stand.stand_comp_interactions[comp_ind]["lower"]
				prev_stand_upper_point = prev_stand.stand_comp_interactions[comp_ind]["upper"]
				next_stand_lower_point = next_stand.stand_comp_interactions[comp_ind]["lower"]
				next_stand_upper_point = next_stand.stand_comp_interactions[comp_ind]["upper"]

				# Get direction vectors	between the lines connecting the two adjacent stand beams
				left_to_right_dir_vector = Geometry.findDirectionVector(prev_stand_upper_point, next_stand_upper_point)
				right_to_left_dir_vector = Geometry.findDirectionVector(next_stand_upper_point, prev_stand_upper_point)

				# Lower and Upper Stand Radii derived from the segment diameter at given point
				lower_stand_radius = Geometry.findStandRadiusAtPoint(prev_stand, prev_stand_lower_point)
				upper_stand_radius = Geometry.findStandRadiusAtPoint(prev_stand, prev_stand_upper_point)

				# Set radii as 0 if calculation invalid
				if lower_stand_radius == -1.0:
					lower_stand_radius = 0
				if upper_stand_radius == -1.0:
					upper_stand_radius = 0

				# Updates the border points (compartment beam edge points) to move by an offset value
				prev_stand_lower_point = Geometry.findPointOnVector(prev_stand_lower_point, left_to_right_dir_vector, lower_stand_radius)
				prev_stand_upper_point = Geometry.findPointOnVector(prev_stand_upper_point, left_to_right_dir_vector, upper_stand_radius)
				next_stand_lower_point = Geometry.findPointOnVector(next_stand_lower_point, right_to_left_dir_vector, lower_stand_radius)
				next_stand_upper_point = Geometry.findPointOnVector(next_stand_upper_point, right_to_left_dir_vector, upper_stand_radius)

				# Gets the intersection point of the two stand lines
				stand_intersection_point = Geometry.findThreeDimIntersection(prev_stand_lower_point, prev_stand_upper_point, next_stand_lower_point, next_stand_upper_point)

				# Define all Compartments Beams for current side
				prev_lower_beam = Beam(self.curr_beam_id, "Compartment", curr_lower_n_segments, curr_lower_n_elems)
				prev_upper_beam = Beam(self.curr_beam_id+1, "Compartment", curr_upper_n_segments, curr_upper_n_elems)
				next_lower_beam = Beam(self.curr_beam_id+2, "Compartment", curr_lower_n_segments, curr_lower_n_elems)
				next_upper_beam = Beam(self.curr_beam_id+3, "Compartment", curr_upper_n_segments, curr_upper_n_elems)

				# Setting the beam classes of the respective beams
				prev_lower_beam.setBeamClass(f"C{comp_ind+1}L")
				prev_upper_beam.setBeamClass(f"C{comp_ind+1}U")
				next_lower_beam.setBeamClass(f"C{comp_ind+1}L")
				next_upper_beam.setBeamClass(f"C{comp_ind+1}U")

				# Set Start And End Points of the respective Beams
				prev_lower_beam.setStartAndEndPoints(prev_stand_lower_point, stand_intersection_point)
				prev_upper_beam.setStartAndEndPoints(stand_intersection_point, prev_stand_upper_point)
				next_lower_beam.setStartAndEndPoints(next_stand_lower_point, stand_intersection_point)
				next_upper_beam.setStartAndEndPoints(stand_intersection_point, next_stand_upper_point)

				# Stores the lengths of all lower and upper compartment beams for each compartment
				if prev_stand_ind == 0:
					self.lower_comp_beam_lengths.append(prev_lower_beam.getBeamLength())
					self.upper_comp_beam_lengths.append(prev_upper_beam.getBeamLength())

				# Create Segments for Prev Lower Beam
				prev_lower_segments = list()
				for segment_ind in range(len(curr_lower_segments)):
					# Update the segment ids of the current segments
					prev_lower_segments.append(Segment.getSegmentFromSettings(self.curr_segment_id, curr_lower_segments[segment_ind]))
					# Increment Segment ID
					self.curr_segment_id += 1
				# Create Segments for Prev Upper Beam
				prev_upper_segments = list()
				for segment_ind in range(len(curr_upper_segments)):
					# Update the segment ids of the current segments
					prev_upper_segments.append(Segment.getSegmentFromSettings(self.curr_segment_id, curr_upper_segments[segment_ind]))
					# Increment Segment ID
					self.curr_segment_id += 1
				# Create Segments for Next Lower Beam
				next_lower_segments = list()
				for segment_ind in range(len(curr_lower_segments)):
					# Update the segment ids of the current segments
					next_lower_segments.append(Segment.getSegmentFromSettings(self.curr_segment_id, curr_lower_segments[segment_ind]))
					# Increment Segment ID
					self.curr_segment_id += 1
				# Create Segments for Next Upper Beam
				next_upper_segments = list()
				for segment_ind in range(len(curr_upper_segments)):
					# Update the segment ids of the current segments
					next_upper_segments.append(Segment.getSegmentFromSettings(self.curr_segment_id, curr_upper_segments[segment_ind]))
					# Increment Segment ID
					self.curr_segment_id += 1

				# Add updated Segments to respective Beams
				prev_lower_beam.setSegments(prev_lower_segments)
				prev_upper_beam.setSegments(prev_upper_segments)
				next_lower_beam.setSegments(next_lower_segments)
				next_upper_beam.setSegments(next_upper_segments)
				# Generate Start And End Points of all Beam Segments
				prev_lower_beam.generateSegmentsStartAndEndPoints()
				prev_upper_beam.generateSegmentsStartAndEndPoints()
				next_lower_beam.generateSegmentsStartAndEndPoints()
				next_upper_beam.generateSegmentsStartAndEndPoints()
				# Append Beams to respective lists of Compartment Beams
				curr_comp_beams.append(prev_lower_beam)
				curr_comp_beams.append(prev_upper_beam)
				curr_comp_beams.append(next_lower_beam)
				curr_comp_beams.append(next_upper_beam)

				self.curr_beam_id += 4

			# Set Beams to Compartment
			self.compartments[comp_ind].setBeams(curr_comp_beams)

	def generateCrossSectionProperties(self):
		try:
			cross_section_properties = dict()	# initializing the properties dictionary
			prop_id = 1	# the property id which is to be incremented

			# Generate Cross-section properties for all stand beam classes
			for beam_ind in range(len(self.stand_beam_classes)):
				# Stand beam settings
				beam_class = self.stand_beam_classes[beam_ind][0]
				beam_length = self.stand_beam_lengths[beam_ind]
				n_elems = self.stand_beam_n_elems[beam_ind]
				segments = self.stand_beam_segments_settings[beam_ind]

				# Get list of cross section properties
				cross_section_properties[beam_class] = CrossSectionProperty.getCrossSectionProperties(prop_id, 
							beam_length, n_elems, segments)
				# Increment property id by the size of the current cross section properties list
				prop_id += len(cross_section_properties[beam_class])

			# Iterating through compartments
			for comp_ind in range(self.num_bays):
				# Creating cross section properties for the Lower Compartment beam class
				beam_class = f"C{comp_ind+1}L"
				beam_length = self.lower_comp_beam_lengths[comp_ind]
				n_elems = self.lower_comp_beam_n_elems[comp_ind]
				segments = self.lower_comp_beam_segments_settings[comp_ind]

				# Get list of cross section properties
				cross_section_properties[beam_class] = CrossSectionProperty.getCrossSectionProperties(prop_id, 
							beam_length, n_elems, segments)
				# Increment property id by the size of the current cross section properties list
				prop_id += len(cross_section_properties[beam_class])

				# Creating cross section properties for the Upper Compartment beam class
				beam_class = f"C{comp_ind+1}U"
				beam_length = self.upper_comp_beam_lengths[comp_ind]
				n_elems = self.upper_comp_beam_n_elems[comp_ind]
				segments = self.upper_comp_beam_segments_settings[comp_ind]

				# Get list of cross section properties
				cross_section_properties[beam_class] = CrossSectionProperty.getCrossSectionProperties(prop_id, 
							beam_length, n_elems, segments)
				# Increment property id by the size of the current cross section properties list
				prop_id += len(cross_section_properties[beam_class])

			# Store the dictionary containing the cross section properties
			self.cross_section_properties = cross_section_properties
		except Exception:
			print(str(traceback.print_exc()))	# DEBUG_TEST

	def generateNodesAndElements(self):
		coordinates = list()
		line_end_points = list()

		# Create Nodes And Elements for all Stand Beams
		for stand_ind in range(self.n_stands):
			stand = self.stands[stand_ind]

			# Iterate through all Stand beams
			for beam in stand.beams:
				# Stand beam settings
				start_point = beam.start_point
				end_point = beam.end_point
				n_elements = beam.n_elements 
				n_nodes = n_elements + 1
				beam_class = beam.getBeamClass()
				properties = self.cross_section_properties[beam_class]

				# Initialize Lists for Nodes and Elements
				nodes = list()
				elements = list()

				# Calculate and retrieve coordinates for each stand beam
				curr_coordinates = Geometry.getCoordinates(start_point, end_point, n_elements)
				line_end_points.append([start_point, end_point])

				# Generate nodes
				for node_local_num in range(0, n_nodes):
					coordinate = curr_coordinates[node_local_num]
					self.node_global_num += 1

					node = Node(node_local_num + 1, self.node_global_num, coordinate)
					nodes.append(node)

				# Generate elements
				for elem_local_num in range(0, n_elements):
					self.elem_global_num += 1

					# Create Element
					elem = Element(elem_local_num + 1, self.elem_global_num)
					# Set Start And End Nodes for each Element
					start_node = nodes[elem_local_num]
					end_node = nodes[elem_local_num+1]
					elem.setStartAndEndNodes(start_node, end_node)
					# Set Cross-Section Property to each Element
					elem.setCrossSectionProperty(properties[elem_local_num])
					# Append Element to list of elements
					elements.append(elem)

				# Set nodes and elements to jacket stand beams
				beam.setNodes(nodes)
				beam.setElements(elements)

				# Concatenate coordinates
				coordinates += curr_coordinates

		# Create Nodes And Elements for all Compartment Beams
		for comp_ind in range(self.num_bays):
			comp = self.compartments[comp_ind]
			beams = comp.beams

			# Iterate through all Beams
			for beam in beams:
				start_point = beam.start_point
				end_point = beam.end_point
				n_elements = beam.n_elements 
				n_nodes = n_elements + 1
				beam_class = beam.getBeamClass()
				properties = self.cross_section_properties[beam_class]

				# Initialize Lists for Nodes and Elements
				nodes = list()
				elements = list()

				# Calculate and retrieve coordinates for each compartment beam
				curr_coordinates = Geometry.getCoordinates(start_point, end_point, n_elements)
				line_end_points.append([start_point, end_point])

				# Generate nodes for each beam
				for node_local_num in range(0, n_nodes):
					coordinate = curr_coordinates[node_local_num]
					self.node_global_num += 1

					node = Node(node_local_num + 1, self.node_global_num, coordinate)
					nodes.append(node)

				# Generate elements for each beam
				for elem_local_num in range(0, n_elements):
					self.elem_global_num += 1

					# Create Element
					elem = Element(elem_local_num + 1, self.elem_global_num)
					# Set Start And End Nodes for each Element
					start_node = nodes[elem_local_num]
					end_node = nodes[elem_local_num+1]
					elem.setStartAndEndNodes(start_node, end_node)
					# Set Cross-Section Property to each Element
					elem.setCrossSectionProperty(properties[elem_local_num])
					# Append Element to list of elements
					elements.append(elem)

				# Set nodes and elements to current jacket compartment beam
				beam.setNodes(nodes)
				beam.setElements(elements)

				# Concatenate coordinates
				coordinates += curr_coordinates

		# Store the lists of all coordinates and line end points in the Jacket object
		self.all_coordinates = coordinates
		self.all_line_end_points = line_end_points

	def findJacketHeight(self):
		# Calculates the perpendicular height of the jacket (i.e. parallel to Z-Axis)
		# TODO: again dividing it by 2.0. There shouldn't be 2.0 if it was radius
		tilt = (float(self.stands_Rfoot) - float(self.stands_Rhead))/2.0
		if tilt < 0: 
			return -1.0
		jacket_height = (float(self.stand_length**2) - float(tilt**2))**0.5
		return jacket_height

	# Check if the sum of compartment heights and gap from below is compatible with jacket height
	def checkHeightsCompatibility(self):
		jacketHeight = float(self.findJacketHeight()) # calculated jacket height
		gapBelow = float(self.LOSG)	# minimum gap from below, from which the compartments begin
		compHeights = 0.0	# sum of all compartment heights
		for comp_id in range(self.num_bays):
			comp = self.compartments[comp_id]
			compHeights += float(comp.height)

		if (compHeights + gapBelow) <= jacketHeight:
			return True
		return False

	def checkLengthRatiosValidity(self):
		# Check if the length ratios sum upto 1 for all segment settings
		for stand_segment_setting in self.stand_beam_segments_settings:
			length_ratio_sum = 0
			for stand_segment in stand_segment_setting:
				length_ratio_sum += float(stand_segment.length_ratio)
			length_ratio_valid = (length_ratio_sum > 0.999) and (length_ratio_sum <= 1.001)
			if not length_ratio_valid:
				return False
		for lower_comp_segment_setting in self.lower_comp_beam_segments_settings:
			length_ratio_sum = 0
			for lower_comp_segment in lower_comp_segment_setting:
				length_ratio_sum += float(lower_comp_segment.length_ratio)
			length_ratio_valid = (length_ratio_sum > 0.999) and (length_ratio_sum <= 1.001)
			if not length_ratio_valid:
				return False
		for upper_comp_segment_setting in self.upper_comp_beam_segments_settings:
			length_ratio_sum = 0
			for upper_comp_segment in upper_comp_segment_setting:
				length_ratio_sum += float(upper_comp_segment.length_ratio)
			length_ratio_valid = (length_ratio_sum > 0.999) and (length_ratio_sum <= 1.001)
			if not length_ratio_valid:
				return False
		# Return True if all length ratios valid
		return True

	def getNumBeams(self):
		try:
			num_beams = int((self.n_stands * (self.num_bays + 2)) + (self.num_bays * self.n_stands * Jacket.getNumBeamsPerComp()))
			return num_beams
		except Exception:
			return -1
	
	# Gets list of all beams sorted by beam id 
	def getBeams(self):
		beams = list()
		for stand_ind in range(self.n_stands):
			beams += self.stands[stand_ind].beams
		for comp_ind in range(self.num_bays):
			beams += self.compartments[comp_ind].beams
		beams.sort(key=lambda x: x.beam_id)
		return beams

	# Finds beam by beam id from list of all beams
	def findBeamById(self, beam_id):
		beams = self.getBeams()
		for beam in beams:
			if str(beam.beam_id) == str(beam_id):
				return beam
		return None

	# Gets a dictionary of all coordinates (value) mapped to their respective jacket beams (key)
	def getBeamCoordinates(self):
		beam_coord_dict = dict()

		# Get list of all beams sorted by beam id 
		beams = self.getBeams()

		# Get list of coordinates for each beam and add them to the dictionary
		for beam in beams:
			curr_coord = list()
			for node in beam.nodes:
				curr_coord.append(node.coordinates)
			beam_coord_dict[str(beam.beam_id)] = curr_coord

		return beam_coord_dict
	
	# Gets list of all cross section properties sorted by prop id 
	def getCrossSectionProperties(self):
		cross_section_properties = list()
		for beam_class in self.cross_section_properties:
			for property in self.cross_section_properties[beam_class]:
				cross_section_properties.append(property)
		# Sort all cross section properties by prop id
		cross_section_properties.sort(key=lambda x: x.prop_id)
		# Return sorted list of cross section properties
		return cross_section_properties

	def writeBeamInput(self):
		# Data to be written
		data = f""

		data += f"!! beam input for NREL 15 MW reference wind turbine (based on HAWC2 input files)\n!! \n"
		data += f"!! number of beams (1), number of cross-section properties (2)\n"
		
		# Calculate total number of jacket beams - If invalid, raise Exception
		num_beams = self.getNumBeams()
		if num_beams == -1:
			raise Exception(f"Error: Total number of beams invalid.")

		# Calculate total number of cross section properties
		num_csps = sum([len(self.cross_section_properties[beam_class]) for beam_class in self.cross_section_properties])
		
		data += f"{num_beams} {num_csps}\n!! \n!! \n"

		# Gets beam coordinates dictionary
		beam_coord_dict = self.getBeamCoordinates()

		# Iterates the beam coordinates dictionary and gets the coordinates as lines for each beam
		for beam_id in beam_coord_dict:
			# Gets beam coordinates
			coordinates = beam_coord_dict[beam_id]
			# Finds beam class corresponding to current beam id
			beam_class = ""
			beam = self.findBeamById(beam_id)
			if beam != None:
				beam_class = beam.getBeamClass()
			beam_name = f"jacket beam {beam_id} ({beam_class})"
			# Fill node data
			data += f"!! {beam_name}: number of nodes (1), number of elements (2)\n"
			data += f"{len(coordinates)} {len(coordinates)-1}\n!! \n!! \n"
			data += f"!! {beam_name}: nodal phi (1, 2, 3), nodal d1 (4, 5, 6), nodal d2 (7, 8, 9), nodal d3 (10, 11, 12)\n"
			# Iterates through the coordinates for the given beam and gets them as lines
			for coordinate in coordinates:
				data += Utilities.getCoordsAsLine(coordinate)
			data += f"!! \n!! \n"

			# Fill elements data
			data += f"!! {beam_name}: connectivities (1, 2), cross-section property (3)\n"
			for elem in beam.elements:
				data += f"{elem.start_node.local_num} {elem.end_node.local_num} {elem.cross_section_property.prop_id}\n"
			data += "!! \n!! \n"
		
		# Iterates through all cross section properties to extract their data
		cross_section_properties = self.getCrossSectionProperties()
		for i in range(len(cross_section_properties)):
			csp = cross_section_properties[i]
			data += f"!!cross-section property #{csp.prop_id} (jacket element {csp.prop_id})\n"
			data += csp.getPropertiesAsLine()
			if i < len(cross_section_properties)-1:	# add footer, if not last element
				data += "!! \n!! \n"

		# Write Beam Input File
		Utilities.writeBeamInput(data)

	def writeLogFile(self):
		# Data to be written
		data = f"==================\nJACKET LOG:\n"
		# General Information
		data += f"\n==================\nGENERAL INFORMATION:\n"
		data += f"\nNo. of stands: {self.n_stands}"
		data += f"\nNo. of compartments: {self.num_bays}"
		data += f"\nStand length: {self.stand_length}"
		data += f"\nDistance above: {self.stands_Rhead}"
		data += f"\nDistance below: {self.stands_Rfoot}"
		data += f"\nGap from below: {self.LOSG}\n"
		# Stands and Compartments
		data += f"\n==================\nSTANDS:\n\n{Stand.getStandsAsStr(self.stands)}\n"
		data += f"\n==================\nCOMPARTMENTS:\n\n{Compartment.getCompartmentsAsStr(self.compartments)}\n"
		data += f"\n=================="
		
		# Write Log File
		Utilities.writeLogFile(data)

	@staticmethod
	def getNumBeamsPerComp():
		return 4

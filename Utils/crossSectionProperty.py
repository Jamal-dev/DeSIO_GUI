import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from Utils.geometry import Geometry
import math
import traceback

class CrossSectionProperty():

	def __init__(self, prop_id):
		self.prop_id = prop_id
		self.properties = []
	
	def setProperties(self, properties):
		self.properties = properties
	
	def getPropertyAsStr(self):
		return f"[Property ID: {self.prop_id}, Properties: {str(self.properties)}]"

	def getPropertiesAsLine(self):
		line = f""
		if self.properties != [] and self.properties != None:
			for i in range(len(self.properties)):
				for j in range(len(self.properties[i])):
					curr_item = str(self.properties[i][j])
					line += f"{float(curr_item):.3f}"	# string formatting
					if j < len(self.properties[i])-1: # ending if not last item
						line += f" "
					else:	# ending if last item
						line += f"\n"
						break
		return line
	
	@staticmethod
	def getCrossSectionProperties(prop_id, beam_length, n_elems, segments):
		try:
			# List of cross section properties to be output
			cross_section_properties = list()

			# Select first segment
			curr_segment_id = 0
			current_segment = segments[curr_segment_id]
			curr_segments_length = current_segment.length_ratio*beam_length
			segments_length_sum = curr_segments_length

			# Select element length settings
			elem_length = beam_length/n_elems
			elem_length_sum = 0.0

			for elem_ind in range(n_elems):
				# Creating Cross Section Property
				csp = CrossSectionProperty(prop_id)
				# Element diameter same for all elements of a beam
				elem_diameter = 0.0
				if (elem_length_sum + elem_length) <= segments_length_sum:	# Element in current segment
					# Set start diameter of segment at element start position
					curr_diameter_start = Geometry.findDiameterAtSegmentPosition(current_segment.diameter_start,
							current_segment.diameter_end, elem_length_sum/curr_segments_length)
					# Increment current sum of element lengths
					elem_length_sum += elem_length
					# Set end diameter of segment at element end position
					curr_diameter_end = Geometry.findDiameterAtSegmentPosition(current_segment.diameter_start,
							current_segment.diameter_end, elem_length_sum/curr_segments_length)
				elif elem_length_sum < segments_length_sum:	# Element in between two segments
					# Getting the first and second length parts of the element across the two segments
					elem_first_part_length = segments_length_sum - elem_length_sum
					elem_second_part_length = elem_length - elem_first_part_length

					# Set start diameter of segment at element start position
					curr_diameter_start = Geometry.findDiameterAtSegmentPosition(current_segment.diameter_start,
							current_segment.diameter_end, elem_length_sum/curr_segments_length)
					# Increment current sum of element lengths
					elem_length_sum += elem_length

					# Select next segment
					if curr_segment_id < len(segments)-1:	# if next segment available
						curr_segment_id += 1
						current_segment = segments[curr_segment_id]
						curr_segments_length = current_segment.length_ratio*beam_length
						segments_length_sum += curr_segments_length

						# Set end diameter of segment at element end position
						curr_diameter_end = Geometry.findDiameterAtSegmentPosition(current_segment.diameter_start,
								current_segment.diameter_end, elem_second_part_length/curr_segments_length)
					else:	# else if next segment not available
						# Set end diameter of segment at element end position
						curr_diameter_end = current_segment.diameter_end
				else:	# Element in next segment (when current elements- and segments-sum match)
					if curr_segment_id < len(segments)-1:	# if next segment available
						# Select next segment
						curr_segment_id += 1
						current_segment = segments[curr_segment_id]
						curr_segments_length = current_segment.length_ratio*beam_length
						segments_length_sum += curr_segments_length

						# Set start diameter of segment at element start position
						curr_diameter_start = current_segment.diameter_start
						# Increment current sum of element lengths
						elem_length_sum += elem_length
						# Set end diameter of segment at element end position
						curr_diameter_end = Geometry.findDiameterAtSegmentPosition(current_segment.diameter_start,
								current_segment.diameter_end, elem_length_sum/curr_segments_length)
					else:	# else if next segment not available, break loop
						break
				
				# Current element diameters are obtained from the center of element (avg of start and end diameters)
				elem_diameter = (curr_diameter_start + curr_diameter_end)/2.0

				# Calculate inner diameters of the element at the start and end positions
				curr_inner_diameter_start = curr_diameter_start - 2*current_segment.thickness_start
				curr_inner_diameter_end = curr_diameter_end - 2*current_segment.thickness_end

				# Inner element diameter from the center of element (avg of start and end inner diameters)
				elem_inner_diameter = (curr_inner_diameter_start + curr_inner_diameter_end)/2.0

				# Element thickness at center
				elem_thickness = (elem_diameter - elem_inner_diameter)/2.0

				# Create current property from acquired parameter values
				curr_elem_property = [elem_length, elem_diameter, elem_diameter, elem_thickness, 
						elem_thickness, current_segment.density, current_segment.e, 
						current_segment.g, current_segment.alpha_s, current_segment.alpha_v,
						 current_segment.scf_start, current_segment.scf_end]
				# Generate the cross section property from element data and append it to the list of properties
				csp.setProperties(CrossSectionProperty.deriveCrossSectionProperty(curr_elem_property))
				cross_section_properties.append(csp)
				# Increment cross section property ID for each element
				prop_id += 1

			# Return the list of cross section properties
			return cross_section_properties
		except Exception:
			print(str(traceback.print_exc()))	# DEBUG
	
	@staticmethod
	def deriveCrossSectionProperty(elem_prop):
		try:
			rho = elem_prop[5]	# Mass Density
			E = elem_prop[6]	# Young's Modulus of Elasticity
			G = elem_prop[7]	# Shear Modulus
			alpha_s = elem_prop[8] if (elem_prop[8] != None) else 0	# Stress induced dissipation
			alpha_v = elem_prop[9] if (elem_prop[9] != None) else 0	# Velocity induced dissipation

			outer_r = elem_prop[1]/2.0	# Outer radius
			inner_r = outer_r - elem_prop[3]	# Inner radius: outer radius - thickness (th.: start = end)

			# Formulae for the cylinder coordinates of a ring
			area = (math.pi)*(outer_r**2 - inner_r**2)
			i_1 = (math.pi/4)*(outer_r**4 - inner_r**4)
			i_2 = i_1
			s_1 = 0
			s_2 = 0
			i_12 = 0

			# Matrices
			c_gg = [
				[2*G*area, 0, 0],
				[0, 2*G*area, 0],
				[0,	0, E*area]
			]

			c_kk = [
				[E*i_1, -E*i_12, 0],
				[-E*i_12, E*i_2, 0],
				[0,	0, 2*G*i_1 + 2*G*i_2]
			]

			c_gk = [
				[0, 0, -2*G*s_1],
				[0, 0, 2*G*s_2],
				[E*s_1, -E*s_2,	0]
			]

			c_kg = [
				[0, 0, E*s_1],
				[0, 0, -E*s_2],
				[-2*G*s_1, 2*G*s_2,	0]
			]

			m = [
				[rho*area, rho*s_2, rho*s_1],
				[rho*s_2, rho*i_2, rho*i_12],
				[rho*s_1, rho*i_12, rho*i_1]
			]

			# Property shown as output
			updated_property = [
				[c_gg[0][0], c_gg[1][1], c_gg[2][2], c_gg[1][2], c_gg[0][2], c_gg[0][1]],
				[c_kk[0][0], c_kk[1][1], c_kk[2][2], c_kk[1][2], c_kk[0][2], c_kk[0][1]],
				[c_gk[0][0], c_gk[1][1], c_gk[2][2], c_gk[1][2], c_gk[0][2], c_gk[0][1]],
				[c_kg[0][0], c_kg[1][1], c_kg[2][2], c_kg[1][2], c_kg[0][2], c_kg[0][1]],
				[m[0][0], m[1][1], m[2][2], m[1][2], m[0][2], m[0][1]],
				[alpha_s, alpha_v]
			]

			# Value correction: Convert all negative zeroes to positive zeroes
			for line_ind in range(len(updated_property)):
				for val_ind in range(len(updated_property[line_ind])):
					if float(updated_property[line_ind][val_ind]) == -0.0:
						updated_property[line_ind][val_ind] = 0.0

			# Return the updated cross section property
			return updated_property
		except Exception:
			print(str(traceback.print_exc()))	# DEBUG

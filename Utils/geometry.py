import math
import numpy as np
import matplotlib.pyplot as plt

class Geometry():

	def __init__(self):
		pass

	@staticmethod
	def getMagnitude(vector):
		return math.sqrt(float(vector[0])**2 + float(vector[1])**2 + float(vector[2])**2)

	@staticmethod
	def normalizeVector(vector):
		try:
			magnitude = Geometry.getMagnitude(vector)
			norm_vector = [float(vector[0])/magnitude, float(vector[1])/magnitude, float(vector[2])/magnitude]
			return norm_vector
		except Exception:
			return [0, 0, 0]

	@staticmethod
	def findOrthonormalBases(dir_vector):
		try:
			# The z-base is the direction vector of the beam
			z_base = [float(dir_vector[0]), float(dir_vector[1]), float(dir_vector[2])]
			z_base = Geometry.normalizeVector(z_base)

			# Find x-base orthogonal to z-base
			x_base = list()
			if z_base[0] == 0 and z_base[1] == 0:
				x_base = [1, 0, 0]
			else:
				z_base_complement = [z_base[0], z_base[1], 0]
				z_base_complement = Geometry.normalizeVector(z_base_complement)
				x_base = np.cross(z_base, z_base_complement)
			x_base = Geometry.normalizeVector(x_base)

			# Find y-base from the cross-product of the z- and x-bases (right hand rule)
			y_base = np.cross(z_base, x_base)
			y_base = Geometry.normalizeVector(y_base)

			# Return orthonormal bases as a tuple
			return (x_base, y_base, z_base)
		except Exception:
			return ([0, 0, 0], [0, 0, 0], [0, 0, 0])

	@staticmethod
	def findRadius(n_points, dist_bw_points):
		if n_points == 3:	# Jacket 3 Stand
			return (float(dist_bw_points)/(3.0**0.5))
		elif n_points == 4:	# Jacket 4 Stand
			return (float(dist_bw_points)/(2.0**0.5))
		else:
			return -1.0
	
	@staticmethod
	def findStandRadiusAtPoint(stand, point):
		stand_radius = None

		# TODO
		# Distance of the point from the start point of the stand
		stand_point_length = Geometry.findPointsDistance(stand.start_point, point)

		# Iterates through all stand beam portions to find the corresponding segment radius
		curr_beam_length_sum = 0.0
		for beam in stand.beams:
			# Current Stand beam settings
			beam_length = beam.getBeamLength()
			segments = beam.segments
			start_point = beam.start_point

			curr_beam_length_sum += beam_length
			if stand_point_length == curr_beam_length_sum:
				# The last segment diameter end of current beam is selected if point is exactly on current beam end
				found_segment = segments[len(segments)-1]
				found_diameter = found_segment.diameter_end
				stand_radius = found_diameter/2.0
				break
			elif stand_point_length < curr_beam_length_sum:
				# Iterates through the segments of current beam to find corresponding diameter
				beam_point_length = Geometry.findPointsDistance(start_point, point)
				prev_length_sum = 0.0
				temp_length_sum = segments[0].length_ratio*beam_length
				for i in range(len(segments)):
					if temp_length_sum >= beam_point_length:
						found_segment = segments[i]

						# Finds the diameter at the given position on the found segment 
						pos_ratio = (beam_point_length - prev_length_sum)/(found_segment.length_ratio*beam_length)
						found_diameter = Geometry.findDiameterAtSegmentPosition(found_segment.diameter_start,
							found_segment.diameter_end, pos_ratio)

						stand_radius = found_diameter/2.0
						break
					if i < len(segments)-1:
						prev_length_sum = temp_length_sum
						temp_length_sum += segments[i+1].length_ratio*beam_length
			
			# Returns stand radius if found
			if stand_radius != None:
				return stand_radius
		return -1.0
	
	@staticmethod
	def findDiameterAtSegmentPosition(diameter_start, diameter_end, pos_ratio):
		diameter_pos = -1.0
		if diameter_start != None and diameter_end != None:
			if diameter_start == diameter_end:
				diameter_pos = diameter_start
			elif diameter_start > diameter_end:
				diameter_pos = diameter_start - ((diameter_start - diameter_end) * pos_ratio)
			else:
				diameter_pos = diameter_start + ((diameter_end - diameter_start) * pos_ratio)
		return diameter_pos

	@staticmethod
	def findPointsDistance(start_point, end_point):
		dir_vector = Geometry.findDirectionVector(start_point, end_point)
		return Geometry.getMagnitude(dir_vector)

	@staticmethod
	def findPointOnVectorFromRatio(start_point, dir_vector, dist_ratio):
		return [start_point[0] + dist_ratio*dir_vector[0], start_point[1] + dist_ratio*dir_vector[1], 
			start_point[2] + dist_ratio*dir_vector[2]]
	
	@staticmethod
	def findPointOnVector(start_point, dir_vector, dist):
		dir_vector = Geometry.normalizeVector(dir_vector)
		return [start_point[0] + dist*dir_vector[0], start_point[1] + dist*dir_vector[1], 
			start_point[2] + dist*dir_vector[2]]

	@staticmethod
	def findDirectionVector(start_point, end_point):
		return [end_point[0]-start_point[0], end_point[1]-start_point[1], end_point[2]-start_point[2]]

	@staticmethod
	def findThreeDimIntersection(prev_lower, prev_upper, next_lower, next_upper):
		line1 = Geometry.findDirectionVector(prev_upper, next_lower)
		line2 = Geometry.findDirectionVector(next_upper, prev_lower)

		x_line1 = [next_lower[0], line1[0]]
		y_line1 = [next_lower[1], line1[1]]
		z_line1 = [next_lower[2], line1[2]]

		x_line2 = [prev_lower[0], line2[0]]
		y_line2 = [prev_lower[1], line2[1]]
		z_line2 = [prev_lower[2], line2[2]]

		eq_xy_vars = [[x_line1[1], -x_line2[1]], [y_line1[1], -y_line2[1]]] 
		eq_xy_vals = [(x_line2[0] - x_line1[0]), (y_line2[0] - y_line1[0])]

		eq_xy_coeffs = np.linalg.solve(np.array(eq_xy_vars), np.array(eq_xy_vals))

		#eq_yz_vars = [[y_line1[1], -y_line2[1]], [z_line1[1], -z_line2[1]]] 
		#eq_yz_vals = [(y_line2[0] - y_line1[0]), (z_line2[0] - z_line1[0])]
		#eq_yz_coeffs = np.linalg.solve(np.array(eq_yz_vars), np.array(eq_yz_vals))

		x = x_line1[0] + (eq_xy_coeffs[0]*x_line1[1])
		y = y_line1[0] + (eq_xy_coeffs[0]*y_line1[1])
		z = z_line1[0] + (eq_xy_coeffs[0]*z_line1[1])

		intersection = [x, y, z]
		return intersection

	@staticmethod
	def findAdjacentCirclePoints(n_points, center, radius):
		# Finds adjacent equidistant circle points
		adjPoints = list()
		for i in range(n_points):
			point = [None, None, center[2]]

			point[0] = center[0] + radius*math.cos(2*math.pi*i/n_points)
			point[1] = center[1] + radius*math.sin(2*math.pi*i/n_points)

			adjPoints.append(point)
		return adjPoints

	@staticmethod
	def getCoordinates(start_point, end_point, n_elements):
		# Determine direction vector from start and end points
		dir_vector = Geometry.findDirectionVector(start_point, end_point)
		# Find orthonormal bases
		orthonormal_bases = Geometry.findOrthonormalBases(dir_vector)

		# Determine Increment from Direction Vector
		x_increment = (dir_vector[0])/n_elements
		y_increment = (dir_vector[1])/n_elements
		z_increment = (dir_vector[2])/n_elements
		n_nodes = n_elements + 1	# No. of nodes

		# Calculate the coordinates and enter them into a list
		coordinates = list()
		x = start_point[0]
		y = start_point[1]
		z = start_point[2]
		i = 0
		while True:
			i += 1
			if i == n_nodes:
				# Append point vector and orthonormal bases to get coordinate
				coordinate = [end_point[0], end_point[1], end_point[2]]
				coordinate += orthonormal_bases[0] + orthonormal_bases[1] + orthonormal_bases[2]
				# Append coordinate to list of coordinates
				coordinates.append(coordinate)
				break
			else:
				# Append point vector and orthonormal bases  to get coordinate
				coordinate = [x, y, z]
				coordinate += orthonormal_bases[0] + orthonormal_bases[1] + orthonormal_bases[2]
				# Append coordinate to list of coordinates
				coordinates.append(coordinate)
				# Increment point values
				x += x_increment
				y += y_increment
				z += z_increment

		# Return the list of coordinates
		return coordinates

	@staticmethod
	def plotDataInAxes(ax, scatter_coordinates, line_end_points, windowTitle='3D-Simulation', orthoBaseLength=0.5, scatterEnabled=True, axesOn=True, orthonormalBaseOn=True):
		# Display scatter points
		if scatterEnabled:
			x_scatter, y_scatter, z_scatter = list(), list(), list()
			for coordinate in scatter_coordinates:
				x_scatter.append(coordinate[0])
				y_scatter.append(coordinate[1])
				z_scatter.append(coordinate[2])
			ax.scatter(x_scatter, y_scatter, z_scatter, c='b',s=50)
		
		# Display Lines connecting edge points
		for coord_pair in line_end_points:
			start_point = coord_pair[0]
			end_point = coord_pair[1]

			# Display 3D-line connecting start and end points
			ax.plot3D([start_point[0], end_point[0]], [start_point[1], end_point[1]], 
				[start_point[2], end_point[2]], 'blue')
		
		# Display the orthonormal bases
		if orthonormalBaseOn:
			for coordinate in scatter_coordinates:
				# Slice the coordinates into their constituents
				point = coordinate[:3]
				x_base = coordinate[3:6]
				y_base = coordinate[6:9]
				z_base = coordinate[9:]

				# Finds the orthonormal base points from the given point, direction vector and base length
				x_base_point = Geometry.findPointOnVector(point, x_base, orthoBaseLength)
				y_base_point = Geometry.findPointOnVector(point, y_base, orthoBaseLength)
				z_base_point = Geometry.findPointOnVector(point, z_base, orthoBaseLength)

				# Display the orthonormal bases as 3D-lines
				ax.plot3D([point[0], x_base_point[0]], [point[1], x_base_point[1]], 
					[point[2], x_base_point[2]], 'red')
				ax.plot3D([point[0], y_base_point[0]], [point[1], y_base_point[1]], 
					[point[2], y_base_point[2]], 'red')
				ax.plot3D([point[0], z_base_point[0]], [point[1], z_base_point[1]], 
					[point[2], z_base_point[2]], 'red')

		# Axes can be disabled
		if not axesOn:
			ax.set_axis_off()
	
	@staticmethod
	def displayPlot(scatter_coordinates, line_end_points, windowTitle='3D-Simulation', orthoBaseLength=0.5, scatterEnabled=True, axesOn=True, orthonormalBaseOn=True):
		# Create figure for displaying 3D graph
		fig = plt.figure()
		fig.canvas.set_window_title(windowTitle)
		# Plot 3D data in axes
		ax = plt.axes(projection='3d')
		Geometry.plotDataInAxes(ax, scatter_coordinates, line_end_points, windowTitle, orthoBaseLength, scatterEnabled, axesOn, orthonormalBaseOn)
		# Display final plot
		plt.show()

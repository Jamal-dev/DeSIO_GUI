import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities

class Segment():

	def __init__(self, segment_id, length_ratio="0.5", diameter_start="2", 
					diameter_end="2", thickness_start="0.1", 
					thickness_end="0.1", density="2", e="10", g="5", alpha_s="0.0", 
					alpha_v="0.0", 
					scf_start="0.0", scf_end="0.0"):
		self.segment_id = segment_id
		self.length_ratio = length_ratio
		self.diameter_start = diameter_start
		self.diameter_end = diameter_end
		self.thickness_start = thickness_start
		self.thickness_end = thickness_end
		self.density = density
		self.e = e
		self.g = g
		self.alpha_s = alpha_s
		self.alpha_v = alpha_v
		self.scf_start = scf_start
		self.scf_end = scf_end
		self.start_point = None
		self.end_point = None
	def __repr__(self) -> str:
		return f"Segment ID: {self.segment_id}"

	def getSegmentAsStr(self):
		start_point_str = "[]"
		end_point_str = "[]"
		if self.start_point != None:
			start_point_str = str(Utilities.getFormattedCoordList(self.start_point))
		if self.end_point != None:
			end_point_str = str(Utilities.getFormattedCoordList(self.end_point))
		return f"[Segment ID: {int(float(self.segment_id))} (Start Point: {start_point_str}, End Point: {end_point_str}), Length Ratio: {float(self.length_ratio):.3f}, Diameter Start: {float(self.diameter_start):.3f}, Diameter End: {float(self.diameter_end):.3f}, Thickness Start: {float(self.thickness_start):.3f}, Thickness End: {float(self.thickness_end):.3f}, Density: {float(self.density):.3f}, E: {float(self.e):.3f}, G: {float(self.g):.3f},  alpha_s: {float(self.alpha_s):.3f}, alpha_v: {float(self.alpha_v):.3f}, SCF Start: {float(self.scf_start):.3f}, SCF End: {float(self.scf_end):.3f}]"

	def setStartAndEndPoints(self, start_point, end_point):
		self.start_point = start_point
		self.end_point = end_point

	

	def convertFields2numeric(self):
		self.segment_id = int(self.segment_id)
		self.length_ratio = self.c2f(self.length_ratio)
		self.diameter_start = self.c2f(self.diameter_start)
		self.diameter_end = self.c2f(self.diameter_end)
		self.thickness_start = self.c2f(self.thickness_start)
		self.thickness_end = self.c2f(self.thickness_end)
		self.density = self.c2f(self.density)
		self.e = self.c2f(self.e)
		self.g = self.c2f(self.g)
		self.alpha_s = self.c2f(self.alpha_s)
		self.alpha_v = self.c2f(self.alpha_v)
		self.scf_start = self.c2f(self.scf_start)
		self.scf_end = self.c2f(self.scf_end)


	def c2f(self,field:str):
		# it converts empty str to 0.0 or other numeric value to float
		if not field:
			return 0.0
		else:
			return float(field)

	@classmethod
	def fromRow(cls, row):
		return cls(int(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6]), float(row[7]), float(row[8]), float(row[9]), float(row[10]), float(row[11]), float(row[12]))

	@staticmethod
	def getSegmentFromSettings(segment_id, segment):
		return Segment(segment_id, segment.length_ratio, segment.diameter_start, segment.diameter_end, segment.thickness_start, segment.thickness_end, segment.density, segment.e, segment.g, segment.alpha_s, segment.alpha_v, segment.scf_start, segment.scf_end)

	@staticmethod
	def getSegmentsAsStr(segments):
		segmentsStr = f"["
		if segments != None:
			i = 0
			for segment in segments:
				segmentsStr += f"\n{segment.getSegmentAsStr()}"
				if (i < len(segments)-1):
					segmentsStr += f","
					i += 1
				else:
					segmentsStr += f"\n"
					break
		segmentsStr += f"]"
		return segmentsStr

	@staticmethod
	def getTableHeaders():
		return ["Segment", "Length Ratio", "Diameter-Start", "Diameter-End", "Thickness-Start", "Thickness-End", "Density", "E", "G", "alpha_s", "alpha_v", "SCF-Start", "SCF-End"]


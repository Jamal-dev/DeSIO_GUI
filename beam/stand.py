import sys
from pathlib import Path
import os
file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir)/ Path("..")
sys.path.append(str(path_main))
from beam.beam import Beam
from Utils.geometry import Geometry

class Stand():

	def __init__(self, stand_id, height):
		self.stand_id = stand_id
		self.height = height
		self.beams = None
		self.stand_comp_interactions = None
		self.start_point = None
		self.end_point = None

	def getStandAsStr(self):
		beamStr = "[]"
		if self.beams != None:
			beamStr = Beam.getBeamsAsStr(self.beams)
		return f"[Stand ID: {int(float(self.stand_id))}, Height: {float(self.height):.3f},\nBeams: {beamStr}\n]"

	@staticmethod
	def getStandsAsStr(stands):
		standsStr = f"["
		if stands != None:
			i = 0
			for stand in stands:
				standsStr += f"\n{stand.getStandAsStr()}"
				if (i < len(stands)-1):
					standsStr += f","
					i += 1
				else:
					standsStr += f"\n"
					break
		standsStr += f"]"
		return standsStr
	
	def getStandLength(self):
		if (self.start_point != None) and (self.end_point != None):
			return Geometry.findPointsDistance(self.start_point, self.end_point)
		return -1.0
	
	def setBeams(self, beams):
		self.beams = beams

	def setStandCompInteractions(self, stand_comp_interactions):
		self.stand_comp_interactions = stand_comp_interactions
	
	def setStartAndEndPoints(self, start_point, end_point):
		self.start_point = start_point
		self.end_point = end_point

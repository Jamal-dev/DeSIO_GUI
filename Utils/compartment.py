from beam.beam import Beam

class Compartment():

	def __init__(self, compartment_id, height):
		self.compartment_id = compartment_id
		self.height = height
		self.beams = None

	def setBeams(self, beams):
		self.beams = beams

	def getCompartmentAsStr(self):
		beamsStr = "[]"
		if self.beams != None:
			beamsStr = Beam.getBeamsAsStr(self.beams)
		return f"[Compartment ID: {int(float(self.compartment_id))}, Height: {float(self.height):.3f},\nBeams: {beamsStr}\n]"

	@staticmethod
	def getCompartmentsAsStr(comps):
		compsStr = f"["
		if comps != None:
			i = 0
			for comp in comps:
				compsStr += f"\n{comp.getCompartmentAsStr()}"
				if (i < len(comps)-1):
					compsStr += f","
					i += 1
				else:
					compsStr += f"\n"
					break
		compsStr += f"]"
		return compsStr

	@staticmethod
	def getTableHeaders():
		return ["ID", "Height"]

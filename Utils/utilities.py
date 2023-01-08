import sys
from pathlib import Path
import os
cur_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(cur_dir)/ Path("..")
sys.path.append(str(path_main))
import time
from PyQt5.QtWidgets import  QMessageBox
from PyQt5.QtWidgets import QDialog, QWidget
from segment_table import Ui_Dialog
# from beam.segment import Segment

class Utilities():

	def __init__(self):
		pass

	# check if string is empty
	@staticmethod
	def isEmpty(strVal):
		if strVal == '':
			return True
		else:
			return False

	# Checks if a string is an integer or not
	@staticmethod
	def isInt(strVal):
		try:
			int(strVal)
			return True
		except ValueError:
			return False

	# Checks if a string is a float or not
	@staticmethod
	def isFloat(strVal):
		try:
			float(strVal)
			return True
		except ValueError:
			return False
	
	# Checks if a string is a positive number or not
	@staticmethod
	def isPos(val):
		try:
			if Utilities.isFloat(val):
				val = float(val)
				if val >= 0:
					return True
			elif Utilities.isInt(val):
				val = int(val)
				if val >= 0:
					return True
			else:
				return False
		except ValueError:
			return False
	# Checks if a string is a positive number  and greater than 0
	@staticmethod
	def isPosG0(val):
		try:
			if Utilities.isFloat(val):
				val = float(val)
				if val > 0:
					return True
			elif Utilities.isInt(val):
				val = int(val)
				if val > 0:
					return True
			else:
				return False
		except ValueError:
			return False
	# Checks if a string is a positive int  and greater than 0
	@staticmethod
	def isPosG0andInt(val):
		try:
			if Utilities.isInt(val):
				val = int(val)
				if val > 0:
					return True
			else:
				return False
		except ValueError:
			return False
	# Checks if a string is a positive float  and greater than 0
	@staticmethod
	def isPosG0andFloat(val):
		try:
			if Utilities.isFloat(val):
				val = float(val)
				if val > 0:
					return True
			else:
				return False
		except ValueError:
			return False
	# Checks if a string is a positive float  and greater or equal to 0
	@staticmethod
	def isPosandFloat(val):
		try:
			if Utilities.isFloat(val):
				val = float(val)
				if val >= 0:
					return True
			else:
				return False
		except ValueError:
			return False
	
	# Convert empty values to zero
	@staticmethod
	def convertEmpty2zero(val):
		for k,v in val.items():
			if not v:
				# empty entry
				val[k] = 0.0
		return val
	
	# Convert dictonary to float
	@staticmethod
	def convertDic2float(val):
		for k,v in val.items():
			val[k] = float(v)
		return val
	
	# To show Information message box
	@staticmethod
	def showInfoMsg(message,tit=None):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Information)
		msg.setWindowTitle(tit)
		msg.setStandardButtons(QMessageBox.Ok )
		msg.setText(message)	
		retval = msg.exec_()
	# To show Warning message box
	@staticmethod
	def showWarningMsg(message,tit=None):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Warning)
		msg.setWindowTitle(tit)
		msg.setStandardButtons(QMessageBox.Ok )
		msg.setText(message)	
		retval = msg.exec_()
	
	# To show error message box
	@staticmethod
	def showErrorMsg(field_name,suggestion=None):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Critical)
		message = field_name + " is not correct."
		msg.setWindowTitle("Check inputs")
		msg.setStandardButtons(QMessageBox.Ok )
		msg.setText(message)
		msg.setInformativeText(suggestion)
		retval = msg.exec_()
	# To show error message box
	@staticmethod
	def showErrorCustomMsg(message,suggestion=None):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Critical)
		msg.setWindowTitle("Check inputs")
		msg.setStandardButtons(QMessageBox.Ok )
		msg.setText(message)
		msg.setInformativeText(suggestion)
		retval = msg.exec_()

	# error message if field is empty
	@staticmethod
	def errorMsg_empty(field_value,field_name):
		if not field_value:
			Utilities.showErrorMsg(field_name,"It should not be empty")
			return False
		return True
	# error message if field is empty path
	@staticmethod
	def errorMsg_emptyPath(field_value,field_name):
		if not field_value:
			Utilities.showErrorMsg(field_name,"Path cannot be empty")
			return False
		if not os.path.exists(field_value):
			Utilities.showErrorMsg(field_name,"It should be a valid path")
			return False
		return True
	# error message if field is between 0 and 1 and float
	@staticmethod
	def errorMsg_0to1_float(field_value,field_name):
		if not field_value:
			field_value = 0.0
		if not Utilities.isFloat(field_value):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		if field_value < 0 or field_value > 1:
			Utilities.showErrorMsg(field_name,"It should be between 0 and 1")
			return False
		return True
	# error message if vector enteries are float or not
	@staticmethod
	def errorMsg_vectorFloat(field_value_1,field_value_2,field_value_3,field_name):
		if not field_value_1:
			field_value_1 = 0.0
		if not field_value_2:
			field_value_2 = 0.0
		if not field_value_3:
			field_value_3 = 0.0
		if not Utilities.isFloat(field_value_1):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		if not Utilities.isFloat(field_value_2):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		if not Utilities.isFloat(field_value_3):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		return True
	# error message if field is between 0 and 1 and the vector sum is 1
	@staticmethod
	def errorMsg_0to1_vectorSum1(field_value_1,field_value_2,field_value_3,field_name):
		if not field_value_1:
			field_value_1 = 0.0
		if not field_value_2:
			field_value_2 = 0.0
		if not field_value_3:
			field_value_3 = 0.0
		if not Utilities.isFloat(field_value_1):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		if not Utilities.isFloat(field_value_2):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		if not Utilities.isFloat(field_value_3):
			Utilities.showErrorMsg(field_name,"It should be a float")
			return False
		field_value_1 = float(field_value_1)
		field_value_2 = float(field_value_2)
		field_value_3 = float(field_value_3)
		if field_value_1 < 0 or field_value_1 > 1:
			Utilities.showErrorMsg(field_name,"It should be between 0 and 1")
			return False
		if field_value_2 < 0 or field_value_2 > 1:
			Utilities.showErrorMsg(field_name,"It should be between 0 and 1")
			return False
		if field_value_3 < 0 or field_value_3 > 1:
			Utilities.showErrorMsg(field_name,"It should be between 0 and 1")
			return False
		if field_value_1 + field_value_2 + field_value_3 != 1:
			Utilities.showErrorMsg(field_name,"The vector sum should be 1")
			return False
		return True

	# error message if int value grater than 0
	@staticmethod
	def errorMsg_greaterThan0_int(field_value,field_name):
		if not field_value:
				Utilities.showErrorMsg(field_name,"It should not be empty")
				return False
		if not Utilities.isPosG0andInt(field_value):
			Utilities.showErrorMsg(field_name,"It should be positive, greater than 0 and int")
			return False
		return True

    # error message if float value grater than or equal to 0
	@staticmethod
	def errorMsg_greaterOrequal0_float(field_value,field_name):
		if not field_value:
			field_value = 0.0
		if not Utilities.isPosandFloat(field_value):
			Utilities.showErrorMsg(field_name,"It should be positive, greater or equal to 0 and float")
			return False
		return True

    # error message if float value grater than 0
	@staticmethod
	def errorMsg_greaterthan0_float(field_value,field_name):
		if not field_value:
			field_value = 0.0
		if not Utilities.isPosG0andFloat(field_value):
			Utilities.showErrorMsg(field_name,"It should be positive, greater than 0 and float")
			return False
		return True
	
	@staticmethod
	def getZeroesList(num_zeroes):
		try:
			num_zeroes = int(num_zeroes)
			if num_zeroes <= 0:
				return list()
			return [0 for i in range(num_zeroes)]
		except Exception:
			return list()

	# Returns a list of strings as a single line string separated by spaces
	@staticmethod
	def getCoordsAsLine(coordItems):
		line = f""
		for i in range(len(coordItems)):
			curr_item = str(coordItems[i])
			line += f"{float(curr_item):8.6f}"	# string formatting
			if i < len(coordItems)-1: # ending if not last item
				line += f" "
			else:	# ending if last item
				line += f"\n"
				break
		return line

	# Converts a coordinate list to a list of formatted coordinate strings
	@staticmethod
	def getFormattedCoordList(coordItems):
		formattedList = list()
		for i in range(len(coordItems)):
			curr_item = str(coordItems[i])
			curr_item = f"{float(curr_item):.6f}"	# string formatting
			formattedList.append(curr_item)
		return formattedList

	# Returns the support structure modes and their respective image paths
	@staticmethod
	def getSupportStructureModes():
		return [
			(f"Monopile", f"support_struct_type_1_monopile.png"),
			(f"Jacket 3-Stand", f"support_struct_type_2_jacket3stand.png"),
			(f"Jacket 4-Stand", f"support_struct_type_3_jacket4stand.png")
		]

	@staticmethod
	def areSegmentsValid(sheet):
		# Get segments data from corresponding segments sheet
		sheet_data = sheet.get_sheet_data(return_copy = True, get_header = False, get_index = False)
		# Iterate through the sheet data to determine validity of the segments values
		for r in range(len(sheet_data)):
			for c in range(len(sheet_data[0])):
				# If any of the values can not be converted to a float value, then return False
				if not Utilities.isFloat(sheet_data[r][c]):
					return False
		# Return True if all values are valid
		return True

	@staticmethod
	def writeBeamInput(data, filename="beaminput.txt"):
		# Create an 'io' input-output folder if it does not already exist
		if not os.path.exists(f"io"):
			os.makedirs(f"io")

		# Define beam input file path
		file_path = path_main / Path(f"io/{filename}")

		# Write the data into the beam input file
		beam_input_file = open(file_path, "w")
		data = f"{data}"
		beam_input_file.write(data)

		# Close the file
		beam_input_file.close()

	@staticmethod
	def writeLogFile(data, filename_start="support-structure-log-data"):
		# Create an 'io\log' input-output folder if it does not exist
		path_log = path_main / Path("io/log")
		if not path_log.exists():
			os.makedirs(path_log)

		# Define log file path
		path_log_file = path_main / Path(f"io/log/{filename_start}_%Y.%m.%d_%H.%M.%S.txt")
		log_file_path = time.strftime(str(path_log_file))
		# Write the data into the log file
		log_file = open(log_file_path, "w")

		data = f"{data}"
		log_file.write(data)

		# Close the file
		log_file.close()

	@staticmethod
	def dispDialogSegmentTable(cur_segment,title:str,id:int):


		# current segment
		Dialog = QDialog()
		ui = Ui_Dialog()
		ui.setupUi(Dialog)
		# Set the labels txts
		ui.lblTitle.setText( title)
		ui.lineLengthRatio.setText(str(cur_segment.length_ratio))
		ui.lineDStar.setText(str(cur_segment.diameter_start))
		ui.lineDEnd.setText(str(cur_segment.diameter_end))
		ui.lineTStar.setText(str(cur_segment.thickness_start))
		ui.lineTend.setText(str(cur_segment.thickness_end))
		ui.lineRho.setText(str(cur_segment.density))
		ui.lineE.setText(str(cur_segment.e))
		ui.lineG.setText(str(cur_segment.g))
		ui.lineAlphaS.setText(str(cur_segment.alpha_s))
		ui.lineAlphaV.setText(str(cur_segment.alpha_v))
		ui.lineScfStart.setText(str(cur_segment.scf_start))
		ui.lineScfEnd.setText(str(cur_segment.scf_end))


		Dialog.show()
		resp = Dialog.exec_()

		values = {}
		values["id"] = int(id)
		
		if resp == QDialog.Accepted:
			values["length_ratio"] = ui.lineLengthRatio.text()
			values["Dstar"]= ui.lineDStar.text()
			values["Dend"]= ui.lineDEnd.text()
			values["Tstar"]= ui.lineTStar.text()
			values["Tend"]= ui.lineTend.text()
			values["rho"]= ui.lineRho.text()
			values["E"] = ui.lineE.text()
			values["G"] = ui.lineG.text()
			values["alpha_s"] = ui.lineAlphaS.text()
			values["alpha_v"] = ui.lineAlphaV.text()       
			values["scfStart"]= ui.lineScfStart.text()
			values["scfEnd"]= ui.lineScfEnd.text()

			cur_segment.segment_id=values["id"] 
			cur_segment.length_ratio=values["length_ratio"]
			cur_segment.diameter_start=values["Dstar"]
			cur_segment.diameter_end = values["Dend"] 
			cur_segment.thickness_start = values["Tstar"] 
			cur_segment.thickness_end = values["Tend"]
			cur_segment.density = values["rho"]
			cur_segment.e = values["E"]
			cur_segment.g = values["G"] 
			cur_segment.alpha_s = values["alpha_s"]
			cur_segment.alpha_v = values["alpha_v"] 
			cur_segment.scf_start = values["scfStart"] 
			cur_segment.scf_end = values["scfEnd"]

			if  Utilities.checkDialogInput(values):
				# everything is good now
				# values = Utilities.convertEmpty2zero(values)
				# values = Utilities.convertDic2float(values)
				cur_segment.convertFields2numeric()
				Dialog.close()
				return True
				
			else:
				# give use another chance for correct inputs
				Utilities.dispDialogSegmentTable(cur_segment,title,id)
				

		else:
			return False
			print("Cancel is pressed")
	@staticmethod
	def checkDialogInput(input_fields):
		if isinstance(input_fields,dict):
			values = input_fields
		else:
			values = {}
			# assuming input_fields is the Segment instance
			values["id"] = input_fields.segment_id 
			values["length_ratio"] = input_fields.length_ratio
			values["Dstar"] = input_fields.diameter_start
			values["Dend"]  = input_fields.diameter_end
			values["Tstar"]  = input_fields.thickness_start
			values["Tend"] = input_fields.thickness_end
			values["rho"] = input_fields.density
			values["E"] = input_fields.e
			values["G"]  = input_fields.g
			values["alpha_s"] = input_fields.alpha_s
			values["alpha_v"]  = input_fields.alpha_v
			values["scfStart"]  = input_fields.scf_start
			values["scfEnd"] = input_fields.scf_end
		# ToDO: Segment class imports gives circular import error;
		# otherwise it was a better solution than using else:
		# if isinstance(input_fields,Segment):
		# 	values = {}
		# 	values["id"] = input_fields.segment_id 
		# 	values["length_ratio"] = input_fields.length_ratio
		# 	values["Dstar"] = input_fields.diameter_start
		# 	values["Dend"]  = input_fields.diameter_end
		# 	values["Tstar"]  = input_fields.thickness_start
		# 	values["Tend"] = input_fields.thickness_end
		# 	values["rho"] = input_fields.density
		# 	values["E"] = input_fields.e
		# 	values["G"]  = input_fields.g
		# 	values["alpha_s"] = input_fields.alpha_s
		# 	values["alpha_v"]  = input_fields.alpha_v
		# 	values["scfStart"]  = input_fields.scf_start
		# 	values["scfEnd"] = input_fields.scf_end
		
		# length_ratio
		v = Utilities.errorMsg_greaterOrequal0_float(values["length_ratio"],"Length ratio")
		if not v:
			return False
		if not values["length_ratio"]:
			values["length_ratio"] = 0.0
		else:
			values["length_ratio"] = float(values["length_ratio"])
		if values["length_ratio"]>1:
			Utilities.showErrorMsg("Length ratio","It should not be greater than 1. The sum of all length ratios must be 1")
			return False
		# E
		v = Utilities.errorMsg_greaterthan0_float(values["E"],"E")
		if not v:
			return False
		
		# G
		v = Utilities.errorMsg_greaterthan0_float(values["G"],"G")
		if not v:
			return False
		
		# relationship between E and G
		# calculating poision ratio E = 2G( 1 + nu) 
		nu = float(values["E"])/(2*float(values["G"])) - 1
		if nu<=-1 or nu>=0.5:
			Utilities.showErrorCustomMsg("E and G values are not in correct range",f"Poison ratio is {nu} which is not between -1 and 0.5")
			return False

		# Dstar
		v = Utilities.errorMsg_greaterOrequal0_float(values["Dstar"],"D_start")
		if not v:
			return False
		
		# Dend
		v = Utilities.errorMsg_greaterOrequal0_float(values["Dend"],"D_end")
		if not v:
			return False

		# alpha_s
		v = Utilities.errorMsg_greaterOrequal0_float(values["alpha_s"],"alpha_s")
		if not values["alpha_s"]:
			values["alpha_s"] = 0.0
		if not v:
			return False
		elif float(values["alpha_s"])<0 or float(values["alpha_s"])>1:
			Utilities.showErrorMsg("alpha_s","It should be greater or equal to 0 and less or equal to 1")
			return False
		
		# alpha_v
		
		v = Utilities.errorMsg_greaterOrequal0_float(values["alpha_v"],"alpha_v")
		if not values["alpha_v"]:
			values["alpha_v"] = 0.0
		if not v:
			return False
		elif float(values["alpha_v"])<0 or float(values["alpha_v"])>1:
			Utilities.showErrorMsg("alpha_v","It should be greater or equal to 0 and less or equal to 1")
			return False

		# Tstar: thickness start
		v = Utilities.errorMsg_greaterOrequal0_float(values["Tstar"],"t_start")
		if not v:
			return False
		
		
		# Tend: thickness end
		v = Utilities.errorMsg_greaterOrequal0_float(values["Tend"],"t_end")
		if not v:
			return False

		# relationship between thickness and diameter
		rad_start = float(values["Dstar"])/2
		rad_end = float(values["Dend"])/2
		if rad_start<= float(values["Tstar"]):
			Utilities.showErrorCustomMsg("D_start/2 can not be less than or equal to thickness_start",f"This can make the area 0 or negative. Please change D_start or thickness_start")
			return False
		
		if rad_end<= float(values["Tend"]):
			Utilities.showErrorCustomMsg("D_end/2 can not be less than or equal to thickness_end",f"This can make the area 0 or negative. Please change D_end or thickness_end")
			return False
		
		# scfStart
		v = Utilities.errorMsg_greaterOrequal0_float(values["scfStart"],"SCF_start")
		if not v:
			return False
		
		# scfEnd
		v = Utilities.errorMsg_greaterOrequal0_float(values["scfEnd"],"SCF_end")
		if not v:
			return False

		# rho
		v = Utilities.errorMsg_greaterOrequal0_float(values["rho"],"Density")
		if not v:
			return False
		
		return True
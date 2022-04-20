import sys
from pathlib import Path
import os
cur_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(cur_dir)/ Path("..")
sys.path.append(str(path_main))
import time
from PyQt5.QtWidgets import  QMessageBox

class Utilities():

	def __init__(self):
		pass

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

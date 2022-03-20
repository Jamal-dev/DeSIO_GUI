import sys
from pathlib import Path
import os

file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir) / Path("..")
sys.path.append(str(path_main))
from Utils.utilities import Utilities as util
from structure.monopile.monopile import Monopile
from structure.jacket.jacket import Jacket
from beam.beam import Beam
from beam.segment import Segment
from beam.stand import Stand
from Utils.segments_userInterface import segments_ui
import traceback

from PyQt5.QtWidgets import QDialog
from segment_table import *


class MonopilePage:
    def __init__(self, parent=None):
        # ,stand_length,no_segments=None,no_elements=None,distance_above=None,distance_below=None,
        # self.tower_page_fields = {}
        # self.tower_page_fields["stand_length"] = stand_length
        # self.tower_page_fields["no_segments"] = no_segments
        # self.tower_page_fields["no_elements"] = no_elements
        # self.tower_page_fields["distance_above"] = distance_above
        # self.tower_page_fields["distance_below"] = distance_below
        self.ui = parent
        self.ui.lineStructureMono_StandLength.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureMono_NoOfSegments.textChanged.connect(self.disableGenbtn)
        self.ui.lineStructureMono_NoOfElements.textChanged.connect(self.disableGenbtn)

        # self.ui.lineStructureJ3_StandLength.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ3_NoOfCompartments.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ3_gapFromBelow.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ3_DistanceAbove.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ3_DistanceBelow.textChanged.connect(self.disableGenbtn)
        #
        # self.ui.lineStructureJ4_StandLength.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ4_NoOfCompartments.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ4_gapFromBelow.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ4_DistanceAbove.textChanged.connect(self.disableGenbtn)
        # self.ui.lineStructureJ4_DistanceBelow.textChanged.connect(self.disableGenbtn)

    def getValuesFromParent(self):
        self.mono_page_fields = {}
        # self.j3_page_fields = {}
        # self.j4_page_field = {}
        self.mono_page_fields["stand_length"] = self.ui.lineStructureMono_StandLength.text()
        self.mono_page_fields["no_segments"] = self.ui.lineStructureMono_NoOfSegments.text()
        self.mono_page_fields["no_elements"] = self.ui.lineStructureMono_NoOfElements.text()
        # self.tower_page_fields["distance_above"] = self.ui.lineStructureTower_DistanceAbove.text()
        # self.tower_page_fields["distance_below"] = self.ui.lineStructureTower_DistanceBelow.text()
        # Grab J3 text values
        # self.j3_page_fields["stand_length"] = self.ui.lineStructureJ3_StandLength.text()
        # self.j3_page_fields["no_of_compartments"] = self.ui.lineStructureJ3_NoOfCompartments.text()
        # self.j3_page_fields["gaps_from_below"] = self.ui.lineStructureJ3_gapFromBelow.text()
        # self.j3_page_fields["distance_above"] = self.ui.lineStructureJ3_DistanceAbove.text()
        # self.j3_page_fields["distance_below"] = self.ui.lineStructureJ3_DistanceBelow.text()
        # # Grab J4 text values
        # self.j4_page_field["stand_length"] = self.ui.lineStructureJ4_StandLength.text()
        # self.j4_page_field["no_of_compartments"] = self.ui.lineStructureJ4_NoOfCompartments.text()
        # self.j4_page_field["gaps_from_below"] = self.ui.lineStructureJ4_gapFromBelow.text()
        # self.j4_page_field["distance_above"] = self.ui.lineStructureJ4_DistanceAbove.text()
        # self.j4_page_field["distance_below"] = self.ui.lineStructureJ4_DistanceBelow.text()

    def disableGenbtn(self):
        self.ui.btnStructureTowerGGenrateFile_2.setEnabled(False)

    def checkMonoPageInputs(self):
        # stand length
        v = util.errorMsg_greaterThan0_int(self.mono_page_fields["stand_length"], "Stnad length")
        if not v:
            return False
        # no of segments
        v = util.errorMsg_greaterThan0_int(self.mono_page_fields["no_segments"], "Number of segments")
        if not v:
            return False

        # no of elements
        v = util.errorMsg_greaterThan0_int(self.mono_page_fields["no_elements"], "Number of elements")
        if not v:
            return False

        # # stand length J3
        # v = util.errorMsg_greaterThan0_int(self.j3_page_fields["stand_length"], "Stnad length")
        # if not v:
        #     return False
        # # no of segments
        # v = util.errorMsg_greaterThan0_int(self.j3_page_fields["no_of_compartments"], "Number of compartments")
        # if not v:
        #     return False
        #
        # # no of elements
        # v = util.errorMsg_greaterThan0_int(self.j3_page_fields["gaps_from_below"], "Gap from below")
        # if not v:
        #     return False
        #
        # # distance above
        # v = util.errorMsg_greaterOrequal0_float(self.j3_page_fields["distance_above"],"Distance above")
        # if not v:
        #     return False
        #
        # # # distance below
        # v = util.errorMsg_greaterOrequal0_float(self.j3_page_fields["distance_below"],"Distance below")
        # if not v:
        #     return False
        #
        # # stand length J4
        # v = util.errorMsg_greaterThan0_int(self.j4_page_field["stand_length"], "Stnad length")
        # if not v:
        #     return False
        # # no of segments
        # v = util.errorMsg_greaterThan0_int(self.j4_page_field["no_of_compartments"], "Number of compartments")
        # if not v:
        #     return False
        #
        # # no of elements
        # v = util.errorMsg_greaterThan0_int(self.j4_page_field["gaps_from_below"], "Gap from below")
        # if not v:
        #     return False
        #
        # # distance above
        # v = util.errorMsg_greaterOrequal0_float(self.j4_page_field["distance_above"], "Distance above")
        # if not v:
        #     return False
        #
        # # # distance below
        # v = util.errorMsg_greaterOrequal0_float(self.j4_page_field["distance_below"], "Distance below")
        # if not v:
        #     return False

        # convert tower fields
        self.mono_page_fields = util.convertEmpty2zero(self.mono_page_fields)
        self.mono_page_fields["stand_length"] = int(self.mono_page_fields["stand_length"])
        self.mono_page_fields["no_segments"] = int(self.mono_page_fields["no_segments"])
        self.mono_page_fields["no_elements"] = int(self.mono_page_fields["no_elements"])
        # self.tower_page_fields["distance_above"] = float(self.tower_page_fields["distance_above"])
        # self.tower_page_fields["distance_below"] = float(self.tower_page_fields["distance_below"])

        return True

    def checkDialogInput(self, input_fields):
        if isinstance(input_fields, dict):
            values = input_fields
        if isinstance(input_fields, Segment):
            values = {}
            values["id"] = input_fields.segment_id
            values["length_ratio"] = input_fields.length_ratio
            values["Dstar"] = input_fields.diameter_start
            values["Dend"] = input_fields.diameter_end
            values["Tstar"] = input_fields.thickness_start
            values["Tend"] = input_fields.thickness_end
            values["rho"] = input_fields.density
            values["E"] = input_fields.e
            values["G"] = input_fields.g
            values["alpha_s"] = input_fields.alpha_s
            values["alpha_v"] = input_fields.alpha_v
            values["scfStart"] = input_fields.scf_start
            values["scfEnd"] = input_fields.scf_end

        # length_ratio
        v = util.errorMsg_greaterOrequal0_float(values["length_ratio"], "Length ratio")
        if not v:
            return False
        if not values["length_ratio"]:
            values["length_ratio"] = 0.0
        else:
            values["length_ratio"] = float(values["length_ratio"])
        if values["length_ratio"] > 1:
            util.showErrorMsg("Length ratio", "It should not be greater than 1. The sum of all length ratios must be 1")
            return False
        # E
        v = util.errorMsg_greaterthan0_float(values["E"], "E")
        if not v:
            return False

        # G
        v = util.errorMsg_greaterthan0_float(values["G"], "G")
        if not v:
            return False
        # calculating poision ratio E = 2G( 1 + nu)
        nu = float(values["E"]) / (2 * float(values["G"])) - 1
        if nu <= -1 or nu >= 0.5:
            util.showErrorCustomMsg("E and G values are not in correct range",
                                    f"Poison ratio is {nu} which is not between -1 and 0.5")
            return False

        # Dstar
        v = util.errorMsg_greaterOrequal0_float(values["Dstar"], "D_start")
        if not v:
            return False

        # Dend
        v = util.errorMsg_greaterOrequal0_float(values["Dend"], "D_end")
        if not v:
            return False

        # alpha_s
        v = util.errorMsg_greaterOrequal0_float(values["alpha_s"], "alpha_s")
        if not values["alpha_s"]:
            values["alpha_s"] = 0.0
        if not v:
            return False
        elif float(values["alpha_s"]) < 0 or float(values["alpha_s"]) > 1:
            util.showErrorMsg("alpha_s", "It should be greater or equal to 0 and less or equal to 1")
            return False

        # alpha_v

        v = util.errorMsg_greaterOrequal0_float(values["alpha_v"], "alpha_v")
        if not values["alpha_v"]:
            values["alpha_v"] = 0.0
        if not v:
            return False
        elif float(values["alpha_v"]) < 0 or float(values["alpha_v"]) > 1:
            util.showErrorMsg("alpha_v", "It should be greater or equal to 0 and less or equal to 1")
            return False

        # Tstar
        v = util.errorMsg_greaterOrequal0_float(values["Tstar"], "t_start")
        if not v:
            return False

        # Tend
        v = util.errorMsg_greaterOrequal0_float(values["Tend"], "t_end")
        if not v:
            return False

        # scfStart
        v = util.errorMsg_greaterOrequal0_float(values["scfStart"], "SCF_start")
        if not v:
            return False

        # scfEnd
        v = util.errorMsg_greaterOrequal0_float(values["scfEnd"], "SCF_end")
        if not v:
            return False

        # rho
        v = util.errorMsg_greaterOrequal0_float(values["rho"], "Density")
        if not v:
            return False

        return True

    def checkAllSegments(self):
        for id, seg in enumerate(self.segments):
            if not self.checkDialogInput(seg):
                util.showErrorCustomMsg(f"Segment {id + 1} field values were not correct",
                                        "Please correct them as per suggestions. Output file is not created yet")
                return False
        return True

    def dispDialogSegmentTable(self, id=None):
        # current segment

        if id is None or type(id) == bool:
            id = -1 * (self.segment_btns.btn_grp.checkedId() + 2)

        # id = self.cur_segment_id
        cur_segment = self.segments[id]
        Dialog = QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        # Set the labels txts
        ui.lblTitle.setText("Enter data for segment " + str(id + 1))
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
        self.segments_valid_output = False
        if resp == QDialog.Accepted:
            values["length_ratio"] = ui.lineLengthRatio.text()
            values["Dstar"] = ui.lineDStar.text()
            values["Dend"] = ui.lineDEnd.text()
            values["Tstar"] = ui.lineTStar.text()
            values["Tend"] = ui.lineTend.text()
            values["rho"] = ui.lineRho.text()
            values["E"] = ui.lineE.text()
            values["G"] = ui.lineG.text()
            values["alpha_s"] = ui.lineAlphaS.text()
            values["alpha_v"] = ui.lineAlphaV.text()
            values["scfStart"] = ui.lineScfStart.text()
            values["scfEnd"] = ui.lineScfEnd.text()

            cur_segment.segment_id = values["id"]
            cur_segment.length_ratio = values["length_ratio"]
            cur_segment.diameter_start = values["Dstar"]
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

            if self.checkDialogInput(values):
                # everything is good now
                # values = util.convertEmpty2zero(values)
                # values = util.convertDic2float(values)
                cur_segment.convertFields2numeric()
                Dialog.close()
                self.segments_valid_output = True

            else:
                # give use another chance for correct inputs
                self.dispDialogSegmentTable(id)

        else:
            print("Cancel is pressed")

    def display(self, text):
        cur_txt = text
        if cur_txt == 'Please select Input':
            self.ui.stackedWidget.setCurrentWidget(self.ui.Main_page)
        elif cur_txt == 'Monopile':
            self.ui.stackedWidget.setCurrentWidget(self.ui.Mono_page)
        elif cur_txt == 'Jacket 3-Stand':
            self.ui.stackedWidget.setCurrentWidget(self.ui.J3_page)
        else:
            self.ui.stackedWidget.setCurrentWidget(self.ui.J4_page)

    def generate_input_files(self):
        self.monopile_stand.setBeams(self.beams)
        self.monopile_stand.beams[0].setSegments(self.segments)
        # Generate log file containing entered segments input data
        # if self.ui.comboBox_Modal.activated[str].connect(self.display) == "Monopile":  # Monopile

        # Create Monopile object
        self.monopile = Monopile(self.monopile_stand)

        # Raise exception if length ratios not valid
        lengthRatiosValidity = self.monopile.checkLengthRatiosValidity()
        if not lengthRatiosValidity:
            util.showErrorCustomMsg(
                "SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.",
                "Output file is not generated")
            # raise Exception(f'SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.')
            return
        # Generate Beam Input Data (Beam, Segments, Nodes, Elements, etc.) of Monopile
        beamInputGenerated = self.monopile.generateBeamInputData()
        if beamInputGenerated:
            # Write Beam Input File
            self.monopile.writeBeamInput()
            # Write Log File
            self.monopile.writeLogFile()
            # Info Message box stating that input files have successfully been generated
            util.showInfoMsg(tit="Input Files Generated",
                                message="The Monopile input files (Beam Input and Log) have successfully been generated.")
        else:
            # Throw exception
            # raise Exception(f'Error: Beam Input Data (Monopile) not generated.')	# DEBUG_TEST
            pass  # DEBUG_TEST
        # except Exception:
        #     # Error Message box stating that writing to log file has been failed
        #     exc_type, exc_value, exc_tb = sys.exc_info()
        #     error_message = f"Error: Input Files (Monopile) could not be generated.\n{exc_type}: {exc_value}\n"
        #     print(str(traceback.print_exc()))  # DEBUG_TEST
        #     util.showInfoMsg(tit="Input Files Error", message=error_message)

    def main_bts(self):
        self.getValuesFromParent()
        if not self.checkMonoPageInputs():
            return
        if not self.checkAllSegments():
            return
        self.generate_input_files()

    def load_monopile(self):
        self.getValuesFromParent()
        if not self.checkMonoPageInputs():
            return
        # Initialization
        self.monopile_stand = Stand(1, self.mono_page_fields["stand_length"])
        self.monopile_beam = Beam(1, "Stand", self.mono_page_fields["no_segments"], self.mono_page_fields["no_elements"])
        # self.segments = list()
        # Append monopile beam to the list of beams (single item list)
        self.beams = list()
        self.beams.append(self.monopile_beam)

        # taking values from segment table
        self.segments = []
        self.segment_btns = segments_ui(self.ui.tabSupportStructure, self.mono_page_fields["no_segments"])
        self.segment_btns.ui_prop()

        for id in range(self.mono_page_fields["no_segments"]):
            self.segments.append(Segment(segment_id=id))
            self.cur_segment_id = id
            self.segment_btns.btns[id].clicked.connect(self.dispDialogSegmentTable)
        self.ui.btnStructureTowerGGenrateFile_2.setEnabled(True)
        for id in range(self.mono_page_fields["no_segments"]):
            print("Data entering for segment =", id)
            self.dispDialogSegmentTable(id)
            if not self.segments_valid_output:
                util.showWarningMsg("Beam file is not generated!", 'Generated File')
                return
        self.generate_input_files()




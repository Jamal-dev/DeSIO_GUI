from PyQt5 import QtWidget, QCheckBox, QLabel, QLineEdit, QPushButton, QVBoxLayout




class ImageSettingsWindow(QtWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Image Settings")
        self.setFixedSize(250, 300)

        allImg_checkBox = QCheckBox("All Images")
        lastImg_checkBox = QCheckBox("Last Image")
        timedImgs_checkBox = QCheckBox("Timed Images")
        start_label = QLabel("Start Duration:")
        self.txt_start = QLineEdit()
        end_label = QLabel("End Duration:")
        self.txt_end = QLineEdit()
        timeInt_label = QLabel("Time Interval:")
        self.txt_timeInt = QLineEdit()
        maxDiff_label = QLabel("Max Difference:")
        self.txt_maxDiff = QLineEdit()
        buttonSave = QPushButton("Save")

        layout = QVBoxLayout(self)
        layout.addWidget(allImg_checkBox)
        layout.addWidget(lastImg_checkBox)
        layout.addWidget(timedImgs_checkBox)
        layout.addWidget(start_label)
        layout.addWidget(self.txt_start)
        layout.addWidget(end_label)
        layout.addWidget(self.txt_end)
        layout.addWidget(timeInt_label)
        layout.addWidget(self.txt_timeInt)
        layout.addWidget(maxDiff_label)
        layout.addWidget(self.txt_maxDiff)
        layout.addWidget(buttonSave)

        buttonSave.clicked.connect(self.clickedSave)
"""
    def clickedSave(self):
        try:
            param = dict()
            with open("config/image.json", "r") as json_image:
                param = json.load(json_image)

            param["start_duration"] = int(self.txt_start.text())
            param["end_duration"] = int(self.txt_end.text())
            param["time_interval"] = float(self.txt_timeInt.text())
            param["max_diff"] = int(self.txt_maxDiff.text())

            with open("config/image.json", "w") as json_image:
                json.dump(param, json_image)

        except IOError as err:
            print("Error writing the image.json file , please double check", err)
"""
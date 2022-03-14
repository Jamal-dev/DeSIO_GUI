import sys
from pathlib import Path
import os

file_runing_dir = os.path.dirname(os.path.abspath(__file__))
path_main = Path(file_runing_dir) / Path("..")
sys.path.append(str(path_main))

from structure.monopile.monopile import Monopile
from structure.jacket.jacket import Jacket
from beam.stand import Stand
from Utils.compartment import Compartment
from beam.beam import Beam
from beam.segment import Segment
from Utils.utilities import Utilities
from Utils.geometry import Geometry
from form import Ui_MainWindow

import traceback
import matplotlib

# matplotlib.use("TkAgg")
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# from matplotlib.figure import Figure
import re


class SupportStructurePage():

    def __init__(self, parent=None):

        self.ui = parent

        # Support Structure Tab Modes
        self.support_struct_modes = Utilities.getSupportStructureModes()

        # Beam Types
        self.beam_types = Beam.getBeamTypes()  # 0: Stand | 1: Compartment

        # Table Headers
        self.comp_table_headers = Compartment.getTableHeaders()
        self.beams_info_table_headers = Beam.getInfoTableHeaders()
        self.segments_sheet_headers = Segment.getTableHeaders()

        self.selected_mode_index = 0
        self.selected_mode_name = self.support_struct_modes[0][0]
        # self.selected_mode_img_path = cur_dir / Path("img/desio/" + self.support_struct_modes[0][1])

        # Info Section (Structure Tab)
        # self.info_section = tk.Frame(self)
        # self.info_section.pack(fill="both", side="left", expand=1, pady=10)
        #
        # self.topology_frame = tk.Frame(self.info_section, highlightbackground="black", highlightthickness=1)
        # self.topology_frame.grid(row=0, column=0, padx=5, pady=5, sticky='news')
        #
        # self.jacket_info_frame = tk.Frame(self.info_section, highlightbackground="black", highlightthickness=1)
        # self.jacket_info_frame.grid(row=0, column=1, padx=5, pady=5, sticky='news')
        #
        # self.beams_segments_frame = tk.Frame(self.info_section)
        # self.beams_segments_frame.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky='news')
        #
        # self.support_struct_mode_label = tk.Label(self.topology_frame, text="Mode: ", font=("Helvetica", 12, "bold"))
        # self.support_struct_mode_label.grid(row=0, column=0, padx=2, pady=2, sticky='news')

        # Combobox displaying the different support structure modes
        # self.support_struct_mode_combo = ttk.Combobox(self.topology_frame, state="readonly",
        #                                               font=('Helvetica', '12', 'bold'),
        #                                               value=[item[0] for item in self.support_struct_modes])
        # self.support_struct_mode_combo.current(self.selected_mode_index)
        # self.support_struct_mode_combo.grid(row=0, column=1, padx=2, pady=2, sticky='news')
        # self.support_struct_mode_combo.bind("<<ComboboxSelected>>", self.change_support_struct_mode)

        # Image Section (Structure Tab)
        # self.image_section = tk.Frame(self)
        # self.image_section.pack(fill="both", side="right", expand=1, pady=10)
        #
        # self.support_struct_img = Image.open(str(self.selected_mode_img_path))
        # self.support_struct_img = self.support_struct_img.resize((self.photo_size["x"], self.photo_size["y"]),
        #                                                          Image.ANTIALIAS)
        # self.support_struct_img = ImageTk.PhotoImage(self.support_struct_img)
        # self.img_label = tk.Label(self.image_section, image=self.support_struct_img, borderwidth="2", relief="solid")
        # self.img_label.pack()

        # Display Monopile tables as the default option
        # self.generate_monopile_frames()

    def change_support_struct_mode(self, parent = None):
        self.ui = parent
        # Update the image on the right side of the screen
        # self.update_image()

        # Hiding the border of the beams frame in the beginning
        # self.beams_segments_frame.config(highlightbackground=None, highlightthickness=0)

        self.selected_comboBox_Mode = self.comboBox_Modal.activated[str].connect(self.text_var)
        if self.selected_comboBox_Mode == self.support_struct_modes[0][0]:  # Monopile
            self.generate_monopile_frames()
        elif self.selected_comboBox_Mode == self.support_struct_modes[1][0]:  # Jacket 3-Stand
            self.generate_jacket_frames()
        elif self.selected_comboBox_Mode == self.support_struct_modes[2][0]:  # Jacket 4-Stand
            self.generate_jacket_frames()

    def text_var(self, text):
        cur_txt = text

        if cur_txt == 'Please select Input':
            self.Mono_groupBox.hide()
            self.J3_groupBox.hide()
            self.Mono_pic.hide()
            self.J3_pic.hide()
            self.J4_pic.hide()

        elif cur_txt == 'Monopile':
            self.Mono_groupBox.show()
            self.J3_groupBox.hide()
            self.Mono_pic.show()
            self.J3_pic.hide()
            self.J4_pic.hide()

        elif cur_txt == 'Jacket 3-Stand':
            self.J3_groupBox.show()
            self.Mono_groupBox.hide()
            self.Mono_pic.hide()
            self.J3_pic.show()
            self.J4_pic.hide()

        else:
            self.J3_groupBox.show()
            self.Mono_groupBox.hide()
            self.Mono_pic.hide()
            self.J3_pic.hide()
            self.J4_pic.show()

    # def update_image(self):
    #     # Refresh the image section frame
    #     for child_widget in self.image_section.winfo_children():
    #         child_widget.destroy()
    #
    #     # Set frame border to zero by default
    #     self.image_section.config(borderwidth="0")
    #
    #     self.selected_mode_name = self.support_struct_mode_combo.get()
    #     self.selected_mode_img_path = ""
    #
    #     for support_struct_mode in self.support_struct_modes:
    #         if support_struct_mode[0] == self.selected_mode_name:
    #             self.selected_mode_img_path = cur_dir / Path(f"img/desio/{support_struct_mode[1]}")
    #             break
    #
    #     # Update mode name and image in the homepage
    #     if self.selected_mode_img_path != "":
    #         self.support_struct_img = Image.open(str(self.selected_mode_img_path))
    #         self.support_struct_img = self.support_struct_img.resize((self.photo_size["x"], self.photo_size["y"]),
    #                                                                  Image.ANTIALIAS)
    #         self.support_struct_img = ImageTk.PhotoImage(self.support_struct_img)
    #
    #         self.img_label = tk.Label(self.image_section, image=self.support_struct_img, borderwidth="2",
    #                                   relief="solid")
    #         self.img_label.pack()
    #
    # def update_graph(self, scatter_coordinates, line_end_points, windowTitle='3D-Simulation', orthoBaseLength=0.5,
    #                  scatterEnabled=True, axesOn=True, orthonormalBaseOn=True):
    #     # Description: Display the generated graph (monopile/jacket) in the image section frame
    #
    #     # Refresh the image section frame
    #     for child_widget in self.image_section.winfo_children():
    #         child_widget.destroy()
    #
    #     # Set solid frame border for graph
    #     self.image_section.config(borderwidth="2", relief="solid")
    #
    #     # Saving the input parameters for later use
    #     self.scatter_coordinates = scatter_coordinates
    #     self.line_end_points = line_end_points
    #     self.windowTitle = windowTitle
    #     self.orthoBaseLength = orthoBaseLength
    #     self.scatterEnabled = scatterEnabled
    #     self.axesOn = axesOn
    #     self.orthonormalBaseOn = orthonormalBaseOn
    #
    #     # Create a new figure
    #     self.fig = Figure(figsize=(5.5, 5.5), dpi=100)
    #
    #     # Display the figure in a canvas
    #     self.canvas = FigureCanvasTkAgg(self.fig, self.image_section)
    #     self.canvas.draw()
    #
    #     # Plot data into the figure
    #     self.ax = self.fig.add_subplot(111, projection="3d")
    #     Geometry.plotDataInAxes(self.ax, scatter_coordinates, line_end_points, windowTitle, orthoBaseLength,
    #                             scatterEnabled, axesOn, orthonormalBaseOn)
    #
    #     # Add the default Matplotlib Navigation Toolbar to Figure
    #     self.toolbar = NavigationToolbar2Tk(self.canvas, self.image_section)
    #     self.toolbar.update()
    #
    #     # Creates maximize (new window) button
    #     self.maximize_button = tk.Button(master=self.toolbar, command=self.maximizeWindow)
    #     self.maximize_button.pack(side="left")
    #
    #     # Sets icon for the maximize (new window) button
    #     self.maximize_icon = tk.PhotoImage(file=cur_dir / Path('img/logos/maximize-2.png'))
    #     self.maximize_button.config(image=self.maximize_icon)
    #
    #     # Tooltip for maximize (new window) button
    #     self.maximize_button_tooltip = CustomToolTip(self.maximize_button, "Maximize - open graph in new window")
    #
    #     # Pack the graph and toolbar into the canvas
    #     self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
    #     self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def maximizeWindow(self):
        # Displays the generated graph in a new window
        Geometry.displayPlot(self.scatter_coordinates, self.line_end_points, self.windowTitle,
                             self.orthoBaseLength, self.scatterEnabled, self.axesOn, self.orthonormalBaseOn)

    def generate_monopile_frames(self):
        for child_widget in self.jacket_info_frame.winfo_children():
            child_widget.destroy()
        for child_widget in self.beams_segments_frame.winfo_children():
            child_widget.destroy()

        # Define the selected mode and the default number of stands for the mode
        global selected_mode_name
        self.selected_mode_name = self.support_struct_mode_combo.get()
        global num_stands
        self.num_stands = 1
        global monopile
        self.monopile = None

        # Define remaining global variables
        global num_compartments
        self.num_compartments = 0  # fixed value: 0
        global num_beams
        self.num_beams = self.num_compartments + self.num_stands
        global beam_name
        self.beam_name = "Beam 1"
        global beam_type
        self.beam_type = self.beam_types[0]  # "Stand"

        # Set Text Widgets
        self.stand_length_label = tk.Label(self.jacket_info_frame, text="Stand length: ",
                                           font=("Helvetica", 12, "bold"))
        self.stand_length_label.grid(row=0, column=0, padx=2, pady=2)
        self.stand_length_field = tk.Text(self.jacket_info_frame, height=1, width=5, font=("Helvetica", 12))
        self.stand_length_field.grid(row=0, column=1, padx=4, pady=2)
        self.num_segments_label = tk.Label(self.jacket_info_frame, text="No. of Segments: ",
                                           font=("Helvetica", 12, "bold"))
        self.num_segments_label.grid(row=1, column=0, padx=2, pady=2)

        self.num_segments_field = tk.Text(self.jacket_info_frame, height=1, width=5, font=("Helvetica", 12))
        self.num_segments_field.grid(row=1, column=1, padx=4, pady=2)

        self.num_elements_label = tk.Label(self.jacket_info_frame, text="No. of Elements: ",
                                           font=("Helvetica", 12, "bold"))
        self.num_elements_label.grid(row=2, column=0, padx=2, pady=2)

        self.num_elements_field = tk.Text(self.jacket_info_frame, height=1, width=5, font=("Helvetica", 12))
        self.num_elements_field.grid(row=2, column=1, padx=4, pady=2)

        # Button to create segments table
        self.create_segments_table_button = tk.Button(self.jacket_info_frame, text="Segments Table",
                                                      font=("Helvetica", 10, "bold"),
                                                      bg="#e3b129",
                                                      command=self.create_segments_table_monopile)
        self.create_segments_table_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

    def create_segments_table_monopile(self):
        # Logic block to accept only positive integers as valid input
        stand_len = self.stand_length_field.get("1.0", 'end-1c')
        stand_len_valid = Utilities.isInt(stand_len) and Utilities.isPos(stand_len)
        segments_num = self.num_segments_field.get("1.0", 'end-1c')
        segments_num_valid = Utilities.isInt(segments_num) and Utilities.isPos(segments_num)
        elements_num = self.num_elements_field.get("1.0", 'end-1c')
        elements_num_valid = Utilities.isInt(elements_num) and Utilities.isPos(elements_num)
        if not stand_len_valid or not segments_num_valid or not elements_num_valid:
            # Error Message Box displayed upon invalid input
            messagebox.showerror(title=f"Invalid Input",
                                 message=f"Error: Invalid Input. Positive Integer values expected as input.")
            return

        # Store the monopile values if all above conditions met
        global stand_length
        self.stand_length = int(stand_len)
        global num_segments
        self.num_segments = int(segments_num)
        global num_elements
        self.num_elements = int(elements_num)
        # Define global variables for the segments sheet
        global segments_custom_sheet
        self.segments_custom_sheet = None
        global segments_sheet
        self.segments_sheet = None

        # Showing the border of the beams frame
        self.beams_segments_frame.config(highlightbackground="black", highlightthickness=1)

        # Set Table
        # Segments Sheet Parameters
        n_rows = self.num_segments
        n_cols = len(self.segments_sheet_headers)

        # Create, Store and Display Sheet
        self.segments_custom_sheet = CustomSheet(self.beams_segments_frame,
                                                 self.segments_sheet_headers, n_rows, n_cols, 100, 600)
        self.segments_sheet = self.segments_custom_sheet.create_sheet()

        # Initialize all cells of 'Segment' column with respective Segment Serial Nos. (read only)
        for r in range(n_rows):
            self.segments_sheet.set_cell_data(r, 0, str(r + 1))
            self.segments_sheet.readonly_cells(r, 0)

        # Frame for 'Generate input files' Button
        self.monopile_gen_input_button_frame = tk.Frame(self.beams_segments_frame)
        self.monopile_gen_input_button_frame.grid(row=1, column=0, columnspan=2, padx=2, pady=2)

        # Add 'Generate input files' button right below the segments table
        self.monopile_gen_input_button = tk.Button(self.monopile_gen_input_button_frame, text="Generate input files",
                                                   font=("Helvetica", 10, "bold"), bg="#2bd659",
                                                   command=self.generate_input_files)
        self.monopile_gen_input_button.pack(padx=5, pady=5)

    # def generate_jacket_frames(self):
    #     # Delete previous children of the frames
    #     for child_widget in self.jacket_info_frame.winfo_children():
    #         child_widget.destroy()
    #     for child_widget in self.beams_segments_frame.winfo_children():
    #         child_widget.destroy()
    #
    #     # Define the selected mode and the default number of stands for the mode
    #     global selected_mode_name
    #     self.selected_mode_name = self.support_struct_mode_combo.get()
    #     global num_stands
    #     if self.selected_mode_name == self.support_struct_modes[1][0]:  # Jacket 3-Stand
    #         self.num_stands = 3
    #     else:  # Jacket 4-Stand
    #         self.num_stands = 4
    #     global jacket
    #     self.jacket = None
    #
    #     # Input Field frame
    #     self.input_field_frame = tk.Frame(self.jacket_info_frame)
    #     self.input_field_frame.pack()
    #
    #     self.num_compartments_label = tk.Label(self.input_field_frame, text="No. of Compartments: ",
    #                                            font=("Helvetica", 12, "bold"))
    #     self.num_compartments_label.grid(row=0, column=0, padx=2, pady=2)
    #
    #     self.stand_length_label = tk.Label(self.input_field_frame, text="Stand length: ",
    #                                        font=("Helvetica", 12, "bold"))
    #     self.stand_length_label.grid(row=1, column=0, padx=2, pady=2)
    #
    #     self.distance_above_label = tk.Label(self.input_field_frame, text="Distance above: ",
    #                                          font=("Helvetica", 12, "bold"))
    #     self.distance_above_label.grid(row=0, column=2, padx=2, pady=2)
    #
    #     self.distance_below_label = tk.Label(self.input_field_frame, text="Distance below: ",
    #                                          font=("Helvetica", 12, "bold"))
    #     self.distance_below_label.grid(row=1, column=2, padx=2, pady=2)
    #
    #     self.gap_from_below_label = tk.Label(self.input_field_frame, text="Gap from below: ",
    #                                          font=("Helvetica", 12, "bold"))
    #     self.gap_from_below_label.grid(row=2, column=0, padx=2, pady=2)
    #
    #     self.num_compartments_field = tk.Text(self.input_field_frame, height=1, width=5, font=("Helvetica", 10))
    #     self.num_compartments_field.grid(row=0, column=1, padx=4, pady=2)
    #
    #     self.stand_length_field = tk.Text(self.input_field_frame, height=1, width=5, font=("Helvetica", 10))
    #     self.stand_length_field.grid(row=1, column=1, padx=4, pady=2)
    #
    #     self.distance_above_field = tk.Text(self.input_field_frame, height=1, width=5, font=("Helvetica", 10))
    #     self.distance_above_field.grid(row=0, column=3, padx=4, pady=2)
    #
    #     self.distance_below_field = tk.Text(self.input_field_frame, height=1, width=5, font=("Helvetica", 10))
    #     self.distance_below_field.grid(row=1, column=3, padx=4, pady=2)
    #
    #     self.gap_from_below_field = tk.Text(self.input_field_frame, height=1, width=5, font=("Helvetica", 10))
    #     self.gap_from_below_field.grid(row=2, column=1, padx=4, pady=2)
    #
    #     # Button to create beam and compartment info tables
    #     self.create_comp_beam_table_button = tk.Button(self.input_field_frame,
    #                                                    text="Compartment/Beam Tables", bg="#e3b129",
    #                                                    font=("Helvetica", 10, "bold"),
    #                                                    command=self.create_comp_and_beam_tables)
    #     self.create_comp_beam_table_button.grid(row=3, column=0, columnspan=4, padx=5, pady=5)
    #
    #     # Create table frames for the compartments, beams and segments
    #     self.compartments_table_frame = tk.Frame(self.beams_segments_frame)
    #     self.compartments_table_frame.grid(row=0, column=0, padx=2, pady=2)
    #
    #     # Frame for Beams Info Table
    #     self.beams_info_table_frame = tk.Frame(self.beams_segments_frame)
    #     self.beams_info_table_frame.grid(row=0, column=1, padx=2, pady=2)
    #
    #     # Frame for Button to create segment tables
    #     self.create_segments_tables_button_frame = tk.Frame(self.beams_segments_frame)
    #     self.create_segments_tables_button_frame.grid(row=1, column=0, columnspan=2, padx=2, pady=2)
    #
    #     # Frame for Segments Tables
    #     self.segments_tables_frame = tk.Frame(self.beams_segments_frame)
    #     self.segments_tables_frame.grid(row=2, column=0, columnspan=2, padx=2, pady=2)
    #
    #     # Frame for 'Generate input files' Button
    #     self.jacket_gen_input_button_frame = tk.Frame(self.beams_segments_frame)
    #     self.jacket_gen_input_button_frame.grid(row=3, column=0, columnspan=2)

    def create_comp_and_beam_tables(self):
        # Logic block to accept only positive values (within bounds) as valid input
        global comp_max_bound
        self.comp_max_bound = 5
        comp_num = self.num_compartments_field.get("1.0", 'end-1c')
        comp_num_valid = Utilities.isInt(comp_num) and Utilities.isPos(comp_num)
        stand_len = self.stand_length_field.get("1.0", 'end-1c')
        stand_len_valid = Utilities.isFloat(stand_len) and Utilities.isPos(stand_len)
        dist_above = self.distance_above_field.get("1.0", 'end-1c')
        dist_above_valid = Utilities.isFloat(dist_above) and Utilities.isPos(dist_above)
        dist_below = self.distance_below_field.get("1.0", 'end-1c')
        dist_below_valid = Utilities.isFloat(dist_below) and Utilities.isPos(dist_below)
        gap_fr_below = self.gap_from_below_field.get("1.0", 'end-1c')
        gap_fr_below_valid = Utilities.isFloat(gap_fr_below) and Utilities.isPos(gap_fr_below)

        if (not gap_fr_below_valid or not dist_below_valid or not dist_above_valid):
            # Error Message Box displayed upon invalid input
            messagebox.showerror(title=f"Invalid Input",
                                 message=f"Error: Invalid Input. Positive values expected as general input.")
            return
        elif (not comp_num_valid or not stand_len_valid):
            # Error Message Box displayed upon invalid input (expected integer, not float - for specific entries)
            messagebox.showerror(title=f"Invalid Input",
                                 message=f"Error: Invalid Input. Positive integer values expected as input for 'No. of compartments' and 'Stand length'.")
            return
        elif int(comp_num) > self.comp_max_bound:
            # Error Message Box displayed upon exceeding Max Bound of Beams
            messagebox.showerror(title=f"Max Bound Reached",
                                 message=f"Error: Input (No. of compartments) exceeds Max Bound (max 5 allowed).")
            return

        # Store the jacket parameter values if all conditions met
        global num_compartments
        self.num_compartments = int(comp_num)
        global stand_length
        self.stand_length = float(stand_len)
        global distance_above
        self.distance_above = float(dist_above)
        global distance_below
        self.distance_below = float(dist_below)
        global gap_from_below
        self.gap_from_below = float(gap_fr_below)
        # Get Beam classes from the no. of compartments
        global beam_classes
        self.beam_classes = Beam.getBeamClasses(self.num_compartments)
        global stand_beam_classes
        self.stand_beam_classes = Beam.getStandBeamClasses(self.num_compartments)
        global comp_beam_classes
        self.comp_beam_classes = Beam.getCompBeamClasses(self.num_compartments)

        # Delete previous children of all the table frames
        for child_widget in self.compartments_table_frame.winfo_children():
            child_widget.destroy()
        for child_widget in self.beams_info_table_frame.winfo_children():
            child_widget.destroy()
        for child_widget in self.create_segments_tables_button_frame.winfo_children():
            child_widget.destroy()
        for child_widget in self.segments_tables_frame.winfo_children():
            child_widget.destroy()
        for child_widget in self.jacket_gen_input_button_frame.winfo_children():
            child_widget.destroy()

        # Task 1: Creating Compartments Table

        self.comp_table_title = tk.Label(self.compartments_table_frame, text="Compartments Table",
                                         font=("Helvetica", 12, "bold"))
        self.comp_table_title.pack(padx=2, pady=2)

        # Compartment Table Frame
        self.comp_table_container = tk.Frame(self.compartments_table_frame)
        self.comp_table_container.pack(fill="both", expand=1, padx=2, pady=2)
        # Compartment Table Parameters
        n_rows = self.num_compartments
        n_cols = len(self.comp_table_headers)

        # Create, Store and Display Sheet
        self.comp_custom_sheet = CustomSheet(self.comp_table_container, self.comp_table_headers,
                                             n_rows, n_cols, 100, 300)
        global comp_sheet
        self.comp_sheet = self.comp_custom_sheet.create_sheet()

        # Initialize all cells of 'Compartment' column with respective Compartment Serial Nos. (read only)
        for r in range(n_rows):
            self.comp_sheet.set_cell_data(r, 0, str(r + 1))
            self.comp_sheet.readonly_cells(r, 0)

        # Task 2: Creating Beams Info Table

        self.beams_info_table_title = tk.Label(self.beams_info_table_frame, text="Beams Info Table",
                                               font=("Helvetica", 12, "bold"))
        self.beams_info_table_title.pack(padx=2, pady=2)

        # Beams Table Frame
        self.beams_info_table_container = tk.Frame(self.beams_info_table_frame)
        self.beams_info_table_container.pack(fill="both", expand=1, padx=2, pady=2)
        # Beams Table Parameters
        global num_beams
        self.num_beams = (self.num_compartments * self.num_stands * Jacket.getNumBeamsPerComp()) + (
                    self.num_stands * (self.num_compartments + 2))
        global num_beam_classes
        self.num_beam_classes = len(self.beam_classes)
        n_rows = self.num_beam_classes
        n_cols = len(self.beams_info_table_headers)

        # Create, Store and Display Sheet
        self.beams_info_custom_sheet = CustomSheet(self.beams_info_table_container, self.beams_info_table_headers,
                                                   n_rows, n_cols, 100, 400)
        global beams_info_sheet
        self.beams_info_sheet = self.beams_info_custom_sheet.create_sheet()

        # Initialize all cells of 'Beam Class' column with respective Beam Class IDs (Read only)
        # Initialize all cells of 'Description' column with descriptions for respective Beam Classes (Read only)
        for r in range(n_rows):
            curr_beam_class = str(self.beam_classes[r][0])
            curr_descr = str(self.beam_classes[r][1])
            self.beams_info_sheet.set_cell_data(r, 0, curr_beam_class)
            self.beams_info_sheet.set_cell_data(r, 1, curr_descr)
            self.beams_info_sheet.readonly_cells(r, 0)
            self.beams_info_sheet.readonly_cells(r, 1)

        # Task 3: Show Border and add Beams-Segments Button

        # Showing the border of the beams frame
        self.beams_segments_frame.config(highlightbackground="black", highlightthickness=1)

        # Adding the Beams_Segments Button
        self.create_segments_tables_button = tk.Button(self.create_segments_tables_button_frame,
                                                       text="Segments Tables", font=("Helvetica", 10, "bold"),
                                                       bg="#e3b129",
                                                       command=self.create_segments_tables_jackets)
        self.create_segments_tables_button.pack(padx=5, pady=5)

    # def create_segments_tables_jackets(self):
    #     # Check if the input no. of segments values are valid for all beams
    #     num_segments_list = []
    #     for r in range(self.beams_info_sheet.get_total_rows()):
    #         curr_num_segments = self.beams_info_sheet.get_cell_data(r, 2)
    #         num_segments_input_valid = Utilities.isInt(curr_num_segments) and Utilities.isPos(curr_num_segments)
    #         if not num_segments_input_valid:
    #             # Error Message Box displayed upon invalid input
    #             messagebox.showerror(title=f"Invalid Input",
    #                                  message=f"Error: Invalid Input. All 'No. of Segments' values in the 'Beams Info' Table must be Positive Integers.")
    #             return
    #         num_segments_list.append(curr_num_segments)
    #
    #     # Delete previous children of the segments tables and 'Generate input files' button frames
    #     for child_widget in self.segments_tables_frame.winfo_children():
    #         child_widget.destroy()
    #     for child_widget in self.jacket_gen_input_button_frame.winfo_children():
    #         child_widget.destroy()
    #
    #     # Create a Tab View (a.k.a 'Beams-Segments' Notebook)
    #     global beams_segments_notebook
    #     self.beams_segments_notebook = ttk.Notebook(self.segments_tables_frame)
    #     self.beams_segments_notebook.grid(row=4, column=0, columnspan=2, padx=2, pady=2)
    #
    #     # Define global variables for the beam tables (dictionaries use beam title as key)
    #     global beam_class_names
    #     self.beam_class_names = [str(beam_class_info[0]) for beam_class_info in self.beam_classes]
    #     global beam_class_tabs
    #     self.beam_class_tabs = dict()
    #     global num_segments
    #     self.num_segments = dict()
    #     global segments_sheet_containers
    #     self.segments_sheet_containers = dict()
    #     global segments_custom_sheets
    #     self.segments_custom_sheets = dict()
    #     global segments_sheets
    #     self.segments_sheets = dict()
    #
    #     for beam_class_ind in range(self.num_beam_classes):
    #         # Initialize Tab
    #         curr_header = self.beam_class_names[beam_class_ind]
    #         self.beam_class_tabs[curr_header] = tk.Frame(self.beams_segments_notebook)
    #
    #         # Create and Add Segments Sheet to Tab
    #         # Segments Sheet Frame
    #         self.segments_sheet_containers[curr_header] = tk.Frame(self.beam_class_tabs[curr_header])
    #         self.segments_sheet_containers[curr_header].pack(fill="both", expand=1, padx=2, pady=2)
    #
    #         # Segments Sheet Parameters
    #         self.num_segments[curr_header] = int(num_segments_list[beam_class_ind])
    #         n_rows = self.num_segments[curr_header]
    #         n_cols = len(self.segments_sheet_headers)
    #
    #         # Create, Store and Display Sheet
    #         self.segments_custom_sheets[curr_header] = CustomSheet(self.segments_sheet_containers[curr_header],
    #                                                                self.segments_sheet_headers, n_rows, n_cols, 100,
    #                                                                800)
    #         self.segments_sheets[curr_header] = self.segments_custom_sheets[curr_header].create_sheet()
    #
    #         # Initialize all cells of 'Segment' column with respective Segment Serial Nos. (read only)
    #         for r in range(n_rows):
    #             self.segments_sheets[curr_header].set_cell_data(r, 0, str(r + 1))
    #             self.segments_sheets[curr_header].readonly_cells(r, 0)
    #
    #         # Pack tab
    #         self.beam_class_tabs[curr_header].pack(fill="both", expand=1)
    #         # Add tab to the 'Beams-Segments' notebook
    #         self.beams_segments_notebook.add(self.beam_class_tabs[curr_header], text=curr_header)
    #     # Adding the 'Generate input files' button right below the segments table
    #     self.jacket_gen_input_button = tk.Button(self.jacket_gen_input_button_frame, text=f"Generate input files",
    #                                              font=("Helvetica", 10, "bold"), bg="#2bd659",
    #                                              command=self.generate_input_files)
    #     self.jacket_gen_input_button.pack(padx=5, pady=5)

    def generate_input_files(self):
        # Generate log file containing entered segments input data
        if self.selected_mode_name == self.support_struct_modes[0][0]:  # Monopile
            try:
                # Check if all segment sheet values are valid and complete - If False, throw exception
                if not Utilities.areSegmentsValid(self.segments_sheet):
                    raise Exception(f'Segments data invalid/incomplete. All values must be positive.')

                global monopile
                self.monopile = self.load_monopile()

                # Generate Beam Input Data (Beam, Segments, Nodes, Elements, etc.) of Monopile
                beamInputGenerated = self.monopile.generateBeamInputData()
                if beamInputGenerated:
                    # Write Beam Input File
                    self.monopile.writeBeamInput()
                    # Write Log File
                    self.monopile.writeLogFile()
                    # Info Message box stating that input files have successfully been generated
                    messagebox.showinfo(title=f"Input Files Generated",
                                        message=f"The Monopile input files (Beam Input and Log) have successfully been generated.")
                    # Display Monopile graph in image frame
                    self.update_graph(self.monopile.coordinates, self.monopile.line_end_points,
                                      '3D-Simulation (Monopile)')
                else:
                    # Throw exception
                    # raise Exception(f'Error: Beam Input Data (Monopile) not generated.')	# DEBUG_TEST
                    pass  # DEBUG_TEST
            except Exception:
                # Error Message box stating that writing to log file has been failed
                exc_type, exc_value, exc_tb = sys.exc_info()
                error_message = f"Error: Input Files (Monopile) could not be generated.\n{exc_type}: {exc_value}\n"
                print(str(traceback.print_exc()))  # DEBUG_TEST
                messagebox.showerror(title=f"Input Files Error", message=error_message)
        else:  # Jackets
            try:
                # Check if all segment sheet values are valid and complete - If False, throw exception
                # Iterate through the beam classes and perform the test for each segments sheet
                for beam_class_ind in range(self.num_beam_classes):
                    # Get header (beam class) to access the segments sheets corresponding to the respective header
                    curr_header = self.beam_class_names[beam_class_ind]
                    # Check if the segments sheet values corresponding to the current beam class are valid
                    if not Utilities.areSegmentsValid(self.segments_sheets[curr_header]):
                        raise Exception(f'Segments data invalid/incomplete. All values must be positive.')

                global jacket
                self.jacket = self.load_jacket()

                # Generate Beam Input Data (Beam, Segments, Nodes, Elements, etc.) of Jacket
                beamInputGenerated = self.jacket.generateBeamInputData()
                if beamInputGenerated:
                    # Write Beam Input File
                    self.jacket.writeBeamInput()
                    # Write Log File
                    self.jacket.writeLogFile()
                    # Info Message box stating that input files have successfully been generated
                    messagebox.showinfo(title=f"Input Files Generated",
                                        message=f"The Jacket input files (Beam Input and Log) have successfully been generated.")
                    # Display Jacket graph in image frame
                    self.update_graph(self.jacket.all_coordinates, self.jacket.all_line_end_points,
                                      '3D-Simulation (Jacket)')
                else:
                    # Throw exception
                    # raise Exception(f'Error: Beam Input Data (Jacket) not generated.')	# DEBUG_TEST
                    pass  # DEBUG_TEST
            except Exception:
                # Error Message box stating that writing to log file has been failed
                exc_type, exc_value, exc_tb = sys.exc_info()
                error_message = f"Error: Input Files (Jacket) could not be generated.\n{exc_type}: {exc_value}\n"
                print(str(traceback.print_exc()))  # DEBUG_TEST
                messagebox.showerror(title=f"Input Files Error", message=error_message)

    def load_monopile(self):
        # Initialization
        monopile_stand = Stand(1, self.stand_length)
        monopile_beam = Beam(1, "Stand", self.num_segments, self.num_elements)
        segments = list()
        # Append monopile beam to the list of beams (single item list)
        beams = list()
        beams.append(monopile_beam)

        # Get segments data from corresponding segments sheet
        segments_sheet_data = self.segments_sheet.get_sheet_data(return_copy=True, get_header=False, get_index=False)
        # Iterate through segments data and get list of segments
        for segment_ind in range(self.num_segments):
            segments.append(Segment.fromRow(segments_sheet_data[segment_ind]))

        # Setting intermediate parameters
        monopile_stand.setBeams(beams)
        beams[0].setSegments(segments)

        # Create Monopile object
        monopile = Monopile(monopile_stand)

        # Raise exception if length ratios not valid
        lengthRatiosValidity = monopile.checkLengthRatiosValidity()
        if not lengthRatiosValidity:
            raise Exception(f'SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.')

        # Return monopile object for later use
        return monopile

    # def load_jacket(self):
    #     # Define Jacket object with initial parameters
    #     jacket = Jacket(self.num_stands, self.num_compartments, self.stand_length, self.distance_above,
    #                     self.distance_below, self.gap_from_below)
    #
    #     # Get the data from the compartments sheet
    #     comp_sheet_data = self.comp_sheet.get_sheet_data(return_copy=True, get_header=False, get_index=False)
    #
    #     # Iterate through stands and compartments and append items to respective lists
    #     jacket_stands = list()
    #     jacket_comps = list()
    #     for stand_ind in range(self.num_stands):
    #         jacket_stands.append(Stand((stand_ind + 1), jacket.stand_length))
    #     for comp_ind in range(self.num_compartments):
    #         curr_comp_height = float(comp_sheet_data[comp_ind][1])
    #         jacket_comps.append(Compartment((comp_ind + 1), curr_comp_height))
    #
    #     # Setting Stands and Compartments of Jacket
    #     jacket.setStands(jacket_stands)
    #     jacket.setCompartments(jacket_comps)
    #
    #     # Raises exception if compartment heights (and gap from below value) not compatible
    #     heightsCompatibility = jacket.checkHeightsCompatibility()
    #     if not heightsCompatibility:
    #         raise Exception(
    #             f'Compartment heights (and "Gap from below" value) not compatible for Jacket.\nExceeds the max Jacket Height bound.\nMax Jacket Height: {jacket.findJacketHeight():.3f}')
    #
    #     # Initializing additional parameters (beam and segment settings)
    #     stand_beam_n_elems = list()
    #     stand_beam_segments_settings = list()
    #     upper_comp_beam_n_elems = list()
    #     upper_comp_beam_segments_settings = list()
    #     lower_comp_beam_n_elems = list()
    #     lower_comp_beam_segments_settings = list()
    #
    #     # Get the data from the beams info sheet
    #     beams_info_sheet_data = self.beams_info_sheet.get_sheet_data(return_copy=True, get_header=False,
    #                                                                  get_index=False)
    #
    #     # Iterate through the beam info data and get no. of elements for each beam class
    #     for beam_class_ind in range(self.num_beam_classes):
    #         # Get header (beam class) to access the segments sheets corresponding to the respective header
    #         curr_header = self.beam_class_names[beam_class_ind]
    #         # Get the data from the segments sheet corresponding to the current beam class
    #         curr_segments_sheet_data = self.segments_sheets[curr_header].get_sheet_data(return_copy=True,
    #                                                                                     get_header=False,
    #                                                                                     get_index=False)
    #
    #         # Append number of elements and list of segments for each umbrella beam class
    #         if beam_class_ind < len(self.stand_beam_classes):  # Stand
    #             stand_beam_n_elems.append(int(beams_info_sheet_data[beam_class_ind][3]))
    #             curr_stand_segments_settings = list()
    #             for segment_ind in range(len(curr_segments_sheet_data)):
    #                 curr_stand_segments_settings.append(Segment.fromRow(curr_segments_sheet_data[segment_ind]))
    #             stand_beam_segments_settings.append(curr_stand_segments_settings)
    #         elif re.match("^C.*L$", self.beam_class_names[beam_class_ind]):  # Compartment (Lower)
    #             lower_comp_beam_n_elems.append(int(beams_info_sheet_data[beam_class_ind][3]))
    #             curr_lower_comp_segments_settings = list()
    #             for segment_ind in range(len(curr_segments_sheet_data)):
    #                 curr_lower_comp_segments_settings.append(Segment.fromRow(curr_segments_sheet_data[segment_ind]))
    #             lower_comp_beam_segments_settings.append(curr_lower_comp_segments_settings)
    #         else:  # Compartment (Upper)
    #             upper_comp_beam_n_elems.append(int(beams_info_sheet_data[beam_class_ind][3]))
    #             curr_upper_comp_segments_settings = list()
    #             for segment_ind in range(len(curr_segments_sheet_data)):
    #                 curr_upper_comp_segments_settings.append(Segment.fromRow(curr_segments_sheet_data[segment_ind]))
    #             upper_comp_beam_segments_settings.append(curr_upper_comp_segments_settings)
    #
    #     # Setting additional parameters (beam and segment settings) of Jacket
    #     jacket.setCompBeamsData(stand_beam_n_elems, stand_beam_segments_settings, upper_comp_beam_n_elems,
    #                             upper_comp_beam_segments_settings, lower_comp_beam_n_elems,
    #                             lower_comp_beam_segments_settings)
    #
    #     # Raise exception if length ratios not valid
    #     lengthRatiosValidity = jacket.checkLengthRatiosValidity()
    #     if not lengthRatiosValidity:
    #         raise Exception(f'SegmentsLengthRatioSumError: Length Ratios do not sum upto 1 for certain segments.')
    #
    #     # Return jacket object for later use
    #     return jacket

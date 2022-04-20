
# DeSIO_GUI

It's a gui for windmill simulation data. Build via user Data and gives numerial output of simulation data.



## Video Presentation

![Demo GUI Representation](https://github.com/Jamal-dev/DeSIO_GUI/blob/main/desio/DeSIO-2022-04-18-21-20-07.gif)


## Screen Shorts

![Tower and monopile beam](https://github.com/Jamal-dev/DeSIO_GUI/blob/main/desio/tower_beam.png)
![Jackert 3](https://github.com/Jamal-dev/DeSIO_GUI/blob/main/desio/support_struct_type_2_jacket3stand.png)
![Jacket 4](https://github.com/Jamal-dev/DeSIO_GUI/blob/main/desio/support_struct_type_3_jacket4stand.png)


## PyQt Designer Currently available Components

- Slider (attachable to parameters)
- QGroupBox (Group all widget)
- QPushButton (Push button)
- QLabel (Structur text via label) 
- QLineEdit (User data input)
- QComboBox (choose available option)
- QTabWidget (Tab widget)


## Design content of GUI 
- forum.py
All geometrical data regarding Mainwindow with Label and structure are contained in this.py file.

- segment_table.py
This .py file contains all geometrical information on the segment table, including Label and structure.

- mplwidget.py
This class file contains all geometrical data about 3-dimensional and 2-dimensional maps created with the help of libraries such as matplotlib and toolkits.



## Documentation

[Documentation](https://github.com/Jamal-dev/DeSIO_GUI/blob/main/DeSIO.pdf)

The user may find all documentation for all classes and functions in the provided pdf file, and the rest of the information is delivered as a comment inside the code.


## Classes stored in different folders

[beam](https://github.com/Jamal-dev/DeSIO_GUI/tree/main/beam)
- beam.py
- element.py
- node.py
- segment.py
- stand.py

[gui_fun_file](https://github.com/Jamal-dev/DeSIO_GUI/tree/main/gui_fun_files)
- monopile_page.py
- tower_page.py

[io](https://github.com/Jamal-dev/DeSIO_GUI/tree/main/io)
- genrated output file

[structure](https://github.com/Jamal-dev/DeSIO_GUI/tree/main/structure)
- jacket.py
- monopile.py
- tower.py

[Utils](https://github.com/Jamal-dev/DeSIO_GUI/tree/main/Utils)
- beamsInfo_userInterface.py
- compartment.py
- compartments_userInterface.py
- crossSectionProperty.py
- geometry.py
- segments_userInterface.py
- Support_Structure.py
- utilities.py



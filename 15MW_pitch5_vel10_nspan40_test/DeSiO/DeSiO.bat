set path_DeSiO=C:\Users\hente\Desktop\DeSiO_compile\DeSiO\DeSiO\x64\Release
set path_mtee=.
REM
set path_fsi=/home/jamal/qt_project_learning/DeSIO_GUI/15MW_pitch5_vel10_nspan40_test/DeSiO
set path_struc=/home/jamal/qt_project_learning/DeSIO_GUI/15MW_pitch5_vel10_nspan40_test/DeSiO-Structure
REM Start fsi simulation with DeSiO
cd %path_fsi%
REM Boost thread priority
SET desio_exe=C:\Users\hente\Desktop\DeSiO_compile\DeSiO\DeSiO\x64\Release\DeSiO.exe
start "" /REALTIME /B /W  %desio_exe% | ".\mtee.exe" "%path_fsi%\output.txt"
pause

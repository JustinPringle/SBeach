SBEACH python module

Currently this version of SBEACH can only be run on UNIX and Windows operating systems.
The build handles metric units only. Will need to adjust source code and compile accordingly.
Paths to Folders will need to be adjusted accordingly.

To run the SBEACH module:
	1. Run createInput.py - This generates the input params used in SBEACH.
	2. Run Test.py to execute SBEACH in python.
	  	- Test.py runs an example from the demo SBEACH code downloaded from Veritech.
		- It also plots an animation.
		- Demo units are given in imperial units and are adjusted in Test.py to metric.
	3. Additional variables used in SBEACH are stored in globVars.py
	4. All docs describing algorithm and variables can be found in Docs folder.

To compile SBEACH:
	UNIX
	- f2py3 was used to compile SBEACH. This module is shipped anaconda. User will need to add /path/to/anacoda/bin to PATH.
	- Open Terminal change dir to folder containing fortran modules and pass f2py3 -c fortranfiles[AS_Processing.F90 AS_Memory.F90 AS_Errors.F90 Model.f90] -m model into the command line.
	- gfortran was used to compile the fortran files. Will need to change the random module in fortran if IFORT is used.
	WINDOWS
	- I compiled with winPython.
	- Copy all files: AS_Processing.F90 AS_Memory.F90 AS_Errors.F90 Model.f90 to C:\path\to\winpython\python-Version.amd64.
	- Open the WinPython Command prompt.
	- type: python Scripts\f2py.py -c fortranfiles[AS_Processing.F90 AS_Memory.F90 AS_Errors.F90 Model.f90] -m model into the command line and compile.
	- a model.pyd is then created in the C:\path\to\winpython\python-Version.amd64 directory.
	
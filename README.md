# SolidAngleViaGaussBonnet
c++ implemention of a method to compute solid angles from boundary curve data, via an immersed Gauss Bonnet theorem, and a dual cone construction

How to use me:

(1) go into Constants.h, and give a filename to read the boundary curve in. it expects a file named YOURNAMEHERE.txt - only insert YOURNAMEHERE into the code. If you want, change the grid settings. 
(2) recompile by typing "make", which runs the makefile
(3) its parallelised with open mp, so set a number of threads by typing "export OMP_NUM_THREADS=8" , for 8 threads say.
(4) run the executable SolidAngleGaussBonnet

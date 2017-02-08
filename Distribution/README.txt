Installing and running WEBMOD (version 1.0.0 - February 8, 2017):

Ensure first that you have a java runtime environment (JRE) current through version 8.
To check if you have java installed, enter 'java -version' at the Command Prompt. A recent version
would show 1.8.0_91. If you need to install or update your java go to 

http://java.sun.com/javase/downloads/index.jsp and download Java Runtime Environment (JRE) for Java SE 8u102

The WEBMOD Users manual is located in the .\doc\ directory

To run WEBMOD on the example problem for Andrews Creek, first test a hydrology only run. Change the 
directory to .\projects\Andrews_tutorial\, then double click on webmod.bat. The .\output directory should have
new output files after a short time.

The user can then run the GUI version by double clicking on .\projects\Andrews_tutorial\webmod_gui.bat and then
following the instructions in the User's manual “Quick Start Guide > Andrews Creek simulation and calibration”,

WEBMOD Directory Structure

/lib (All JARS)
/bin (webmod.exe - This version was compiled for a PC, let us know if you need a Unix version)
/doc (PDFs of manuals and documentation)
/projects (Directories for each watershed model)
  /Andrews (Example for Andrews Creek model described in Example Problems in the WEBMOD User's manual)
    /Input (Data and parameter files)
    /Output (Model output)
    /Control (Control files and name files created here) 
    Andrews.xlsm (Excel workbook with inputs, outputs, and source data
    webmod.bat (Runs WEBMOD for Andrews Creek model in batch mode)
    webmod_gui.bat (Runs WEBMOD using MMS Tool GUI for Andrews Creek in MMS GUI)
    webmod_paramtool.bat (Open Paramtool GUI populated with Andrews Creek dimensions and parameters)
    webmod_print.bat (Writes the name files to the /control directory that provide self documentation)
  /Andrews_tutorial (Interactive project included in WEBMOD User's Manual)
    same structure as /Andrews but lacking the Excel workbook.
  /DR2 (Example for DR2 watershed described in Example Problems in the WEBMOD User's manual)
    same structure as /Andrews but Excel workbook is named DR2.xlsm
    
  
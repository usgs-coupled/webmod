Installing and running WEBMOD (version 1.0.0 - May 1, 2011):

Ensure first that you have a java runtime environment (JRE) current through version 6, update 5.
To check if you have java installed, enter 'java -version' at the Command Prompt. A recent version
would show 1.6.0_21. If you need to install or update your java go to 

http://java.sun.com/javase/downloads/index.jsp and download Java Runtime Environment (JRE) 6 Update 5

To run WEBMOD on the example problem for Loch Vale, go to .\projects\LochVale\ and double click on
Loch.bat or Loch_gui.bat. Before you can run the GUI version, you will need File descriptors in the 
control directory. These are produced whenever the non-gui version (Loch.bat) is run or by running 
Loch_print.bat.


WEBMOD Directory Structure

/lib (All JARS)
/bin (PRMS2012.exe - This version was compiled for a PC, let us know if you need a Unix version)
/doc (PDFs of manuals and documentation)
/projects (Directories for each example or application)
  /LochVale (Example for Loch Vale project)
    /Input (Data and parameter files)
    /Output (Model output)
    /Control (Control files and name files created here) 
  Loch.bat (Runs PRMS for East River in batch mode)
  Loch_gui.bat (Runs PRMS for East River in MMS GUI)
  Loch_print.bat (Writes the name files to directory containing control file)
  Luca.bat (Runs Luca on the WEBMOD LochVale Model) 



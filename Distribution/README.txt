Installing and running WEBMOD (version 1.0.0 - August 18, 2016):

Ensure first that you have a java runtime environment (JRE) current through version 8.
To check if you have java installed, enter 'java -version' at the Command Prompt. A recent version
would show 1.8.0_91. If you need to install or update your java go to 

http://java.sun.com/javase/downloads/index.jsp and download Java Runtime Environment (JRE) for Java SE 8u102

To run WEBMOD on the example problem for Andrews Creek, go to .\projects\Andrews\ and double click on
webmod.bat or webmod_gui.bat. Before you can run the GUI version, you will need File descriptors in the 
control directory. These are produced whenever the non-gui version (webmod.bat) is run or by running 
webmod_print.bat.


WEBMOD Directory Structure

/lib (All JARS)
/bin (webmod.exe - This version was compiled for a PC, let us know if you need a Unix version)
/doc (PDFs of manuals and documentation)
/projects (Directories for each example or application)
  /Andrews (Example for Andrews Creek model described in Example Problems in the WEBMOD User's manual)
    /Input (Data and parameter files)
    /Output (Model output)
    /Control (Control files and name files created here) 
    webmod.bat (Runs WEBMOD for Andrews Creek model in batch mode)
    webmod_gui.bat (Runs WEBMOD using MMS Tool GUI for Andrews Creek in MMS GUI)
    webmod_paramtool.bat (Open Paramtool GUI populated with Andrews Creek dimensions and parameters)
    webmod_print.bat (Writes the name files to the /control directory that provide self documentation)
  /Andrews_tutorial (Interactive project included in WEBMOD User's Manual)
    same structure as /Andrews
  /DR2 (Example for DR2 watershed described in Example Problems in the WEBMOD User's manual)
    same structure as /Andrews
  
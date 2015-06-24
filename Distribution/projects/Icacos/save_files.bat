@echo off
IF "%1" == "" GOTO Usage
SET RunID=%1
set src_folder=.
set dst_folder=.\Results\%RunID%
set file_list=save_files.txt

if not exist "%dst_folder%" mkdir "%dst_folder%"

for /f "delims=" %%f in (%file_list%) do (
    xcopy "%src_folder%\%%f" "%dst_folder%\"
    )
GOTO Finish
:Usage
  ECHO   save_files RunID
  ECHO   saves files listed in save_file.txt to Results folder in directory RunID 
  pause
:Finish

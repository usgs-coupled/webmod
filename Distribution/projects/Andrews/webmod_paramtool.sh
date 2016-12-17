#!/bin/sh
JARDIR=../../lib
../../bin/webmod -C./control/webmod.control -print
java -cp $JARDIR/oui4.jar oui.paramtool.ParamTool ./input/webmod.params

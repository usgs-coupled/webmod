#
# Makefile for mmf
# markstro
# Tue Jul 26 08:43:10 MST 2005
#
# $Id: Makefile 4627 2008-10-01 16:48:11Z markstro $
#

MLIBDIR  = /raid1/models/lib

CC       = gcc
RM       = rm -f
AR       = ar ruv
RANLIB   = ranlib

#CFLAGS   = -g -Wall -Werror-implicit-function-declaration
CFLAGS   = -O -Wall -Werror-implicit-function-declaration

SRCS = 	mmf.c parse_args.c alloc_space.c build_lists.c \
	setup_cont.c decl_control.c control_addr.c \
	control_var.c read_params.c sort_dims.c sort_params.c sort_vars.c \
	var_addr.c declvar.c str_to_vals.c\
	declparam.c param_addr.c getdim.c timing.c getparam.c umalloc_etc.c \
	julday.c getvar.c dprint.c julconvert.c oprint.c readvar.c decldim.c \
	get_times.c batch_run.c read_control.c dim_addr.c reset_dim.c read_line.c \
	get_elem_add.c read_vars.c \
	save_vars.c save_params.c load_param.c check_vars.c \
	create_vstats.c free_vstats.c write_vstats.c uprint.c \
	stats.c call_modules.c call_setdims.c \
	read_datainfo.c putvar.c print_params.c print_vars.c \
	print_model_info.c batch_run_functions.c graph_single_run.c \
	control_array.c


MMSOBJS = ${SRCS:.c=.o}

.c.o:
	$(CC) $(CFLAGS) -c $<

#
# Private makefile defines
#
MMSLIB = $(MLIBDIR)/libmmf.a

#
# Standard Targets for Users
#

$(MMSLIB): $(MMSOBJS)
	$(AR) $(MMSLIB) $(MMSOBJS)
	$(RANLIB) $(MMSLIB)

clean:
	${RM} $(MMSLIB) *.o *.ln *~ *.a ${OBJS} core lint \#*

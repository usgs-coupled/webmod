# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc

# Include paths...
Debug_Include_Path=
Release_Include_Path=

# Library paths...
Debug_Library_Path=
Release_Library_Path=

# Additional libraries...
Debug_Libraries=
Release_Libraries=

# Preprocessor definitions...
Debug_Preprocessor_Definitions=-D GCC_BUILD -D _DEBUG -D _LIB -D _CRT_SECURE_NO_DEPRECATE -D _CRT_NONSTDC_NO_DEPRECATE 
Release_Preprocessor_Definitions=-D GCC_BUILD -D NDEBUG -D _LIB -D _CRT_SECURE_NO_DEPRECATE -D _CRT_NONSTDC_NO_DEPRECATE 

# Implictly linked object files...
Debug_Implicitly_Linked_Objects=
Release_Implicitly_Linked_Objects=

# Compiler flags...
Debug_Compiler_Flags=-O0 
Debug_Compiler_Flags=-O0 
Release_Compiler_Flags=-O2 
Release_Compiler_Flags=-O2 

# Builds all configurations for this project...
.PHONY: build_all_configurations
build_all_configurations: Debug Debug Release Release 

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders gccDebug/src/alloc_space.o gccDebug/src/batch_run.o gccDebug/src/batch_run_functions.o gccDebug/src/build_lists.o gccDebug/src/call_modules.o gccDebug/src/call_setdims.o gccDebug/src/check_vars.o gccDebug/src/control_addr.o gccDebug/src/control_array.o gccDebug/src/control_var.o gccDebug/src/create_vstats.o gccDebug/src/decl_control.o gccDebug/src/decldim.o gccDebug/src/declparam.o gccDebug/src/declvar.o gccDebug/src/dim_addr.o gccDebug/src/dprint.o gccDebug/src/free_vstats.o gccDebug/src/get_elem_add.o gccDebug/src/get_times.o gccDebug/src/getdim.o gccDebug/src/getdimname.o gccDebug/src/getparam.o gccDebug/src/getvar.o gccDebug/src/graph_single_run.o gccDebug/src/julconvert.o gccDebug/src/julday.o gccDebug/src/load_param.o gccDebug/src/mmf.o gccDebug/src/oprint.o gccDebug/src/param_addr.o gccDebug/src/parse_args.o gccDebug/src/print_model_info.o gccDebug/src/print_params.o gccDebug/src/print_vars.o gccDebug/src/putvar.o gccDebug/src/read_control.o gccDebug/src/read_datainfo.o gccDebug/src/read_line.o gccDebug/src/read_params.o gccDebug/src/read_vars.o gccDebug/src/readvar.o gccDebug/src/reset_dim.o gccDebug/src/save_params.o gccDebug/src/save_vars.o gccDebug/src/setup_cont.o gccDebug/src/sort_dims.o gccDebug/src/sort_params.o gccDebug/src/sort_vars.o gccDebug/src/stats.o gccDebug/src/str_to_vals.o gccDebug/src/timing.o gccDebug/src/umalloc_etc.o gccDebug/src/uprint.o gccDebug/src/var_addr.o gccDebug/src/write_vstats.o 
	ar rcs ../gccDebug/libmmf_c.a gccDebug/src/alloc_space.o gccDebug/src/batch_run.o gccDebug/src/batch_run_functions.o gccDebug/src/build_lists.o gccDebug/src/call_modules.o gccDebug/src/call_setdims.o gccDebug/src/check_vars.o gccDebug/src/control_addr.o gccDebug/src/control_array.o gccDebug/src/control_var.o gccDebug/src/create_vstats.o gccDebug/src/decl_control.o gccDebug/src/decldim.o gccDebug/src/declparam.o gccDebug/src/declvar.o gccDebug/src/dim_addr.o gccDebug/src/dprint.o gccDebug/src/free_vstats.o gccDebug/src/get_elem_add.o gccDebug/src/get_times.o gccDebug/src/getdim.o gccDebug/src/getdimname.o gccDebug/src/getparam.o gccDebug/src/getvar.o gccDebug/src/graph_single_run.o gccDebug/src/julconvert.o gccDebug/src/julday.o gccDebug/src/load_param.o gccDebug/src/mmf.o gccDebug/src/oprint.o gccDebug/src/param_addr.o gccDebug/src/parse_args.o gccDebug/src/print_model_info.o gccDebug/src/print_params.o gccDebug/src/print_vars.o gccDebug/src/putvar.o gccDebug/src/read_control.o gccDebug/src/read_datainfo.o gccDebug/src/read_line.o gccDebug/src/read_params.o gccDebug/src/read_vars.o gccDebug/src/readvar.o gccDebug/src/reset_dim.o gccDebug/src/save_params.o gccDebug/src/save_vars.o gccDebug/src/setup_cont.o gccDebug/src/sort_dims.o gccDebug/src/sort_params.o gccDebug/src/sort_vars.o gccDebug/src/stats.o gccDebug/src/str_to_vals.o gccDebug/src/timing.o gccDebug/src/umalloc_etc.o gccDebug/src/uprint.o gccDebug/src/var_addr.o gccDebug/src/write_vstats.o  $(Debug_Implicitly_Linked_Objects)

# Compiles file src/alloc_space.c for the Debug configuration...
-include gccDebug/src/alloc_space.d
gccDebug/src/alloc_space.o: src/alloc_space.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/alloc_space.c $(Debug_Include_Path) -o gccDebug/src/alloc_space.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/alloc_space.c $(Debug_Include_Path) > gccDebug/src/alloc_space.d

# Compiles file src/batch_run.c for the Debug configuration...
-include gccDebug/src/batch_run.d
gccDebug/src/batch_run.o: src/batch_run.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/batch_run.c $(Debug_Include_Path) -o gccDebug/src/batch_run.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/batch_run.c $(Debug_Include_Path) > gccDebug/src/batch_run.d

# Compiles file src/batch_run_functions.c for the Debug configuration...
-include gccDebug/src/batch_run_functions.d
gccDebug/src/batch_run_functions.o: src/batch_run_functions.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/batch_run_functions.c $(Debug_Include_Path) -o gccDebug/src/batch_run_functions.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/batch_run_functions.c $(Debug_Include_Path) > gccDebug/src/batch_run_functions.d

# Compiles file src/build_lists.c for the Debug configuration...
-include gccDebug/src/build_lists.d
gccDebug/src/build_lists.o: src/build_lists.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/build_lists.c $(Debug_Include_Path) -o gccDebug/src/build_lists.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/build_lists.c $(Debug_Include_Path) > gccDebug/src/build_lists.d

# Compiles file src/call_modules.c for the Debug configuration...
-include gccDebug/src/call_modules.d
gccDebug/src/call_modules.o: src/call_modules.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/call_modules.c $(Debug_Include_Path) -o gccDebug/src/call_modules.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/call_modules.c $(Debug_Include_Path) > gccDebug/src/call_modules.d

# Compiles file src/call_setdims.c for the Debug configuration...
-include gccDebug/src/call_setdims.d
gccDebug/src/call_setdims.o: src/call_setdims.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/call_setdims.c $(Debug_Include_Path) -o gccDebug/src/call_setdims.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/call_setdims.c $(Debug_Include_Path) > gccDebug/src/call_setdims.d

# Compiles file src/check_vars.c for the Debug configuration...
-include gccDebug/src/check_vars.d
gccDebug/src/check_vars.o: src/check_vars.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/check_vars.c $(Debug_Include_Path) -o gccDebug/src/check_vars.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/check_vars.c $(Debug_Include_Path) > gccDebug/src/check_vars.d

# Compiles file src/control_addr.c for the Debug configuration...
-include gccDebug/src/control_addr.d
gccDebug/src/control_addr.o: src/control_addr.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/control_addr.c $(Debug_Include_Path) -o gccDebug/src/control_addr.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/control_addr.c $(Debug_Include_Path) > gccDebug/src/control_addr.d

# Compiles file src/control_array.c for the Debug configuration...
-include gccDebug/src/control_array.d
gccDebug/src/control_array.o: src/control_array.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/control_array.c $(Debug_Include_Path) -o gccDebug/src/control_array.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/control_array.c $(Debug_Include_Path) > gccDebug/src/control_array.d

# Compiles file src/control_var.c for the Debug configuration...
-include gccDebug/src/control_var.d
gccDebug/src/control_var.o: src/control_var.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/control_var.c $(Debug_Include_Path) -o gccDebug/src/control_var.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/control_var.c $(Debug_Include_Path) > gccDebug/src/control_var.d

# Compiles file src/create_vstats.c for the Debug configuration...
-include gccDebug/src/create_vstats.d
gccDebug/src/create_vstats.o: src/create_vstats.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/create_vstats.c $(Debug_Include_Path) -o gccDebug/src/create_vstats.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/create_vstats.c $(Debug_Include_Path) > gccDebug/src/create_vstats.d

# Compiles file src/decl_control.c for the Debug configuration...
-include gccDebug/src/decl_control.d
gccDebug/src/decl_control.o: src/decl_control.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/decl_control.c $(Debug_Include_Path) -o gccDebug/src/decl_control.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/decl_control.c $(Debug_Include_Path) > gccDebug/src/decl_control.d

# Compiles file src/decldim.c for the Debug configuration...
-include gccDebug/src/decldim.d
gccDebug/src/decldim.o: src/decldim.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/decldim.c $(Debug_Include_Path) -o gccDebug/src/decldim.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/decldim.c $(Debug_Include_Path) > gccDebug/src/decldim.d

# Compiles file src/declparam.c for the Debug configuration...
-include gccDebug/src/declparam.d
gccDebug/src/declparam.o: src/declparam.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/declparam.c $(Debug_Include_Path) -o gccDebug/src/declparam.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/declparam.c $(Debug_Include_Path) > gccDebug/src/declparam.d

# Compiles file src/declvar.c for the Debug configuration...
-include gccDebug/src/declvar.d
gccDebug/src/declvar.o: src/declvar.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/declvar.c $(Debug_Include_Path) -o gccDebug/src/declvar.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/declvar.c $(Debug_Include_Path) > gccDebug/src/declvar.d

# Compiles file src/dim_addr.c for the Debug configuration...
-include gccDebug/src/dim_addr.d
gccDebug/src/dim_addr.o: src/dim_addr.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/dim_addr.c $(Debug_Include_Path) -o gccDebug/src/dim_addr.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/dim_addr.c $(Debug_Include_Path) > gccDebug/src/dim_addr.d

# Compiles file src/dprint.c for the Debug configuration...
-include gccDebug/src/dprint.d
gccDebug/src/dprint.o: src/dprint.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/dprint.c $(Debug_Include_Path) -o gccDebug/src/dprint.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/dprint.c $(Debug_Include_Path) > gccDebug/src/dprint.d

# Compiles file src/free_vstats.c for the Debug configuration...
-include gccDebug/src/free_vstats.d
gccDebug/src/free_vstats.o: src/free_vstats.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/free_vstats.c $(Debug_Include_Path) -o gccDebug/src/free_vstats.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/free_vstats.c $(Debug_Include_Path) > gccDebug/src/free_vstats.d

# Compiles file src/get_elem_add.c for the Debug configuration...
-include gccDebug/src/get_elem_add.d
gccDebug/src/get_elem_add.o: src/get_elem_add.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/get_elem_add.c $(Debug_Include_Path) -o gccDebug/src/get_elem_add.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/get_elem_add.c $(Debug_Include_Path) > gccDebug/src/get_elem_add.d

# Compiles file src/get_times.c for the Debug configuration...
-include gccDebug/src/get_times.d
gccDebug/src/get_times.o: src/get_times.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/get_times.c $(Debug_Include_Path) -o gccDebug/src/get_times.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/get_times.c $(Debug_Include_Path) > gccDebug/src/get_times.d

# Compiles file src/getdim.c for the Debug configuration...
-include gccDebug/src/getdim.d
gccDebug/src/getdim.o: src/getdim.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/getdim.c $(Debug_Include_Path) -o gccDebug/src/getdim.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/getdim.c $(Debug_Include_Path) > gccDebug/src/getdim.d

# Compiles file src/getdimname.c for the Debug configuration...
-include gccDebug/src/getdimname.d
gccDebug/src/getdimname.o: src/getdimname.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/getdimname.c $(Debug_Include_Path) -o gccDebug/src/getdimname.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/getdimname.c $(Debug_Include_Path) > gccDebug/src/getdimname.d

# Compiles file src/getparam.c for the Debug configuration...
-include gccDebug/src/getparam.d
gccDebug/src/getparam.o: src/getparam.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/getparam.c $(Debug_Include_Path) -o gccDebug/src/getparam.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/getparam.c $(Debug_Include_Path) > gccDebug/src/getparam.d

# Compiles file src/getvar.c for the Debug configuration...
-include gccDebug/src/getvar.d
gccDebug/src/getvar.o: src/getvar.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/getvar.c $(Debug_Include_Path) -o gccDebug/src/getvar.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/getvar.c $(Debug_Include_Path) > gccDebug/src/getvar.d

# Compiles file src/graph_single_run.c for the Debug configuration...
-include gccDebug/src/graph_single_run.d
gccDebug/src/graph_single_run.o: src/graph_single_run.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/graph_single_run.c $(Debug_Include_Path) -o gccDebug/src/graph_single_run.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/graph_single_run.c $(Debug_Include_Path) > gccDebug/src/graph_single_run.d

# Compiles file src/julconvert.c for the Debug configuration...
-include gccDebug/src/julconvert.d
gccDebug/src/julconvert.o: src/julconvert.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/julconvert.c $(Debug_Include_Path) -o gccDebug/src/julconvert.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/julconvert.c $(Debug_Include_Path) > gccDebug/src/julconvert.d

# Compiles file src/julday.c for the Debug configuration...
-include gccDebug/src/julday.d
gccDebug/src/julday.o: src/julday.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/julday.c $(Debug_Include_Path) -o gccDebug/src/julday.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/julday.c $(Debug_Include_Path) > gccDebug/src/julday.d

# Compiles file src/load_param.c for the Debug configuration...
-include gccDebug/src/load_param.d
gccDebug/src/load_param.o: src/load_param.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/load_param.c $(Debug_Include_Path) -o gccDebug/src/load_param.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/load_param.c $(Debug_Include_Path) > gccDebug/src/load_param.d

# Compiles file src/mmf.c for the Debug configuration...
-include gccDebug/src/mmf.d
gccDebug/src/mmf.o: src/mmf.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/mmf.c $(Debug_Include_Path) -o gccDebug/src/mmf.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/mmf.c $(Debug_Include_Path) > gccDebug/src/mmf.d

# Compiles file src/oprint.c for the Debug configuration...
-include gccDebug/src/oprint.d
gccDebug/src/oprint.o: src/oprint.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/oprint.c $(Debug_Include_Path) -o gccDebug/src/oprint.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/oprint.c $(Debug_Include_Path) > gccDebug/src/oprint.d

# Compiles file src/param_addr.c for the Debug configuration...
-include gccDebug/src/param_addr.d
gccDebug/src/param_addr.o: src/param_addr.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/param_addr.c $(Debug_Include_Path) -o gccDebug/src/param_addr.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/param_addr.c $(Debug_Include_Path) > gccDebug/src/param_addr.d

# Compiles file src/parse_args.c for the Debug configuration...
-include gccDebug/src/parse_args.d
gccDebug/src/parse_args.o: src/parse_args.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/parse_args.c $(Debug_Include_Path) -o gccDebug/src/parse_args.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/parse_args.c $(Debug_Include_Path) > gccDebug/src/parse_args.d

# Compiles file src/print_model_info.c for the Debug configuration...
-include gccDebug/src/print_model_info.d
gccDebug/src/print_model_info.o: src/print_model_info.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/print_model_info.c $(Debug_Include_Path) -o gccDebug/src/print_model_info.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/print_model_info.c $(Debug_Include_Path) > gccDebug/src/print_model_info.d

# Compiles file src/print_params.c for the Debug configuration...
-include gccDebug/src/print_params.d
gccDebug/src/print_params.o: src/print_params.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/print_params.c $(Debug_Include_Path) -o gccDebug/src/print_params.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/print_params.c $(Debug_Include_Path) > gccDebug/src/print_params.d

# Compiles file src/print_vars.c for the Debug configuration...
-include gccDebug/src/print_vars.d
gccDebug/src/print_vars.o: src/print_vars.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/print_vars.c $(Debug_Include_Path) -o gccDebug/src/print_vars.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/print_vars.c $(Debug_Include_Path) > gccDebug/src/print_vars.d

# Compiles file src/putvar.c for the Debug configuration...
-include gccDebug/src/putvar.d
gccDebug/src/putvar.o: src/putvar.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/putvar.c $(Debug_Include_Path) -o gccDebug/src/putvar.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/putvar.c $(Debug_Include_Path) > gccDebug/src/putvar.d

# Compiles file src/read_control.c for the Debug configuration...
-include gccDebug/src/read_control.d
gccDebug/src/read_control.o: src/read_control.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/read_control.c $(Debug_Include_Path) -o gccDebug/src/read_control.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/read_control.c $(Debug_Include_Path) > gccDebug/src/read_control.d

# Compiles file src/read_datainfo.c for the Debug configuration...
-include gccDebug/src/read_datainfo.d
gccDebug/src/read_datainfo.o: src/read_datainfo.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/read_datainfo.c $(Debug_Include_Path) -o gccDebug/src/read_datainfo.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/read_datainfo.c $(Debug_Include_Path) > gccDebug/src/read_datainfo.d

# Compiles file src/read_line.c for the Debug configuration...
-include gccDebug/src/read_line.d
gccDebug/src/read_line.o: src/read_line.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/read_line.c $(Debug_Include_Path) -o gccDebug/src/read_line.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/read_line.c $(Debug_Include_Path) > gccDebug/src/read_line.d

# Compiles file src/read_params.c for the Debug configuration...
-include gccDebug/src/read_params.d
gccDebug/src/read_params.o: src/read_params.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/read_params.c $(Debug_Include_Path) -o gccDebug/src/read_params.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/read_params.c $(Debug_Include_Path) > gccDebug/src/read_params.d

# Compiles file src/read_vars.c for the Debug configuration...
-include gccDebug/src/read_vars.d
gccDebug/src/read_vars.o: src/read_vars.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/read_vars.c $(Debug_Include_Path) -o gccDebug/src/read_vars.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/read_vars.c $(Debug_Include_Path) > gccDebug/src/read_vars.d

# Compiles file src/readvar.c for the Debug configuration...
-include gccDebug/src/readvar.d
gccDebug/src/readvar.o: src/readvar.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/readvar.c $(Debug_Include_Path) -o gccDebug/src/readvar.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/readvar.c $(Debug_Include_Path) > gccDebug/src/readvar.d

# Compiles file src/reset_dim.c for the Debug configuration...
-include gccDebug/src/reset_dim.d
gccDebug/src/reset_dim.o: src/reset_dim.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/reset_dim.c $(Debug_Include_Path) -o gccDebug/src/reset_dim.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/reset_dim.c $(Debug_Include_Path) > gccDebug/src/reset_dim.d

# Compiles file src/save_params.c for the Debug configuration...
-include gccDebug/src/save_params.d
gccDebug/src/save_params.o: src/save_params.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/save_params.c $(Debug_Include_Path) -o gccDebug/src/save_params.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/save_params.c $(Debug_Include_Path) > gccDebug/src/save_params.d

# Compiles file src/save_vars.c for the Debug configuration...
-include gccDebug/src/save_vars.d
gccDebug/src/save_vars.o: src/save_vars.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/save_vars.c $(Debug_Include_Path) -o gccDebug/src/save_vars.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/save_vars.c $(Debug_Include_Path) > gccDebug/src/save_vars.d

# Compiles file src/setup_cont.c for the Debug configuration...
-include gccDebug/src/setup_cont.d
gccDebug/src/setup_cont.o: src/setup_cont.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/setup_cont.c $(Debug_Include_Path) -o gccDebug/src/setup_cont.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/setup_cont.c $(Debug_Include_Path) > gccDebug/src/setup_cont.d

# Compiles file src/sort_dims.c for the Debug configuration...
-include gccDebug/src/sort_dims.d
gccDebug/src/sort_dims.o: src/sort_dims.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/sort_dims.c $(Debug_Include_Path) -o gccDebug/src/sort_dims.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/sort_dims.c $(Debug_Include_Path) > gccDebug/src/sort_dims.d

# Compiles file src/sort_params.c for the Debug configuration...
-include gccDebug/src/sort_params.d
gccDebug/src/sort_params.o: src/sort_params.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/sort_params.c $(Debug_Include_Path) -o gccDebug/src/sort_params.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/sort_params.c $(Debug_Include_Path) > gccDebug/src/sort_params.d

# Compiles file src/sort_vars.c for the Debug configuration...
-include gccDebug/src/sort_vars.d
gccDebug/src/sort_vars.o: src/sort_vars.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/sort_vars.c $(Debug_Include_Path) -o gccDebug/src/sort_vars.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/sort_vars.c $(Debug_Include_Path) > gccDebug/src/sort_vars.d

# Compiles file src/stats.c for the Debug configuration...
-include gccDebug/src/stats.d
gccDebug/src/stats.o: src/stats.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/stats.c $(Debug_Include_Path) -o gccDebug/src/stats.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/stats.c $(Debug_Include_Path) > gccDebug/src/stats.d

# Compiles file src/str_to_vals.c for the Debug configuration...
-include gccDebug/src/str_to_vals.d
gccDebug/src/str_to_vals.o: src/str_to_vals.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/str_to_vals.c $(Debug_Include_Path) -o gccDebug/src/str_to_vals.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/str_to_vals.c $(Debug_Include_Path) > gccDebug/src/str_to_vals.d

# Compiles file src/timing.c for the Debug configuration...
-include gccDebug/src/timing.d
gccDebug/src/timing.o: src/timing.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/timing.c $(Debug_Include_Path) -o gccDebug/src/timing.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/timing.c $(Debug_Include_Path) > gccDebug/src/timing.d

# Compiles file src/umalloc_etc.c for the Debug configuration...
-include gccDebug/src/umalloc_etc.d
gccDebug/src/umalloc_etc.o: src/umalloc_etc.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/umalloc_etc.c $(Debug_Include_Path) -o gccDebug/src/umalloc_etc.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/umalloc_etc.c $(Debug_Include_Path) > gccDebug/src/umalloc_etc.d

# Compiles file src/uprint.c for the Debug configuration...
-include gccDebug/src/uprint.d
gccDebug/src/uprint.o: src/uprint.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/uprint.c $(Debug_Include_Path) -o gccDebug/src/uprint.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/uprint.c $(Debug_Include_Path) > gccDebug/src/uprint.d

# Compiles file src/var_addr.c for the Debug configuration...
-include gccDebug/src/var_addr.d
gccDebug/src/var_addr.o: src/var_addr.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/var_addr.c $(Debug_Include_Path) -o gccDebug/src/var_addr.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/var_addr.c $(Debug_Include_Path) > gccDebug/src/var_addr.d

# Compiles file src/write_vstats.c for the Debug configuration...
-include gccDebug/src/write_vstats.d
gccDebug/src/write_vstats.o: src/write_vstats.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/write_vstats.c $(Debug_Include_Path) -o gccDebug/src/write_vstats.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/write_vstats.c $(Debug_Include_Path) > gccDebug/src/write_vstats.d

# Builds the Release configuration...
.PHONY: Release
Release: create_folders gccRelease/src/alloc_space.o gccRelease/src/batch_run.o gccRelease/src/batch_run_functions.o gccRelease/src/build_lists.o gccRelease/src/call_modules.o gccRelease/src/call_setdims.o gccRelease/src/check_vars.o gccRelease/src/control_addr.o gccRelease/src/control_array.o gccRelease/src/control_var.o gccRelease/src/create_vstats.o gccRelease/src/decl_control.o gccRelease/src/decldim.o gccRelease/src/declparam.o gccRelease/src/declvar.o gccRelease/src/dim_addr.o gccRelease/src/dprint.o gccRelease/src/free_vstats.o gccRelease/src/get_elem_add.o gccRelease/src/get_times.o gccRelease/src/getdim.o gccRelease/src/getdimname.o gccRelease/src/getparam.o gccRelease/src/getvar.o gccRelease/src/graph_single_run.o gccRelease/src/julconvert.o gccRelease/src/julday.o gccRelease/src/load_param.o gccRelease/src/mmf.o gccRelease/src/oprint.o gccRelease/src/param_addr.o gccRelease/src/parse_args.o gccRelease/src/print_model_info.o gccRelease/src/print_params.o gccRelease/src/print_vars.o gccRelease/src/putvar.o gccRelease/src/read_control.o gccRelease/src/read_datainfo.o gccRelease/src/read_line.o gccRelease/src/read_params.o gccRelease/src/read_vars.o gccRelease/src/readvar.o gccRelease/src/reset_dim.o gccRelease/src/save_params.o gccRelease/src/save_vars.o gccRelease/src/setup_cont.o gccRelease/src/sort_dims.o gccRelease/src/sort_params.o gccRelease/src/sort_vars.o gccRelease/src/stats.o gccRelease/src/str_to_vals.o gccRelease/src/timing.o gccRelease/src/umalloc_etc.o gccRelease/src/uprint.o gccRelease/src/var_addr.o gccRelease/src/write_vstats.o 
	ar rcs ../gccRelease/libmmf_c.a gccRelease/src/alloc_space.o gccRelease/src/batch_run.o gccRelease/src/batch_run_functions.o gccRelease/src/build_lists.o gccRelease/src/call_modules.o gccRelease/src/call_setdims.o gccRelease/src/check_vars.o gccRelease/src/control_addr.o gccRelease/src/control_array.o gccRelease/src/control_var.o gccRelease/src/create_vstats.o gccRelease/src/decl_control.o gccRelease/src/decldim.o gccRelease/src/declparam.o gccRelease/src/declvar.o gccRelease/src/dim_addr.o gccRelease/src/dprint.o gccRelease/src/free_vstats.o gccRelease/src/get_elem_add.o gccRelease/src/get_times.o gccRelease/src/getdim.o gccRelease/src/getdimname.o gccRelease/src/getparam.o gccRelease/src/getvar.o gccRelease/src/graph_single_run.o gccRelease/src/julconvert.o gccRelease/src/julday.o gccRelease/src/load_param.o gccRelease/src/mmf.o gccRelease/src/oprint.o gccRelease/src/param_addr.o gccRelease/src/parse_args.o gccRelease/src/print_model_info.o gccRelease/src/print_params.o gccRelease/src/print_vars.o gccRelease/src/putvar.o gccRelease/src/read_control.o gccRelease/src/read_datainfo.o gccRelease/src/read_line.o gccRelease/src/read_params.o gccRelease/src/read_vars.o gccRelease/src/readvar.o gccRelease/src/reset_dim.o gccRelease/src/save_params.o gccRelease/src/save_vars.o gccRelease/src/setup_cont.o gccRelease/src/sort_dims.o gccRelease/src/sort_params.o gccRelease/src/sort_vars.o gccRelease/src/stats.o gccRelease/src/str_to_vals.o gccRelease/src/timing.o gccRelease/src/umalloc_etc.o gccRelease/src/uprint.o gccRelease/src/var_addr.o gccRelease/src/write_vstats.o  $(Release_Implicitly_Linked_Objects)

# Compiles file src/alloc_space.c for the Release configuration...
-include gccRelease/src/alloc_space.d
gccRelease/src/alloc_space.o: src/alloc_space.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/alloc_space.c $(Release_Include_Path) -o gccRelease/src/alloc_space.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/alloc_space.c $(Release_Include_Path) > gccRelease/src/alloc_space.d

# Compiles file src/batch_run.c for the Release configuration...
-include gccRelease/src/batch_run.d
gccRelease/src/batch_run.o: src/batch_run.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/batch_run.c $(Release_Include_Path) -o gccRelease/src/batch_run.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/batch_run.c $(Release_Include_Path) > gccRelease/src/batch_run.d

# Compiles file src/batch_run_functions.c for the Release configuration...
-include gccRelease/src/batch_run_functions.d
gccRelease/src/batch_run_functions.o: src/batch_run_functions.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/batch_run_functions.c $(Release_Include_Path) -o gccRelease/src/batch_run_functions.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/batch_run_functions.c $(Release_Include_Path) > gccRelease/src/batch_run_functions.d

# Compiles file src/build_lists.c for the Release configuration...
-include gccRelease/src/build_lists.d
gccRelease/src/build_lists.o: src/build_lists.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/build_lists.c $(Release_Include_Path) -o gccRelease/src/build_lists.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/build_lists.c $(Release_Include_Path) > gccRelease/src/build_lists.d

# Compiles file src/call_modules.c for the Release configuration...
-include gccRelease/src/call_modules.d
gccRelease/src/call_modules.o: src/call_modules.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/call_modules.c $(Release_Include_Path) -o gccRelease/src/call_modules.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/call_modules.c $(Release_Include_Path) > gccRelease/src/call_modules.d

# Compiles file src/call_setdims.c for the Release configuration...
-include gccRelease/src/call_setdims.d
gccRelease/src/call_setdims.o: src/call_setdims.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/call_setdims.c $(Release_Include_Path) -o gccRelease/src/call_setdims.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/call_setdims.c $(Release_Include_Path) > gccRelease/src/call_setdims.d

# Compiles file src/check_vars.c for the Release configuration...
-include gccRelease/src/check_vars.d
gccRelease/src/check_vars.o: src/check_vars.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/check_vars.c $(Release_Include_Path) -o gccRelease/src/check_vars.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/check_vars.c $(Release_Include_Path) > gccRelease/src/check_vars.d

# Compiles file src/control_addr.c for the Release configuration...
-include gccRelease/src/control_addr.d
gccRelease/src/control_addr.o: src/control_addr.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/control_addr.c $(Release_Include_Path) -o gccRelease/src/control_addr.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/control_addr.c $(Release_Include_Path) > gccRelease/src/control_addr.d

# Compiles file src/control_array.c for the Release configuration...
-include gccRelease/src/control_array.d
gccRelease/src/control_array.o: src/control_array.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/control_array.c $(Release_Include_Path) -o gccRelease/src/control_array.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/control_array.c $(Release_Include_Path) > gccRelease/src/control_array.d

# Compiles file src/control_var.c for the Release configuration...
-include gccRelease/src/control_var.d
gccRelease/src/control_var.o: src/control_var.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/control_var.c $(Release_Include_Path) -o gccRelease/src/control_var.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/control_var.c $(Release_Include_Path) > gccRelease/src/control_var.d

# Compiles file src/create_vstats.c for the Release configuration...
-include gccRelease/src/create_vstats.d
gccRelease/src/create_vstats.o: src/create_vstats.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/create_vstats.c $(Release_Include_Path) -o gccRelease/src/create_vstats.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/create_vstats.c $(Release_Include_Path) > gccRelease/src/create_vstats.d

# Compiles file src/decl_control.c for the Release configuration...
-include gccRelease/src/decl_control.d
gccRelease/src/decl_control.o: src/decl_control.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/decl_control.c $(Release_Include_Path) -o gccRelease/src/decl_control.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/decl_control.c $(Release_Include_Path) > gccRelease/src/decl_control.d

# Compiles file src/decldim.c for the Release configuration...
-include gccRelease/src/decldim.d
gccRelease/src/decldim.o: src/decldim.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/decldim.c $(Release_Include_Path) -o gccRelease/src/decldim.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/decldim.c $(Release_Include_Path) > gccRelease/src/decldim.d

# Compiles file src/declparam.c for the Release configuration...
-include gccRelease/src/declparam.d
gccRelease/src/declparam.o: src/declparam.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/declparam.c $(Release_Include_Path) -o gccRelease/src/declparam.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/declparam.c $(Release_Include_Path) > gccRelease/src/declparam.d

# Compiles file src/declvar.c for the Release configuration...
-include gccRelease/src/declvar.d
gccRelease/src/declvar.o: src/declvar.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/declvar.c $(Release_Include_Path) -o gccRelease/src/declvar.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/declvar.c $(Release_Include_Path) > gccRelease/src/declvar.d

# Compiles file src/dim_addr.c for the Release configuration...
-include gccRelease/src/dim_addr.d
gccRelease/src/dim_addr.o: src/dim_addr.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/dim_addr.c $(Release_Include_Path) -o gccRelease/src/dim_addr.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/dim_addr.c $(Release_Include_Path) > gccRelease/src/dim_addr.d

# Compiles file src/dprint.c for the Release configuration...
-include gccRelease/src/dprint.d
gccRelease/src/dprint.o: src/dprint.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/dprint.c $(Release_Include_Path) -o gccRelease/src/dprint.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/dprint.c $(Release_Include_Path) > gccRelease/src/dprint.d

# Compiles file src/free_vstats.c for the Release configuration...
-include gccRelease/src/free_vstats.d
gccRelease/src/free_vstats.o: src/free_vstats.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/free_vstats.c $(Release_Include_Path) -o gccRelease/src/free_vstats.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/free_vstats.c $(Release_Include_Path) > gccRelease/src/free_vstats.d

# Compiles file src/get_elem_add.c for the Release configuration...
-include gccRelease/src/get_elem_add.d
gccRelease/src/get_elem_add.o: src/get_elem_add.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/get_elem_add.c $(Release_Include_Path) -o gccRelease/src/get_elem_add.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/get_elem_add.c $(Release_Include_Path) > gccRelease/src/get_elem_add.d

# Compiles file src/get_times.c for the Release configuration...
-include gccRelease/src/get_times.d
gccRelease/src/get_times.o: src/get_times.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/get_times.c $(Release_Include_Path) -o gccRelease/src/get_times.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/get_times.c $(Release_Include_Path) > gccRelease/src/get_times.d

# Compiles file src/getdim.c for the Release configuration...
-include gccRelease/src/getdim.d
gccRelease/src/getdim.o: src/getdim.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/getdim.c $(Release_Include_Path) -o gccRelease/src/getdim.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/getdim.c $(Release_Include_Path) > gccRelease/src/getdim.d

# Compiles file src/getdimname.c for the Release configuration...
-include gccRelease/src/getdimname.d
gccRelease/src/getdimname.o: src/getdimname.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/getdimname.c $(Release_Include_Path) -o gccRelease/src/getdimname.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/getdimname.c $(Release_Include_Path) > gccRelease/src/getdimname.d

# Compiles file src/getparam.c for the Release configuration...
-include gccRelease/src/getparam.d
gccRelease/src/getparam.o: src/getparam.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/getparam.c $(Release_Include_Path) -o gccRelease/src/getparam.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/getparam.c $(Release_Include_Path) > gccRelease/src/getparam.d

# Compiles file src/getvar.c for the Release configuration...
-include gccRelease/src/getvar.d
gccRelease/src/getvar.o: src/getvar.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/getvar.c $(Release_Include_Path) -o gccRelease/src/getvar.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/getvar.c $(Release_Include_Path) > gccRelease/src/getvar.d

# Compiles file src/graph_single_run.c for the Release configuration...
-include gccRelease/src/graph_single_run.d
gccRelease/src/graph_single_run.o: src/graph_single_run.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/graph_single_run.c $(Release_Include_Path) -o gccRelease/src/graph_single_run.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/graph_single_run.c $(Release_Include_Path) > gccRelease/src/graph_single_run.d

# Compiles file src/julconvert.c for the Release configuration...
-include gccRelease/src/julconvert.d
gccRelease/src/julconvert.o: src/julconvert.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/julconvert.c $(Release_Include_Path) -o gccRelease/src/julconvert.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/julconvert.c $(Release_Include_Path) > gccRelease/src/julconvert.d

# Compiles file src/julday.c for the Release configuration...
-include gccRelease/src/julday.d
gccRelease/src/julday.o: src/julday.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/julday.c $(Release_Include_Path) -o gccRelease/src/julday.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/julday.c $(Release_Include_Path) > gccRelease/src/julday.d

# Compiles file src/load_param.c for the Release configuration...
-include gccRelease/src/load_param.d
gccRelease/src/load_param.o: src/load_param.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/load_param.c $(Release_Include_Path) -o gccRelease/src/load_param.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/load_param.c $(Release_Include_Path) > gccRelease/src/load_param.d

# Compiles file src/mmf.c for the Release configuration...
-include gccRelease/src/mmf.d
gccRelease/src/mmf.o: src/mmf.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/mmf.c $(Release_Include_Path) -o gccRelease/src/mmf.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/mmf.c $(Release_Include_Path) > gccRelease/src/mmf.d

# Compiles file src/oprint.c for the Release configuration...
-include gccRelease/src/oprint.d
gccRelease/src/oprint.o: src/oprint.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/oprint.c $(Release_Include_Path) -o gccRelease/src/oprint.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/oprint.c $(Release_Include_Path) > gccRelease/src/oprint.d

# Compiles file src/param_addr.c for the Release configuration...
-include gccRelease/src/param_addr.d
gccRelease/src/param_addr.o: src/param_addr.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/param_addr.c $(Release_Include_Path) -o gccRelease/src/param_addr.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/param_addr.c $(Release_Include_Path) > gccRelease/src/param_addr.d

# Compiles file src/parse_args.c for the Release configuration...
-include gccRelease/src/parse_args.d
gccRelease/src/parse_args.o: src/parse_args.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/parse_args.c $(Release_Include_Path) -o gccRelease/src/parse_args.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/parse_args.c $(Release_Include_Path) > gccRelease/src/parse_args.d

# Compiles file src/print_model_info.c for the Release configuration...
-include gccRelease/src/print_model_info.d
gccRelease/src/print_model_info.o: src/print_model_info.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/print_model_info.c $(Release_Include_Path) -o gccRelease/src/print_model_info.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/print_model_info.c $(Release_Include_Path) > gccRelease/src/print_model_info.d

# Compiles file src/print_params.c for the Release configuration...
-include gccRelease/src/print_params.d
gccRelease/src/print_params.o: src/print_params.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/print_params.c $(Release_Include_Path) -o gccRelease/src/print_params.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/print_params.c $(Release_Include_Path) > gccRelease/src/print_params.d

# Compiles file src/print_vars.c for the Release configuration...
-include gccRelease/src/print_vars.d
gccRelease/src/print_vars.o: src/print_vars.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/print_vars.c $(Release_Include_Path) -o gccRelease/src/print_vars.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/print_vars.c $(Release_Include_Path) > gccRelease/src/print_vars.d

# Compiles file src/putvar.c for the Release configuration...
-include gccRelease/src/putvar.d
gccRelease/src/putvar.o: src/putvar.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/putvar.c $(Release_Include_Path) -o gccRelease/src/putvar.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/putvar.c $(Release_Include_Path) > gccRelease/src/putvar.d

# Compiles file src/read_control.c for the Release configuration...
-include gccRelease/src/read_control.d
gccRelease/src/read_control.o: src/read_control.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/read_control.c $(Release_Include_Path) -o gccRelease/src/read_control.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/read_control.c $(Release_Include_Path) > gccRelease/src/read_control.d

# Compiles file src/read_datainfo.c for the Release configuration...
-include gccRelease/src/read_datainfo.d
gccRelease/src/read_datainfo.o: src/read_datainfo.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/read_datainfo.c $(Release_Include_Path) -o gccRelease/src/read_datainfo.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/read_datainfo.c $(Release_Include_Path) > gccRelease/src/read_datainfo.d

# Compiles file src/read_line.c for the Release configuration...
-include gccRelease/src/read_line.d
gccRelease/src/read_line.o: src/read_line.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/read_line.c $(Release_Include_Path) -o gccRelease/src/read_line.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/read_line.c $(Release_Include_Path) > gccRelease/src/read_line.d

# Compiles file src/read_params.c for the Release configuration...
-include gccRelease/src/read_params.d
gccRelease/src/read_params.o: src/read_params.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/read_params.c $(Release_Include_Path) -o gccRelease/src/read_params.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/read_params.c $(Release_Include_Path) > gccRelease/src/read_params.d

# Compiles file src/read_vars.c for the Release configuration...
-include gccRelease/src/read_vars.d
gccRelease/src/read_vars.o: src/read_vars.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/read_vars.c $(Release_Include_Path) -o gccRelease/src/read_vars.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/read_vars.c $(Release_Include_Path) > gccRelease/src/read_vars.d

# Compiles file src/readvar.c for the Release configuration...
-include gccRelease/src/readvar.d
gccRelease/src/readvar.o: src/readvar.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/readvar.c $(Release_Include_Path) -o gccRelease/src/readvar.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/readvar.c $(Release_Include_Path) > gccRelease/src/readvar.d

# Compiles file src/reset_dim.c for the Release configuration...
-include gccRelease/src/reset_dim.d
gccRelease/src/reset_dim.o: src/reset_dim.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/reset_dim.c $(Release_Include_Path) -o gccRelease/src/reset_dim.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/reset_dim.c $(Release_Include_Path) > gccRelease/src/reset_dim.d

# Compiles file src/save_params.c for the Release configuration...
-include gccRelease/src/save_params.d
gccRelease/src/save_params.o: src/save_params.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/save_params.c $(Release_Include_Path) -o gccRelease/src/save_params.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/save_params.c $(Release_Include_Path) > gccRelease/src/save_params.d

# Compiles file src/save_vars.c for the Release configuration...
-include gccRelease/src/save_vars.d
gccRelease/src/save_vars.o: src/save_vars.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/save_vars.c $(Release_Include_Path) -o gccRelease/src/save_vars.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/save_vars.c $(Release_Include_Path) > gccRelease/src/save_vars.d

# Compiles file src/setup_cont.c for the Release configuration...
-include gccRelease/src/setup_cont.d
gccRelease/src/setup_cont.o: src/setup_cont.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/setup_cont.c $(Release_Include_Path) -o gccRelease/src/setup_cont.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/setup_cont.c $(Release_Include_Path) > gccRelease/src/setup_cont.d

# Compiles file src/sort_dims.c for the Release configuration...
-include gccRelease/src/sort_dims.d
gccRelease/src/sort_dims.o: src/sort_dims.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/sort_dims.c $(Release_Include_Path) -o gccRelease/src/sort_dims.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/sort_dims.c $(Release_Include_Path) > gccRelease/src/sort_dims.d

# Compiles file src/sort_params.c for the Release configuration...
-include gccRelease/src/sort_params.d
gccRelease/src/sort_params.o: src/sort_params.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/sort_params.c $(Release_Include_Path) -o gccRelease/src/sort_params.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/sort_params.c $(Release_Include_Path) > gccRelease/src/sort_params.d

# Compiles file src/sort_vars.c for the Release configuration...
-include gccRelease/src/sort_vars.d
gccRelease/src/sort_vars.o: src/sort_vars.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/sort_vars.c $(Release_Include_Path) -o gccRelease/src/sort_vars.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/sort_vars.c $(Release_Include_Path) > gccRelease/src/sort_vars.d

# Compiles file src/stats.c for the Release configuration...
-include gccRelease/src/stats.d
gccRelease/src/stats.o: src/stats.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/stats.c $(Release_Include_Path) -o gccRelease/src/stats.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/stats.c $(Release_Include_Path) > gccRelease/src/stats.d

# Compiles file src/str_to_vals.c for the Release configuration...
-include gccRelease/src/str_to_vals.d
gccRelease/src/str_to_vals.o: src/str_to_vals.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/str_to_vals.c $(Release_Include_Path) -o gccRelease/src/str_to_vals.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/str_to_vals.c $(Release_Include_Path) > gccRelease/src/str_to_vals.d

# Compiles file src/timing.c for the Release configuration...
-include gccRelease/src/timing.d
gccRelease/src/timing.o: src/timing.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/timing.c $(Release_Include_Path) -o gccRelease/src/timing.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/timing.c $(Release_Include_Path) > gccRelease/src/timing.d

# Compiles file src/umalloc_etc.c for the Release configuration...
-include gccRelease/src/umalloc_etc.d
gccRelease/src/umalloc_etc.o: src/umalloc_etc.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/umalloc_etc.c $(Release_Include_Path) -o gccRelease/src/umalloc_etc.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/umalloc_etc.c $(Release_Include_Path) > gccRelease/src/umalloc_etc.d

# Compiles file src/uprint.c for the Release configuration...
-include gccRelease/src/uprint.d
gccRelease/src/uprint.o: src/uprint.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/uprint.c $(Release_Include_Path) -o gccRelease/src/uprint.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/uprint.c $(Release_Include_Path) > gccRelease/src/uprint.d

# Compiles file src/var_addr.c for the Release configuration...
-include gccRelease/src/var_addr.d
gccRelease/src/var_addr.o: src/var_addr.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/var_addr.c $(Release_Include_Path) -o gccRelease/src/var_addr.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/var_addr.c $(Release_Include_Path) > gccRelease/src/var_addr.d

# Compiles file src/write_vstats.c for the Release configuration...
-include gccRelease/src/write_vstats.d
gccRelease/src/write_vstats.o: src/write_vstats.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/write_vstats.c $(Release_Include_Path) -o gccRelease/src/write_vstats.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/write_vstats.c $(Release_Include_Path) > gccRelease/src/write_vstats.d

# Creates the intermediate and output folders for each configuration...
.PHONY: create_folders
create_folders:
	mkdir -p gccDebug/src
	mkdir -p ../gccDebug
	mkdir -p gccRelease/src
	mkdir -p ../gccRelease

# Cleans intermediate and output files (objects, libraries, executables)...
.PHONY: clean
clean:
	rm -f gccDebug/src/*.o
	rm -f gccDebug/src/*.d
	rm -f ../gccDebug/*.a
	rm -f ../gccDebug/*.so
	rm -f ../gccDebug/*.dll
	rm -f ../gccDebug/*.exe
	rm -f gccRelease/src/*.o
	rm -f gccRelease/src/*.d
	rm -f ../gccRelease/*.a
	rm -f ../gccRelease/*.so
	rm -f ../gccRelease/*.dll
	rm -f ../gccRelease/*.exe

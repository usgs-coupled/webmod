#define EXTERNAL
#define MAIN
#define PHREEQC_IDENT
#include "global.h"
#include "phqalloc.h"
#include "message.h"
#undef PHREEQC_IDENT


/*     $Date: 2004/02/17 22:18:49 $ */
static char const rcsid[] = "$RCSfile: main.c,v $  $Revision: 2.29 $";

int main(int argc, char *argv[]);

#ifdef SKIP
static int echo_file_handler(const int type, const char *err_str, const int stop, void *cookie);
static int error_file_handler(const int type, const char *err_str, const int stop, void *cookie);
static int log_file_handler(const int type, const char *err_str, const int stop, void *cookie);
static int output_handler(const int type, const char *err_str, const int stop, void *cookie);
static int stop_handler(const int type, const char *err_str, const int stop, void *cookie);
#endif
static int process_file_names(int argc, char *argv[]);
#ifdef DOS
static int write_banner(void);
#endif
extern int advection(void);
extern int backspace (FILE *file, int spaces);
extern int clean_up(void);
extern int copy_entities(void);
extern int copy_token (char *token_ptr, char **ptr, int *length);
extern int dup_print(const char *ptr, int emphasis);
extern int error_msg (const char *err_str, const int stop);
extern FILE *file_open (char *query, char *default_name, const char *status, int batch);
void *free_check_null(void *ptr);
extern int get_line( FILE *fp);
extern void initialize(void);
extern int initial_exchangers(int print);
extern int initial_gas_phases(int print);
extern int initial_solutions(int print);
extern int initial_surfaces(int print);
extern int inverse_models(void);
extern int reactions(void);
extern int read_input(void);
extern int replace(const char *str1, const char *str2, char *str);
extern int status(int count, const char *str);
extern int strcmp_nocase(const char *str1, const char *str2);
extern char * string_duplicate (const char *token);
extern int string_trim (char *str);
extern int tidy_model(void);
extern int transport(void);
extern int warning_msg (const char *err_str);

/* ----------------------------------------------------------------------
 *   MAIN
 * ---------------------------------------------------------------------- */
int main(int argc, char *argv[])
/*
 *   Main program for PHREEQC
 */
{
	char token[MAX_LENGTH];
	if (rcsid == NULL) fprintf(stderr," ");
/*
 *   Open files
 */
	if (argc > 4) {
		error_file = fopen(argv[4], "w");
	} else {
		error_file = stderr;	/*  fopen("screen.out", "w"); */
	}
#       ifdef DOS
               write_banner(); 
#       endif
	state = INITIALIZE;
	initialize();
	phast = FALSE;
/*
 *   Add callbacks for error_msg and warning_msg
 */
	add_message_callback(default_handler, NULL); 

	process_file_names(argc, argv);
/*
 *   Use to cause output to be completely unbuffered
 */
/*	setbuf(output, NULL); */
/*
 *   Use to cause log_file to be completely unbuffered
 */
/*	setbuf(log_file, NULL); */
/*
 *   Initialize arrays
 */
	simulation = 0;
/*
 *   Read data base
 */
#ifdef PHREEQCI_GUI
	_ASSERT(pr.headings == TRUE);
#endif
	dup_print("Reading data base.", TRUE);
	input = database_file;
	read_input();
	tidy_model();
	fclose(database_file);
	database_file = NULL;
	status(0, NULL);
/*
 *   Read input data for simulation
 */
	input = input_file;
	for (simulation = 1; ; simulation++) {

		sprintf(token, "Reading input data for simulation %d.", simulation);

		output_msg(OUTPUT_GUI_ERROR, "\nSimulation %d\n", simulation);

		dup_print(token, TRUE);
		if (read_input() == EOF) break;
		if (title_x != NULL) {
			sprintf(token, "TITLE");
			dup_print(token, TRUE);
			if (pr.headings == TRUE) output_msg(OUTPUT_MESSAGE,"%s\n\n", title_x);
		}
		tidy_model();
/*
 *   Calculate distribution of species for initial solutions
 */
		if (new_solution) initial_solutions(TRUE);
/*
 *   Calculate distribution for exchangers
 */
		if (new_exchange) initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
		if (new_surface) initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
		if (new_gas_phase) initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
		reactions();
/*
 *   Calculate inverse models
 */
		inverse_models();
/*
 *   Calculate advection
 */
		if (use.advect_in == TRUE) {
			dup_print ("Beginning of advection calculations.", TRUE);
			advection(); 
		}
/*
 *   Calculate transport
 */
		if (use.trans_in == TRUE) {
			dup_print ("Beginning of transport calculations.", TRUE);
			transport(); 
		}
/*
 *   Copy
 */
		if (new_copy) copy_entities();
/*
 *   End of simulation
 */
		dup_print( "End of simulation.", TRUE);
	}
	if (pr.status == TRUE) {
#if defined(PHREEQCI_GUI)
		state = -1;
		status(0, "\r\nDone.");
#else
		status(0, "\nDone.");
#endif
		output_msg(OUTPUT_SCREEN,"\n");
#if defined(PHREEQCI_GUI)
		SendMessage(g_status.hText, EM_REPLACESEL, (WPARAM)0, (LPARAM)"\r\n");
#endif
	}
	dup_print ("End of run.", TRUE);
	output_msg(OUTPUT_SCREEN,"\nEnd of Run.\n");
#if defined(PHREEEQCI_GUI)
	sprintf(g_szLineBuf, "\r\nEnd of Run.\r\n"); 
	SendMessage(g_status.hText, EM_REPLACESEL, (WPARAM)0, (LPARAM)g_szLineBuf);
#endif
	clean_up();
#if !defined(_WINDOWS)
	exit(0);
#endif
	return(0);
}
/* ---------------------------------------------------------------------- */
int process_file_names(int argc, char *argv[])
/* ---------------------------------------------------------------------- */
{
	int l;
	char token[2*MAX_LENGTH], default_name[2*MAX_LENGTH];
	char query[2*MAX_LENGTH];
	char in_file[2*MAX_LENGTH], out_file[2*MAX_LENGTH];
	char *env_ptr;
	char *ptr;
/*
 *   Open user-input file
 */
	strcpy(query,"Name of input file?");
	if (argc <= 1) {
		default_name[0]='\0';
		input_file = file_open(query, default_name, "r", FALSE);
	} else {
		strcpy(default_name, argv[1]);
		input_file = file_open(query, default_name, "r", TRUE);
	}
	output_msg(OUTPUT_SCREEN, "Input file: %s\n\n", default_name);
#if defined(PHREEQCI_GUI)
	sprintf(g_szLineBuf, "Input file: %s\r\n\r\n", default_name); 
	SendMessage(g_status.hText, EM_REPLACESEL, (WPARAM)0, (LPARAM)g_szLineBuf);
#endif
	strcpy(in_file, default_name);
/*
 *   Open file for output
 */
	strcpy(query,"Name of output file?");
#ifdef DOS
	replace("."," ",default_name);
#endif
	ptr = default_name;
	copy_token(token, &ptr, &l);
	strcat(token,".out");
	if (argc <= 1) {
		output = file_open(query, token, "w", FALSE);
	} else if (argc == 2) {
		output = file_open(query, token, "w", TRUE);
	} else if (argc >= 3) {
		strcpy(token, argv[2]);
		output = file_open(query, token, "w", TRUE);
	}		
	output_msg(OUTPUT_SCREEN, "Output file: %s\n\n", token);
#if defined(PHREEQCI_GUI)
	sprintf(g_szLineBuf, "Output file: %s\r\n\r\n", token); 
	SendMessage(g_status.hText, EM_REPLACESEL, (WPARAM)0, (LPARAM)g_szLineBuf);
#endif
	strcpy(out_file, token);
/*
 *   Open file for errors
 */
	if ((log_file = fopen("phreeqc.log","w")) == NULL) {
		error_msg ("Can't open log file, phreeqc.log.", STOP);
	}
	/*
	 *  Read input file for DATABASE keyword
	 */
	user_database = NULL;
	if (get_line(input_file) == KEYWORD) {
		ptr = line;
		copy_token(token, &ptr, &l);
		if (strcmp_nocase(token, "database") == 0) {
			user_database = string_duplicate(ptr);
			if (string_trim(user_database) == EMPTY) {
				warning_msg("DATABASE file name is missing; default database will be used.");
					user_database = NULL;
			}
		}
	}
	fclose(input_file);
	if ((input_file = fopen(in_file,"r")) == NULL) {;
		error_msg ("Can't reopen input file.", STOP);
	}
/*
 *   Open data base
 */
	strcpy(query,"Name of database file?");
	env_ptr = getenv("PHREEQC_DATABASE");
	if (user_database != NULL) {
		strcpy(token, user_database);
	} else if (env_ptr != NULL) {
		strcpy(token, env_ptr);
	} else {
		strcpy(token, default_data_base);
	} 
	if (argc <= 1) {
		database_file = file_open(query, token, "r", FALSE);
	} else if (argc < 4) {
		database_file = file_open(query, token, "r", TRUE);
	} else if (argc >= 4) {
		if (user_database == NULL) {
			strcpy(token, argv[3]);
		} else {
#ifndef PHREEQCI_GUI
			warning_msg("Database file from DATABASE keyword is used; command line argument ignored.");
#endif
		}
		database_file = file_open(query, token, "r", TRUE);
	}
	output_msg(OUTPUT_SCREEN, "Database file: %s\n\n", token);
#if defined(PHREEQCI_GUI)
	sprintf(g_szLineBuf, "Database file: %s\r\n\r\n", token); 
	SendMessage(g_status.hText, EM_REPLACESEL, (WPARAM)0, (LPARAM)g_szLineBuf);
#endif


	output_msg(OUTPUT_MESSAGE, "   Input file: %s\n", in_file);
	output_msg(OUTPUT_MESSAGE, "  Output file: %s\n", out_file);
	output_msg(OUTPUT_MESSAGE, "Database file: %s\n\n", token);

	user_database = free_check_null(user_database);
	return(OK);
}
/* ---------------------------------------------------------------------- */
int write_banner(void)
/* ---------------------------------------------------------------------- */
{
output_msg(OUTPUT_SCREEN, "              €ﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂ€\n");
output_msg(OUTPUT_SCREEN, "              ∫                                            ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫              * PHREEQC-2.8 *               ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫                                            ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫      A hydrogeochemical transport model    ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫                                            ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫                     by                     ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫       D.L. Parkhurst and C.A.J. Appelo     ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫                                            ∫\n");
output_msg(OUTPUT_SCREEN, "              ∫               April 15, 2003               ∫\n");
output_msg(OUTPUT_SCREEN, "              €‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹€\n\n");

return 0;
}

// phreeqc namespace
//

#if !defined(_INC_GLOBALNS)
#define _INC_GLOBALNS

#if defined(_DEBUG)
#pragma warning(disable : 4786)   // disable truncation warning
#endif

#include <math.h>

namespace phreeqc {
extern "C" {

#define EXTERNAL extern
#include "phreeqc/global.h"
#include "phreeqc/phqalloc.h"
#include "phreeqc/message.h"

void malloc_error(void);
int check_key (char *str);
int copy_entities(void);
int copy_token (char *token_ptr, char **ptr, int *l);

int clean_up(void);
void initialize(void);
int initial_exchangers(int print);
int initial_gas_phases(int print);
int initial_solutions(int print);
int initial_surfaces(int print);
int inverse_models(void);
int reactions(void);
int read_input(void);
int tidy_model(void);
int string_trim_right(char *str);
int tidy_punch(void);
int dup_print(const char *ptr, int emphasis);

extern int n_user_punch_index;
}

}

#endif /* _INC_GLOBALNS */


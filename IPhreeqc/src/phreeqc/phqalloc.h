#define USE_PHRQ_ALLOC
#if !defined (INCLUDE_PHRQALLOC_H)
#define INCLUDE_PHRQALLOC_H

#if defined (USE_PHRQ_ALLOC)

/*     $Date: 2002/04/25 19:01:01 $ */
#ifdef PHREEQC_IDENT
static char const rcsidphqalloc[] = "$RCSfile: phqalloc.h,v $  $Revision: 2.2 $";
#endif

#if !defined(NDEBUG)
extern void *PHRQ_malloc(size_t, const char *, int);
extern void *PHRQ_calloc(size_t, size_t, const char *, int);
extern void *PHRQ_realloc(void *, size_t, const char *, int);
#else
extern void *PHRQ_malloc(size_t);
extern void *PHRQ_calloc(size_t, size_t);
extern void *PHRQ_realloc(void *, size_t);
#endif

extern void PHRQ_free(void *);
extern void PHRQ_free_all(void);

#if !defined(NDEBUG)
#define   PHRQ_malloc(s)         PHRQ_malloc(s, __FILE__, __LINE__)
#define   PHRQ_calloc(c, s)      PHRQ_calloc(c, s, __FILE__, __LINE__)
#define   PHRQ_realloc(p, s)     PHRQ_realloc(p, s, __FILE__, __LINE__)
#endif

#define   free(p)                PHRQ_free(p)

#else  /* defined (USE_PHRQ_ALLOC) */

#define PHRQ_malloc malloc
#define PHRQ_realloc realloc
#define PHRQ_calloc calloc

#endif /* defined (USE_PHRQ_ALLOC) */

#endif /* !defined (INCLUDE_PHRQALLOC_H) */

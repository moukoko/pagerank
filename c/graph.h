#ifndef __GRAPH__H__
#define __GRAPH__H__

#include <stdint.h>


#define DEFAULT_ALPHA             0.85
#define DEFAULT_CONVERGENCE       1e-05
#define DEFAULT_MAX_ITERATIONS    10000UL
#define DEFAULT_NUMERIC           0
#define DEFAULT_DELIM             " "

struct links {
   size_t num_links;
   size_t size;
   size_t *p;
};

struct graph {
  double alpha;
  double convergence;
  uint64_t max_iterations;
  char *delim;
  int numeric;
  size_t size;
  double *pr;
  size_t *num_outgoing;
  struct links *pages;
};

void graph_init(struct graph **graph, size_t size);
void graph_read_file(const char *filename, struct graph *g);
void graph_pagerank(struct graph *graph);
void graph_print(const char *filename, struct graph *graph);
void graph_print_pagerank(const char *filename, struct graph *graph);


#endif /* __GRAPH__H__ */

#include <numa.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "graph.h"

/* Initialize graph. Allocate memory for internal 
 * arrays */
void graph_init(struct graph **graph, size_t size)
{
  if (!size)
    *graph = NULL;

  *graph = malloc(sizeof(struct graph));
  if (!*graph) {
    perror("Could not allocate memory for graph\n");
  }

  (*graph)->size = size;
  (*graph)->delim = " ";
  (*graph)->max_iterations = DEFAULT_MAX_ITERATIONS;
  (*graph)->convergence = DEFAULT_CONVERGENCE;
  (*graph)->alpha = DEFAULT_ALPHA;
  (*graph)->pr = numa_alloc_onnode(size * sizeof(double), 0);
  (*graph)->pages = numa_alloc_onnode(size * sizeof(struct links), 0);
  memset((*graph)->pages, 0, sizeof(struct links) * size);
  (*graph)->num_outgoing = numa_alloc_onnode(size * sizeof(size_t), 0);
}

/* Compute the pagerank values for 'graph' */
void graph_pagerank(struct graph *graph)
{
  double diff, sum_pr, dangling_pr, *old_pr, *pr, cpr;
  size_t i, k, j, size, out;
  unsigned long num_iterations = 0;
  double one_Av, one_Iv, h, h_v;


  pr = graph->pr;
  size = graph->size;
  old_pr = malloc(size * sizeof(double));
  if (!old_pr) {
    perror("Could not allocate data structure");
    exit(1);
  }
  
  pr[0] = 1;
  
  printf("Running PageRank. Up to %lu iterations, until: %e error\n", graph->max_iterations, graph->convergence);
  do {
    sum_pr = 0.0;
    dangling_pr = 0.0;

    for (k = 0; k < size; ++k) {
      cpr = pr[k];
      sum_pr += cpr;
      if (graph->num_outgoing[k] == 0)
        dangling_pr += cpr;
    }

    if (num_iterations == 0)
      memcpy(old_pr, pr, graph->size * sizeof(double));
    else {
      /* Normalize so that we start with sum equal to one */
      for (i = 0; i < size; ++i)
        old_pr[i] = pr[i] / sum_pr;
    }

    /*
     * After normalization the elements of the pagerank vector sum
     * to one
     */
    sum_pr = 1;

    /* An element of the A x I vector; all elements are identical */
    one_Av = graph->alpha * dangling_pr / size;

    /* An element of the 1 x I vector; all elements are identical */
    one_Iv = (1 - graph->alpha) * sum_pr / size;

    /* The difference to be checked for convergence */
    diff = 0;
    for (i = 0; i < size; ++i) {
      /* The corresponding element of the H multiplication */
      h = 0;
      for (j = 0; j < graph->pages[i].num_links; ++j) {
        out = graph->pages[i].p[j];
        h_v = (graph->num_outgoing[out]) ? 1.0 / graph->num_outgoing[out] : 0.0;
        h += h_v * old_pr[out];
      }
      h *= graph->alpha;
      pr[i] = h + one_Av + one_Iv;
      diff += fabs(pr[i] - old_pr[i]);
    }
    num_iterations++;
  } while (diff > graph->convergence && num_iterations < graph->max_iterations);
  printf("Ran %lu iterations. Achieved %e error\n", num_iterations, diff);
}

/* Add a link in the graph 'g' from 'from' to 'to' */
static void _graph_add_link(struct graph *g, size_t from, size_t to)
{
  size_t to_out;
  struct links *pages = g->pages;
  
  g->num_outgoing[from]++;
  to_out = g->pages[to].num_links;
  
  if (pages[to].size == 0) {
    pages[to].p = numa_alloc_onnode(10, 0);
    if(!pages[to].p) {
      perror("Could not allocate array");
      exit(1);
    }
    pages[to].size = 10;
  }

  if (pages[to].num_links + 1 > pages[to].size) {
    pages[to].p = numa_realloc(pages[to].p, pages[to].size, 2*pages[to].size);
    if (!pages[to].p) {
      perror("Could not grow array");
      exit(1);
    }
    pages[to].size *= 2;
  }

  pages[to].p[to_out] = from;
  pages[to].num_links++;
}

/* parse a single line from file to get a link.
 * File line is of form:
 * "'from''delim''to'" */
/* TODO: Currently input file is assumed to contain numerical ids.
 * Change this to also handle alpharithmetical */
static void _graph_parse_link(struct graph *g, char *line, size_t len)
{
  char *ret;
  size_t from, to;
  
  line[len-3] = '\0';
  ret = strtok(line, g->delim);
  if (!ret) {
    fprintf(stderr, "Maeormed input file\n");
    exit(1);
  }

  from = strtoul(ret, NULL, 10);

  ret = strtok(NULL, g->delim);
  if (!ret) {
    fprintf(stderr, "Maeormed input file\n");
    exit(1);
  }

  to = strtoul(ret, NULL, 10);
  _graph_add_link(g, from, to);
}

void graph_read_file(const char *filename, struct graph *g)
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(filename, "r");
  if (!fp) {
    perror("Could not read file");
    exit(1);
  }

  while ((read = getline(&line, &len, fp)) != -1) {
    _graph_parse_link(g, line, len);
  }
}

void graph_print(const char *filename, struct graph *g)
{
  size_t i, j, num_in;
  FILE *fp;

  if (!filename)
    fp = stdout;
  else {
    fp = fopen(filename, "w");
    if (!fp) {
      perror("Error opening file for printing graph");
      exit(1);
    }
  }

  fprintf(fp, "alpha: %e\n", g->alpha);
  fprintf(fp, "delimiter: \"%s\"\n", g->delim);
  fprintf(fp, "Number of nodes: %lu\n", g->size);
  for (i = 0; i < g->size; ++i) {
    fprintf(fp, "NODE %lu\n", i);
    fprintf(fp, "NUM_OUTGOING %lu\n", g->num_outgoing[i]);
    num_in = g->pages[i].num_links;
    fprintf(fp, "INCOMING NODES: %lu\n", num_in);
    for (j = 0; j < num_in; ++j)
      fprintf(fp, "%lu%c", g->pages[i].p[j], (j == num_in - 1) ? '\n' : ' ');
  }
}

void graph_print_pagerank(const char *filename, struct graph *g)
{
  size_t i;
  FILE *fp;

  if (!filename)
    fp = stdout;
  else {
    fp = fopen(filename, "w");
    if (!fp) {
      perror("Error opening file for printing graph");
      exit(1);
    }
  }

  for (i = 0; i < g->size; ++i)
    fprintf(fp, "%lu: %e\n", i, g->pr[i]);

  if (filename)
    fclose(fp);
}

#include <stdio.h>
#include <stdlib.h>

#include "graph.h"

void print_usage(void)
{
  fprintf(stderr, "usage: pagerank infile size [outfile]\n");
}



int main(int argc, char *argv[])
{
  struct graph *g;
  size_t size;

  if (argc < 3) {
    print_usage();
    exit(1);
  }
  
  size = atol(argv[2]);

  graph_init(&g, size);
//  graph_print(NULL, g);
  graph_read_file(argv[1], g);
  graph_print(NULL, g);
  graph_pagerank(g);
  graph_print_pagerank(NULL, g);

  return 0;
}

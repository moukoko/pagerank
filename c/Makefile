CC=gcc
CFLAGS=-O3 -g -ggdb -Wall
LFLAGS=
LIBS=-lnuma -lm
SRC=graph.c pagerank.c
OBJ=${SRC:.c=.o}
BIN=pagerank

all: ${BIN}

pagerank: ${OBJ}
	mkdir -p bin
	${CC} ${LFLAGS} $^ -o ./bin/$@ ${LIBS}

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	${RM} ${OBJ} bin/*

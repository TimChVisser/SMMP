CURR_DIR=$(notdir $(basename $(shell pwd)))
PRJ=$(CURR_DIR)
SRC=$(wildcard *.c)
OBJ=$(patsubst %.c,%.o,$(SRC))

CC=gcc
INCLUDES=-I../../include
ifndef DEBUG
CFLAGS=-O2 -g -std=c99 -ffast-math -ftree-vectorize  -ffast-math -mavx
LIB=
else
CFLAGS=-O0 -g3 -std=c99
LIB=
endif
#  -fopt-info-vec-missed
all: $(PRJ)

$(PRJ): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJ) -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@ $(LIB)

clean:
	-rm -f $(OBJ) $(PRJ)

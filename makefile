#!/bin/bash

# Copyright 2014 Georgios Karagiannis
#
# Georgios Karagiannis
#
# Postdoctoral research associate
# Department of Mathematics, Purdue University
# 150 N. University Street
# West Lafayette, IN 47907-2067, USA
#
# Telephone: +1 765 494-3405
#
# Email: gkaragia@purdue.edu
#
# Contact email: georgios.stats@gmail.com

CC=icc
CFLAGS=-O2
LDFLAGS=
CPPFLAGS=

FUN=cost_protein3D.c

SOURCES=pisaa.c \
	Crossover_operations.c \
	Mutation_operations.c \
	HitAndRun_update.c \
	Self_adjastment_prosedure.c \
	$(FUN) \
	uniformdirectionrng.c \
	normalrng.c \
	uniformrng.c \
	nrutil.c 
	
OBJECTS=$(SOURCES:.c=.o)

EXECUTABLE=exe

# ACTIONS ###############

build: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# CLEAR

clean:
	rm -rf *o exe

# DETAILS

details:
	@echo CC       : $(CC)
	@echo CFLAGS   : $(CFLAGS)
	@echo LDFLAGS  : $(LDFLAGS)
	@echo FUN      : $(FUN)
	@echo CPPFLAGS : $(CPPFLAGS)

# RUN

run:
	./exe

test_run:
	./exe -ID 3 \
		-Nmonomer 55 \
		-Data ./data.ABseq.55 \
		-Niter 20000000 \
		-Npop 30 \
		-Nsam 2000 \
		-Gwarm 1000000 \
		-Ghigh 1.0 \
		-Gpow 0.55 \
		-Hlow -170.0 \
		-Hhigh 0.0 \
		-Hsize $(shell echo "scale=0; ( 0.0-(-165.0))/0.1+1;" | bc) \
		-Hzeta 0.0 \
		-Hconst 100.0 \
		-Twarm 1 \
		-Tlow 0.01 \
		-Thigh 5.0 \
		-Tpow 0.5 \
		-Tini 100.0 \
		-Tref 0.00001 \
		-SMO0 0.5 \
		-SMO1 $(shell echo "scale=6;sqrt(2.0);" | bc) \
		-SMO2 1.0 \
		-SMO3 1.0 \
		-SCO1 0.1 \
		-SCO2 0.5 \
		-SCO3 0.5 \
		-Sini 1.0 \
		-Sref 0.0001



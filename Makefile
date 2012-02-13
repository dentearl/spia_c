SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
##############################
# spia makefile
#  dent earl, dearl soe ucsc edu
#  1 may 2009
##############################
binPath = ./bin
srcPath = ./src

CC = gcc 
CFLAGS = -Wall -Werror -std=c99 -pedantic -g

.PHONY: all clean test

all: ${binPath}/spia

${binPath}/spia: ${srcPath}/spia.c ${srcPath}/spia.h ${srcPath}/fileHandling.o ${srcPath}/hashListFunctions.o ${srcPath}/linearAlg.o ${srcPath}/probFunctions.o
	mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@.tmp $^ -I external/uthash-1.5/src/ -I77 -llapack -lblas
	mv $@.tmp $@

${srcPath}/%.o: ${srcPath}/%.c ${srcPath}/spia.h
	${CC} ${CFLAGS} -c -o $@.tmp $< -I external/uthash-1.5/src/
	mv $@.tmp $@

##############################
# system tests
test: test/out/small.txt test/out/testBasic.txt test/out/testACGT.txt test/out/testSPath.txt test/out/testBeta.txt test/out/testD.txt 

test/out/small.txt: ${binPath}/spia
	mkdir -p $(dir $@)
	${binPath}/spia --pathFiles test/altTest/paths/ --de test/altTest/de.txt --nBoots 1000 \
		--array test/altTest/all.txt --verbose 2&> $@.tmp
	mv $@.tmp $@

test/out/testBasic.txt: ${binPath}/spia
	mkdir -p $(dir $@)
	${binPath}/spia --dir test/pathways/ --de test/data/DE_Colorectal.tab \
		--nBoots 2000 --array test/data/ALL_Colorectal.tab 2&> $@.tmp
	mv $@.tmp $@

test/out/testACGT.txt: ${binPath}/spia
	mkdir -p $(dir $@)
	${binPath}/spia --pathFiles test/GBM/pathways/ --de test/GBM/sampleDE.txt \
		--nBoots 2000 --array test/GBM/allDE.txt 2&> $@.tmp
	mv $@.tmp $@

test/out/testSPath.txt: ${binPath}/spia
	mkdir -p $(dir $@)
	${binPath}/spia --pathFiles test/alternateFormatPathways/ --de test/data/DE_Colorectal.tab \
		--nBoots 2000 --array test/data/ALL_Colorectal.tab 2&> $@.tmp
	mv $@.tmp $@

test/out/testBeta.txt: ${binPath}/spia
	mkdir -p $(dir $@)
	${binPath}/spia --betaCoFile test/data/beta.txt --dir test/pathways/ \
		--de test/data/DE_Colorectal.tab --nBoots 2000 \
		--array test/data/ALL_Colorectal.tab --verbose 2&> $@.tmp
	mv $@.tmp $@

test/out/testD.txt: ${binPath}/spia
	mkdir -p $(dir $@)
	${binPath}/spia --dir test/pathways/ --de test/data/DE_Colorectal.tab \
		--nBoots 2000 --array test/data/ALL_Colorectal.tab --debug --verbose 2&> $@.tmp
	mv $@.tmp $@

##############################
clean:
	rm -rf ${srcPath}/*.o ${binPath} test/out/

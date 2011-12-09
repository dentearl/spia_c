SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
##############################
# spia makefile
#  dent earl, dearl soe ucsc edu
#  1 may 2009
##############################
binPath = bin
srcPath = src

.PHONY: all clean test

all: ${binPath}/spia

${binPath}/spia: ${srcPath}/spia.c ${srcPath}/spia.h ${srcPath}/fileHandling.o ${srcPath}/hashListFunctions.o ${srcPath}/linearAlg.o ${srcPath}/probFunctions.o
	mkdir -p $(dir $@)
	gcc -g -o $@.tmp $^ -I external/uthash-1.5/src/ -F77 -I77 -llapack -lblas
	mv $@.tmp $@

${srcPath}/%.o: ${srcPath}/%.c
	gcc -c -g -o $@.tmp $< -I external/uthash-1.5/src/
	mv $@.tmp $@

##############################
# system tests
test: test/testBasic.txt test/testTCGA.txt test/testSPath.txt test/testBeta.txt test/testD.txt 

test/testBasic.txt: ${binPath}/spia
	${binPath}/spia --dir test/testPathways/ --de test/testData/DE_Colorectal.tab \
		--nBoots 2000 --array test/testData/ALL_Colorectal.tab 2&> $@.tmp
	mv $@.tmp $@

test/testTCGA.txt: ${binPath}/spia
	${binPath}/spia --pathFiles test/tcgaGBM/pathways/ --de test/tcgaGBM/sample_291.txt \
		--nBoots 2000 --array test/tcgaGBM/tcga_all.txt 2&> $@.tmp
	mv $@.tmp $@

test/testSPath.txt: ${binPath}/spia
	${binPath}/spia --pathFiles test/alternateFormatPathways/ --de test/testData/DE_Colorectal.tab \
		--nBoots 2000 --array test/testData/ALL_Colorectal.tab 2&> $@.tmp
	mv $@.tmp $@

test/testBeta.txt: ${binPath}/spia
	${binPath}/spia --betaCoFile test/testData/beta.txt --dir test/testPathways/ \
		--de test/testData/DE_Colorectal.tab --nBoots 2000 \
		--array test/testData/ALL_Colorectal.tab --verbose 2&> $@.tmp
	mv $@.tmp $@

test/testD.txt: ${binPath}/spia
	${binPath}/spia --dir test/testPathways/ --de test/testData/DE_Colorectal.tab \
		--nBoots 2000 --array test/testData/ALL_Colorectal.tab --debug --verbose 2&> $@.tmp
	mv $@.tmp $@

##############################
clean:
	rm -rf ${srcPath}/*.o ${binPath} test/*.txt

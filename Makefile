########################################
# spia makefile
#  dent earl, dearl soe ucsc edu
#  1 may 2009
########################################
make:
	cat Makefile

spia: spia.c spia.h fileHandling.o hashListFunctions.o linearAlg.o probFunctions.o
	gcc -g spia.c -o spia fileHandling.o hashListFunctions.o linearAlg.o probFunctions.o -I uthash-1.5/src/ -F77 -I77 -llapack -lblas

fileHandling.o: fileHandling.c
	gcc -c -g fileHandling.c -I uthash-1.5/src/

hashListFunctions.o: hashListFunctions.c
	gcc -c -g hashListFunctions.c -I uthash-1.5/src/
linearAlg.o: linearAlg.c
#	gcc -c linearAlg.c -F77 -I77 -I ../uthash-1.5/src/
	gcc -c -g linearAlg.c  -I uthash-1.5/src/
probFunctions.o: probFunctions.c
	gcc -c -g probFunctions.c  -I uthash-1.5/src/


##################################################
# test functions
test:
	./spia --dir testPathways/ --de testData/DE_Colorectal.tab --nBoots 2000 --array testData/ALL_Colorectal.tab
testSPath:
	./spia --pathFiles alternateFormatPathways/ --de testData/DE_Colorectal.tab --nBoots 2000 --array testData/ALL_Colorectal.tab 
testBeta:
	./spia --betaCoFile testData/beta.txt --dir testPathways/ --de testData/DE_Colorectal.tab --nBoots 2000 --array testData/ALL_Colorectal.tab --verbose
testd:
	./spia --dir testPathways/ --de testData/DE_Colorectal.tab --nBoots 2000 --array testData/ALL_Colorectal.tab --debug --verbose
charlie:
	./spia --dir charlieParse/ --de testData/DE_Colorectal.tab --nBoots 2000 --array testData/ALL_Colorectal.tab

##################################################
clean:
	\rm -rf *.o

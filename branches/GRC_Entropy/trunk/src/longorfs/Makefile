CLEANME = *.o *~

all: long-orfs

long-orfs: long-orfs.o delcher.o gene.o
	g++ -o long-orfs long-orfs.o delcher.o gene.o -lm
long-orfs.o: long-orfs.cc
	g++ -c long-orfs.cc
delcher.o: delcher.cc
	g++ -c delcher.cc -I ./
gene.o: gene.cc
	g++ -c gene.cc -I ./
clean:
	rm -rf $(CLEANME)
 

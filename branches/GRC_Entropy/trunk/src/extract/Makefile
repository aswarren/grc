CLEANME = *.o *~

all: extract

extract: extract.o delcher.o
	g++ -o extract extract.o delcher.o -lm
extract.o: extract.cc
	g++ -c extract.cc
delcher.o: delcher.cc
	g++ -c delcher.cc -I ./
clean:
	rm -rf $(CLEANME)
 

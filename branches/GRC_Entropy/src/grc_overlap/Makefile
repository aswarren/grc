all: grc_overlap

grc_overlap: main.o
	g++ -o grc_overlap main.o
main.o: main.cpp
	g++ -c main.cpp -I ./Includes/
clean:
	rm -rf *.o grc_overlap

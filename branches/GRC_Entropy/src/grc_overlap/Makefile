all: grc_overlap

grc_overlap: grc_overlap.o GO.o
	g++ -o grc_overlap grc_overlap.o GO.o
grc_overlap.o: grc_overlap.cpp
	g++ -c -g grc_overlap.cpp -I ./Includes/
GO.o: GO.cpp
	g++ -c -g GO.cpp -I ./Includes/
clean:
	rm -rf *.o grc_overlap

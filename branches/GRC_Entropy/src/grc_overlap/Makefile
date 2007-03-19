all: grc_overlap

grc_overlap: main.o
	g++ -o grc_overlap main.o
main.o: main.cpp
	g++ -c -g main.cpp -I ./Includes/
GO.o: GO.cpp
	g++ -c -g GO.cpp -I ./Includes/
clean:
	rm -rf *.o grc_overlap

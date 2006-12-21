all: grc_compare

grc_compare: main.o GO.o
	g++ -o grc_compare main.o GO.o
main.o: main.cpp
	g++ -c main.cpp -I ./Includes/
GO.o: GO.cpp
	g++ -c GO.cpp -I ./Includes/
clean:
	rm -rf *.o grc_compare

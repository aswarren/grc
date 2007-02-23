all: grc_compare

grc_compare: main.o GO.o
	g++ -o grc_compare main.o GO.o
main.o: main.cpp
	g++ -c -g main.cpp -I ./Includes/
GO.o: GO.cpp
	g++ -c -g GO.cpp -I ./Includes/
clean:
	rm -rf *.o grc_compare

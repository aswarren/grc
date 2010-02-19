SHARE=../shared

all: grc_compare

grc_compare: grc_compare.o GO.o
	g++ -o grc_compare grc_compare.o GO.o
grc_compare.o: grc_compare.cpp
	g++ -c -g grc_compare.cpp -I ./Includes/ -I $(SHARE)/
GO.o: $(SHARE)/GO.cpp
	g++ -c -g $(SHARE)/GO.cpp -I $(SHARE)/
clean:
	rm -rf *.o grc_compare

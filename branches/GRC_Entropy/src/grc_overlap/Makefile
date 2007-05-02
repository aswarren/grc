SHARE=../shared

all: grc_overlap

grc_overlap: grc_overlap.o GO.o CalcPack.o
	g++ -o grc_overlap grc_overlap.o GO.o CalcPack.o
grc_overlap.o: grc_overlap.cpp
	g++ -c -g grc_overlap.cpp -I ./Includes/ -I $(SHARE)/
GO.o: $(SHARE)/GO.cpp
	g++ -c -g $(SHARE)/GO.cpp -I $(SHARE)/
CalcPack.o: $(SHARE)/CalcPack.cpp
	g++ -c -g $(SHARE)/CalcPack.cpp -I $(SHARE)/
clean:
	rm -rf *.o grc_overlap

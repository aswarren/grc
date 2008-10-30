SHARE=../shared

all: grc_annotate

grc_annotate: grc_annotate.o GO.o CalcPack.o
	g++ -o grc_annotate grc_annotate.o GO.o CalcPack.o
grc_annotate.o: grc_annotate.cpp
	g++ -c -g grc_annotate.cpp -I ./Includes/ -I $(SHARE)/
GO.o: $(SHARE)/GO.cpp
	g++ -c -g $(SHARE)/GO.cpp -I $(SHARE)/
CalcPack.o: $(SHARE)/CalcPack.cpp
	g++ -c -g $(SHARE)/CalcPack.cpp -I $(SHARE)/
clean:
	rm -rf *.o grc_annotate

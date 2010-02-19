SHARE=../shared

all: grc_orfs

grc_orfs: grc_orfs.o CalcPack.o
	g++ -o grc_orfs grc_orfs.o CalcPack.o
grc_orfs.o: grc_orfs.cpp
	g++ -c -g grc_orfs.cpp -I $(SHARE)/ -I Includes/
CalcPack.o: $(SHARE)/CalcPack.cpp
	g++ -c -g $(SHARE)/CalcPack.cpp -I $(SHARE)/
clean:
	rm -rf *.o grc_orfs

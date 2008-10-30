SHARE=../shared

all: grc_translate

grc_translate: grc_translate.o
	g++ -o grc_translate grc_translate.o
grc_translate.o: grc_translate.cpp
		g++ -c grc_translate.cpp -I $(SHARE)/
clean:
	rm -rf *.o grc_translate

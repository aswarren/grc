all: grc_translate

grc_translate: main.o
	g++ -o grc_translate main.o
main.o: main.cpp
	g++ -c main.cpp 
clean:
	rm -rf *.o grc_translate

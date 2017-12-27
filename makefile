SI= -I sumgra/
INCLUDES= -I include/
CC=g++ -std=c++11
#CFLAGS=-c -O3 -fopenmp
CFLAGS=-c -O3

all: Mugram.o Extension.o File.o Frequency.o GraphIso.o SubgraphMiner.o SubgraphIso.o Index.o Match.o Trie.o
		$(CC) -O3 -o mugram Mugram.o Extension.o File.o Frequency.o GraphIso.o SubgraphMiner.o SubgraphIso.o Index.o Match.o Trie.o

Mugram.o: Mugram.cpp
	$(CC) $(CFLAGS) $(INCLUDES) Mugram.cpp

Extension.o: src/Extension.cpp include/Extension.h
	$(CC) $(CFLAGS) $(INCLUDES) src/Extension.cpp

File.o: src/File.cpp include/File.h 
	$(CC) $(CFLAGS) $(INCLUDES) src/File.cpp

Frequency.o: src/Frequency.cpp include/Frequency.h
	$(CC) $(CFLAGS) $(INCLUDES) src/Frequency.cpp

GraphIso.o: src/GraphIso.cpp include/GraphIso.h
	$(CC) $(CFLAGS) $(INCLUDES) src/GraphIso.cpp

SubgraphMiner.o: src/SubgraphMiner.cpp include/SubgraphMiner.h
	$(CC) $(CFLAGS) $(INCLUDES) src/SubgraphMiner.cpp

SubgraphIso.o: src/SubgraphIso.cpp include/SubgraphIso.h
	$(CC) $(CFLAGS) $(INCLUDES) src/SubgraphIso.cpp

Index.o: sumgra/Index.cpp sumgra/Index.h
	$(CC) $(CFLAGS) $(SI) sumgra/Index.cpp

Match.o: sumgra/Match.cpp sumgra/Match.h
	$(CC) $(CFLAGS) $(SI) sumgra/Match.cpp

Trie.o: sumgra/Trie.cpp sumgra/Trie.h
	$(CC) $(CFLAGS) $(SI) sumgra/Trie.cpp


clean:
	rm mugram
	rm -f *.o

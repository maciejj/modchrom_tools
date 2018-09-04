#
#g++ -pthread -o mc2 mc2.c pdb.o
g++ -pthread  -O3 -mtune=native -march=native -o fitdna  fitdna_updated.c pdb.o

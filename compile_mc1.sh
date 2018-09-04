#
#g++ -pthread -o mc2 mc2.c pdb.o
g++ -pthread -O3 -mtune=native -march=native -o mc1_branching_noempty_test  mc1.c pdb.o

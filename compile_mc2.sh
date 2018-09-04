#
#g++ -pthread -o mc2 mc2.c pdb.o
g++ -pthread -O3 -mtune=native -march=native -o mc2_branching_noempty_test  mc2.c pdb.o

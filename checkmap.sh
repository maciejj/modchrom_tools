i=1
branchlen=60
domain=60
restraints=500
branchingpoints=350


#with branching points set
# for branchlen in  40 50 60 70 80 90 100; do #7
#  for domain in 20 40 60 80 100 120 140; do #7
#   for restraints in 500 600 700 800; do #4
#    for branchingpoints in 250 300 350 400 450; do #5
#     for i in `seq 1 1 10`; do #10
     if [ -e m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb ] && [ ! -f m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb ]; then
     ./onlyrest_fixed_offset.pl m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb ./${restraints}_harmonic.restlist  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.checkmap.out m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb &     
     fi
 #    done
     wait

#done
#done
#done
#done
 

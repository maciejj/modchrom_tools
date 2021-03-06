#!/bin/bash

show_help() {
cat << EOF
Usage: readInput.sh [OPTIONS]
       -i, --input file   provide previous generated inputfile
EOF
}
input=""
while [[ $# -gt 0 ]]
do
  case $1 in
   -i|--input)
    shift
    input=$1
    ;;
  *)
    echo -e "Error: $0 invalid option '$1'"
    show_help
    exit 1
  esac
shift
done


#tools location
tools_loc="/home/maciej/work/chromosome/webserver/modchrom_tools"


#input="./input/1-60-20-500-250.out"
#----------------get parameters from input file-----------------
copy=`grep "copy" $input | awk '{print substr($0,7,1)}'`
#echo $copy
i="$(echo -e "${copy}" | sed -e 's/^[[:space:]]*//')"
branchingpoints=`grep "branchingpoints" $input | awk '{print substr($0,17,4)}'`
branchingpoints="$(echo -e "${branchingpoints}" | sed -e 's/^[[:space:]]*//')"
restraints=`grep "restraints" $input | awk '{print substr($0,12,4)}'`
restraints="$(echo -e "${restraints}" | sed -e 's/^[[:space:]]*//')"
domain=`grep "domain" $input | awk '{print substr($0,8,3)}'`
domain="$(echo -e "${domain}" | sed -e 's/^[[:space:]]*//')"
branchlen=`grep "branchlength" $input | awk '{print substr($0,14,3)}'`
branchlen="$(echo -e "${branchlen}" | sed -e 's/^[[:space:]]*//')"
chromosomelen=816394 #basepairs
cg15bp_beads=$(($chromosomelen/15))

mc1Steps=20000 #2000000
MDminSteps=100 #800
MD1Steps=2000  #20000
MD2Steps=2000  #20000
MD3Steps=10000 #100000
MD4Steps=25000 #2500000
mc2Steps=100000 #2000000

#-----------------------------monte carlo part 1----------------
if [ ! -f m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.out ]; then
      touch m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.out
      echo m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} >> $(hostname).mc1.startedjobs

#MC first part
      ${tools_loc}/mc1_branching_noempty_test ${mc1Steps} m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} ${restraints}_{harmonic,lower}.restlist $domain $branchlen $branchingpoints $chromosomelen >&  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.out

      echo m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} >> $(hostname).mc1.finishedjobs
fi
echo finished with mc1
#-----------------------molecular dynamics--------------------------------

export CHARMMEXEC="/apps/openmpi-x86_64/bin/mpirun -np 1 /home/asli/c40a2_red/exec/em64t_M/charmm"

 if [  ! -f m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md.pdb ]; then
      touch m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md.pdb
      echo m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} >> $(hostname).md.startedjobs

#MD heat
  blennm=`echo $branchlen | awk '{print $1*0.1}'`
  sed -e "s/BRANCHLEN/$blennm/g" ${tools_loc}/dnacg.TEMPLATE.prm > dnacg.m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.prm

  convpdb.pl -center m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.pdb > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.center.pdb
  convpdb.pl -crdext m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.center.pdb > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.crd
  minCHARMM.pl -custom m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.rest -log m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.min.log -par param=x,xpar=dnacg.m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.prm,xtop=${tools_loc}/dnacg.rtf,nblisttype=bycb,minsteps=0,sdsteps=${MDminSteps},cuton=9,cutoff=10,cutnb=11 -psf m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.psf  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.crd -crdout > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.min.crd

  mdCHARMM.pl -custom m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.rest -log m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md1.log -par param=x,xpar=dnacg.m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.prm,xtop=${tools_loc}/dnacg.rtf,nblisttype=bycb,dynsteps=${MD1Steps},dyntemp=100,cuton=9,cutoff=10,cutnb=11,lang,langfbeta=5,dyntstep=0.001,echeck=1000000  -psf m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.psf  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.min.crd -crdout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md1.crd

  mdCHARMM.pl -custom m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.rest -log m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md2.log -par param=x,xpar=dnacg.m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.prm,xtop=${tools_loc}/dnacg.rtf,nblisttype=bycb,dynsteps=${MD2Steps},dyntemp=200,cuton=9,cutoff=10,cutnb=11,lang,langfbeta=5,dyntstep=0.002,dynoutfrq=1000,echeck=40000 -trajout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md2.dcd -psf m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.psf  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md1.crd -crdout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md2.crd

  mdCHARMM.pl -custom m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.rest -log m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md3.log -par param=x,xpar=dnacg.m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.prm,xtop=${tools_loc}/dnacg.rtf,nblisttype=bycb,dynsteps=${MD3Steps},dyntemp=300,cuton=9,cutoff=10,cutnb=11,lang,langfbeta=5,dyntstep=0.002,dynoutfrq=1000,echeck=40000 -trajout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md3.dcd -psf m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.psf  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md2.crd -crdout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md3.crd
#done

#MD production
  mdCHARMM.pl -custom m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.rest -log m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md4.log -par param=x,xpar=dnacg.m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.prm,xtop=${tools_loc}/dnacg.rtf,nblisttype=bycb,dynsteps=${MD4Steps},dyntemp=300,cuton=9,cutoff=10,cutnb=11,lang,langfbeta=5,dyntstep=0.002,dynoutfrq=5000,echeck=10000 -trajout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md4.dcd -psf m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.psf  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md3.crd -crdout m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md4.crd

 ${tools_loc}/crd2pdb.pl m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md4.crd > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md.pdb

 echo m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} >> $(hostname).md.finishedjobs

fi
 echo finished with md


#----------------------------monte carlo part 2---------------------------
export CHARMMEXEC="/apps/openmpi-x86_64/bin/mpirun -np 1 /home/asli/c40a2_red/exec/em64t_M/charmm"

if [  ! -f m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.pdb ]; then
      touch m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.pdb
      echo m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} >> $(hostname).mc2.startedjobs

#MC2

 ${tools_loc}/crd2pdb.pl m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md4.crd > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.pdb
 grep CONE m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.pdb >> m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.pdb
 egrep '(CDN|PDN)' m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.md4.crd | awk '{printf("%d %18.10f %18.10f %18.10f\n",NR,$5*10.0,$6*10.0,$7*10.0);}' > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.coord

 ${tools_loc}/mc2_branching_noempty_test ${mc2Steps} m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x  ${restraints}_{harmonic,lower}.restlist $domain $branchlen ${chromosomelen} >&  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.x.mcnuc3sx.out

 echo m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints} >> $(hostname).mc2.finishedjobs

fi

#------------------------fitdna----------------------------

#stage1=/home/alexandra/scripts
stage1=./

if [ -e m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.pdb ] && [ ! -f m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb ]; then
	${tools_loc}/fitdna_backup ${stage1}/m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}x.pdb ${branchlen} >  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb
        fi

echo finished fitdna

#--------------------checkmap----------------------------------
if [ -e m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb ] && [ ! -f m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb ]; then
     ${tools_loc}/onlyrest_fixed_offset.pl m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb ./${restraints}_harmonic.restlist  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.checkmap.out m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb
     
#	${tools_loc}/onlyrest.pl m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.pdb ./${restraints}_harmonic.restlist  m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}f.checkmap.out m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb
     fi

#numoflines= `tail -1 m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb | awk '{print substr($0,7,6)}'`
#head $numoflines m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb
echo finished checkmap


#------------------------psf to crd-------------------------------

if [ -f "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb" ] && [ ! -f  "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.crd" ]; then
        genPSF.pl -par param=x,xpar=${tools_loc}/dnacgnew.prm,xtop=${tools_loc}/dnacgnew.rtf,cuton=20,cutoff=25,cutnb=30,patch=CYC:PROA.${cg15bp_beads}:PROA.1,nobuildall m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.psf
        convpdb.pl -scale 10 -crdext m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.pdb > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.crd
        fi
echo finished psf to crd

#------------------------------orient.sh----------------

 if [ -f "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.crd" ] && [ ! -f "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}_oriented.pdb" ]; then

cat > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.inp << EOL
open unit 10 read form name "${tools_loc}/dnacgnew.rtf"
read rtf card unit 10
close unit 10
open unit 10 read form name "${tools_loc}/dnacgnew.prm"
bomlev -1
read para card unit 10
close unit 10
faster on
open unit 10 read form name "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.psf"
read psf card unit 10
open unit 10 read form name "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}fo.crd"
read coor card unit 10
close unit 10
bomlev -2
coor orient sele all end
open unit 10 write card name "m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}_oriented.pdb"
write coor pdb unit 10
*
close unit 10
stop
EOL

/home/asli/charmmc41/exec/em64t_M/2016son/charmm < m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.inp > m${i}_${domain}_${branchlen}_${restraints}_${branchingpoints}.orient.log 2> /dev/null 
fi
echo finished orient

exit 0

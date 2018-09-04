#!/bin/bash
location=$(pwd)

machines=(green1 green3 green4 green5 green6)
for machine in "${machines[@]}"; do
for r in $(seq 1 1 1);do
ssh ${machine} -T -o StrictHostKeyChecking=no << EOF
cd ${location}
./getcpus.pl
#nohup ./run-for-single-all.sh >& ${location}/\$(hostname).\$\$.log & echo \$!
mac=\$(hostname)
#echo "\${mac}    \$!" >> ${location}/$$.machines.pids
exit
EOF
echo "Wait for next run: "
for sec in 3 2 1; do
 echo -ne "${sec}s "
 sleep 1
done
done
echo " "
done


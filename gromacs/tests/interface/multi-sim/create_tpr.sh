
for i in a b c d
do
  cd $i
  gmx_d grompp -f system.mdp -p ../../Common/system.top -c ../../Common/system.gro -o test.tpr
  cd ..
done


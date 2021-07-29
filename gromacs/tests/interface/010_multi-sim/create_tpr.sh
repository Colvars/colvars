gmx_d=gmx_d

cd a
$gmx_d grompp -f system.mdp -p ../../Common/system.top -c init.gro -o test.tpr
cd ..

for i in b c d
do
  cd $i
  $gmx_d grompp -f system.mdp -p ../../Common/system.top -c ../../Common/system.gro -o test.tpr
  cd ..
done


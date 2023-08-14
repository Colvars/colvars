#!/bin/csh -e

echo "Select task: \n"
echo "1: build system, create box, add water,  add ion\n"
echo "2: energy minimize, NVT, and NPT run"
echo "3: prepare MD run"
echo "4: post-processing: extract and center protein"
echo "5: clean the directory"

set task = $< 

setenv pdbfile ../../pdb/ad.pdb

setenv flag '-nt 4 -pin on -v'

setenv gmx gmx

if ($task == 1) then

  echo "Build system:"

  echo "6" | $gmx pdb2gmx -f $pdbfile -water tip3p -o top.gro 

  echo "create box:"
# Create box
  ${gmx} editconf -f top.gro -o box_gro -c -d 1.2 -bt cubic 

# Add Water
  echo "Add water:" 
  ${gmx} solvate -cp box_gro -cs spc216.gro -o solv.gro -p topol.top 

# ION
  echo "Add ion:" 
  ${gmx} grompp -f minim.mdp -maxwarn 4 -c solv.gro -p topol.top -o ions.tpr 
  echo "13" | ${gmx} genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral 

endif

if ($task == 2) then

# Energy minimization
  ${gmx} grompp -f minim.mdp -c ions.gro -p topol.top -o em.tpr
  ${gmx} mdrun $flag -deffnm em

# NVT 
  ${gmx} grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
  ${gmx} mdrun $flag -deffnm nvt

# NPT 
  ${gmx} grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
  ${gmx} mdrun $flag -deffnm npt

endif

if ($task == 3) then 
  ${gmx} grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
endif

if ($task == 4) then 
  ${gmx} trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol -ur compact
endif

if ($task == 5) then
  echo "clean the directory"
  set files=`find . -type f ! -name '*.mdp' ! -name '*.sh' ! -name '*.dat'` 
  set mdoutfile=`find . -type f -name 'mdout.mdp'` 
  echo "files to delete:"
  foreach file ($files $mdoutfile)
    echo $file
  end
  echo "[yes/no]?"
  set tmp = $<
  if ($tmp == 'yes') then
	rm -f $files $mdoutfile
	rm -r figs
  endif
endif


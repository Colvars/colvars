
# MD SETUP
timestep          0.5

colvars		  on

# INPUT
structure               ../Common/da.psf 
parameters              ../Common/par_all22_prot.inp
paraTypeCharmm          on
coordinates             ../Common/da.min.pdb

if { ${outName} == "test.${replicaID}" } {
    temperature             250
} else {
    set inName "test.${replicaID}"
    bincoordinates       $inName.coor
    binvelocities        $inName.vel
}

# OUTPUT

outputname              ${outName}
outputenergies          2000
outputtiming            20000
binaryoutput            yes
restartname             ${outName}r
restartfreq             40000
binaryrestart           yes


seed		87654321

# CONSTANT-T
langevin                on
langevinTemp		    300.0
langevinDamping         10.0
langevinHydrogen	    off

# SPACE PARTITIONING
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  1.0

# CUT-OFFS
switching               on
switchdist              9.0
cutoff                  11.0
pairlistdist            12.0

# RESPA 
fullElectFrequency      1
nonbondedFreq           1

# 1-4 NON-BONDED
exclude                 scaled1-4
1-4scaling              1.0

# COM
commotion               no

# SHAKE
rigidBonds              none

cv config "
colvarsTrajFrequency    1

indexFile ../Common/da.ndx


colvar {

    name d

    outputAppliedForce on
    width 0.5

    lowerBoundary 0.0
    upperBoundary 25.0

    upperWall 25.0
    upperWallConstant 100.0

    distance {
        group1 {
            indexGroup Protein_C-alpha_1
        }
        group2 {
            indexGroup Protein_C-alpha_10
        }
    }
} 


metadynamics {
    colvars d

    hillWeight 0.01
    hillWidth      1.2533141373155001  # Old default
    newHillFrequency 10

    multipleReplicas  yes
    replicaID         ${replicaID}
    replicasRegistry  replicas.registry.txt
    replicaUpdateFrequency 10
}
"

if { [info exists inName] > 0 } {
    cv load $inName.colvars.state
}

run 100

colvar {

    name one

    outputAppliedForce on

    width 0.5

    distanceZ {ifdef(`axis',`
        `axis' (0.3, -0.4, 0.5)',`')
        main {
            indexGroup group5`'ifdef(`fitgroup',`
            centerToReference yes
            rotateToReference yes
            fittingGroup {
                indexGroup heavy_atoms
            }
            refPositionsFile heavy_atoms_refpos.xyz')
        }
        ref {
            indexGroup group1`'ifdef(`fitgroup',`
            centerToReference yes
            rotateToReference yes
            fittingGroup {
                indexGroup heavy_atoms
            }
            refPositionsFile heavy_atoms_refpos.xyz')
        }`'ifdef(`axis',`',`
        ref2 {
            indexGroup group10`'ifdef(`fitgroup',`
            centerToReference yes
            rotateToReference yes
            fittingGroup {
                indexGroup heavy_atoms
            }
            refPositionsFile heavy_atoms_refpos.xyz')
        }')
    }
} 

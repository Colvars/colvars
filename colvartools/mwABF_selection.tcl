# Scripts for multiple-walker ABF in NAMD with selection rules

# NOTE: this script does not perform shared ABF - sharing must be enabled
# separately in the ABF configuration (shared on, sharedFreq xxx)

# Source in NAMD config, and call:
# run_mwABF $cycleNum $selectFreq ($sampleRad $percentStop $biasName)

# cycleNum      Number of cycles of sharing ABF data to be performed
# selectFreq    Number of timesteps between  between rounds of replica selection
#               NOTE: Total simulation time is cycleNum * selectFreq
# sampleRad     Half-length (in bins) of the interval around the current position of a replica
#               in which to count samples for selection purposes (default: 0)
# percentStop   Relative difference (in percent) between lowest and highest sample counts
#               below which no selection is performed (recommended range: 20 to 90, default: 40)
# biasName      Name of the ABF bias to perform selection on (default: "abf1", Colvars default)

# Author: Jeff Comer <jeffcomer at gmail>
# Comer et al. JCTC 2014 https://doi.org/10.1021/ct500874p
# With improvements by Jérôme Hénin <jerome.henin at cnrs dot fr>

# The following is needed for exchanging coordinates between NAMD replicas
# but cannot be "set from script", ie after a scripting command has been run
# So this file should be sourced before any script command (cv, run etc.)
replicaUniformPatchGrids        on


proc run_mwABF_selection { cycleNum selectFreq {sampleRad 0} {percentStop 40} {biasName "abf1"}} {

    set repNum [numReplicas]
    if { $repNum < 2 } {
        print "Error: cannot perform shared ABF with fewer than 2 replicas.\nUse the +replicas command-line flag of NAMD (MPI build).\n\n"
        exit 1
    }
    set rep [myReplica]
    set ncolvars 1
    # Parse ABF config to find out number of colvars
    set config [cv bias $biasName getconfig]
    if { [regexp -nocase {^\s*colvars\s+([^\s{}#]+)} $config match colvars] } {
        set ncolvars [llength $colvars]
        print "Shared ABF $biasName running on $ncolvars colvars: $colvars"
    }

    set s "REPLICA_SETUP"
    foreach var {cycleNum selectFreq repNum rep sampleRad percentStop} {
        set s "$s $var [set $var]"
    }
    print "$s"

    for {set i 0} {$i < $cycleNum} {incr i} {
        print "CYCLE $i"
        # Run the steps until next selection
        # Bias sharing is controlled by the ABF bias at the C++ level
        run $selectFreq

        set count [cv bias $biasName local_sample_count $sampleRad]

        replicaBarrier
        if {$rep > 0} {
            ## Send the count to Replica 0.
            replicaSend $count 0
            ## Receive the source and destination replicas.
            set srcDestList [replicaRecv 0]
        } else {
            ## Build the weight list.
            # The weight is the inverse count + 1 to avoid divide by 0
            set weightList [list [expr {1.0/($count + 1)}]]
            print "REPLICA_COUNT $rep $count"
            set countMin $count
            set countMax $count
            for {set r 1} {$r < $repNum} {incr r} {
                set repCount [replicaRecv $r]
                print "REPLICA_COUNT $r $repCount"
                if {$repCount < $countMin} {
                    set countMin $repCount
                }
                if {$repCount > $countMax} {
                    set countMax $repCount
                }
                # The weight is the inverse count + 1 to avoid divide by 0
                lappend weightList [expr {1.0/($repCount + 1)}]
            }

            ## Normalize the weight list.
            set weightList [normalizeWeights $weightList]
            set s "REPLICA_WEIGHT_LIST"
            foreach w $weightList {
                set s "$s [format " %.3g" $w]"
            }
            print $s

            print "REPLICA_MINMAX $countMin $countMax"
            if {$countMin < 1} {set countMin 1}
            set percentDif [expr {(100*($countMax - $countMin))/$countMax}]

            ## Generate the list of exchanges "srcDestList"
            if {$percentDif < $percentStop} {
                print "REPLICA_SELECTION_DISABLED $percentDif"
                ## If the relative difference between the min and max counts
                ##is less than the threshold, we don't do exchanges.
                set srcDestList {}
            } else {
                print "REPLICA_SELECTION_ENABLED $percentDif"
                set cloneList [resampleWalkers $weightList]
                print "REPLICA_CLONE_LIST $cloneList"

                set srcDestList [resampleExchanges $cloneList]
            }
            ## Replica 0 sends the srcDestList to all other replicas,
            ## so they know who to receive from and who to send to.
            for {set r 1} {$r < $repNum} {incr r} {
                replicaSend $srcDestList $r
            }
        } ;# End Replica 0 work

        ## Everyone should have an identical copy of the srcDestList.
        if {[llength $srcDestList] > 0} {
            print "REPLICA_SRC_DEST_LIST $srcDestList"
            print "REPLICA_BARRIER"
            replicaBarrier

            # Do the coordinate exchanges if this replica is the source or destination.
            foreach srcDest $srcDestList {
                set src [lindex $srcDest 0]
                set dest [lindex $srcDest 1]

                if {$src == $rep} {
                    print "REPLICA_ATOM_SEND $dest"
                    replicaAtomSend $dest
                } elseif {$dest == $rep} {
                    print "REPLICA_ATOM_RECV $src"
                    replicaAtomRecv $src
                }
            }
        } else {
            print "REPLICA_SRC_DEST_EMPTY"
        }
    }
    return $i ;# Number of cycles completed
}

proc normalizeWeights {weightList} {
   ## Normalize the weights.
    set weightSum 0.0
    foreach weight $weightList {
        set weightSum [expr {$weightSum + $weight}]
    }

    set r 0
    set wnList {}
    foreach weight $weightList {
        lappend wnList [expr {double($weight)/$weightSum}]
        incr r
    }
    return $wnList
}


# Determine the number of clones for each walker.
# Weights must be normalized.
proc resampleWalkers {weightList} {
    set num [llength $weightList]

    ## Get the number of clones for each walker.
    set wbar(0) [lindex $weightList 0]
    set u [expr {rand()}]
    set cloneList [list [expr { int($num*$wbar(0)+$u) }] ]
    for {set r 1} {$r < $num} {incr r} {
        set r0 [expr {$r-1}]
        set wbar($r) [expr {$wbar($r0) + [lindex $weightList $r]}]
        lappend cloneList [expr {int($num*$wbar($r)+$u) - int($num*$wbar($r0)+$u) }]
    }
    return $cloneList
}

# Determine the minimal exchanges that must be made to resample.
proc resampleExchanges {cloneList} {
    ## Make a list of exchanges.
    set cloneZeroList {}
    set cloneMultList {}
    set r 0
    foreach cloneNum $cloneList {
        if {$cloneNum == 0} { lappend cloneZeroList $r }
        if {$cloneNum > 1} { lappend cloneMultList $r }
        incr r
    }

    # Is nothing cloned?
    if {[llength $cloneZeroList] == 0} {
        return {}
    }

    ## Walkers cloned multiple times are copied to walkers cloned zero times.
    ## Make the list of exchanges srcDestList
    set srcDestList {}
    set zeroInd 0
    foreach mult $cloneMultList {
        # We get one clone just by leaving the walker where it is.
        set extraNum [expr {[lindex $cloneList $mult]-1}]
        for {set j 0} {$j < $extraNum} {incr j} {
            set dest [lindex $cloneZeroList $zeroInd]
            lappend srcDestList [list $mult $dest]
            incr zeroInd
        }
    }
    return $srcDestList
}

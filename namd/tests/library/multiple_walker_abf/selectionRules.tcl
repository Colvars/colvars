# Procs for implementing selection rules.
# Author: Jeff Comer <jeffcomer at gmail>

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

proc localCount {bin binNum sampleRad biasName} {
    set b0 [expr {$bin-$sampleRad}]
    set b1 [expr {$bin+$sampleRad}]

    set count 0
    for {set b $b0} {$b <= $b1} {incr b} {
	if {$b < 0} {
	    set i 0
	} elseif {$b >= $binNum} {
	    set i [expr {$binNum-1}]
	} else {
	    set i $b
	}
	incr count [cv bias $biasName count $i]
    }
    return $count
}

# Weights are by convention normalized.
proc weightEntropy {weightList} {
    set entropySum 0.0
    set num [llength $weightList]
    foreach w $weightList {
	set entropySum [expr {$entropySum + $w($r)*log($w($r)*$num)}]
    }
    return $entropySum
}


# Initialize global variables

set boltzmann_kB 0.001987191
set temperature [langevinTemp]

set n_replicas [numReplicas]

set be_simulation_step [cv getstepabsolute]


# Create data structures (if not done already)
if { [info exists bias_names] == 0 } {
    for { set i_replica 0 } { ${i_replica} < ${n_replicas} } { incr i_replica } {
        # A group of biases that act on one replica at a time
        set bias_names(${i_replica}) [list]
        # A group of biases, included in the above, that need state exchange
        set stateful_bias_names(${i_replica}) [list]
    }
    unset i_replica
}


# Use current values of the active flag to detect which set of biases is active
proc set_bias_index {} {

    global be_simulation_step
    global bias_index

    if { [info exists be_simulation_step] == 0 } {
        abort "Variable be_simulation_step undefined"
    }

    # Find out which set of biases is active on this replica
    if { ${be_simulation_step} == 0 } {

        # Initialize the index and the active flags for the first time
        set bias_index [myReplica]
        select_bias ${bias_index}

    } else {

        global bias_names
        set n_replicas [numReplicas]
        set n_active_indices 0
        for { set i 0 } { ${i} < ${n_replicas} } { incr i } {
            foreach bias_name $bias_names(${i}) {
                print "BE: REPLICA [myReplica] STEP ${be_simulation_step} BIAS INDEX ${bias_name} HAS STATE [cv bias ${bias_name} get active]"
                if { [cv bias ${bias_name} get active] == 1 } {
                    # Test the first bias in a group only
                    set bias_index ${i}
                    incr n_active_indices
                    break
                }
            }
        }

        if { ${n_active_indices} == 0 } {
            # This is the unbiased replica
            set bias_index [expr ${n_replicas} - 1]
        }

        if { ${n_active_indices} > 1 } {
            abort "Number of active biases is ${n_active_indices}, should be 1 or 0"
        }
    }

    print "BE: REPLICA [myReplica] STEP ${be_simulation_step} BIAS INDEX ${bias_index}"
}


# Set up data structures for metadynamics biases that act on one variable
proc init_be_metadynamics_1d {} {
    global bias_names
    global stateful_bias_names

    set colvar_names [cv list colvars]
    set all_bias_names [cv list biases]

    set num_colvars [llength ${colvar_names}]

    set mtd_biases [list]
    foreach colvar_name ${colvar_names} {
        if { [lsearch ${all_bias_names} "mtd_${colvar_name}"] >= 0 } {
            lappend mtd_biases "mtd_${colvar_name}"
        }
    }

    set num_biases [llength ${mtd_biases}]

    if { ${num_biases} > [numReplicas] } {
        abort "Number of biases is ${num_biases}, but there are only [numReplicas] replicas."
    }

    if { ${num_biases} < [expr [numReplicas] - 1] } {
        abort "Number of biases is ${num_biases}, but there are [numReplicas] replicas and only one can be without bias."
    }

    # Fill each array with the names of the metadynamics bias
    for { set i 0 } { ${i} < ${num_biases} } { incr i } {
        set bias_names(${i}) [list [lindex ${mtd_biases} ${i}]]
        set stateful_bias_names(${i}) [list [lindex ${mtd_biases} ${i}]]
    }

    # Pad the array with a neutral replica
    for { set i ${num_biases} } { ${i} < [numReplicas] } { incr i } {
        set bias_names(${i}) [list]
        set stateful_bias_names(${i}) [list]
    }
}


# Sum up all the energies of the given biases
proc get_total_bias_energy { bias_names } {
    set result 0.0
    foreach bias_name ${bias_names} {
        set result [expr ${result} + [cv bias ${bias_name} energy]]
    }
    return ${result}
}


proc select_bias { i_bias } {
    global bias_names
    # First disable all exchangeable biases
    for { set i 0 } { ${i} < [numReplicas] } { incr i } {
        foreach bias_name $bias_names(${i}) {
            cv bias ${bias_name} set active 0
        }
    }

    # Now turn on the selected one
    foreach bias_name $bias_names(${i_bias}) {
        cv bias ${bias_name} set active 1
        cv bias ${bias_name} update
    }
}


proc print_biases_active_state {} {
    global be_simulation_step
    global bias_names
    # First disable all exchangeable biases
    for { set i 0 } { ${i} < [numReplicas] } { incr i } {
        foreach bias_name $bias_names(${i}) {
            print "BE: REPLICA [myReplica] STEP ${be_simulation_step}: BIAS ${bias_name} ACTIVE FLAG = [cv bias ${bias_name} get active]"
        }
    }
}


proc get_boltzmann_weight { energy_diff } {
    global temperature boltzmann_kB
    return [expr exp(-1.0 * ${energy_diff} / \
                         (${boltzmann_kB} * ${temperature}))]
}


# Send the state of given bias from replica i to replica j
proc share_bias_state { bias_name this_replica i j } {
    if { ${this_replica} == ${i} } {
        replicaSend [cv bias ${bias_name} savetostring] ${j}
    }
    if { ${this_replica} == ${j} } {
        cv bias ${bias_name} loadfromstring [replicaRecv ${i}]
        cv bias ${bias_name} update
    }
}


proc attempt_bias_exchange { i_swap_left i_swap_right } {

    global be_simulation_step
    global bias_index
    global bias_names
    global stateful_bias_names

    set this_replica [myReplica]

    if { ${bias_index} < 0 } {
        abort "BE: REPLICA ${this_replica} UNINITIALIZED OR INVALID BIAS INDEX"
    }

    set line_prefix "BE: REPLICA ${this_replica} STEP ${be_simulation_step}: "


    set this_energy [get_total_bias_energy $bias_names(${bias_index})]

    # Exchange energies of currently active biases

    if { ${this_replica} == ${i_swap_left} } {
        set other_replica ${i_swap_right}
        set other_energy [replicaSendrecv "${bias_index} ${this_energy}" ${i_swap_right}]
        set other_bias_index [lindex ${other_energy} 0]
        set bias_index_left ${bias_index}
        set bias_index_right ${other_bias_index}
    }

    if { ${this_replica} == ${i_swap_right} } {
        set other_replica ${i_swap_left}
        set other_energy [replicaRecv ${i_swap_left}]
        replicaSend "${bias_index} ${this_energy}" ${i_swap_left}
        set other_bias_index [lindex ${other_energy} 0]
        set bias_index_left ${other_bias_index}
        set bias_index_right ${bias_index}
    }

    set other_energy [lindex ${other_energy} 1]

    # Report current values of the bias energy
    print [format "${line_prefix}replica ${this_replica} with bias ${bias_index}, energy: %.6f" ${this_energy}]
    print [format "${line_prefix}replica ${other_replica} with bias ${other_bias_index}, energy: %.6f" ${other_energy}]

    # Synchronize the state of biases that have one
    foreach bias_name $stateful_bias_names(${bias_index_left}) {
        share_bias_state ${bias_name} ${this_replica} \
            ${i_swap_left} ${i_swap_right}
    }
    foreach bias_name $stateful_bias_names(${bias_index_right}) {
        share_bias_state ${bias_name} ${this_replica} \
            ${i_swap_right} ${i_swap_left}
    }

    # Now activate the other replica's bias and recompute the energy

    select_bias ${other_bias_index}
    set this_swapped_energy [get_total_bias_energy \
                                 $bias_names(${other_bias_index})]

    # Share the swapped energy values between the two replicas

    if { ${this_replica} == ${i_swap_left} } {
        set random_number [expr rand()]
        set other_swapped_energy \
            [replicaSendrecv "${this_swapped_energy} ${random_number}" \
                 ${i_swap_right}]
    }

    if { ${this_replica} == ${i_swap_right} } {
        set other_info [split [replicaRecv ${i_swap_left}]]
        set other_swapped_energy [lindex ${other_info} 0]
        set random_number [lindex ${other_info} 1]
        replicaSend ${this_swapped_energy} ${i_swap_left}
    }

    # Report energy values with swapped biases
    print [format "${line_prefix}replica ${this_replica} with bias ${other_bias_index}, energy: %.6f" ${this_swapped_energy}]
    print [format "${line_prefix}replica ${other_replica} with bias ${bias_index}, energy: %.6f" ${other_swapped_energy}]

    set energy_diff [expr (${this_swapped_energy} + ${other_swapped_energy}) - \
                         (${this_energy} + ${other_energy}) ]
    print [format "${line_prefix}energy difference: %.6f" ${energy_diff}]
    set boltzmann_weight [get_boltzmann_weight ${energy_diff}]
    print [format "${line_prefix}Boltzmann vs. random: %.6f vs. %.6f" \
               ${boltzmann_weight} ${random_number}]

    if { ${boltzmann_weight} > ${random_number} } {
        print "${line_prefix}ACCEPTING EXCHANGE"
        set bias_index ${other_bias_index}
    } else {
        print "${line_prefix}REJECTING EXCHANGE"
    }

    select_bias ${bias_index}

    print "BE: REPLICA [myReplica] STEP ${be_simulation_step} BIAS INDEX ${bias_index}"
}


# https://wiki.tcl-lang.org/page/Shuffle+a+list
proc shuffle10 { list } {
    set len [llength $list]
    set len2 $len
    for {set i 0} {$i < $len-1} {incr i} {
        set n [expr {int($i + $len2 * rand())}]
        incr len2 -1
        # Swap elements at i & n
        set temp [lindex $list $i]
        lset list $i [lindex $list $n]
        lset list $n $temp
    }
    return $list
}


# Generate a permutation and share it between replicas
proc get_random_permutation { all_replicas } {
    if { [myReplica] == 0 } {
        set permutation [shuffle10 ${all_replicas}]
        for { set i_rep 1 } { ${i_rep} < [numReplicas] } { incr i_rep } {
            replicaSend ${permutation} ${i_rep}
        }
    } else {
        set permutation [replicaRecv 0]
    }
    return ${permutation}
}


# Attempt bias exchanges between all replicas
proc attempt_exchanges {} {

    global be_simulation_step

    set this_replica [myReplica]
    set n_replicas [numReplicas]

    set all_replicas [list]
    for { set i 0 } { ${i} < ${n_replicas} } { incr i } {
        lappend all_replicas ${i}
    }
    set permutation [get_random_permutation ${all_replicas}]

    for { set i 0 } { ${i} < ${n_replicas} } { incr i 2 } {

        set i_swap_left [lindex ${permutation} ${i}]
        set i_swap_right [lindex ${permutation} [expr ${i} + 1]]

        if { (${this_replica} == ${i_swap_left}) ||
             (${this_replica} == ${i_swap_right}) } {
            print "BE: REPLICA [myReplica] STEP ${be_simulation_step}: ATTEMPTING EXCHANGE BETWEEN REPLICAS ${i_swap_left} and ${i_swap_right}"
            attempt_bias_exchange ${i_swap_left} ${i_swap_right}
        }
    }

    replicaBarrier
}


# Run a MD segment followed by bias exchange attempts
proc run_be_segment { n_steps_segment } {
    global be_simulation_step
    print_biases_active_state
    run norepeat ${n_steps_segment}
    incr be_simulation_step ${n_steps_segment}
    replicaBarrier
    attempt_exchanges
}


# Initialization

if { [info exists bias_index] == 0 } {
    set bias_index -1
    set be_simulation_step 0
}

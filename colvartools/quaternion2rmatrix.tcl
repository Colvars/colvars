proc quaternion_to_matrix { quaternion_string } {

    # this function returns the rotation matrix R(q) as a function of the
    # quaternion q;
    # q is provided as a string, in the format of the collective
    # variables module: (q0, q1, q2, q3)

    # if set to a non-zero value, the function returns a 4x4 matrix suitable for
    # being used in VMD; otherwise, a regular 3x3 matrix
    set vmd_format 1

    set quat_args [split ${quaternion_string} " ,()"]
    set quat_args_parsed [list]
    foreach arg ${quat_args} {
        if { ${arg} != "" } {
            lappend quat_args_parsed ${arg}
        }
    }

    set q0 [lindex ${quat_args_parsed} 0]
    set q1 [lindex ${quat_args_parsed} 1]
    set q2 [lindex ${quat_args_parsed} 2]
    set q3 [lindex ${quat_args_parsed} 3]

    if { ${vmd_format} } {
        return [list \
                  [list [expr ${q0}*${q0} + ${q1}*${q1} - ${q2}*${q2} - ${q3}*${q3}] [expr 2.0 * (${q1}*${q2} - ${q0}*${q3})] [expr 2.0 * (${q0}*${q2} + ${q1}*${q3})] 0.0] \
                  [list [expr 2.0 * (${q0}*${q3} + ${q1}*${q2})] [expr ${q0}*${q0} - ${q1}*${q1} + ${q2}*${q2} - ${q3}*${q3}] [expr 2.0 * (${q2}*${q3} - ${q0}*${q1})] 0.0] \
                  [list [expr 2.0 * (${q1}*${q3} - ${q0}*${q2})] [expr 2.0 * (${q0}*${q1} + ${q2}*${q3})] [expr ${q0}*${q0} - ${q1}*${q1} - ${q2}*${q2} + ${q3}*${q3}] 0.0] \
                  [list 0.0 0.0 0.0 1.0] \
                 ]
    } else {
        return [list \
                  [list [expr ${q0}*${q0} + ${q1}*${q1} - ${q2}*${q2} - ${q3}*${q3}] [expr 2.0 * (${q1}*${q2} - ${q0}*${q3})] [expr 2.0 * (${q0}*${q2} + ${q1}*${q3})]] \
                  [list [expr 2.0 * (${q0}*${q3} + ${q1}*${q2})] [expr ${q0}*${q0} - ${q1}*${q1} + ${q2}*${q2} - ${q3}*${q3}] [expr 2.0 * (${q2}*${q3} - ${q0}*${q1})]] \
                  [list [expr 2.0 * (${q1}*${q3} - ${q0}*${q2})] [expr 2.0 * (${q0}*${q1} + ${q2}*${q3})] [expr ${q0}*${q0} - ${q1}*${q1} - ${q2}*${q2} + ${q3}*${q3}]] \
                 ]
    }

}

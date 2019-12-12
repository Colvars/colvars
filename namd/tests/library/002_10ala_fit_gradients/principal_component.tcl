# VMD tcl code to compute eigenvector with differenceVector

if { 1 } {

    # vmd -e principal_component.tcl \
    #     -f structure_a.pdb \
    #     -f structure_b.pdb \
    #     -f sim.psf -args *dcd


    set bb_0 [atomselect 0 "Ca"]
    set bb_1 [atomselect 1 "Ca"]
    ${bb_1} move [measure fit ${bb_1} ${bb_0}]

    set mol_0_name [file rootname [molinfo 0 get name]]
    set mol_1_name [file rootname [molinfo 1 get name]]

    set diff_x [vecsub [${bb_1} get "x"] [${bb_0} get "x"]]
    set diff_y [vecsub [${bb_1} get "y"] [${bb_0} get "y"]]
    set diff_z [vecsub [${bb_1} get "z"] [${bb_0} get "z"]]

    puts "diff_x = ${diff_x}"
    puts "diff_y = ${diff_y}"
    puts "diff_z = ${diff_z}"

    set diff_norm2 [expr [veclength2 ${diff_x}] + [veclength2 ${diff_y}] + [veclength2 ${diff_z}]]
    set mode_x [vecscale ${diff_x} [expr 1.0/${diff_norm2}]]
    set mode_y [vecscale ${diff_y} [expr 1.0/${diff_norm2}]]
    set mode_z [vecscale ${diff_z} [expr 1.0/${diff_norm2}]]

    package require pbctools

    set mol_name [molinfo top get name]
    set mol_name_ext [file extension ${mol_name}]
    set mol_name [file rootname ${mol_name}]
    if { ${mol_name_ext} == ".pdb" } {
        animate delete beg 0 end 0 top
        set mol_name_ext [file extension ${mol_name}]
        if { ${mol_name_ext} == ".reduced" } {
            set mol_name [file rootname ${mol_name}]
        }
    }

    set all [atomselect top "all"]
    set bb  [atomselect top "Ca"]

    set num_frames_total 0

    # set output [open "${mol_name}.proj.dat" "w"]

}

# foreach filename ${argv} {

#     animate delete all
#     mol addfile ${filename} waitfor all
#     set filename_base [file rootname ${filename}]
#     set num_frames [molinfo top get numframes]

#     for { set i 0 } { ${i} < ${num_frames} } { incr i } {

#         set frame [expr ${num_frames_total} + ${i} + 1]

#         puts -nonewline stderr "Reading frame ${i}\r"
#         ${all} frame ${i}
#         ${bb}  frame ${i}
#         ${all} move [measure fit ${bb} ${bb_0} weight occupancy]

#         set diff_x [vecsub [${bb} get "x"] [${bb_0} get "x"]]
#         set diff_y [vecsub [${bb} get "y"] [${bb_0} get "y"]]
#         set diff_z [vecsub [${bb} get "z"] [${bb_0} get "z"]]

#         set proj [expr [vecdot ${mode_x} ${diff_x}] + [vecdot ${mode_y} ${diff_y}] + [vecdot ${mode_z} ${diff_z}]]

#         puts -nonewline ${output} [format " %6d  " ${frame} ]
#         puts -nonewline ${output} [format "  %13.9f\n" ${proj}]
#     }

#     set num_frames_total [expr ${num_frames_total} + ${num_frames}]
# }

# close ${output}

# quit

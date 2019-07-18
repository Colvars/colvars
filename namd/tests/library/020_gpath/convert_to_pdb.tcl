mol new trialanine.parm7
set num_images 30
# set num_heavy_atoms 20
for {set i 0} {$i < $num_images} {incr i} {
    set image_filename [format "string-%02d.coor" $i]
    mol addfile $image_filename
}
set num_frames [molinfo top get numframes]
set all_atoms [atomselect top all]
for {set i 0} {$i < $num_frames} {incr i} {
    set pdb_filename [format "string-%02d.pdb" $i]
    $all_atoms frame $i
    $all_atoms update
    $all_atoms writepdb $pdb_filename
}

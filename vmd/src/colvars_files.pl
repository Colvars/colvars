# List of files for the Colvars module

our ( $colvars_defines );

our ( @colvars_cc );
our ( @colvars_cu );
our ( @colvars_ccpp );
our ( @colvars_h );

$colvars_defines = " -DVMDCOLVARS";

@colvars_cc      = ();
@colvars_cu      = ();
@colvars_ccpp    = ('colvaratoms.C',
                    'colvarbias.C',
                    'colvarbias_abf.C',
                    'colvarbias_alb.C',
                    'colvarbias_histogram.C',
                    'colvarbias_meta.C',
                    'colvarbias_restraint.C',
                    'colvar.C',
                    'colvarcomp.C',
                    'colvarcomp_angles.C',
                    'colvarcomp_coordnums.C',
                    'colvarcomp_distances.C',
                    'colvarcomp_protein.C',
                    'colvarcomp_rotations.C',
                    'colvardeps.C',
                    'colvargrid.C',
                    'colvarmodule.C',
                    'colvarparse.C',
                    'colvarproxy.C',
                    'colvarproxy_vmd.C',
                    'colvarscript.C',
                    'colvartypes.C',
                    'colvarvalue.C');
@colvars_h    =    ('colvar_UIestimator.h',
                    'colvaratoms.h',
                    'colvarbias.h',
                    'colvarbias_abf.h',
                    'colvarbias_alb.h',
                    'colvarbias_histogram.h',
                    'colvarbias_meta.h',
                    'colvarbias_restraint.h',
                    'colvarcomp.h',
                    'colvardeps.h',
                    'colvargrid.h',
                    'colvar.h',
                    'colvarmodule.h',
                    'colvarparse.h',
                    'colvarproxy.h',
                    'colvarproxy_vmd.h',
                    'colvarscript.h',
                    'colvars_version.h',
                    'colvartypes.h',
                    'colvarvalue.h');

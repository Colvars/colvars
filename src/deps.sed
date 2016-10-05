
# delete irrelevant lines
/^ *f_.*(f_.*/ !d

# node attributes based on description
s/^ *f_description(f_cvb_\(.*\), "\(.*\)");/  cvb_\1 [label = "\2"] [fontcolor = "red"];/
s/^ *f_description(f_cv_\(.*\), "\(.*\)");/  cv_\1 [label = "\2"] [fontcolor = "black"];/
s/^ *f_description(f_cvc_\(.*\), "\(.*\)");/  cvc_\1 [label = "\2"] [fontcolor = "blue"];/
s/^ *f_description(f_ag_\(.*\), "\(.*\)");/  ag_\1 [label = "\2"] [fontcolor = "darkgreen"];/

# edges for dependencies
s/^ *f_req_self(f_\(.*\), f_\(.*\));.*/  \1 -> \2;/
s/^ *f_req_children(f_\(.*\), f_\(.*\));.*/  \1 -> \2;/
s/^ *f_req_alt2(f_\(.*\), f_\(.*\), f_\(.*\));.*/  \1 -> {\2 \3} [style=dotted];/
s/^ *f_req_alt3(f_\(.*\), f_\(.*\), f_\(.*\), f_\(.*\));.*/  \1 -> {\2 \3 \4} [style=dotted];/

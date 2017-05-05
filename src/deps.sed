
# delete irrelevant lines
/^ *f_.*(f_.*\|^ *init_feature(/ !d

# node attributes based on description
s/^ *init_feature(f_cvb_\(.*\), "\(.*\)".*);/  cvb_\1 [label = "\2"] [fontcolor = "red"];/
s/^ *init_feature(f_cv_\(.*\), "\(.*\)".*);/  cv_\1 [label = "\2"] [fontcolor = "black"];/
s/^ *init_feature(f_cvc_\(.*\), "\(.*\)".*);/  cvc_\1 [label = "\2"] [fontcolor = "blue"];/
s/^ *init_feature(f_ag_\(.*\), "\(.*\)".*);/  ag_\1 [label = "\2"] [fontcolor = "darkgreen"];/

# edges for dependencies
s/^ *f_req_self(f_\(.*\), f_\(.*\));.*/  \1 -> \2;/
s/^ *f_req_exclude(f_\(.*\), f_\(.*\));.*/  \1 -> \2 [dir = "both"] [arrowhead = "tee"] [arrowtail = "tee"] [style = dotted];/
s/^ *f_req_children(f_\(.*\), f_\(.*\));.*/  \1 -> \2;/
s/^ *f_req_alt2(f_\(.*\), f_\(.*\), f_\(.*\));.*/  \1 -> {\2 \3} [style = dashed];/
s/^ *f_req_alt3(f_\(.*\), f_\(.*\), f_\(.*\), f_\(.*\));.*/  \1 -> {\2 \3 \4} [style = dotted];/

set output [open "cvscript-tcl.tex" "w"]

foreach prefix [list "cv" "colvar" "bias"] \
    introstr [list "\\cvsubsec\{Commands to manage the Colvars module\}\{sec:cvscript_tcl_cv\}" \
                  "\\cvsubsec\{Commands to manage individual collective variables\}\{sec:cvscript_tcl_colvar\}" \
                  "\\cvsubsec\{Commands to manage individual biases\}\{sec:cvscript_tcl_bias\}"] \
    {
        puts ${output} ${introstr}
        puts ${output} "\\begin\{itemize\}"

        foreach cmd [cv listcommands] {

            if { (${cmd} == "cv_colvar") || (${cmd} == "cv_bias") } continue

            if { ${prefix} == "cv" } {
                if { [string range ${cmd} 0 2] != "cv_" } { continue }
                set subcmd [string range ${cmd} 3 end]
                set lines [split [cv help ${subcmd}] "\n"]
            }
            if { ${prefix} == "bias" } {
                if { [string range ${cmd} 0 4] != "bias_"} { continue }
                set subcmd [string range ${cmd} 5 end]
                set lines [split [cv bias name help ${subcmd}] "\n"]
            }
            if { ${prefix} == "colvar" } {
                if { [string range ${cmd} 0 6] != "colvar_"} { continue }
                set subcmd [string range ${cmd} 7 end]
                set lines [split [cv colvar name help ${subcmd}] "\n"]
            }

            set line [lindex ${lines} 0]
            # Sanitize for LaTeX
            set line [regsub -all "_" "${line}" "\\_"]
            puts ${output} "\\item \\texttt\{${line}\}"

            foreach line [lrange ${lines} 2 end] {
                set line [regsub -all "_" "${line}" "\\_"]
                if { ${line} == "" || ${line} == "----------" } continue
                set line [regsub -all "Parameters" "${line}" "\\textbf\{Parameters\}"]
                puts ${output} "\\\\"
                puts ${output} "\\texttt\{${line}\}"
            }
        }
        puts ${output} "\\end\{itemize\}"
    }
close ${output}


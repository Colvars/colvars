if { [file exists ../build/libcolvars.so] > 0 } {
    load ../build/libcolvars.so
}

if { [info exists output] == 0 } {
    set output [open "cvscript-tcl.tex" "w"]
}

proc gen_cmdline_latex_reference { output { main_cmd "cv" } } {

    foreach prefix [list "cv" "colvar" "bias"] \
        introstr [list "\\cvsubsec\{Commands to manage the Colvars module\}\{sec:cvscript_cmdline_cv\}" \
                      "\\cvsubsec\{Commands to manage individual collective variables\}\{sec:cvscript_cmdline_colvar\}" \
                      "\\cvsubsec\{Commands to manage individual biases\}\{sec:cvscript_cmdline_bias\}"] {

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
            # Allow overriding the main command (for fix_modify)
            set line [regsub -all "^cv" "${line}" "${main_cmd}"]
            puts ${output} "\\item \\texttt\{${line}\}"

            foreach line [lrange ${lines} 2 end] {
                set line [regsub -all "_" "${line}" "\\_"]
                if { ${line} == "" || ${line} == "-------" } continue
                set line [regsub -all "Returns" "${line}" "\\textbf\{Returns\}"]
                if { ${line} == "" || ${line} == "----------" } continue
                set line [regsub -all "Parameters" "${line}" "\\textbf\{Parameters\}"]
                puts ${output} "\\\\"
                puts ${output} "\\texttt\{${line}\}"
            }
        }
        puts ${output} "\\end\{itemize\}"
    }
}


set output [open "cvscript-tcl.tex" "w"]
gen_cmdline_latex_reference ${output} "cv"
close ${output}

set output [open "cvscript-fix-modify.tex" "w"]
gen_cmdline_latex_reference ${output} "fix\\_modify Colvars"
close ${output}

exit

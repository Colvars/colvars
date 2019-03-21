.SILENT:

VMFILES = cv_dashboard.tcl cv_dashboard_main.tcl cv_dashboard_editor.tcl cv_dashboard_plots.tcl pkgIndex.tcl templates
VMVERSION = 1.0
DIR = $(PLUGINDIR)/noarch/tcl/cv_dashboard$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying cv_dashboard $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp -r $(VMFILES) $(DIR) 


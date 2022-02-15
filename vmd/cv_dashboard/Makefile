.SILENT:

VMFILES := $(wildcard *.tcl) templates
VMVERSION = 1.5
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


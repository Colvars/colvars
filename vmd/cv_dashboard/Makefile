.SILENT:

VMFILES := $(wildcard *.tcl) templates VERSION
cat := $(if $(filter $(OS),Windows_NT),type,cat)
VMVERSION := $(shell $(cat) VERSION)
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


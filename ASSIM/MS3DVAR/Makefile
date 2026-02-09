# MS3DVAR Root Makefile
#
# This Makefile helps build all MS3DVAR variants from the root directory
#
#======================================================================

# All variants
VARIANTS = LR MR FILTER MS

# Executables produced by each variant
EXECUTABLES = $(VARIANTS:%=%/msda_$(subst da_,,$@))

#======================================================================
# Default target
#======================================================================

.PHONY: all help
all: help

#======================================================================
# Help target
#======================================================================

help:
	@echo ""
	@echo "MS3DVAR Multi-Scale Data Assimilation Build System"
	@echo "==================================================="
	@echo ""
	@echo "Available targets:"
	@echo ""
	@echo "  make all_variants    - Build all MS3DVAR variants"
	@echo "  make lr              - Build Low Resolution variant only"
	@echo "  make mr              - Build Medium Resolution variant only"
	@echo "  make filter          - Build Filter variant only"
	@echo "  make ms              - Build Multi-Scale variant only"
	@echo ""
	@echo "  make clean           - Clean all variants"
	@echo "  make clean_lr        - Clean Low Resolution variant"
	@echo "  make clean_mr        - Clean Medium Resolution variant"
	@echo "  make clean_filter    - Clean Filter variant"
	@echo "  make clean_ms        - Clean Multi-Scale variant"
	@echo ""
	@echo "  make clobber         - Remove all executables and build files"
	@echo ""
	@echo "  make status          - Show build status of all variants"
	@echo "  make setup           - Setup Makedefs for all variants"
	@echo ""
	@echo "Structure:"
	@echo "  COMMON/           - Common code (164 files) shared by all variants"
	@echo "  LR/               - Low Resolution variant"
	@echo "  MR/               - Medium Resolution variant"
	@echo "  FILTER/           - Filter variant"
	@echo "  MS/               - Multi-Scale variant"
	@echo ""

#======================================================================
# Build all variants
#======================================================================

.PHONY: all_variants
all_variants: $(VARIANTS)

#======================================================================
# Build individual variants
#======================================================================

.PHONY: lr mr filter ms
lr:
	@echo "Building Low Resolution variant..."
	cd LR && $(MAKE)

mr:
	@echo "Building Medium Resolution variant..."
	cd MR && $(MAKE)

filter:
	@echo "Building Filter variant..."
	cd FILTER && $(MAKE)

ms:
	@echo "Building Multi-Scale variant..."
	cd MS && $(MAKE)

#======================================================================
# Build targets for VARIANTS list
#======================================================================

.PHONY: $(VARIANTS)
$(VARIANTS):
	@echo "Building variant: $@..."
	cd $@ && $(MAKE)

#======================================================================
# Setup targets
#======================================================================

.PHONY: setup setup_lr setup_mr setup_filter setup_ms
setup: setup_lr setup_mr setup_filter setup_ms

setup_lr:
	@if [ ! -f LR/Makedefs ]; then \
		echo "Creating LR/Makedefs from Makedefs.generic"; \
		cp LR/Makedefs.generic LR/Makedefs; \
	else \
		echo "LR/Makedefs already exists"; \
	fi

setup_mr:
	@if [ ! -f MR/Makedefs ]; then \
		echo "Creating MR/Makedefs from Makedefs.generic"; \
		cp MR/Makedefs.generic MR/Makedefs; \
	else \
		echo "MR/Makedefs already exists"; \
	fi

setup_filter:
	@if [ ! -f FILTER/Makedefs ]; then \
		echo "Creating FILTER/Makedefs from Makedefs.generic"; \
		cp FILTER/Makedefs.generic FILTER/Makedefs; \
	else \
		echo "FILTER/Makedefs already exists"; \
	fi

setup_ms:
	@if [ ! -f MS/Makedefs ]; then \
		echo "Creating MS/Makedefs from Makedefs.generic"; \
		cp MS/Makedefs.generic MS/Makedefs; \
	else \
		echo "MS/Makedefs already exists"; \
	fi

#======================================================================
# Status target
#======================================================================

.PHONY: status
status:
	@echo ""
	@echo "MS3DVAR Build Status"
	@echo "===================="
	@echo ""
	@echo "Low Resolution (LR):"
	@if [ -f LR/msLR ]; then \
		echo "  [BUILT] msLR executable exists"; \
		ls -lh LR/msLR | awk '{print "  Size: " $$5 "  Modified: " $$6 " " $$7 " " $$8}'; \
	else \
		echo "  [NOT BUILT] msLR executable not found"; \
	fi
	@echo ""
	@echo "Medium Resolution (MR):"
	@if [ -f MR/msMR ]; then \
		echo "  [BUILT] msMR executable exists"; \
		ls -lh MR/msMR | awk '{print "  Size: " $$5 "  Modified: " $$6 " " $$7 " " $$8}'; \
	else \
		echo "  [NOT BUILT] msMR executable not found"; \
	fi
	@echo ""
	@echo "Filter (FILTER):"
	@if [ -f FILTER/msFILTER ]; then \
		echo "  [BUILT] msFILTER executable exists"; \
		ls -lh FILTER/msFILTER | awk '{print "  Size: " $$5 "  Modified: " $$6 " " $$7 " " $$8}'; \
	else \
		echo "  [NOT BUILT] msFILTER executable not found"; \
	fi
	@echo ""
	@echo "Multi-Scale (MS):"
	@if [ -f MS/msMS ]; then \
		echo "  [BUILT] msMS executable exists"; \
		ls -lh MS/msMS | awk '{print "  Size: " $$5 "  Modified: " $$6 " " $$7 " " $$8}'; \
	else \
		echo "  [NOT BUILT] msMS executable not found"; \
	fi
	@echo ""

#======================================================================
# Clean targets
#======================================================================

.PHONY: clean clean_lr clean_mr clean_filter clean_ms
clean: clean_lr clean_mr clean_filter clean_ms
	@echo "All variants cleaned"

clean_lr:
	@echo "Cleaning Low Resolution variant..."
	cd LR && $(MAKE) clean

clean_mr:
	@echo "Cleaning Medium Resolution variant..."
	cd MR && $(MAKE) clean

clean_filter:
	@echo "Cleaning Filter variant..."
	cd FILTER && $(MAKE) clean

clean_ms:
	@echo "Cleaning Multi-Scale variant..."
	cd MS && $(MAKE) clean

#======================================================================
# Clobber targets (remove executables too)
#======================================================================

.PHONY: clobber clobber_lr clobber_mr clobber_filter clobber_ms
clobber: clobber_lr clobber_mr clobber_filter clobber_ms
	@echo "All variants clobbered (executables removed)"

clobber_lr:
	@echo "Clobbering Low Resolution variant..."
	cd LR && $(MAKE) clobber

clobber_mr:
	@echo "Clobbering Medium Resolution variant..."
	cd MR && $(MAKE) clobber

clobber_filter:
	@echo "Clobbering Filter variant..."
	cd FILTER && $(MAKE) clobber

clobber_ms:
	@echo "Clobbering Multi-Scale variant..."
	cd MS && $(MAKE) clobber

#======================================================================
# Info target
#======================================================================

.PHONY: info
info:
	@echo ""
	@echo "MS3DVAR Code Structure Information"
	@echo "==================================="
	@echo ""
	@echo "Common code:"
	@echo "  Location: COMMON/"
	@echo "  Files: $$(ls COMMON/*.F COMMON/*.h 2>/dev/null | wc -l) source files"
	@echo ""
	@echo "Variant-specific code:"
	@for variant in $(VARIANTS); do \
		count=$$(ls $$variant/*.F $$variant/*.h 2>/dev/null | wc -l); \
		echo "  $$variant: $$count files"; \
	done
	@echo ""
	@echo "CROCO integration:"
	@echo "  CROCO sources: ../../../OCEAN/"
	@echo "  Using VPATH in Makefiles to locate sources"
	@echo ""

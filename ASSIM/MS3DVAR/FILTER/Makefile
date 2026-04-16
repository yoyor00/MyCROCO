# Makefile for MS3DVAR Filter variant
#
# This is a simple variant-specific Makefile that includes
# the common build rules from ../Makefile.inc

# Variant-specific settings
VARIANT_NAME = FILTER
SBIN = ms3dvar-filter

# FILTER-specific sources (not in common Makefile.inc)
VARIANT_SRCS = netcdf_filtered.F das_filter.F

# Include common build rules
include Makefile.inc

# Add variant-specific objects
VARIANT_OBJS = $(VARIANT_SRCS:.F=.o)

# Override the link rule to include variant objects
$(SBIN): cppdefs.h param.h scalars.h ncscrum.h \
         $(OBJS) $(OBJS90) $(VARIANT_OBJS)
	$(LDR) $(FFLAGS) -o $(SBIN) $(OBJS) $(OBJS90) $(VARIANT_OBJS) $(LCDF) $(LDFLAGS)

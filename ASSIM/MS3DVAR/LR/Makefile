# Makefile for MS3DVAR Low Resolution (LR) variant
#
# This is a simple variant-specific Makefile that includes
# the common build rules from ../Makefile.inc

# Variant-specific settings
VARIANT_NAME = LR
SBIN = ms3dvar-lr

# LR-specific sources (not in common Makefile.inc)
VARIANT_SRCS = das_get_initial.F

# Include common build rules
include Makefile.inc

# Add variant-specific objects
VARIANT_RCS = $(VARIANT_SRCS:.F=.f)
VARIANT_OBJS = $(VARIANT_RCS:.f=.o)

# Override the link rule to include variant objects
$(SBIN): $(OBJS) $(OBJS90) $(VARIANT_OBJS)
	$(LDR) $(FFLAGS) -o $(SBIN) $(OBJS) $(OBJS90) $(VARIANT_OBJS) $(LCDF) $(LDFLAGS)

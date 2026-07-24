# Makefile for MS3DVAR Low Resolution (LR) scale
#
# This is a simple scale-specific Makefile that includes
# the common build rules from ../Makefile.inc

# Scale-specific settings
SCALE_NAME = LR
SBIN = ms3dvar-lr

# LR-specific sources (not in common Makefile.inc)
SCALE_SRCS = ms3dvar_get_increment.F

# Include common build rules
include Makefile.inc

# Add scale-specific objects
SCALE_RCS = $(SCALE_SRCS:.F=.f)
SCALE_OBJS = $(SCALE_RCS:.f=.o)

# Override the link rule to include scale objects
$(SBIN): $(OBJS) $(OBJS90) $(SCALE_OBJS)
	$(LDR) $(FFLAGS) -o $(SBIN) $(OBJS) $(OBJS90) $(SCALE_OBJS) $(LCDF) $(LDFLAGS)

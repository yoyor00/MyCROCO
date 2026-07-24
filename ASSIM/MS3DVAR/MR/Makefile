# Makefile for MS3DVAR Medium Resolution (MR) scale
#
# This is a simple scale-specific Makefile that includes
# the common build rules from ../Makefile.inc

# Scale-specific settings
SCALE_NAME = MR
SBIN = ms3dvar-mr

# MR-specific sources (not in common Makefile.inc)
SCALE_SRCS = ms3dvar_add_lr.F ms3dvar_get_increment_lr.F

# Include common build rules
include Makefile.inc

# Add scale-specific objects
SCALE_RCS = $(SCALE_SRCS:.F=.f)
SCALE_OBJS = $(SCALE_RCS:.f=.o)

# Override the link rule to include scale objects
$(SBIN): $(OBJS) $(OBJS90) $(SCALE_OBJS)
	$(LDR) $(FFLAGS) -o $(SBIN) $(OBJS) $(OBJS90) $(SCALE_OBJS) $(LCDF) $(LDFLAGS)

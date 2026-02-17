# $Id$
#
#======================================================================
# CROCO is a branch of ROMS developped at IRD and INRIA, in France
# The two other branches from UCLA (Shchepetkin et al)
# and Rutgers University (Arango et al) are under MIT/X style license.
# CROCO specific routines (nesting) are under CeCILL-C license.
#
# CROCO website : http://www.croco-ocean.org
#======================================================================
#
# Common Makefile for MS3DVAR variants
# This file is included by each variant's Makefile
#
# Required variables to be set by variant Makefile:
#   VARIANT_NAME  - Name of variant (e.g., "lr", "mr", "filter", "ms")
#   SBIN          - Executable name (e.g., "msda_lr")
#
#======================================================================

#======================================================================
# Set machine dependent definitions and rules.
#======================================================================

include Makedefs

#======================================================================
# Set paths to common code and CROCO sources
#======================================================================

# Path to common MS3DVAR code
COMMON_DIR = ../COMMON

# Path to CROCO sources
CROCO_SRC = ../../../OCEAN

# Add paths to VPATH for source file search
VPATH = $(COMMON_DIR):$(CROCO_SRC)

#======================================================================
# Model Configuration:
#======================================================================
# SRCS: source codes files
# Variant-specific files are in current directory
# Common MS3DVAR files are in $(COMMON_DIR)
# CROCO files are in $(CROCO_SRC)
#======================================================================

 SRCS = main.F		read_inp.F\
	timers_roms.F	init_scalars.F	init_arrays.F\
	set_scoord.F	setup_grid1.F	setup_grid2.F\
	ana_initial.F	analytical.F\
	set_depth.F\
	diag.F		checkdims.F	grid_stiffness.F\
	debug.F\
        param.F         ncscrum.F   scalars.F\
	put_global_atts.F\
	nf_fread.F nf_fread_x.F	nf_fread_y.F\
	get_date.F	lenstr.F	closecdf.F	insert_node.F\
	fillvalue.F	nf_add_attribute.F		set_cycle.F\
	def_grid_2d.F   def_grid_3d.F   def_his.F       def_rst.F\
	wrt_grid.F      wrt_rst.F\
	get_grid.F	get_initial.F\
	get_ssh.F\
	das_init_arrays.F  das_init_scalars.F\
	das_setup1.F  das_setup2.F das_spline1d.F das_exchange.F\
	das_MessPass2D.F das_s2z.F das_z2s.F das_analysis.F \
	das_ana_bvar.F das_costb.F das_vec.F das_adj_costb.F\
    das_get_d2c.F\
	das_get_tmi.F  das_get_mc.F das_get_goes.F das_lbfgs.F das_innov_sst.F\
	das_cost_sst.F das_adj_csst.F das_ca2zero.F das_adj_zero.F\
	das_adj_zero_s.F das_read_inp.F das_matvec.F das_adj_matvec.F\
	das_get_corr.F das_varvec.F das_adj_varvec.F\
	das_rho_eos9.F das_fill9.F das_extrapolate9.F das_geos_uv_p.F das_get_bvar.F\
    das_diag_steric_sshb.F\
	das_get_ssh_ref.F \
	das_get_zeta_ubal_from_tot.F das_adj_get_zeta_ubal_from_tot.F\
\
	das_zeta_p_smooth.F das_zeta_p_smooth2.F das_diag_pressure.F\
	das_adj_diag_p.F das_get_hfradar.F das_geo_ratio.F das_adj_geo_ratio.F\
	das_cost_reg_uv.F das_adj_cost_reg_uv.F\
    das_d2c_uv.F\
\
	das_get_js1.F das_get_swot.F das_innov_js1.F das_innov_swot.F das_rho_eos_lin.F\
	das_get_tp.F das_innov_tp.F\
	das_diag_ssh.F das_cost_js1.F das_cost_swot.F das_wrt.F das_adj_wrt.F \
    das_theadj_wrt.F\
	das_adj_rho_eos.F das_adj_diag_ssh.F das_adj_cjs1.F das_adj_cswot.F\
	das_adj_diag_ssh_steric.F\
	das_adj_geos_uv_p.F das_psichi_uv.F das_adj_psichi_uv.F\
	das_cross_ts.F das_adj_cross_ts.F das_set_geo_ratio.F\
	das_zetah_smooth.F das_adj_zetah_smooth.F das_copy_zetah_smooth.F\
    das_zeta_h_smooth.F das_adj_zeta_h_smooth.F\
	das_adj_copy_zeteh_smooth.F das_adj_glider_sio.F das_cost_glider_sio.F\
	das_adj_glider_whoi.F das_cost_glider_whoi.F das_adj_map.F das_cost_map.F\
	das_adj_glider_ptsur.F das_cost_glider_ptsur.F\
	das_pressue_smooth.F das_adj_pressure_smooth.F\
	das_copy_pressure_smooth.F das_adj_copy_pressure_smooth.F\
	das_get_sats.F das_innov_sats.F das_cost_sats.F das_adj_csats.F\
\
	das_get_insitu.F das_innov_sio.F das_cost_glider.F das_adj_cglider.F\
	das_innov_whoi.F das_innov_flight.F das_cost_flight.F das_adj_cflight.F\
	das_innov_ptsur.F das_innov_martn.F das_innov_map.F\
	das_get_auv.F das_innov_calauv.F das_cost_calauv.F das_adj_cauv.F\
	das_innov_dorauv.F das_cost_dorauv.F das_adj_dauv.F\
	das_innov_hfradar.F das_cost_hfradar.F das_adj_chfradar.F\
	das_copy.F das_setup.F das_get_initial.F das_cost_reg.F das_adj_cost_reg.F\
\
	das_nf_fread.F\
\
	das_copy_smooth.F das_innov_smooth.F das_diag_ssh_steric.F\
	das_zeta_smooth.F das_copy_uv_smooth.F das_innov_uv_smooth.F\
	das_zeta_s_smooth.F das_adj_zeta_s_smooth.F\
	das_ts_s_smooth.F das_adj_ts_s_smooth.F\
	das_copy_ts_s_smooth.F das_adj_copy_ts_s_smooth.F das_ts_das_smooth.F\
	das_copy_zeta_smooth.F das_adj_copy_zeta_smooth.F\
	das_zeta_das_smooth.F\
	das_get_TS_profile.F das_innov_TS_profile.F das_cost_TS_profile.F\
	das_adj_TS_profile.F\
\
	MPI_Setup.F	MessPass2D.F	MessPass3D.F	exchange.F\
	autotiling.F



 SRCS90 = \
    buffer.F90\
    toolorigindate.F90\
    tooldatosec.F90 toolsectodat.F90 tooldecompdat.F90 tooldatetosec.F90\
\
\
\


#======================================================================

 RCS = $(SRCS:.F=.f)

 RCS90 = $(SRCS90:.F90=.f90)

 OBJS = $(RCS:.f=.o)
 OBJS90 = $(RCS90:.f90=.o)

# SBIN is set by the variant Makefile

%.mod : %.o

.PHONY: all clean clobber
all: $(SBIN) testkeys

#======================================================================

$(SBIN): cppdefs.h param.h scalars.h ncscrum.h \
         $(OBJS) $(OBJS90)
	$(LDR) $(FFLAGS) -o $(SBIN) $(OBJS) $(OBJS90) $(LCDF) $(LDFLAGS)

$(OBJS): cppdefs.h param.h
$(OBJS90): cppdefs.h param.h

testkeys: $(SBIN)
	@echo ""
	@echo '********************************************'
	@echo '****** MS3DVAR-$(VARIANT_NAME) Compilation Done ******'
	@echo '********************************************'
	@echo ""

#======================================================================
# Cleaning targets
#======================================================================

clean:
	/bin/rm -f core *.o *.f *.f90 *.trace *.mod
	/bin/rm -rf *.dSYM ./rii_files

clobber: clean
	/bin/rm -f $(SBIN) *.exe

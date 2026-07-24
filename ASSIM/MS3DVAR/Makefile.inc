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
# Common Makefile for MS3DVAR scales
# This file is included by each scale's Makefile
#
# Required variables to be set by scale Makefile:
#   SCALE_NAME    - Name of scale (e.g., "lr", "mr", "filter", "ms")
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
# Scale-specific files are in current directory
# Common MS3DVAR files are in $(COMMON_DIR)
# CROCO files are in $(CROCO_SRC)
#======================================================================

 SRCS = ms3dvar_main.F		read_inp.F\
	init_scalars.F	init_arrays.F\
	set_scoord.F	setup_grid1.F	setup_grid2.F\
	ana_initial.F	analytical.F\
	set_depth.F\
	rho_eos.F\
	diag.F		checkdims.F	grid_stiffness.F\
	debug.F\
        param.F         ncscrum.F   scalars.F\
	put_global_atts.F\
	setup_kwds.F	check_kwds.F	check_srcs.F\
	check_switches1.F	check_switches2.F\
	nf_fread.F nf_fread_x.F	nf_fread_y.F\
	lenstr.F	closecdf.F	insert_node.F\
	fillvalue.F	nf_add_attribute.F		set_cycle.F\
	def_grid_2d.F   def_grid_3d.F   def_his.F       def_rst.F\
	wrt_grid.F      wrt_rst.F\
	get_grid.F	get_initial.F\
	get_ssh.F\
	ms3dvar_init_arrays.F  ms3dvar_init_scalars.F\
	ms3dvar_setup1.F  ms3dvar_setup2.F ms3dvar_spline1d.F ms3dvar_exchange.F\
	ms3dvar_MessPass2D.F ms3dvar_s2z.F ms3dvar_z2s.F ms3dvar_analysis.F ms3dvar_analysis_zeta.F\
	ms3dvar_ana_bvar.F ms3dvar_costb.F ms3dvar_vec.F ms3dvar_adj_costb.F\
    ms3dvar_get_d2c.F\
	ms3dvar_get_tmi.F  ms3dvar_get_mc.F ms3dvar_get_goes.F ms3dvar_lbfgs.F ms3dvar_innov_sst.F\
	ms3dvar_cost_sst.F ms3dvar_adj_csst.F ms3dvar_ca2zero.F ms3dvar_adj_zero.F\
	ms3dvar_adj_zero_s.F ms3dvar_read_obs_filelist.F ms3dvar_matvec.F ms3dvar_adj_matvec.F\
	ms3dvar_get_corr.F ms3dvar_varvec.F ms3dvar_adj_varvec.F\
	ms3dvar_rho_eos9.F ms3dvar_fill9.F ms3dvar_extrapolate9.F ms3dvar_geos_uv_p.F ms3dvar_get_bvar.F\
    ms3dvar_diag_steric_sshb.F\
	ms3dvar_get_ssh_ref.F \
	ms3dvar_get_zeta_ubal_from_tot.F ms3dvar_adj_get_zeta_ubal_from_tot.F\
\
	ms3dvar_zeta_p_smooth.F ms3dvar_zeta_p_smooth2.F ms3dvar_diag_pressure.F\
	ms3dvar_adj_diag_p.F ms3dvar_get_hfradar.F ms3dvar_geo_ratio.F ms3dvar_adj_geo_ratio.F\
	ms3dvar_cost_reg_uv.F ms3dvar_adj_cost_reg_uv.F\
    ms3dvar_d2c_uv.F\
\
	ms3dvar_get_js1.F ms3dvar_get_swot.F ms3dvar_innov_js1.F ms3dvar_innov_swot.F ms3dvar_rho_eos_lin.F\
	ms3dvar_get_tp.F ms3dvar_innov_tp.F\
	ms3dvar_diag_ssh.F ms3dvar_cost_js1.F ms3dvar_cost_swot.F ms3dvar_wrt.F ms3dvar_adj_wrt.F \
    ms3dvar_theadj_wrt.F\
	ms3dvar_adj_rho_eos.F ms3dvar_adj_diag_ssh.F ms3dvar_adj_cjs1.F ms3dvar_adj_cswot.F\
	ms3dvar_adj_diag_ssh_steric.F\
	ms3dvar_adj_geos_uv_p.F ms3dvar_psichi_uv.F ms3dvar_adj_psichi_uv.F\
	ms3dvar_cross_ts.F ms3dvar_adj_cross_ts.F ms3dvar_set_geo_ratio.F\
	ms3dvar_zetah_smooth.F ms3dvar_adj_zetah_smooth.F ms3dvar_copy_zetah_smooth.F\
    ms3dvar_zeta_h_smooth.F ms3dvar_adj_zeta_h_smooth.F\
	ms3dvar_adj_copy_zeteh_smooth.F ms3dvar_adj_glider_sio.F ms3dvar_cost_glider_sio.F\
	ms3dvar_adj_glider_whoi.F ms3dvar_cost_glider_whoi.F ms3dvar_adj_map.F ms3dvar_cost_map.F\
	ms3dvar_adj_glider_ptsur.F ms3dvar_cost_glider_ptsur.F\
	ms3dvar_pressue_smooth.F ms3dvar_adj_pressure_smooth.F\
	ms3dvar_copy_pressure_smooth.F ms3dvar_adj_copy_pressure_smooth.F\
	ms3dvar_get_sats.F ms3dvar_innov_sats.F ms3dvar_cost_sats.F ms3dvar_adj_csats.F\
\
	ms3dvar_get_insitu.F ms3dvar_innov_sio.F ms3dvar_cost_glider.F ms3dvar_adj_cglider.F\
	ms3dvar_innov_whoi.F ms3dvar_innov_flight.F ms3dvar_cost_flight.F ms3dvar_adj_cflight.F\
	ms3dvar_innov_ptsur.F ms3dvar_innov_martn.F ms3dvar_innov_map.F\
	ms3dvar_get_auv.F ms3dvar_innov_calauv.F ms3dvar_cost_calauv.F ms3dvar_adj_cauv.F\
	ms3dvar_innov_dorauv.F ms3dvar_cost_dorauv.F ms3dvar_adj_dauv.F\
	ms3dvar_innov_hfradar.F ms3dvar_cost_hfradar.F ms3dvar_adj_chfradar.F\
	ms3dvar_copy.F ms3dvar_setup.F ms3dvar_cost_reg.F ms3dvar_adj_cost_reg.F\
	ms3dvar_add_lr_s.F ms3dvar_add_zeta_lr_s.F ms3dvar_ca2zero_s.F\
	ms3dvar_interp_lr.F ms3dvar_nf_fread_lr.F\
\
	ms3dvar_nf_fread.F\
\
	ms3dvar_copy_smooth.F ms3dvar_innov_smooth.F ms3dvar_diag_ssh_steric.F\
	ms3dvar_zeta_smooth.F ms3dvar_copy_uv_smooth.F ms3dvar_innov_uv_smooth.F\
	ms3dvar_zeta_s_smooth.F ms3dvar_adj_zeta_s_smooth.F\
	ms3dvar_ts_s_smooth.F ms3dvar_adj_ts_s_smooth.F\
	ms3dvar_copy_ts_s_smooth.F ms3dvar_adj_copy_ts_s_smooth.F ms3dvar_ts_analysis_smooth.F\
	ms3dvar_copy_zeta_smooth.F ms3dvar_adj_copy_zeta_smooth.F\
	ms3dvar_zeta_analysis_smooth.F\
	ms3dvar_get_TS_profile.F ms3dvar_innov_TS_profile.F ms3dvar_cost_TS_profile.F\
	ms3dvar_adj_TS_profile.F\
\
	MPI_Setup.F	MessPass2D.F	MessPass3D.F	exchange.F



 SRCS90 = \
    buffer.F90\
    toolorigindate.F90\
    tooldatosec.F90 toolsectodat.F90 tooldecompdat.F90 tooldatetosec.F90\
\
\
\


#======================================================================

 TOOLS = mpc cross_matrix cppcheck srcscheck checkkwds
 FULL_TOOLS = mpc cross_matrix cppcheck srcscheck checkkwds partit ncjoin ncrename


 RCS = $(SRCS:.F=.f)

 RCS90 = $(SRCS90:.F90=.f90)

 OBJS = $(RCS:.f=.o)
 OBJS90 = $(RCS90:.f90=.o)

# SBIN is set by the scale Makefile

%.mod : %.o
	@touch $@

.PHONY: all tools clean clobber

#
# Everything
# ==========
all: tools depend $(SBIN) okmsg


#======================================================================

# Executable rule :
# - removing headers deps (included at precompilation stage)
$(SBIN): $(OBJS) $(OBJS90)
	$(LDR) $(FFLAGS) -o $(SBIN) $(OBJS) $(OBJS90) $(LCDF) $(LDFLAGS)


#======================================================================
# Tools for generating keyword/switch check source files
# (mirrors CROCO OCEAN Makefile approach
# Use FULL_TOOLS instead of TOOLS to compile ncjoin, partit, ncrename
#======================================================================

tools: $(TOOLS)

mpc: mpc.F
	$(CPP) -P $(CPPFLAGS) mpc.F > mpc_.f
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpc mpc_.f

cross_matrix: mpc cross_matrix.F
	$(CPP) -P $(CPPFLAGS) cross_matrix.F | ./mpc > cross_matrix_.f
	$(CFT) -c $(FFLAGS) cross_matrix_.f -o cross_matrix.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cross_matrix cross_matrix.o

cppcheck: mpc cppcheck.F
	$(CPP) -P $(CPPFLAGS) cppcheck.F | ./mpc > cppcheck_.f
	$(CFT) -c $(FFLAGS) cppcheck_.f -o cppcheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cppcheck cppcheck.o

srcscheck: mpc srcscheck.F
	$(CPP) -P $(CPPFLAGS) srcscheck.F | ./mpc > srcscheck_.f
	$(CFT) -c $(FFLAGS) srcscheck_.f -o srcscheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck.o
	rm -f check_srcs.F

checkkwds: mpc checkkwds.F
	$(CPP) -P $(CPPFLAGS) checkkwds.F | ./mpc > checkkwds_.f
	$(CFT) -c $(FFLAGS) checkkwds_.f -o checkkwds.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o checkkwds checkkwds.o
	rm -f setup_kwds.F

checkdefs: check_switches1.F setup_kwds.F

# Adding the newly included cppdefs file type
check_switches1.F: cppcheck cppdefs.h
	cat croco_cppdefs.h cppdefs.h ms3dvar_cppdefs.h cppdefs_dev.h > mergcpp.txt
	./cppcheck

check_srcs.F: srcscheck Makefile
	./srcscheck

setup_kwds.F: checkkwds read_inp.F
	./checkkwds

partit: partit.F insert_node.o lenstr.o
	$(CPP) -P $(CPPFLAGS) partit.F > partit._.f
	$(CFT) -c $(FFLAGS) partit._.f -o partit.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o  lenstr.o $(LCDF)

ncjoin: ncjoin.F lenstr.o
	$(CPP) -P $(CPPFLAGS) ncjoin.F | ./mpc  > ncjoin._.f
	$(CFT) -c $(FFLAGS) ncjoin._.f -o ncjoin.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin ncjoin.o lenstr.o $(LCDF)


ncrename: ncrename.F lenstr.o
	$(CPP) -P $(CPPFLAGS) ncrename.F | ./mpc > ncrename._.f
	$(CFT) -c $(FFLAGS) ncrename._.f -o ncrename.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncrename ncrename.o lenstr.o $(LCDF)

depend: checkdefs cross_matrix
	./cross_matrix  *.F90 *.F
	
okmsg: $(SBIN)
	@echo ""
	@echo '********************************************'
	@echo '****** MS3DVAR-$(SCALE_NAME) Compilation Done ******'
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
-include Make.depend

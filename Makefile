# $Id: Makefile 1483 2014-03-15 17:05:10Z rblod $
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
# Universal machine independent makefile for ROMS model 
#
#======================================================================
# Set machine dependent definitions and rules.
#======================================================================

%.tap.f: %.F
	$(CPP) -P $(CPPFLAGS) -D__TAPENADE__ $*.F | ./mpc > $*.tap.f


include Makedefs

#======================================================================
# Model Configuration:
#======================================================================
# SRCS: source codes files are sorted into groups, separated by
# blanc lines:
#   1) main driving part;		 2) 2D time stepping engine;
#   3) 3D time stepping engine;		 4) non-hydrostatic engine
#   5) sea-water EOS;			 6) vertical mixing schemes;    
#   7) on-fly model diagnostics;	 8) netCDF I/O routines;	   
#   9) main forcing;			10) Wave forcing routines;
#  11) Online surface forcing;		12) Floats routines;
#  13) Station diagnostics;		14) biology/sediment routines;
#  15) PISCES biology;			16) MPI routines;
#  17) AGRIF routines;			18) OASIS coupling interface;
#
# SRC90: additional F90 routines for PISCES code and Non-hydrostatic
#        analysis routines
#======================================================================

 SRCS = 		step.F		read_inp.F\
	timers_roms.F	init_scalars.F	init_arrays.F	set_weights.F\
	set_scoord.F	ana_grid.F	setup_grid1.F	setup_grid2.F\
	set_nudgcof.F	ana_initial.F	analytical.F	zonavg.F\
\
	step2d.F 	u2dbc.F		v2dbc.F		zetabc.F\
	obc_volcons.F\
\
	pre_step3d.F	step3d_t.F	step3d_uv1.F	step3d_uv2.F\
	prsgrd.F	rhs3d.F		set_depth.F	omega.F\
	uv3dmix.F	uv3dmix_spg.F	t3dmix.F	t3dmix_spg.F\
	hmix_coef.F	wetdry.F\
	u3dbc.F		v3dbc.F		t3dbc.F\
\
	step3d_fast.F\
	step3d_w.F	rhs3d_w_nh.F	initial_nbq.F	grid_nbq.F\
	unbq_bc.F	vnbq_bc.F	wnbq_bc.F	rnbq_bc.F\
        w3dbc.F		nbq_bry_store.F\
\
	rho_eos.F	ab_ratio.F	alfabeta.F\
\
	ana_vmix.F	bvf_mix.F	lmd_vmix.F 	gls_mixing.F\
	lmd_skpp.F	lmd_bkpp.F	lmd_swfrac.F	lmd_wscale.F\
\
	diag.F		wvlcty.F	checkdims.F	grid_stiffness.F\
	bio_diag.F	setup_kwds.F    check_kwds.F	check_srcs.F\
	check_switches1.F		check_switches2.F\
	debug.F\
\
        param.F         ncscrum.F   scalars.F\
	output.F	put_global_atts.F\
	nf_fread.F	nf_fread_x.F	nf_fread_y.F	nf_read_bry.F\
	get_date.F	lenstr.F	closecdf.F	insert_node.F\
	fillvalue.F	nf_add_attribute.F		set_cycle.F\
	def_grid_2d.F   def_grid_3d.F   def_his.F       def_rst.F\
	def_diags.F	def_diagsM.F	def_bio_diags.F\
	wrt_grid.F      wrt_his.F       wrt_avg.F	wrt_rst.F\
	wrt_diags.F	wrt_diags_avg.F	wrt_diagsM.F	wrt_diagsM_avg.F\
	wrt_bio_diags.F	wrt_bio_diags_avg.F\
	set_avg.F	set_diags_avg.F	set_diagsM_avg.F\
	set_bio_diags_avg.F\
	def_diags_vrt.F	wrt_diags_vrt.F\
	set_diags_vrt.F	set_diags_vrt_avg.F	wrt_diags_vrt_avg.F\
	def_diags_ek.F	wrt_diags_ek.F\
	set_diags_ek.F	set_diags_ek_avg.F	wrt_diags_ek_avg.F\
	def_diags_pv.F	wrt_diags_pv.F\
	set_diags_pv.F	set_diags_pv_avg.F	wrt_diags_pv_avg.F\
	def_diags_eddy.F\
	set_diags_eddy_avg.F	wrt_diags_eddy_avg.F\
	def_surf.F	wrt_surf.F\
	set_surf_avg.F	wrt_surf_avg.F\
\
	get_grid.F	get_initial.F	get_vbc.F	get_wwave.F\
	get_tclima.F    get_uclima.F    get_ssh.F       get_sss.F\
	get_smflux.F    get_stflux.F    get_srflux.F    get_sst.F\
	mod_tides_mas.F tidedata.F      mas.F         get_tides.F\
        clm_tides.F     get_bulk.F      bulk_flux.F\
	get_bry.F       get_bry_bio.F	sstskin.F\
	get_psource.F   get_psource_ts.F\
        cfb_stress.F\
\
	mrl_wci.F  	wkb_wwave.F 	wkbbc.F		get_bry_wkb.F\
\
	online_bulk_var.F		online_get_bulk.F\
	online_interp.F			online_interpolate_bulk.F\
	online_set_bulk.F\
\
	init_floats.F	wrt_floats.F	step_floats.F	rhs_floats.F\
	interp_rho.F	def_floats.F	init_arrays_floats.F\
	random_walk.F	get_initial_floats.F\
\
	init_sta.F	wrt_sta.F	step_sta.F	interp_sta.F\
	def_sta.F	init_arrays_sta.F\
\
	biology.F	o2sato.F	sediment.F	bbl.F\
\
	MPI_Setup.F	MessPass2D.F	MessPass3D.F	exchange.F\
	autotiling.F\
\
	zoom.F		update2D.F	set_nudgcof_fine.F\
	zoombc_2D.F	zoombc_3D.F	uv3dpremix.F\
	t3dpremix.F     update3D.F	zoombc_3Dfast.F\
	Agrif2Model.F\
\
	send_xios_diags.F\
\
	cpl_prism_define.F	cpl_prism_put.F 	cpl_prism_init.F\
	cpl_prism_get.F 	cpl_prism_getvar.F      cpl_prism_grid.F\

 SRCS90 = \
	par_pisces.F90	ocean2pisces.F90	trc.F90		sms_pisces.F90\
	p4zche.F90	p4zint.F90		p4zlys.F90	p4zflx.F90\
	p4zlim.F90	p4zsink.F90		p4zmicro.F90	p4zmeso.F90\
	p4zmort.F90	p4zopt.F90		p4zprod.F90	p4zrem.F90\
	p4zsed.F90	p4zbio.F90		trcwri_pisces.F90\
	trcsms_pisces.F90			trcini_pisces.F90\
	pisces_ini.F90\
\
	module_parameter_oa.F90\
	module_oa_time.F90	module_oa_space.F90	module_oa_periode.F90\
	module_oa_variables.F90	module_oa_type.F90	module_oa_stock.F90\
	module_oa_level.F90	module_oa_interface.F90\
	module_oa_upd.F90	croco_oa.F90		var_oa.F90\
\
	tooldatosec.F90	toolsectodat.F90 tooldecompdat.F90

AMRDIR = AGRIF/AGRIF_YOURFILES

#======================================================================

 TOOLS = mpc cross_matrix cppcheck srcscheck checkkwds partit ncjoin ncrename


 RCS = $(SRCS:.F=.f)

 RCS90 = $(SRCS90:.F90=.f90)

 OBJS = $(RCS:.f=.o)
 OBJS90 = $(RCS90:.f90=.o)

 SBIN = croco

 AMRRCS=$(AMRSRCS:.F=.f)

 AMROBJS=$(AMRRCS:.f=.o)


%.mod : %.o
	@touch $@

#======================================================================
#
# Tapenade:
# =========
#

ADJ_SRCS=cost_fun.F get_vbc.F step.F step2d.F pre_step3d.F set_depth.F grid_stiffness.F rho_eos.F ana_vmix.F omega.F prsgrd.F rhs3d.F step3d_uv1.F step3d_uv2.F step3d_t.F t3dbc.F u3dbc.F v3dbc.F v2dbc.F u2dbc.F exchange.F analytical.F MessPass2D.F zetabc.F set_avg.F debug.F dummy.F get_ij.F distance.F xtime.F
ADJ_PSRCS=$(ADJ_SRCS:.F=.tap.f)
TAP_TARGET=autodiff
ADJ_OBJS=$(TAP_TARGET)_b.o m1qn3.o treeverse.o adBinomial.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o adj_driver.o cost_fun.o get_ij.o distance.o xtime.o code_insertion.o

TGT_SRCS= $(ADJ_SRCS) ana_initial.F ana_grid.F wrt_his.F def_his.F setup_grid1.F setup_grid2.F insert_node.F checkdims.F put_global_atts.F def_grid_3d.F wrt_grid.F lenstr.F nf_fread.F fillvalue.F wvlcty.F nf_add_attribute.F init_scalars.F timers_roms.F init_arrays.F get_date.F nf_fread_y.F nf_fread_x.F set_scoord.F set_weights.F get_initial.F closecdf.F cost_driver.F analytical.F
TGT_PSRCS=$(TGT_SRCS:.F=.tap.f)
TGT_PSRCS_DBG= main_cost.tap.f $(TGT_PSRCS)
TGT_OBJS=$(TAP_TARGET)_d.o m1qn3.o treeverse.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o tgt_driver.o cost_fun.o get_ij.o distance.o xtime.o code_insertion.o
TGT_OBJS_DBG=$(TAP_TARGET)_context_d.o m1qn3.o adContext.o treeverse.o adBufferC.o adBuffer.o  adStack.o code_insertion.o  read_obs.o analytical.o t3dbc.o v3dbc.o exchange.o zetabc.o u2dbc.o v2dbc.o pre_step3d.o u3dbc.o get_vbc.o diag.o ana_vmix.o output.o  check_kwds.o  read_inp.o init_scalars.o timers_roms.o init_arrays.o ana_grid.o setup_grid1.o setup_grid2.o set_scoord.o set_weights.o set_depth.o grid_stiffness.o get_initial.o ana_initial.o set_depth.o rho_eos.o omega.o wrt_his.o step.o closecdf.o wrt_rst.o lenstr.o

ADJ_OBJS_DBG1=$(TAP_TARGET)_context1_d.o m1qn3.o adDebug.o treeverse.o adBufferC.o adBuffer.o  adStack.o code_insertion.o  read_obs.o analytical.o t3dbc.o v3dbc.o exchange.o zetabc.o u2dbc.o v2dbc.o pre_step3d.o u3dbc.o get_vbc.o diag.o ana_vmix.o output.o  check_kwds.o  read_inp.o init_scalars.o timers_roms.o init_arrays.o ana_grid.o setup_grid1.o setup_grid2.o set_scoord.o set_weights.o set_depth.o grid_stiffness.o get_initial.o ana_initial.o set_depth.o rho_eos.o omega.o wrt_his.o step.o closecdf.o wrt_rst.o lenstr.o  cost_fun.o

ADJ_OBJS_DBG2=$(TAP_TARGET)_context2_b.o m1qn3.o adDebug.o treeverse.o adBufferC.o adBuffer.o  adStack.o code_insertion.o  read_obs.o analytical.o t3dbc.o v3dbc.o exchange.o zetabc.o u2dbc.o v2dbc.o pre_step3d.o u3dbc.o get_vbc.o diag.o ana_vmix.o output.o  check_kwds.o  read_inp.o init_scalars.o timers_roms.o init_arrays.o ana_grid.o setup_grid1.o setup_grid2.o set_scoord.o set_weights.o set_depth.o grid_stiffness.o get_initial.o ana_initial.o set_depth.o rho_eos.o omega.o wrt_his.o step.o closecdf.o wrt_rst.o lenstr.o  cost_fun.o



DIV_OBJS=m1qn3.o treeverse.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o div_driver.o cost_fun.o get_ij.o distance.o xtime.o

TGT_CONTEXT_OBJS=$(TAP_TARGET)_d.o cost_fun.o contextAD.o

#======================================================================

#
# Everything
# ==========
all: tools depend $(SBIN) $(SBIN)_tgt_dbg $(SBIN)_adj_dbg1 $(SBIN)_adj_dbg2

#
# Executables files.
# =========== =====
#
$(SBIN): $(OBJS90) $(OBJS) main.o $(MPI_DIR_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) #-lampiPlainC

$(SBIN)_adj:  $(ADJ_OBJS) $(OBJS90) $(OBJS) main_adj.o $(MPI_ADJ_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape  -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_adc:  $(ADJ_OBJS) $(TAP_TARGET)_d.o check_driver.o $(OBJS90) $(OBJS) main_adc.o $(MPI_ADJ_OBJS) 
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape   -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_tgt: $(TGT_OBJS) $(OBJS90) $(OBJS) main_adj.o $(MPI_TGT_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape   -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_tgt_dbg: $(TGT_OBJS_DBG) $(OBJS90) $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape   -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_adj_dbg1: $(ADJ_OBJS_DBG1) $(OBJS90) $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape   -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_adj_dbg2: $(ADJ_OBJS_DBG2) $(OBJS90) $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape   -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_div: $(DIV_OBJS) $(OBJS90) $(OBJS) main_adj.o $(MPI_TGT_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lblas #-lampiCommon  -lampiTape   -lampiBookkeeping -lblas -lampiPlainC

# $Id: Makefile 3922 2011-05-19 08:54:39Z llh $

profile.o: profile.c
	$(CC) $(FLAGS) -c profile.c

testMemSizef.o : testMemSizef.f
	$(CFT) $(FFLAGS) -c testMemSizef.f

testMemSizec.o : testMemSizec.c
	$(CC) $(FFLAGS) -c testMemSizec.c

testMemSize : testMemSizef.o testMemSizec.o
	$(CFT) $(FFLAGS) testMemSizef.o testMemSizec.o -o testMemSize

treeverse.f: treeverse.F
	/bin/cp $^ $@

adBuffer.f: adBuffer.F
	/bin/cp $^ $@

treeverse.o: treeverse.f
	$(CFT) $(FFLAGS) -c $^

adBinomial.o: adBinomial.c
	$(CC) $(FFLAGS) -c $^

adBuffer.o : adBuffer.f
	$(CFT) $(FFLAGS) -c $^

adBufferC.o: adBufferC.c
	$(CC) $(FFLAGS) -c $^

adStack.o : adStack.c
	$(CC) $(FFLAGS) -c adStack.c

adDebug.o : adDebug.c
	$(CC) $(FFLAGS) -c adDebug.c

dpStack.o : dpStack.c
	$(CC) $(FFLAGS) -c dpStack.c

dpTest.o : dpTest.f
	$(CFT) $(FFLAGS) -c dpTest.f

validityTest.o : validityTest.f
	$(CFT) $(FFLAGS) -c validityTest.f

treeverseFtest : treeverseFtest.f treeverse.f
	$(CFT) $(FFLAGS) -c treeverseFtest.f
	$(CFT) $(FFLAGS) -c treeverse.f
	$(CFT) $(FFLAGS) $(LDFLAGS) treeverseFtest.o treeverse.o -o treeverseFtest

treeverseCtest : treeverseCtest.c treeverse.c treeverse.h
	$(CC) $(FFLAGS) -c treeverse.c
	$(CC) $(FFLAGS) -c treeverseCtest.c
	$(CC) $(FFLAGS) $(LDFLAGS) treeverseCtest.o treeverse.o -o treeverseCtest

adBufferFtest : adStack.c adBuffer.f adBufferFtest.f
	$(CFT) $(FFLAGS) $^ -o $@

adBufferCtest : adStack.c adBuffer.c adBufferCtest.c
	$(CFT) $(FFLAGS) $^ -o $@

m1qn3.o: m1qn3.F
	$(CFT) $(FFLAGS) -c $^ -o $@


pre_step3d.tap.f: pre_step3d.F
	$(CPP) -P $(CPPFLAGS) -D__TAPENADE__ pre_step3d.F > pre_step3d.1.f
	emacs -Q --script fsplit.el pre_step3d.1.f pre_step3d.tap.f	

pre_step3d.f: pre_step3d.F
	$(CPP) -P $(CPPFLAGS) pre_step3d.F > pre_step3d.2.f
	emacs -Q --script fsplit.el pre_step3d.2.f pre_step3d.f

pre_step3d.o: pre_step3d.f
	$(CFT) -c $(FFLAGS) pre_step3d.f -o $@


#
# Auxiliary utility programs and List of Dependecies:
# ========= ======= ======== === ==== == ============
#

tools: $(TOOLS)

mpc: mpc.F
	$(CPP) -P $(CPPFLAGS) mpc.F > mpc_.f
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpc mpc_.f

cross_matrix: mpc cross_matrix.F 
	$(CPP) -P $(CPPFLAGS) cross_matrix.F | ./mpc > cross_matrix_.f
	$(CFT) -c $(FFLAGS) cross_matrix_.f -o cross_matrix.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cross_matrix cross_matrix.o

cppcheck:  mpc cppcheck.F 
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

check_switches1.F: cppcheck cppdefs.h
	cat cppdefs.h cppdefs_dev.h > mergcpp.txt
	./cppcheck
check_srcs.F: srcscheck Makefile
	./srcscheck
setup_kwds.F: checkkwds read_inp.F
	./checkkwds
partit: partit.o insert_node.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o lenstr.o $(LCDF)

ncjoin: ncjoin.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin ncjoin.o lenstr.o $(LCDF)

ncrename: ncrename.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncrename ncrename.o lenstr.o $(LCDF)

depend: checkdefs cross_matrix
	./cross_matrix  *.F90 *.F 

mymodules: $(MOBJS) $(AMROBJS)
#
# Target to create tar file.
# ====== == ====== === =====
#
tarfile:
	tar cvf croco.tar Make* *.h *.F *.F90 *.in *.in.*
#
# Cleaning:
# =========
#
rmtools:
	/bin/rm -f $(TOOLS)
clean:
	/bin/rm -rf core *.o $(AMRDIR)/*.o *.i *.s *.f *.f90 *.trace *.mod ${COMP_FILES}

clobber: clean
	/bin/rm -rf $(SBIN) $(TOOLS) ./rii_files


plotter: plotter.F
	f77 -n32 -o plotter plotter.F $(LIBNCAR)

$(TAP_TARGET)_context_d.o: $(TAP_TARGET)_context_d.f
	ln -sf empty_code_insertion.h code_insertion.h
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

$(TAP_TARGET)_context1_d.o: $(TAP_TARGET)_context1_d.f
	ln -sf empty_code_insertion.h code_insertion.h
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

$(TAP_TARGET)_context2_b.o: $(TAP_TARGET)_context2_b.f
	ln -sf empty_code_insertion.h code_insertion.h
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

$(TAP_TARGET)_d.o: $(TAP_TARGET)_d.f
	ln -sf empty_code_insertion.h code_insertion.h
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

$(TAP_TARGET)_b.o: $(TAP_TARGET)_b.f
	ln -sf adtool_ampi_turn_code_insertion.h code_insertion.h
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

$(TAP_TARGET)_b.f: $(ADJ_PSRCS)
	ln -sf empty_code_insertion.h code_insertion.h
	${TAPENADE} $^ -noisize -noisize77 -msglevel 10 -msginfile -nocheckpoint "step3d_uv_thread step3d_t_thread omega_tile rho_eos rho_eos_tile set_vbc prsgrd rhs3d pre_step3d pre_step3d_tile step3d_uv2_tile set_depth set_depth_tile set_huv set_huv_tile set_huv2 set_huv2_tile exchange_r2d_tile exchange_u2d_tile exchange_v2d_tile exchange_u3d_tile exchange_v3d_tile exchange_r3d_tile exchange_r2d_tile exchange_w3d_tile prsgrd_tile rhs3d_tile" -head "cost_fun(ad_x)\(cost)" -r8 -reverse -output $(TAP_TARGET) $(AMPIINC)
	ln -sf adtool_ampi_turn_code_insertion.h code_insertion.h  
	sed -i 's/REAL, DIMENSION(\*, \*, \*)/REAL, DIMENSION(:, :, :)/g' $(TAP_TARGET)_b.f 
	sed -i 's/REAL, DIMENSION(:, :, :), POINTER a/REAL, DIMENSION(:, :, :), POINTER :: a/g' $(TAP_TARGET)_b.f
	sed -i 's/PARAMETER n/PARAMETER :: n/g' $(TAP_TARGET)_b.f 
	sed -i 's/PARAMETER m/PARAMETER :: m/g' $(TAP_TARGET)_b.f

cmaker.f: cmaker.F
	$(CPP) -P $(CPPFLAGS) $^ | ./mpc > $@

cmaker: cmaker.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lampiPlainC

main_tgt.f: main.F
	$(CPP) -P $(CPPFLAGS) -DTANGENT_CHECK $^ | ./mpc > $@

main_adj.f: main.F
	$(CPP) -P $(CPPFLAGS) -DSTATE_CONTROL $^ | ./mpc > $@

main_adc.f: main.F
	$(CPP) -P $(CPPFLAGS) -DSTATE_CONTROL -DAD_CHECK $^ | ./mpc > $@

main_cost.tap.f: main.F
	$(CPP) -P $(CPPFLAGS) -DSTATE_CONTROL -DAD_COST $^ | ./mpc > $@


$(TAP_TARGET)_d.f: $(TGT_PSRCS) #main_tgt.f
	ln -sf empty_code_insertion.h code_insertion.h
	${TAPENADE} $^ -noisize -noisize77 -tracelevel 10 -msglevel 20 -msginfile -head "cost_fun(cost)/(ad_x)" -r8 -output $(TAP_TARGET) $(AMPIINC)
	ln -sf adtool_ampi_turn_code_insertion.h code_insertion.h
	sed -i 's/REAL, DIMENSION(\*, \*, \*)/REAL, DIMENSION(:, :, :)/g' $(TAP_TARGET)_d.f 
	sed -i 's/REAL, DIMENSION(:, :, :), POINTER a/REAL, DIMENSION(:, :, :), POINTER :: a/g' $(TAP_TARGET)_d.f

$(TAP_TARGET)_context_d.f: $(TGT_PSRCS_DBG)
	ln -sf empty_code_insertion.h code_insertion.h
	${TAPENADE} -d -fixinterface -context -noisize -noisize77 -tracelevel 10 -msglevel 20 -msginfile -head "cost_fun(cost)/(ad_x)" -r8  -output $(TAP_TARGET)_context $(AMPIINC) $^
	ln -sf adtool_ampi_turn_code_insertion.h code_insertion.h
	sed -i 's/REAL, DIMENSION(\*, \*, \*)/REAL, DIMENSION(:, :, :)/g' $(TAP_TARGET)_context_d.f 
	sed -i 's/REAL, DIMENSION(:, :, :), POINTER a/REAL, DIMENSION(:, :, :), POINTER :: a/g' $(TAP_TARGET)_context_d.f
	sed -i 's/CHARACTER\*(/CHARACTER(/g' $(TAP_TARGET)_context_d.f	
	sed -i 's/\(ADDEBUGTGT_CALL([^,]*\), 0, 0/\1, 10, 10/g' $(TAP_TARGET)_context_d.f

$(TAP_TARGET)_context1_d.f: $(TGT_PSRCS_DBG)
	ln -sf empty_code_insertion.h code_insertion.h
	${TAPENADE} -d -fixinterface -context -debugADJ -noisize -noisize77 -msginfile -head "cost_fun(cost)/(ad_x)" -r8  -output $(TAP_TARGET)_context1 $(AMPIINC) $^
	ln -sf adtool_ampi_turn_code_insertion.h code_insertion.h
	sed -i 's/REAL, DIMENSION(\*, \*, \*)/REAL, DIMENSION(:, :, :)/g' $(TAP_TARGET)_context1_d.f 
	sed -i 's/REAL, DIMENSION(:, :, :), POINTER a/REAL, DIMENSION(:, :, :), POINTER :: a/g' $(TAP_TARGET)_context1_d.f
	sed -i 's/CHARACTER\*(/CHARACTER(/g' $(TAP_TARGET)_context_d.f	
	sed -i 's/\(ADDEBUGTGT_CALL([^,]*\), 0, 0/\1, 10, 10/g' $(TAP_TARGET)_context1_d.f

$(TAP_TARGET)_context2_b.f: $(TGT_PSRCS_DBG)
	ln -sf empty_code_insertion.h code_insertion.h
	${TAPENADE} -b -fixinterface -context -debugADJ -noisize -noisize77 -msginfile -nocheckpoint "step3d_uv_thread step3d_t_thread omega_tile rho_eos rho_eos_tile set_vbc prsgrd rhs3d pre_step3d pre_step3d_tile set_depth set_depth_tile set_huv set_huv_tile set_huv2 set_huv2_tile exchange_r2d_tile exchange_u2d_tile exchange_v2d_tile exchange_u3d_tile exchange_v3d_tile exchange_r3d_tile exchange_r2d_tile exchange_w3d_tile prsgrd_tile rhs3d_tile" -head "cost_fun(cost)/(ad_x)" -r8  -output $(TAP_TARGET)_context2 $(AMPIINC) $^
	ln -sf adtool_ampi_turn_code_insertion.h code_insertion.h
	sed -i 's/REAL, DIMENSION(\*, \*, \*)/REAL, DIMENSION(:, :, :)/g' $(TAP_TARGET)_context2_b.f 
	sed -i 's/REAL, DIMENSION(:, :, :), POINTER a/REAL, DIMENSION(:, :, :), POINTER :: a/g' $(TAP_TARGET)_context2_b.f
	sed -i 's/CHARACTER\*(/CHARACTER(/g' $(TAP_TARGET)_context_d.f	
	sed -i 's/\(ADDEBUGBWD_CALL([^,]*\), 0/\1, 10/g' $(TAP_TARGET)_context2_b.f
	sed -i 's/IF (.FALSE./IF (.TRUE./g' $(TAP_TARGET)_context2_b.f

fortranSupport.o : fortranSupport.F
	$(CFT) $(FFLAGS) $(CPPFLAGS) -I /usr/include -I /usr/local/include -c $^ -o $@

ampiSupport.o : ampiSupport.c
	$(CC) $(CPPFLAGS) -I /usr/include -I /usr/local/include -c $^ -o $@


# Special treatment for barrier function:
# THERE SHALL BE NO OPTIMIZATION HERE!!!!
#
my_barrier.o: my_barrier.f
	$(CFT) -c -O0 my_barrier.f
#
# Include automatically generated dependency list:
#

-include Make.depend


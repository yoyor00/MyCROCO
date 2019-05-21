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
#   1) main driving part;	   2) 2D time stepping engine;
#   3) 3D time stepping engine;	   4) sea-water EOS routines;
#   5) vertical mixing schemes;    6) on-fly model diagnostic routines;
#   7) netCDF I/O routines;	   8) main forcing routines;
#   9) Wave forcing routines;	  10) Online surface forcing routines;
#  11) Floats routines;		  12) Stations routines;
#  13) biology/sediment routines; 14) PISCES routines
#  15) MPI routines;		  16) AGRIF routines;           
#  17) non-hydrostatic engine     18) OASIS coupling interface;
#
# SRC90: additional F90 routines for the non-hydrostatic engine
#======================================================================

 SRCS = 		step.F		read_inp.F\
	timers_roms.F	init_scalars.F	init_arrays.F	set_weights.F\
	set_scoord.F	ana_grid.F	setup_grid1.F	setup_grid2.F\
	set_nudgcof.F	ana_initial.F	analytical.F	zonavg.F\
\
	step2d.F	u2dbc.F		v2dbc.F		zetabc.F\
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
	gls_prestep.F	gls_mix2017.F\
\
	diag.F		wvlcty.F	checkdims.F	grid_stiffness.F\
	bio_diag.F	setup_kwds.F    check_kwds.F	check_srcs.F\
	check_switches1.F		check_switches2.F\
	debug.F \
\
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
\
	get_grid.F	get_initial.F	get_vbc.F	get_wwave.F\
	get_tclima.F    get_uclima.F    get_ssh.F       get_sss.F\
	get_smflux.F    get_stflux.F    get_srflux.F    get_sst.F\
	mod_tides_mas.F tidedata.F      mas.F         get_tides.F\
        clm_tides.F     get_bulk.F      bulk_flux.F\
	get_bry.F       get_bry_bio.F	sstskin.F\
	get_psource.F   get_psource_ts.F\
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
	autotiling.F	MessPass3D_nbq.F\
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
	module_oa_parameter.F90\
	module_oa_time.F90	module_oa_space.F90	module_oa_periode.F90\
	module_oa_variables.F90	module_oa_type.F90	module_oa_stock.F90\
	module_oa_level.F90	module_oa_interface.F90\
	module_oa_upd.F90	croco_oa.F90		var_oa.F90\
\
	tooldatosec.F90	toolsectodat.F90 tooldecompdat.F90

AMRDIR = AGRIF/AGRIF_YOURFILES

#======================================================================

 RCS = $(SRCS:.F=.f)
 RCS90 = $(SRCS90:.F90=.f90)

 OBJS = $(RCS:.f=.o)
 OBJS90 = $(RCS90:.f90=.o)

 SBIN = croco

 AMRRCS=$(AMRSRCS:.F=.f)

 AMROBJS=$(AMRRCS:.f=.o)

#======================================================================
#
# Tapenade:
# =========
#

ADJ_SRCS=cost_fun.F step.F step2d.F  v2dbc.F u2dbc.F exchange.F analytical.F MessPass2D.F zetabc.F set_avg.F debug.F dummy.F get_ij.F distance.F xtime.F
ADJ_PSRCS=$(ADJ_SRCS:.F=.tap.f)
TAP_TARGET=autodiff
ADJ_OBJS=$(TAP_TARGET)_b.o m1qn3.o treeverse.o adBinomial.o adBufferC.o adBuffer.o adStack.o read_obs.o check_driver.o optim_driver.o adj_driver.o cost_fun.o get_ij.o distance.o xtime.o

TGT_SRCS=$(ADJ_SRCS)
TGT_PSRCS=$(TGT_SRCS:.F=.tap.f)
TGT_OBJS=$(TAP_TARGET)_d.o m1qn3.o treeverse.o adBufferC.o adBuffer.o adStack.o read_obs.o optim_driver.o tgt_driver.o cost_fun.o get_ij.o distance.o xtime.o

DIV_OBJS=m1qn3.o treeverse.o adBufferC.o adBuffer.o adStack.o read_obs.o optim_driver.o div_driver.o cost_fun.o get_ij.o distance.o xtime.o

TGT_CONTEXT_OBJS=$(TAP_TARGET)_d.o cost_fun.o contextAD.o

#======================================================================

#
# Everything
# ==========
all: tools depend $(SBIN) $(SBIN)_adj $(SBIN)_adc $(SBIN)_tgt $(SBIN)_div

#
# Executables files.
# =========== =====
#
$(SBIN): $(OBJS90) $(OBJS) main.o $(MPI_DIR_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lampiPlainC

$(SBIN)_adj:  $(ADJ_OBJS) $(OBJS90) $(OBJS) main_adj.o $(MPI_ADJ_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lampiCommon  -lampiTape -lampiADtoolStubsOO -lampiADtoolStubsST -lampiBookkeeping -lblas -lampiPlainC 

$(SBIN)_adc:  $(ADJ_OBJS) $(OBJS90) $(OBJS) main_adc.o $(MPI_ADJ_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lampiCommon  -lampiTape -lampiADtoolStubsOO -lampiADtoolStubsST -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_tgt: $(TGT_OBJS) $(OBJS90) $(OBJS) main_adj.o $(MPI_TGT_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lampiCommon  -lampiTape -lampiADtoolStubsOO -lampiADtoolStubsST -lampiBookkeeping -lblas -lampiPlainC

$(SBIN)_div: $(DIV_OBJS) $(OBJS90) $(OBJS) main_adj.o $(MPI_TGT_OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LCDF) $(LMPI) -lampiCommon  -lampiTape -lampiADtoolStubsOO -lampiADtoolStubsST -lampiBookkeeping -lblas -lampiPlainC

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


#
# Auxiliary utility programs and List of Dependecies:
# ========= ======= ======== === ==== == ============
#
  TOOLS = mpc cross_matrix cppcheck srcscheck checkkwds partit ncjoin ncrename

tools: $(TOOLS)

mpc: mpc.F
	$(CPP) -P $(CPPFLAGS) mpc.F > mpc_.f
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpc mpc_.f
cross_matrix: cross_matrix.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cross_matrix cross_matrix.o
cppcheck: cppcheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cppcheck cppcheck.o
srcscheck: srcscheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck.o
	rm -f check_srcs.F
checkkwds: checkkwds.o
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
	./cross_matrix module*.F90 *.F90 *.F

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

$(TAP_TARGET)_b.f: $(ADJ_PSRCS)
	${TAPENADE} $^ -noisize -noisize77 -tracelevel 10 -msglevel 20 -msginfile -head "cost_fun(ad_x)\(cost)" -r8 -reverse -output $(TAP_TARGET) -I /usr/include -I /usr/local/include

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

$(TAP_TARGET)_d.f: $(TGT_PSRCS) #main_tgt.f
	${TAPENADE} $^ -noisize -noisize77 -tracelevel 10 -msglevel 20 -msginfile -head "cost_fun(cost)/(ad_x)" -r8 -output $(TAP_TARGET) -I /usr/include -I /usr/local/include

$(TAP_TARGET)_context_d.f: $(TGT_PSRCS) main_tgt.f
	${TAPENADE} $^ -noisize -noisize77 -tracelevel 10 -msglevel 20 -msginfile -head "cost_fun(cost)/(ad_x)" -r8 -context -output $(TAP_TARGET) -I /usr/include -I /usr/local/include


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

include Make.depend


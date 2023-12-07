#======================================================================
# CROCO is a branch of ROMS developped at IRD and INRIA, in France
# The two other branches from UCLA (Shchepetkin et al) 
# and Rutgers University (Arango et al) are under MIT/X style license.
# CROCO specific routines (nesting) are under CeCILL-C license.
# 
# CROCO website : http://www.croco-ocean.org
#======================================================================

#======================================================================
#
# Automatic Differentiation Makefile
#
#======================================================================

%.tap.f: %.F
	$(CPP) -P $(CPPFLAGS) -D__TAPENADE__ $*.F | ./mpc > $*.tap.f

#======================================================================
#
# Tapenade:
# =========
#

ADJ_SRCS=cost_fun.F init_control.F set_state.F save_restore.F t3dmix.F uv3dmix.F get_vbc.F step.F step2d.F pre_step3d.F set_depth.F grid_stiffness.F rho_eos.F ana_vmix.F omega.F prsgrd.F rhs3d.F step3d_uv1.F step3d_uv2.F step3d_t.F t3dbc.F u3dbc.F v3dbc.F v2dbc.F u2dbc.F exchange.F analytical.F MessPass2D.F MessPass3D.F zetabc.F set_avg.F debug.F dummy.F get_ij.F distance.F xtime.F get_date.F
ADJ_PSRCS=$(ADJ_SRCS:.F=.tap.f)
TAP_TARGET=autodiff
ADJ_OBJS=$(TAP_TARGET)_b.o m1qn3.o treeverse.o adBinomial.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o adj_driver.o cost_fun.o init_control.o set_state.o save_restore.o get_ij.o distance.o xtime.o code_insertion.o

TGT_SRCS= $(ADJ_SRCS) ana_initial.F ana_grid.F wrt_his.F def_his.F setup_grid1.F setup_grid2.F insert_node.F checkdims.F put_global_atts.F def_grid_3d.F wrt_grid.F lenstr.F nf_fread.F fillvalue.F wvlcty.F nf_add_attribute.F init_scalars.F timers_roms.F init_arrays.F nf_fread_y.F nf_fread_x.F set_scoord.F set_weights.F get_initial.F closecdf.F cost_driver.F analytical.F
TGT_PSRCS=$(TGT_SRCS:.F=.tap.f)
TGT_PSRCS_DBG= main_cost.tap.f $(TGT_PSRCS)
TGT_OBJS=$(TAP_TARGET)_d.o m1qn3.o treeverse.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o tgt_driver.o cost_fun.o init_control.o set_state.o save_restore.o get_ij.o distance.o xtime.o code_insertion.o
ADJ_TGT_OBJS=$(TAP_TARGET)_b.o $(TAP_TARGET)_d.o m1qn3.o treeverse.o adBinomial.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o adj_tgt_driver.o cost_fun.o init_control.o set_state.o save_restore.o get_ij.o distance.o xtime.o code_insertion.o 
TGT_OBJS_DBG=$(TAP_TARGET)_context_d.o m1qn3.o adContext.o treeverse.o adBufferC.o adBuffer.o  adStack.o code_insertion.o  read_obs.o analytical.o t3dbc.o v3dbc.o exchange.o zetabc.o u2dbc.o v2dbc.o pre_step3d.o u3dbc.o get_vbc.o diag.o ana_vmix.o output.o  check_kwds.o  read_inp.o init_scalars.o timers_roms.o init_arrays.o ana_grid.o setup_grid1.o setup_grid2.o set_scoord.o set_weights.o set_depth.o grid_stiffness.o get_initial.o ana_initial.o set_depth.o rho_eos.o omega.o wrt_his.o step.o closecdf.o wrt_rst.o lenstr.o

ADJ_OBJS_DBG1=$(TAP_TARGET)_context1_d.o m1qn3.o adDebug.o treeverse.o adBufferC.o adBuffer.o adStack.o code_insertion.o  read_obs.o analytical.o t3dbc.o v3dbc.o exchange.o zetabc.o u2dbc.o v2dbc.o pre_step3d.o u3dbc.o get_vbc.o diag.o ana_vmix.o output.o  check_kwds.o  read_inp.o init_scalars.o timers_roms.o init_arrays.o ana_grid.o setup_grid1.o setup_grid2.o set_scoord.o set_weights.o set_depth.o grid_stiffness.o get_initial.o ana_initial.o set_depth.o rho_eos.o omega.o wrt_his.o step.o closecdf.o wrt_rst.o lenstr.o  cost_fun.o init_control.o set_state.o save_restore.o

ADJ_OBJS_DBG2=$(TAP_TARGET)_context2_b.o m1qn3.o adDebug.o treeverse.o adBufferC.o adBuffer.o  adStack.o adBinomial.o code_insertion.o  read_obs.o analytical.o t3dbc.o v3dbc.o exchange.o zetabc.o u2dbc.o v2dbc.o pre_step3d.o u3dbc.o get_vbc.o diag.o ana_vmix.o output.o  check_kwds.o  read_inp.o init_scalars.o timers_roms.o init_arrays.o ana_grid.o setup_grid1.o setup_grid2.o set_scoord.o set_weights.o set_depth.o grid_stiffness.o get_initial.o ana_initial.o set_depth.o rho_eos.o omega.o wrt_his.o step.o closecdf.o wrt_rst.o lenstr.o  cost_fun.o init_control.o set_state.o save_restore.o



DIV_OBJS=m1qn3.o treeverse.o adBufferC.o adBuffer.o adDebug.o adStack.o read_obs.o optim_driver.o div_driver.o cost_fun.o init_control.o set_state.o save_restore.o get_ij.o distance.o xtime.o

TGT_CONTEXT_OBJS=$(TAP_TARGET)_d.o cost_fun.o init_control.o set_state.o save_restore.o contextAD.o

LAMPI_PLAINC= -lampiPlainC
LAMPI= -lampiCommon  -lampiTape  -lampiBookkeeping -lblas -lampiPlainC

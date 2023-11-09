###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Add the source files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})

###########################################################
# generate config.h
configure_file(config.h.in config.h)
configure_file(config_post.h.in config_post.h)

###########################################################
# -I directives
include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${NETCDFF_INCLUDE_DIR})
set(CROCO_FORTRAN_CPP_FLAGS ${CROCO_FORTRAN_CPP_FLAGS}
                            -DHAVE_CMAKE_CONFIG
                            -I${NETCDFF_INCLUDE_DIR}
                            -I${CMAKE_CURRENT_SOURCE_DIR}
                            -I${CMAKE_CURRENT_BINARY_DIR}
                            -I${CMAKE_BINARY_DIR}
)
if (MPI_FOUND)
	set(CROCO_FORTRAN_CPP_FLAGS ${CROCO_FORTRAN_CPP_FLAGS} -I${MPI_INCLUDE_DIR})
endif()

###########################################################
# Apply croco flags
add_compile_options(${CROCO_FORTRAN_FLAGS})
add_link_options   (${CROCO_FORTRAN_FLAGS})

###########################################################
# path to "mpc", the ocean internal source preprocessor
set(CROCO_MPC        $<TARGET_FILE:mpc>)
set(CROCO_AGRIF_CONV $<TARGET_FILE:conv>)

###########################################################
file(GLOB CROCO_CPP_H
    *.h
    *.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
    ${CMAKE_CURRENT_BINARY_DIR}/config_post.h
    ${CPPDEF_OVERRIDE}
    ${CPPDEF_DEV_OVERRIDE}
    ${PARAM_H_OVERRIDE}
)

###########################################################
# List files
#======================================================================
# SRCS: source codes files are sorted into groups, separated by
# blanc lines:
#   1) main driving part;         2) 2D time stepping engine;
#   3) 3D time stepping engine;   4) non-hydrostatic engine
#   5) sea-water EOS;             6) vertical mixing schemes;
#   7) on-fly model diagnostics;  8) netCDF I/O routines;
#   9) main forcing;              10) Wave forcing routines;
#  11) Online surface forcing;    12) Floats routines;
#  13) Station diagnostics;       14) biology/sediment routines;
#  15) PISCES biology;            16) MPI routines;
#  17) AGRIF routines;            18) OASIS coupling interface;
#
# SRC90: additional F90 routines for PISCES code and Non-hydrostatic
#        analysis routines
#======================================================================
# 1) Main driving part
list(APPEND CROCO_SRCS_F77 main.F step.F read_inp.F timers_roms.F init_scalars.F
                           init_arrays.F set_weights.F set_scoord.F ana_grid.F
                           setup_grid1.F setup_grid2.F set_nudgcof.F
                           ana_initial.F analytical.F zonavg.F)
# 2) 2D time stepping engine
list(APPEND CROCO_SRCS_F77 ts_2d/step2d.F ts_2d/u2dbc.F ts_2d/v2dbc.F ts_2d/zetabc.F
                           ts_2d/obc_volcons.F)
# 3) 3D time stepping engine
list(APPEND CROCO_SRCS_F77 ts_3d/pre_step3d.F ts_3d/step3d_t.F ts_3d/step3d_uv1.F
                           ts_3d/step3d_uv2.F ts_3d/prsgrd.F ts_3d/rhs3d.F
                           ts_3d/set_depth.F ts_3d/omega.F ts_3d/uv3dmix.F
                           ts_3d/uv3dmix_spg.F ts_3d/t3dmix.F ts_3d/t3dmix_spg.F
                           ts_3d/hmix_coef.F ts_3d/wetdry.F ts_3d/wavedry.F
                           ts_3d/u3dbc.F ts_3d/v3dbc.F ts_3d/t3dbc.F)
# 4) non-hydrostatic engine
list(APPEND CROCO_SRCS_F77 step3d_fast.F step3d_w.F rhs3d_w_nh.F initial_nbq.F
                           grid_nbq.F unbq_bc.F vnbq_bc.F wnbq_bc.F rnbq_bc.F
                           w3dbc.F nbq_bry_store.F)
# 5) sea-water EOS
list(APPEND CROCO_SRCS_F77 seawater_eos/rho_eos.F
                           seawater_eos/ab_ratio.F
                           seawater_eos/alfabeta.F
                           seawater_eos/alfabeta_k.F)
# 6) vertical mixing schemes
list(APPEND CROCO_SRCS_F77 wave_forcing/ana_vmix.F wave_forcing/bvf_mix.F
                           wave_forcing/lmd_vmix.F wave_forcing/gls_mixing.F
                           wave_forcing/lmd_skpp.F wave_forcing/lmd_bkpp.F
                           wave_forcing/lmd_swfrac.F wave_forcing/lmd_wscale.F)
# 7) on-fly model diagnostics
list(APPEND CROCO_SRCS_F77 diag.F wvlcty.F checkdims.F grid_stiffness.F
                           bio_diag.F check_kwds.F
                           check_switches2.F debug.F)
# 8) netCDF I/O routines
list(APPEND CROCO_SRCS_F77 param.F ncscrum.F scalars.F output.F
                           put_global_atts.F nf_fread.F nf_fread_x.F
                           nf_fread_y.F nf_read_bry.F get_date.F
                           closecdf.F fillvalue.F
                           nf_add_attribute.F set_cycle.F def_grid_2d.F
                           def_grid_3d.F def_his.F def_rst.F def_diags.F
                           def_diagsM.F def_bio_diags.F wrt_grid.F wrt_his.F
                           wrt_avg.F wrt_rst.F wrt_diags.F wrt_diags_avg.F
                           wrt_diagsM.F wrt_diagsM_avg.F wrt_bio_diags.F
                           wrt_bio_diags_avg.F set_avg.F set_diags_avg.F
                           set_diagsM_avg.F set_bio_diags_avg.F def_diags_vrt.F
                           wrt_diags_vrt.F set_diags_vrt.F set_diags_vrt_avg.F
                           wrt_diags_vrt_avg.F def_diags_ek.F wrt_diags_ek.F
                           set_diags_ek.F set_diags_ek_avg.F wrt_diags_ek_avg.F
                           def_diags_pv.F wrt_diags_pv.F set_diags_pv.F
                           set_diags_pv_avg.F wrt_diags_pv_avg.F
                           def_diags_eddy.F set_diags_eddy_avg.F
                           wrt_diags_eddy_avg.F def_surf.F wrt_surf.F
                           set_surf_avg.F wrt_surf_avg.F)
# 9) main forcing
list(APPEND CROCO_SRCS_F77 get_grid.F get_initial.F get_vbc.F get_wwave.F
                           get_tclima.F get_uclima.F get_ssh.F get_sss.F
                           get_smflux.F get_stflux.F get_srflux.F get_sst.F
                           mod_tides_mas.F tidedata.F mas.F get_tides.F
                           clm_tides.F get_bulk.F bulk_flux.F get_bry.F
                           get_bry_bio.F sstskin.F get_psource.F
                           get_psource_ts.F get_btflux.F)
# 10) Wave forcing routines
list(APPEND CROCO_SRCS_F77 wave_forcing/mrl_wci.F wave_forcing/wkb_wwave.F wave_forcing/wkbbc.F wave_forcing/get_bry_wkb.F)
# 11) Online surface forcing
list(APPEND CROCO_SRCS_F77 surface_forcing/online_bulk_var.F
                           surface_forcing/online_get_bulk.F
                           surface_forcing/online_interp.F
                           surface_forcing/online_interpolate_bulk.F
                           surface_forcing/online_set_bulk.F)
# 12) Floats routines
list(APPEND CROCO_SRCS_F77 init_floats.F wrt_floats.F step_floats.F rhs_floats.F
                           interp_rho.F def_floats.F init_arrays_floats.F
                           random_walk.F get_initial_floats.F)
# 13) Station diagnostics
list(APPEND CROCO_SRCS_F77 station_diagnostics/init_sta.F
                           station_diagnostics/wrt_sta.F
                           station_diagnostics/step_sta.F
                           station_diagnostics/interp_sta.F
                           station_diagnostics/def_sta.F
                           station_diagnostics/init_arrays_sta.F)
# 14) biology/sediment routines
list(APPEND CROCO_SRCS_F77 biology.F o2sato.F sediment.F bbl.F)
# 16) MPI routines
list(APPEND CROCO_SRCS_F77 MPI_Setup.F MessPass2D.F MessPass3D.F exchange.F
                           autotiling.F)
# 16) AGRIF routines
list(APPEND CROCO_SRCS_F77 zoom.F update2D.F set_nudgcof_fine.F zoombc_2D.F
                           zoombc_3D.F uv3dpremix.F t3dpremix.F update3D.F
                           zoombc_3Dfast.F Agrif2Model.F)
# 17) XIOS routines
#list(APPEND CROCO_SRCS_F77 send_xios_diags.F)
# 18) OASIS coupling interface
list(APPEND CROCO_SRCS_F77 cpl_prism_define.F cpl_prism_put.F cpl_prism_init.F
                           cpl_prism_get.F cpl_prism_getvar.F cpl_prism_grid.F)
#list(APPEND CROCO_SRCS_F77 module_substance.F module_MUSTANG.F)
list(APPEND CROCO_SRCS_F77 online_spectral_diags.F var3d_oa_out.F var2d_oa_out.F
                           scal0d_oa_out_loc.F scal0d_oa_out_full.F output_oa.F
                           copy_to_devices.F exchange_device_host.F init_acc.F)
# SRC90: additional F90 routines for PISCES code and Non-hydrostatic
#        analysis routines
list(APPEND CROCO_SRCS_F90 online_analysis/module_oa_space.F90
                           online_analysis/module_oa_periode.F90
                           online_analysis/module_oa_variables.F90
                           online_analysis/module_oa_type.F90
                           online_analysis/module_tile_oa.F90
                           online_analysis/module_oa_level.F90
                           online_analysis/module_oa_interface.F90
                           online_analysis/var_oa.F90
                           )
# libs
list(APPEND CROCO_SRC_F90 lib/tooldatosec.F90
                          lib/toolsectodat.F90
                          lib/tooldecompdat.F90
                          lib/tooldatetosec.F90
                          )
# Pices files
#list(APPEND CROCO_SRCS_F90 par_pisces.F90 ocean2pisces.F90 trc.F90
#                           sms_pisces.F90 p4zche.F90 p4zint.F90 p4zlys.F90
#                           p4zflx.F90 p4zlim.F90 p4zsink.F90 p4zmicro.F90
#                           p4zmeso.F90 p4zmort.F90 p4zopt.F90  p4zprod.F90
#                           p4zrem.F90 p4zsbc.F90 p4zsed.F90 trcwri_pisces.F90
#                           p4zagg.F90 p4zfechem.F90 p4zligand.F90 p4zpoc.F90
#                           p4zsms.F90 p5zlim.F90 p5zmeso.F90 p5zmicro.F90
#                           p5zmort.F90 p5zprod.F90 p4zbio.F90 trcsms_pisces.F90

#                           sedarr.F90 sedbtb.F90 sedchem.F90 sedco3.F90
#                           seddiff.F90 seddsr.F90 seddta.F90 sed.F90
#                           sedini.F90 sedinitrc.F90 sedinorg.F90 sedmat.F90
#                           sedmbc.F90 sedmodel.F90 sed_oce.F90 sedorg.F90
#                           sedrst.F90 sedsfc.F90 sedstp.F90 sedwri.F90
#                           setavg_sed.F90 online_analysis/module_parameter_oa.F90
#
#                           flocmod.F90 comsubstance.F90 substance.F90
#
#                           comMUSTANG.F90 coupler_MUSTANG.F90
#                           sed_MUSTANG_CROCO.F90 sed_MUSTANG.F90 initMUSTANG.F90
#                           plug_MUSTANG_CROCO.F90

###########################################################
# files for copy to device
set(COPY_TO_DEVICE_VAR_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean2d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/scalars.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean3d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/grid.h
    ${CMAKE_CURRENT_SOURCE_DIR}/coupling.h
    ${CMAKE_CURRENT_SOURCE_DIR}/private_scratch.h
    ${CMAKE_CURRENT_SOURCE_DIR}/mixing.h
    ${CMAKE_CURRENT_SOURCE_DIR}/forces.h
    ${CMAKE_CURRENT_SOURCE_DIR}/work.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ncscrum.h
    ${CMAKE_CURRENT_SOURCE_DIR}/averages.h
    ${CMAKE_CURRENT_SOURCE_DIR}/lmd_kpp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/climat.h
    ${CMAKE_CURRENT_SOURCE_DIR}/nbq.h
    ${CMAKE_CURRENT_SOURCE_DIR}/sources.h
    ${CMAKE_CURRENT_SOURCE_DIR}/wkb_wwave.h
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary.h
)

###########################################################
# Regenerate the copy_to_device.F file to auto-inject the vars to upload
# to the GPU
if (OPENACC)
    # build command to gen copy_to_device
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/copy_to_devices_generated.h
        COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/common2device.py
                ${COPY_TO_DEVICE_VAR_FILES}
                -o ${CMAKE_CURRENT_BINARY_DIR}/copy_to_devices_generated.h
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/common2device.py
        DEPENDS ${COPY_TO_DEVICE_VAR_FILES}
        VERBATIM
    )
    
    # add it in deps loop
    list(APPEND CROCO_CPP_H ${CMAKE_CURRENT_BINARY_DIR}/copy_to_devices_generated.h)
endif()

###########################################################
# Make absolute paths for source prepare pass
croco_make_absolute_paths(CROCO_SRCS_F77 ${CMAKE_CURRENT_SOURCE_DIR})
croco_make_absolute_paths(CROCO_SRCS_F90 ${CMAKE_CURRENT_SOURCE_DIR})

###########################################################
# generate some sources files byb using some binary programs

# setup_kwds.F
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/setup_kwds.F
	COMMAND checkkwds && mv setup_kwds.F ${CMAKE_CURRENT_BINARY_DIR}/setup_kwds.F
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS checkkwds
	VERBATIM
)

# check_srcs.F
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/check_srcs.F
	COMMAND srcscheck
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	DEPENDS srcscheck
	VERBATIM
)

# check_switches1.F
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/check_switches1.F
	COMMAND cat cppdefs.h cppdefs_dev.h > mergcpp.txt && ${CMAKE_CURRENT_BINARY_DIR}/cppcheck && mv check_switches1.F ${CMAKE_CURRENT_BINARY_DIR}/check_switches1.F
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS cppcheck cppdefs.h cppdefs_dev.h
	VERBATIM
)

# Append the generated file to ocean croco
list(APPEND CROCO_SRCS_F77 ${CMAKE_CURRENT_BINARY_DIR}/setup_kwds.F)
list(APPEND CROCO_SRCS_F77 ${CMAKE_CURRENT_BINARY_DIR}/check_srcs.F)
list(APPEND CROCO_SRCS_F77 ${CMAKE_CURRENT_BINARY_DIR}/check_switches1.F)

###########################################################
# ISSUE: remark, this is to solve an issue due to the trick
# of pre-applying MPC for "ocean" which is we let lenstr.F
# in the list, make add_executable(partit) skipping it.
add_library(croco_common STATIC lenstr.F insert_node.F)

###########################################################
# build ncjoin utility
set(NCJOIN_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncjoin.F)
croco_cpp_and_mpc_preprocess(NCJOIN_SRC_F)
add_executable(ncjoin ${NCJOIN_SRC_F})
target_link_libraries(ncjoin ${NETCDFF_LIBRARY})
target_link_libraries(ncjoin croco_common)
set_property(TARGET ncjoin PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build ncrename utility
set(NCRENAME_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncrename.F)
croco_cpp_and_mpc_preprocess(CRENAME_SRC_F)
add_executable(ncrename ${NCRENAME_SRC_F})
target_link_libraries(ncrename ${NETCDFF_LIBRARY})
target_link_libraries(ncrename croco_common)
set_property(TARGET ncrename PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build mpc utility
add_executable(mpc mpc.F)
set_property(TARGET mpc PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build checkkwds utility (used to produce setup_kwds.F)
add_executable(checkkwds checkkwds.F)
set_property(TARGET checkkwds PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build checkkwds utility (used to produce setup_kwds.F)
add_executable(cppcheck cppcheck.F)
set_property(TARGET cppcheck PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build checkkwds utility (used to produce check_srcs.F)
add_executable(srcscheck srcscheck.F)
set_property(TARGET srcscheck PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build cross_matrix utility
add_executable(cross_matrix cross_matrix.F)
set_property(TARGET cross_matrix PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build partit utility
add_executable(partit partit.F)
target_link_libraries(partit ${NETCDFF_LIBRARY})
target_link_libraries(partit croco_common)
set_property(TARGET partit PROPERTY LINKER_LANGUAGE Fortran)

###########################################################
# build croco ocean simulation

# TODO: remove when finished
if (WITH_PARALLEL STREQUAL "openacc-psyclone")
	croco_psyclone_pre_filter_acc(CROCO_SRCS_F77)
endif ()

# prepare the sources with filters
if (ENABLE_AGRIF)
	croco_agrif_preprocess(CROCO_SRCS_F77)
	croco_agrif_preprocess(CROCO_SRCS_F90)
else()
	croco_cpp_and_mpc_preprocess(CROCO_SRCS_F77)
endif()
if (OPENACC)
	croco_change_loop_preprocess(CROCO_SRCS_F77)
endif ()

# build executable
set(CMAKE_Fortran_PREPROCESS ON)
add_executable(croco ${CROCO_SRCS_F90} ${CROCO_SRCS_F77})
target_link_libraries(croco ${NETCDFF_LIBRARY})
target_link_libraries(croco croco_common)
set_property(TARGET croco PROPERTY LINKER_LANGUAGE Fortran)

# Inject psyclone in the loop if enabled
if (WITH_PARALLEL STREQUAL "openacc-psyclone")
	set_target_properties(croco PROPERTIES Fortran_COMPILER_LAUNCHER ${CMAKE_SOURCE_DIR}/PSYCLONE/psyclone.compiler.wrapper.py)
endif ()

# enable AGRIF
if (ENABLE_AGRIF)
	target_link_libraries(croco agrif)
endif()

###########################################################
# Make some of the executables landing directly in root dir
set_property(TARGET partit croco cross_matrix PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

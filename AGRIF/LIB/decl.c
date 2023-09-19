/******************************************************************************/
/*                                                                            */
/*     CONV (converter) for Agrif (Adaptive Grid Refinement In Fortran)       */
/*                                                                            */
/* Copyright or   or Copr. Laurent Debreu (Laurent.Debreu@imag.fr)            */
/*                        Cyril Mazauric (Cyril_Mazauric@yahoo.fr)            */
/* This software is governed by the CeCILL-C license under French law and     */
/* abiding by the rules of distribution of free software.  You can  use,      */
/* modify and/ or redistribute the software under the terms of the CeCILL-C   */
/* license as circulated by CEA, CNRS and INRIA at the following URL          */
/* "http://www.cecill.info".                                                  */
/*                                                                            */
/* As a counterpart to the access to the source code and  rights to copy,     */
/* modify and redistribute granted by the license, users are provided only    */
/* with a limited warranty  and the software's author,  the holder of the     */
/* economic rights,  and the successive licensors  have only  limited         */
/* liability.                                                                 */
/*                                                                            */
/* In this respect, the user's attention is drawn to the risks associated     */
/* with loading,  using,  modifying and/or developing or reproducing the      */
/* software by the user in light of its specific status of free software,     */
/* that may mean  that it is complicated to manipulate,  and  that  also      */
/* therefore means  that it is reserved for developers  and  experienced      */
/* professionals having in-depth computer knowledge. Users are therefore      */
/* encouraged to load and test the software's suitability as regards their    */
/* requirements in conditions enabling the security of their systems and/or   */
/* data to be ensured and,  more generally, to use and operate it in the      */
/* same conditions as regards security.                                       */
/*                                                                            */
/* The fact that you are presently reading this means that you have had       */
/* knowledge of the CeCILL-C license and that you accept its terms.           */
/******************************************************************************/
/* version 1.7                                                                */
/******************************************************************************/

#include "decl.h"


 variable *curvar = NULL;

 listvar *List_ModuleUsedInModuleUsed_Var = NULL;
 listvar *List_ModuleUsed_Var = NULL;
 listvar *listduplicated = NULL;

 listvar *List_GlobalParameter_Var = NULL;
 listvar *List_Global_Var = NULL;
 listvar *List_Data_Var = NULL;
 listvar *List_Data_Var_Cur = NULL;
 listvar *List_Save_Var = NULL;
 listvar *List_SubroutineArgument_Var = NULL;
 listvar *List_SubroutineDeclaration_Var = NULL;
 listvar *List_UsedInSubroutine_Var = NULL;
 listvar *List_Parameter_Var = NULL;
 listvar *List_Dimension_Var = NULL;
 listvar *List_FunctionType_Var = NULL;
 listvar *List_NotGridDepend_Var = NULL;
 listvar *List_Common_Var = NULL;


 listname *List_Pointer_Var = NULL;
 listname *List_ImplicitNoneSubroutine = NULL;
 
 listname *List_Do_labels = NULL; 
 /* A list that contains the do labels if any */

 listusemodule *List_NameOfModuleUsed = NULL;
 listusemodule *List_Include = NULL;
 listusemodule *listofmoduletmp = NULL;
 listusemodule *tmpuselocallist = NULL;

 listparameter *List_GlobParamModuleUsedInModuleUsed_Var = NULL;
 listparameter *List_GlobParamModuleUsed_Var = NULL;

 listnom *List_ContainsSubroutine = NULL;
 listnom *List_Subroutine_For_Alloc = NULL;
 listnom *listofmodules = NULL;
 listnom *listofkind = NULL;
 listnom *List_NameOfModule = NULL;
 listnom *List_NameOfCommon = NULL;
 listnom *List_SubroutineWhereAgrifUsed = NULL;

 listallocate *List_Allocate_Var = NULL;

 listvarpointtovar *List_CouplePointed_Var = NULL;
                           /*  variables which are pointed to an other one    */

 listindice *Listofavailableindices = NULL;
                           /* List of available indices in the tabvars table  */
 listindice **Listofavailableindices_glob = NULL;

 listdim *curdim = NULL;
 listdim *commondim = NULL;

/******************************************************************************/
/****************   *** COMMON Variables ***  *********************************/
/******************************************************************************/

 int positioninblock = 0;
 char commonvar[LONG_VNAME];
 char commonblockname[LONG_VNAME];

/******************************************************************************/
/****************   *** AGRIF Variables ***   *********************************/
/******************************************************************************/
 int inagrifcallargument = 0;
 int afterpercent = 0;
 int sameagrifargument = 0;
 int InAgrifParentDef = 0;
 char sameagrifname[LONG_VNAME];
/******************************************************************************/
/****************   *** VAR DEF Variables ***   *******************************/
/******************************************************************************/
 int indicemaxtabvars[NB_CAT_VARIABLES];     /* Number of variables in the model i.e. last      */
                           /*    indice used in  the tabvars table            */
 int PublicDeclare = 0;        /* Variable has been declared as PUBLIC */
 int PrivateDeclare = 0;       /* Variable has been declared as PRIVATE */
 int ExternalDeclare = 0;      /* Variable has been declared as EXTERNAL */
 int InitialValueGiven = 0;    /* An initial value has been given */
 int Allocatabledeclare = 0;
 int Targetdeclare = 0;
 int SaveDeclare = 0;
 int functiondeclarationisdone = 0;
 int pointerdeclare = 0;
 int optionaldeclare = 0;
 int inside_type_declare = 0;
 int VariableIsParameter = 0;
 int dimsgiven = 0;
 int shouldincludempif = 0;
 int c_star = 0;
 char DeclType[LONG_VNAME];
 char nameinttypename[LONG_VNAME];
 char nameinttypenameback[LONG_VNAME];
 int GlobalDeclaration = 0;
 int GlobalDeclarationType = 0;
 char InitValue[LONG_M];
 char IntentSpec[LONG_M];
 char NamePrecision[LONG_C];
 char CharacterSize[LONG_VNAME];
 char vallengspec[LONG_VNAME];
 int isrecursive = 0;
 int is_result_present = 0;

/******************************************************************************/
/****************   *** CONV Variables ***   **********************************/
/******************************************************************************/
 int dimprob = 0;             /* dimension of the problem : 1 for 1D,2 for 2D,   */
                           /*    3 for 3D                                     */
 int onlyfixedgrids = 0;       /* = 1 if onlyfixedgrids is true                   */
 int todebug = 0;
 int fixedgrids = 0;           /* = 1 if fixedgrids is true                       */
 char nbmaillesX[LONG_VNAME];	// number of cells in the x direction
 char nbmaillesY[LONG_VNAME];	// number of cells in the y direction
 char nbmaillesZ[LONG_VNAME];	// number of cells in the z direction
 int IndicenbmaillesX = 0;
 int IndicenbmaillesY = 0;
 int IndicenbmaillesZ = 0;

 int inmodulemeet = 0;
 int incalldeclare = 0;
 int aftercontainsdeclare = 0; /* Signale si l'on vient d'un contains ou non */
 int retour77 = 0;
 int callagrifinitgrids = 0;
 int callmpiinit = 0;
 int firstpass = 0;
 int pointedvar = 0;
 int NbMailleXDefined = 0;
 int agrif_parentcall = 0;
 int didvariableadded = 0;
 int SubloopScalar = 0;        /* = 1 we should put in argument of sub_loop       */
                           /*    only                                         */
                           /*    scalar and not table u(1,1,1) in place of u  */
 int inprogramdeclare = 0;
 int insubroutinedeclare = 0;
 int inmoduledeclare = 0;
 int dimsempty = 0;
 int created_dimensionlist = 0;
 int incontainssubroutine = 0;

 char meetagrifinitgrids[LONG_M];
 char mpiinitvar[LONG_M];
 char toprintglob[LONG_M];
 char tmpvargridname[LONG_M];
 char dependfilename[LONG_FNAME];
 char charusemodule[LONG_VNAME];
 char subofagrifinitgrids[LONG_M];
 char curmodulename[LONG_VNAME];
 char subroutinename[LONG_VNAME];
 char old_subroutinename[LONG_VNAME]; // For internal subprogramm
 char cur_filename[LONG_FNAME];		// Name of the current parsed Fortran file
 char config_file[LONG_FNAME];		// Name of conv configuration file (ex: amr.in)
 char work_dir[LONG_FNAME];			// Work directory         (default: './')
 char include_dir[LONG_FNAME];		// Include directory      (default: './AGRIF_INC')
 char output_dir[LONG_FNAME];		// output directory       (default: './AGRIF_MODELFILES')
 char input_dir[LONG_FNAME];		// source input directory (default: './')

 FILE *fortran_out = NULL;          /* Output File                                    */
 FILE *fortran_in = NULL;           /* Input File                                     */
 FILE *oldfortran_out = NULL;
 FILE *old_oldfortran_out = NULL; // For internal subprogramm
 FILE *subloop = NULL;
 FILE *module_declar = NULL;
 FILE *allocationagrif = NULL;

 long int pos_cur = 0;         /* current position in the output file             */
 long int pos_curagrifparent = 0;
                           /* current position in the output file             */
 long int pos_curcall = 0;     /* current position in the output file             */
 long int pos_curuse = 0;      /* current position in the output file             */
 long int pos_curuseold = 0;   /* current position in the output file             */
 long int pos_curfunction = 0; /* current position in the output file             */
 long int pos_cur_decl = 0;    /* current position in the output file             */
 long int pos_curdata = 0;     /* current position in the output file             */
 long int pos_curparameter = 0;/* current position in the output file             */
 long int pos_curcommon = 0;   /* current position in the output file             */
 long int pos_cursave = 0;     /* current position in the output file             */
 long int pos_curdimension = 0;/* current position in the output file             */
 long int pos_curinclude = 0;  /* final position of a line in file                */
 long int pos_end = 0;         /* final position of a line in file                */
 long int pos_endsubroutine = 0;
                           /* final position of a line in file                */

size_t length_last = 0;
size_t length_first = 0;
size_t length_v_vallengspec = 0;
size_t length_v_commoninfile = 0;
size_t length_v_precision = 0;
size_t length_v_IntentSpec = 0;
size_t length_v_initialvalue = 0;
size_t length_v_readedlistdimension = 0;
size_t length_a_nomvar = 0;
size_t length_toprintglob = 0;
size_t length_tmpvargridname = 0;
size_t length_ligne_Subloop = 0;
size_t length_toprint_utilagrif = 0;
size_t length_toprinttmp_utilchar = 0;
size_t length_ligne_writedecl = 0;
size_t length_newname_toamr = 0;
size_t length_newname_writedecl = 0;
size_t length_ligne_toamr = 0;
size_t length_tmpligne_writedecl = 0;
 int value_char_size = 0;
 int value_char_size1 = 0;
 int value_char_size2 = 0;
 int value_char_size3 = 0;


 int inallocate = 0;
 int infixed = 0;
 int infree = 0;

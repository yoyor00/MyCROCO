
#if defined OBC_WEST || defined AGRIF_OBC_WEST
# ifdef Z_FRC_BRY
!$acc& ,zetabry_west
# endif
# ifdef M2_FRC_BRY
!$acc& ,ubarbry_west, vbarbry_west 
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc& ,ubry_west, vbry_west 
#  endif
#  ifdef T_FRC_BRY
!$acc& ,tbry_west
#  endif
# endif
#endif

#if defined OBC_EAST || defined AGRIF_OBC_EAST
# ifdef Z_FRC_BRY
!$acc& ,zetabry_east       
# endif
# ifdef M2_FRC_BRY
!$acc& ,ubarbry_east, vbarbry_east       
# endif
# ifdef SOLVE3D 
#  ifdef M3_FRC_BRY
!$acc& ,ubry_east, vbry_east       
#  endif
#  ifdef T_FRC_BRY
!$acc& ,tbry_east      
#  endif
# endif
#endif

#if defined OBC_SOUTH || defined AGRIF_OBC_SOUTH
# ifdef Z_FRC_BRY 
!$acc& ,zetabry_south       
# endif
# ifdef M2_FRC_BRY
!$acc& ,ubarbry_south, vbarbry_south       
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc& ,ubry_south, vbry_south       
#  endif
#  ifdef T_FRC_BRY
!$acc& ,tbry_south       
#  endif
# endif
#endif

#if defined OBC_NORTH || defined AGRIF_OBC_NORTH
# ifdef Z_FRC_BRY
!$acc& ,zetabry_north       
# endif
# ifdef M2_FRC_BRY
!$acc& ,ubarbry_north, vbarbry_north        
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc& ,ubry_north, vbry_north        
#  endif
#  ifdef T_FRC_BRY
!$acc& ,tbry_north        
#  endif
# endif
#endif

#if defined M3FAST && defined NBQ_FRC_BRY
# if defined OBC_WEST || defined AGRIF_OBC_WEST
!$acc& ,unbqbry_west, vnbqbry_west      
#  ifdef NBQ
!$acc& ,wnbqbry_west, rnbqbry_west      
#  endif
# endif
# if defined OBC_EAST || defined AGRIF_OBC_EAST
!$acc& ,unbqbry_east, vnbqbry_east      
#  ifdef NBQ
!$acc& ,wnbqbry_east, rnbqbry_east 
#  endif
# endif
# if defined OBC_SOUTH || defined AGRIF_OBC_SOUTH
!$acc& ,unbqbry_south, vnbqbry_south      
#  ifdef NBQ
!$acc& ,wnbqbry_south, rnbqbry_south      
#  endif
# endif
# if defined OBC_NORTH || defined AGRIF_OBC_NORTH
!$acc& ,unbqbry_north, vnbqbry_north      
#  ifdef NBQ
!$acc& ,wnbqbry_north, rnbqbry_north      
#  endif
# endif
#endif /* M3FAST */


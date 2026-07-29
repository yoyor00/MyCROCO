# ifdef BLOOM

#  define key_oxygen
#  undef key_BLOOM_opt2
#  ifndef key_BLOOM_opt2
#   if defined MUSTANG
#    define key_BLOOM_insed
#   else
#    undef key_BLOOM_insed
#   endif
                      /*   Biology diagnostics    */
                      /* ifdef BLOOM you must have defined DIAGNOSTICS_BIO */
#   define DIAGNOSTICS_BIO
#  endif

#endif



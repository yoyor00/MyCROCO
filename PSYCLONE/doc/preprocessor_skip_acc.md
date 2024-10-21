psyclone.preprocessor.skip.acc
==============================

I'm currently based on the manually made OpenACC branch of CROCO which means:

1. I need to enable the `OPENACC` key to get the files I didn't yet patched
(IOs, GPU-upload-download).
1. But for the files I'm transforming with PSyClone I need to get their CPU
version so I need to disable the `OPENACC` key for those files only.

I currently made this by duplicating `cppdef` and by disabling the `OPENACC`
key.

Todo:
-----

I made it in a dirty way I currently don't like. We might want to play more with
some just `-DDISABLE_OPENACC` option to the compiler. As it was originaly
temporary I currently kept the hacky verison I made.

Fill free to match as soon as possible.
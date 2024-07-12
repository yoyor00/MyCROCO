(use-modules (gnu packages mpi)
             (gnu packages autotools)
             (gnu packages base)
             (gnu packages gcc)
             (guix packages)
             (guix git-download)
             (guix build-system gnu)
             (guix licenses))

(define-public adjoinable-mpi
  (package
    (name "adjoinable-mpi")
    (version "0.1")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/bremond/adjoinablempi")
             (commit"30c3a638807eafbb73d4df17e327c2225cf6bf0f")))
       (sha256
        (base32
         "1snk8dd3k8shbihzg3yq1kdqjimd11biypzgx9sc351wxa9si2gf"))))
    (build-system gnu-build-system)
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)
       ("libtool" ,libtool)
       ("gnu-make" ,gnu-make)
       ("gcc" ,gcc)
       ("gfortran" ,gfortran)))
    (inputs
     `(("mpi" ,mpich)))
    (outputs '("out" "debug"))
    (arguments
     `(#:configure-flags
       '("--enable-fortranCompatible" "--with-gnu-ld" "--enable-debug")
       #:tests? #f))
    (synopsis "The Adjoinable MPI (AMPI) library for automatic
differentiation of MPI program.")
    (home-page "https://www.mcs.anl.gov/~utke/AdjoinableMPI/AdjoinableMPIDox/index.html")
    (description "The Adjoinable MPI (AMPI) library provides a
modified set of MPI subroutines that are constructed such that an
adjoint in the context of algorithmic differentiation (AD) can be
computed. The library is designed to be supported by a variety of AD
tools and to enable also the computation of (higher-order) forward
derivatives.")
    (license (x11-style "" "See file headers."))))

(packages->manifest (list adjoinable-mpi))

Programming style and guidelines
================================

.. _Programming_style :

Code organization
-----------------

The code is divided into three types of files. The comBIOLink contains the declaration of common variables so that BIOLink shares variables with CROCO, along with coupleur_define_BIOLink that defines the equivalence of variables between the biogeochemical models and CROCO. The coupleur_BIOLink contains the subroutines and functions of BIOLink. And finally, module_BIOLink is just a convenient structure that groups headers that are used by all BIOLink files and tooljulien is a function that is used by some subroutines of BIOLink.

The comBIOLink and coupleur_BIOLink are divided in three files: The standard one, the "_physics" one and the "_helping" one. The standard ones are used by BIOLink to do only its core mission of exchanging tracer concentration and names between biogeochemical models and CROCO. The physics ones are dedicated to the few subroutines of BIOLink that are used for the physics and the helping ones for the subroutines that helps the biogeochemical models with missing routines.

Version of Fortran
------------------

The code is written in fortran 90 but uses legacy elements from FORTRAN 77. The code still uses mostly subroutine instead of functions and somewhere a fixed format is still used. The fortran 90 part is mostly used for the free-formatting and the use of derived type in PEPTIC.

Code guidelines
---------------

The code is still under reorganization to follow a few guidelines.

* Commentary should be put under the code if it concerns one line and above if it is longer.
                                                                   
* if's should be indented and endif should contain which if they close in comment.

* Spacing should be used between blocks of instructions ( such as declarations or allocation of variables)
                                                                   
* Commentaries should be in english for global accessibility
                                                                   
* For allocation/declaration of variables, the variables are grouped by theme ( hydrodynamical, biological, other). The theme should be indicated by the following structure :

.. f:program:: Section comment

     !*************************************************************************!
     !*************************************************************************!
     !******************** Variables from the hydro model *********************!
     !*************************************************************************!
     !*************************************************************************!

     
     Subtheme can also be indicated by the following structure :
     
.. f:program:: Subsection comment    
     
      !============================================
      !  Variables related to the verification of the conservation routine
      !============================================
     
     


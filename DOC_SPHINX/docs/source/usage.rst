Usage
=====

.. _usage:

General description of a coupler
--------------------------------

The general objective of a coupler is to transfer information between one model and another. Here those informations are relatives to the biological/chemical variables that are considered as passive variables in the physical model and reactive in biological model. 

.. figure:: /General_interest_coupler.png


Usage
-----

To use it one needs to activate the cppkeys BIOLink in the cppdefs file, as well as the cppkeys linked to the model being used with CROCO ( BLOOM/PEPTIC/METEOR). Specific files linked with those models can be required.

Internal functioning
--------------------

BIOLink uses the tracer variables of CROCO and adds sources and sink terms on their evolution equation. At each baroclinic timestep, BIOLink is called and transfers the values of concentration of tracers to the biogeochemical models. It can transfer other variables such as temperature or salinity if they are required by the biogeochemical model.

The biogeochmical models computes the source and sink terms. Those terms are recovered by BIOLink that updates the concentration of tracers and sends it back to CROCO. It is presented in the scheme behind :

.. figure:: /Coupler_functioning.png 

The physical evolution of tracers is entirely computed by CROCO. BIOLink only handles the source and sink terms. Therefore, each cell is treated independantly in BIOLink.

Supplementary functions
-----------------------

To help interfacing CROCO with the biogeochemical models, some supplementary functions have been added. They include a simplified model for evaluating the PAR, as well as a computation for evaluating the conservation of tracers. Those functions and the related variables have been grouped in a file of the code.



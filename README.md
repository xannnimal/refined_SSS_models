# refined_SSS_models
Functions and examples for investigating and developing refined SSS models: spheroidal and multi-origin expansions, iterative implementations, and examples with different MEG helmet geometries.
Spheroidal Functions: From Tierney et al, https://www.biorxiv.org/content/10.1101/2023.09.11.557150v1
SSS functions: Taulu et al, single-origin SSS expansions. Details on functions can be found in the headers of files.



List of Functions:
"Combine_basis" combines two overlapping interior SSS basis expansions into one matrix
"gen_opm_geometry" generates positions and orientations OPM sensors and "coils" given a sensing direction
"gen_squid_geometry" generates magnetomenter positions from ".fif" file as input, typically uses MNE-Python/C examples of 306-sensor MEG helmets
"single_dipole_sim/noise" uses FieldTrip cfg to simulate dipole source inside spherical head conductor model. Can change sources, locations, noise ect using cfg options
"multiVSHin_singleVSHout" calculates SSS expansion using multi-origin interior and single-origin exterior
"spheroidIN_spheroidOUT/vshOUT" spheroidal interior basis with spheroidal/single-origin exterior
"xi" iterative SSS model

TODO
(1) investigate ways to create noise using cfg.courcemodel and "single_dipole_sim_noise". (2) impliment iterative model for multi-sphere expansion. (3) work on testing everything with real raw data from SQUID sensors

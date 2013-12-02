##ana2pi
* Overview
* Details of various Programs
* Bugs

### Overview
This package is implemented to extract "a set of defined" Observables by analyzing data from the interaction of an electron with a proton that results in the production of 2 pions and 1 proton(double-charged-pion Electroproduction off the proton). The "set of defined" Observables will be used by the JM Model to extract Electrocouplings of the Virtual Photon with the Proton. These Electrocouplings, in turn, will be used as a part of a "larger" research infrastructure to understand how the degrees of freedom inside the proton manifest themselves at various energy scales.

### Details of various Programs:

##### `proc_simstats.C` & `plot_simstats.py`[11-22-13]
These two programs are used to test the hypothesis for __Complete Simulation__ (defined in the text below)

`proc_simstats.C`:

Consolidates from various [top]_[q2w_range].root files, in a CSV format, store the following tabular data:
`Sim`,`Top`,`q2wbinnum`,`q2wbin`,`Varset`,`nFB_ST`,`nEB_SR`,`nFB_SR`,`nEB_SA`

where:

*	`nFB_ST` = Total number of Bins that are Filled (`FB`) by `genev` (`ST`) in the 2pi Reaction Phase Space (`PS`)
*	`nEB_SR` = Out of `nFB_ST`, the number of Bins that DO NOT HAVE Reconstructed events.
*	`nFB_SR` = Out of `nFB_ST`, the number of Bins that HAVE Reconstructed events.
*	`nEB_SA` = Out of `nFB_ER` (number of Bins that are Filled by Reconstructed events in Experiment), the number
         of bins that have no Acceptance.

`plot_simstats.py` is used to test the hypothesis stated below.

__Hypothesis for Complete Simulation__

Before stating the hypothesis, I want to note the following 2 important assumptions:

1. It is assumed that `genev`=Nature and therefore, `genev` throws events within the Natural Phase Space for `2pi` electroproduction
2. It is assumed that `GSIM`= CLAS Detector

Counting `nEB_SR/nFB_SR` is the primary means for keeping track of when the Simulation is Complete. Ideally,
for all `FB_ST` that lie within the Fiducial Volume of the CLAS detector, there should be a `FB_SR`. Any of 
the `nFB_ST` events that lie outside this Fiducial volume, should count for `nEB_SR`.Therefore:

`nFB_SR` = `nFB_ST`-`nEB-SR`

In addition, `nEB_SA` can serve as another (independent?) gauge. For if the Simulation is Complete,
`nEB_SA` should approach its asymptotic limit too. Ideally, if `ER` data has no Background and the Simulation is 
ideal, this asymptotic limit should approach 0. However, the limit will approach a number that may be representative
of the number of bins in `ER` data that is filled with Background events.

Hence, once Simulation is complete:

1. `nFB_ST`,`nFB_SR`,`nEB_SR`,`nEB_SA` should approach their asymptotic limits
2. `nFB_ST` = `nFB_SR` + `nEB_SR`

###Bugs
1. [12-02-13]__ProcCopyH10__: Definitely does not work when using TPROOF. The h10-Tree is Cloned the first time ProcCopyH10:handle() is called. However, each subsequent time a Slave is passed a new TChain, which happens often due to multi-threading, the TSelector::Init() is called, which resets the memory address into which the TTreee Branches are loaded and therefore, the Cloned h10-Tree does not point to the correct memory address anymore.   


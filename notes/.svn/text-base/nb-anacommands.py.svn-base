# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Chain of analysis:
# * __[exp] h10 $\rightarrow$ yields $\rightarrow$ processed yields__
# * __[sim] h10 $\rightarrow$ h10skim $\rightarrow$ yields $\rightarrow$ processed yields $\rightarrow$ process simstats $\rightarrow$ plot simstats __
#        

# <markdowncell>

# ##1. 
# ##[exp] h10 $\rightarrow$ yields
# ##[sim] h10skim $\rightarrow$ yields
# 
# This step uses the script __make.py__ that wraps __runSelector.C__
# 
# * __make.py__ 
#     1. makes h10 based skims (for example, h10skim) from h10
#     2. makes yields from h10 and various h10 based skims (for example, h10skim)
# * __make.py__ uses 
#     1. __make_tools.py__ 
#     2. __addy__ (to add yields from each __sim__ and put it in appropriate __cumsim__)
# * Run __make.py --h__ to get all the options
# 
# ###For Simulation data, make sure that 
# 
# 1. __ongoing/e1fs&lt;x&gt;.lst__ and __ongoing/e1fs&lt;x&gt;_cumulative.lst__ have been updated __AND__
# 2. __h10skim__ have been made for simulation data:
# 
#     * __make.py__ --h10skim=e1fs1
#     * __make.py__ --h10skim=e1fs2

# <markdowncell>

# ## Three Methods to obtain Yields
# __Method 1.__ was written keeping in mind the entire reprocessing of data, while __Method 2__ is more flexible and is to be used, for example, when an additional __sim__ finishes on the farm and it needs to be processed and added to __cumsim__ or if for any other reason, any individual __sim-__ or just the __exp-yield__ needs to be reprocessed.
# 
# __Method 3__ is the most flexible method to be used for debugging and "quick & dirty"/test jobs.
# 
# __NOTE__ __log__ from __make.py__ is located in __$anadir/makelog/__
# 
# 1. __makea &lt;bkdirname&gt;__ can be directly called which makes __sim-,mc-,exp-yields__ from __e1f1 & e1fd2__ by calling __make.py__ as follows (__log__ for __make__ in __$anadir/makelog/__):
# 
#    1. __for e1f1__
#        * make.py --yield=e1fs1   --skim=h10skim --bkdir=bkdirname 
#        * make.py --mcyield=e1fs1 --skim=h10skim --bkdir=bkdirname
#        * make.py --yield=e1fd1   --skim=h10     --bkdir=bkdirname
#    2. __for e1f2__
#        * make.py --yield=e1fs2   --skim=h10skim --bkdir=bkdirname
#        * make.py --mcyield=e1fs2 --skim=h10skim --bkdir=bkdirname
#        * make.py --yield=e1fd2   --skim=h10     --bkdir=bkdirname
#        
# 2. Any of the individual commands from __makea__ can directly be called. Note that __--bkdir__ option can be left blank here, and as a result, wherever yield.root/mcyield.root already exist, they will not be overwritten and the program will do nothing
# 
#     1. __for e1f1__
#         * make.py --yield=e1fs1   --skim=h10skim  
#         * make.py --mcyield=e1fs1 --skim=h10skim 
#         * make.py --yield=e1fd1   --skim=h10     
#     2. __for e1f2__
#         * make.py --yield=e1fs2   --skim=h10skim 
#         * make.py --mcyield=e1fs2 --skim=h10skim 
#         * make.py --yield=e1fd2   --skim=h10     
# 
# 3. For __debug__ jobs, use __runSelector.C__ in local directory

# <markdowncell>

# 
# 
# 
# 

# <markdowncell>

# ##2. yields $\rightarrow$ processed yields
# 
# * __procy__ is the main tcsh script that calls __procYields.C::procYields()__
# * All scripts specfic to __data type (e1fs1/e1fd1/e1fs2/e1fd2)__ and __usage of Unpolarized/Polarized Beam Helicity states__ are are soft linked to this main script
# 
# __NOTE__ that __procy__ has been clumsily updated to be sensitive to picking up if Yield processing needs to be dependent on the state of the Beam Helicity
# 
# ### All Tops (if Simulation, then for ~~each __cumsim__~~ latest cumsim), for example: 
# __NOTE__ if the following tcsh shell scripts are called with $1=1, then all existing ProcessedYields will be backed up under the current time stamp; otherwise only where ProcessedYields do not exists, will the following commands have any effect
# 
# * __e1fs1__: procy_e1fs1 
# * __e1fd1__: procy_e1fd1
# * __e1fs2__: procy_e1fs2
# * __e1fd2__: procy_e1fd2
# * __e1f_all__: procya
# 
# * Each Top individually, for example, __Top 5__:
#     * __e1f1, Top 1__: root -b -q 'procYields.C("yield.root", "1:2:3:4", 1, 1.4, 1.5, 8, 1.6, 1.8)'
#     * __e1f2, Top 1__: root -b -q 'procYields.C("yield.root", "1:2:3:4", 1, 2.0, 2.4, 24, 1.3, 1.9)'
#         
# ##Helicity Dependent   
# * __e1fs1__: procyh_e1fs1 
# * __e1fd1__: procyh_e1fd1
# * __e1fs2__: procyh_e1fs2
# * __e1fd2__: procyh_e1fd2
# * __e1f_all__: procyha
# 
# * Each Top individually, for example, __Top 5__:
#     * __e1f1, Top 1__: root -b -q 'procYields.C("yield.root", "1:2:3:4", 1, 1.4, 1.5, 8, 1.6, 1.8,1)'
#     * __e1f2, Top 1__: root -b -q 'procYields.C("yield.root", "1:2:3:4", 1, 2.0, 2.4, 24, 1.3, 1.9,1)'

# <markdowncell>

# ##3. processed yields $\rightarrow$ process simstats $\rightarrow$ plot simstats
# 1. Make sure that __proc_simstats.h__ is updated with latest __cumsim__ (__I should automate this__)  
# 2. Make sure that __proc_yields__ has been called for the latest __cumsim__
# 3. proc_simstats("e1fs2")
#         
# ### Helicity dependent 
# * Add optional 2nd argument (set to UNPOL by default) i.e. proc_simstats("e1fs2",POS)
# 
# ### Plot Simstats
# * plot_simstats.py --e1fs2 --hel=UNPOL [--help to confirm]
# 
# ### Helicity dependent
# * plot_simstats.py --e1fs2 --hel=POS
#     

# <markdowncell>

# ## Extracting Observables
# 
# ### Polarization Obervables
# 
# Each of the following commands, executed directly in __anadirs__ extracts Beam __Polarization Observable $R^{\theta}_{LT^{'}}$__ using different techniques/methods:
# 
# 1. __xsec-tools::plotphi__ called by __xpo__ $\longrightarrow$ __polobs__: First pass to extract $R^{\theta}$ by making __projphi[thetabin]__ and then performing $\int\sin\phi$
# 1. __xsec-tools::plotr__ called by __xpo2__ $\longrightarrow$ __polobs2__: hY5D*$\sin\phi$ and directly project on to __theta__ to extract $R^{\theta}$; __"hack" method to extract $R^{\theta}$ only for h=POS__
# 1. __xsec-tools::plotr3__ called by __xpo3__ $\longrightarrow$ __polobs3__: hY5D*$\sin\phi$ and directly project on to __theta__ to extract $R^{\theta}$; __ Only in Experimental Phase Space__
# 1. __xsec-tools::plotr4__ called by __xpo4__ $\longrightarrow$ __polobs4__: Basically a modification of #1 in that __projphi[thetabin]__ are made using only using __hY5D[ACC_CORR]__ instead of __hY5D[FULL]__

# <markdowncell>


# <codecell>



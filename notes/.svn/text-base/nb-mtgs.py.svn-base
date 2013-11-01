# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Mtg 10-02-13
# * __Discuss which Topologies so far seem usable given current accumulation of statistics, by referring to:__
#     
#     1. $\Large Y^{F}_{i} = Y^{\text{C-combins}}_{i} + Y^{\text{C-noncombins}_{i}}_{i} +Y^{\text{H-combins}}_{i} + Y^{\text{H-noncombins}_{i}}_{i}$
#     
#         1. $Y^{\text{C-combins}}_{t1} > Y^{\text{C-combins}}_{t2}$. This is because:
#             * $\Large \bf \left(\frac{n^{ER}_{t1}}{n^{ER}_{t2}}\right)_j \approx 1$. 
#                 * In words: Even though the Acceptance may be accurate i.e $(A_{t1})_{j} < (A_{t2})_{j}$, the fact that $({ER}_{t1})_{j} \not< ({ER}_{t2})_{j}$,  the Acceptance Corrected yield in Top 1 is inaccurately calculated to be higher i.e. $({EC}_{t1})_{j} > ({EC}_{t2})_{j}$
#             * __This also causes the Hole-Normalization Factor for Top1 to be biased to higher values__
#             
#         1. $Y^{\text{H}}_{t1} >> Y^{\text{H}}_{t2}$. This is because
#             * $Y^{\text{H}}_{t1}$ are _highly_ overestimated due to A.
#             * I don't think I am quite sensitive to the following reasons, but am listing them for completeness and to investigate in the future
#                 * $\text{nSR=0}_{t1}, \text{nSR=0}_{t2}$ are approaching their stable values. However till their ratio i.e. $\frac{\text{nSR=0}_{t1}}{\text{nSR=0}_{t2}} \neq \text{Constant}$, I could be biased in my Hole-Filling:
#                     * Infact this ratio continues to increase with increased Simulation Statistics. This implies that I am Hole-Filling more bins where ER=0 due to statistics in Top1 than in Top2.
#                 * There is always the possibility that an "Imperfect Simulation" is affecting the Topology with more Physical Holes (in this case Top 1) since Hole-Filling is Model Dependent.
#             
#             > Therefore, $\Large Y^{F}_{t1} >> Y^{F}_{t2}$ not as much due to "Imperfect Simulation" (or that I am not sensitive to that effect yet), but because
#             >
#             > $\Large \bf \left(\frac{n^{ER}_{t1}}{n^{ER}_{t2}}\right)_j \approx 1$
#             
#     1. 1D Observables in each Topology
#         * _Compare with Isupov (and others) and does a particular Top agree more with Current Results?_
#         
# * __Discuss Pol. Obs. thus far extracted__
#     * __ Notable Observations__
#         * Error/5D-bin $\equiv \delta n_{i} \approx \sqrt n_{i}$
#         * This is true even after Acceptance Correction. The error propagation ($C = \frac{R}{A}$) is dominated by the relative errors in the Reconstructed Events per bin (relative errors in the Acceptance per bin is negligible, approaching  0 with increasing number of event Simulated): 
#         
#             * $\bf{n^{R}_{i} \approx 1}$ and therefore, $\Large \bf \frac{\delta n^{R}_{i}}{n^{R}_{i}} \approx 1$
#     * __ Open Issues__
#         * Statistical errors do not appear to be approaching required sensitivity
#             * Am losing events somewhere in my analysis stage or is it that my Acceptance is underestimated?
#             * Increase Q2Wbin size?
#             * Is assumption of Poisson Statistics per 5D-bin appropriate?
#             
#             
# ## Mtg 10-02-13 Action Items to work on
# 
# 1. Do Simulation results show any $\sin\phi$ dependence?
# 1. For Pol. Obs: When to multiply by $\sin\phi$:
#     + Multiply hY5D by $\sin\phi$ and then Project on to $\phi$
#     + Directly Project on to $\phi$ and then multiply by $\sin\phi$
# 1. Does my argument on why $Y^F_{t1} >  Y^F_{t2}$ still hold (given that I had not taken into account just how "sparse" my hY5Ds were i.e. nevt $\leq$ nbins and therefore the reason why $\bf \left(\frac{n^{ER}_{t1}}{n^{ER}_{t2}}\right)_j \approx 1$. On this was premise was based my entire argument.)
# 

# <markdowncell>

# ## Mtg 10-09-13
# 
# 1. Presented results from Action Item 1 (AI) from last meeting (10-02-13)
#     * Simulation results __do__ show a $\sin\phi$ dependence; in other words $\bf \mathcal{R}^{1\theta}_{{LT}^{'}}(\theta) \neq 0$, __for both SR and ST events__ (results in polobs.mtg.100913 and polobs.sim.mc.100913 respectively)
#     * <font color="red"> Need to confirm with Victor why we are observing a $\sin\phi$ dependence in Simulation if the Model is constrained by Unpolarized Cross Section results from Experiment $\longrightarrow$ AI 1. </font>
#     
#  1. Observations on Pol. Obs. Extracted:
#      * Tops. with maximum Holes have a $\bf \mathcal{R}^{1\theta}_{{LT}^{'}}(\theta)$ that is closer to Simulation
#          + This __perhaps__ is due to the fact that these Tops. are __heavily__ influenced by Hole filling from Thrown events
#      * Tops. with least holes:
#          + $\mathcal{R}^{1\theta}_{{LT}^{'}}(\theta)$ for h=0 move closer to 0, which is what one would expect
#          + <font color="red"> __However__, $\mathcal{R}^{1\theta}_{{LT}^{'}}(\theta)$ for h=+/- do not move away significantly from 0 either. </font>
#              * My hypothesis at this time is that even though there are fewer Holes, they are still significant compared to the number of ER bins.
#                  * Should check if this changes with "new" method (AI 2) of extractign Pol. Obs. in which Hole filling will be "smarter" about Hole filling for h=+ and h=- cases.
#                  * To confirm this, in __simstats__, should additional diagnostics to keep track of: (AI 3)
#                      1. Fer
#                      1. Fer/Esr
#        
# ## Mtg 10-09-13 Action Items
# 
# 1. Need to confirm with Victor why we are observing a $\sin\phi$ dependence in Simulation if the Model is constrained by Unpolarized Cross Section results from Experiment
# 1. "New" method to extract Beam Pol. Obs. that is "smarter" about Hole filling depending on if h=+/-.
#     * Alternatively, try Method 3. ?
# 1. In __simstats__, should additional diagnostics to keep track of
#     1. Fer
#     1. Fer/Esr

# <markdowncell>

# ## Mtg 10-16-13
# 
# ### Currently, I am trying to extract $\mathcal{R}^{1\theta}_{{LT}^{'}}$
# 
# * Reminder:
#     1. $ M_{p\pi^{+}}, M_{\pi^{+}\pi^{-}}, \theta_{\pi^{-}}, \phi_{\pi^{-}}, \alpha_{[p^{'}\pi^{+}][p\pi^{-}]}$
#     2. $ M_{p\pi^{+}}, M_{\pi^{+}\pi^{-}}, \theta_{p}, \phi_{p}, \alpha_{[\pi^{+}\pi^{-}][p^{'}p]}$
#     3. $ M_{p\pi^{+}}, M_{p\pi^{-}}, \theta_{\pi^{+}}, \phi_{\pi^{+}}, \alpha_{[p^{'}\pi^{-}][p\pi^{+}]}$
#  
# #### $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}d\phi^{i}} }\right)^{h}  = A^{ij} + B^{ij}\cos\phi^{i} + C^{ij}\cos2\phi^{i} + hPR^{ij}_{LT^{'}}\sin\phi^{i} $$
# 
# * Method: For each $X^{ij}$ bin, project 5D-yield on to $\phi^{i}$ and then use the orthogonality of $sin$ and $cos$ to extract $\mathcal{R}^{ij}_{{LT}^{'}}$
# 
# #### $\int f(\phi)\sin\phi d\phi \longrightarrow hP(\mathcal{R}^{ij}_{k})_{LT^{'}}\int\sin^2\phi d\phi = hP(\mathcal{R}^{ij}_{k})_{LT^{'}}\pi $
# 
# ### Hole Filling:
# * How should we fill Holes?
#     + $\mathcal{R}^{ij}_{{LT}^{'}}$ extracted from Simulation $\neq 0$
#     + $\therefore$ tried Hole-Filling two ways:
#         1. Assuming $\sigma(\text{genev}) = \sigma^{\text{unp}}$
#         2. Assuming $\sigma(\text{genev}) = \sigma^{\text{pos}}$
#             1. __filled Holes__ for only 5D-yields filled with events from __positive__ Beam Polarization $\rightarrow (\text{hY5D_FULL}^{+}$).
#             2. __did not__ fill Holes for hY5D filled with events from __negative__ Beam Polarization $\rightarrow (\text{hY5D_ACC_CORR}^{-}$).
#             3. Then followed "Method" for $\text{hY5D_FULL}^{+}$, but for $\text{hY5D_ACC_CORR}^{-}$, instead of using $\int f(\phi)\sin\phi d\phi$, I used:
#             
#             __-__$\int f(\phi)\sin\phi d\phi$
#             
#             $\therefore$ extracting a measure "proportional" to the actual $\mathcal{R}^{ij}_{{LT}^{'}}$
#             
#             __NOTE__ In both the ways, 5D-yields for unpolarized yields was treated the same way and the Method was directly applied to it
# 
# ### Current Observations:
# * $\mathcal{R}^{ij}_{{LT}^{'}}$ (h=unp/pos/neg) extracted from fully exclusive topology (Top 1) seems to be strongly "biased" (to Model, since it has maximum Holes?)
#     + Additionally $\mathcal{R}^{ij}_{{LT}^{'}}$ (h=unp/pos/neg) are similar for both Q2W regions
# * <font color="red"> $\mathcal{R}^{ij}_{{LT}^{'}}$ (h=unp/pos/neg) extracted from maximal topology (Top 5) does not show enough statistical sensitivity </font>
# ### Accessing e16 Beam Helicity information
# * .bos
#     + helvar = ?
# * h10 (bos --> h10 = nt10maker:h2root)
#     + no beam helicity information
# * h21 (bos --> h21 using ?:h2root)
#     + helvar = irun [run number * helicity status (+1/-1)]
# * * *
# 
# ### CLAS code
# CLAS code is all located in /group/clas/builds ?

# <markdowncell>

# ## General Notes
# * __e1f__ and __e16__ datasets to use (as Evan found out on 10-16-13 that e16 pass2/v1 had "timing problems")
#     + e1f: pass2/v1
#     + e16: pass1/v2 (has Timing Correction, but not EC energy corrections) 
#         + e16: pass2/v1 (has no Timing Correction, but EC energy loss corrections)
# 
# ## General Action Items
# * __Fix the technical matters of obtaining hY5D-HOLE__ (noted in xsectools::plotxsec_CommonBins())
# 
# ## Long Term Action Items
# * __Cuts & Corrections__
#     + e1f:
#         - Implement Hadron Fiducial Cuts
#         - Revise all other Corrections and Cuts
#     + e16:
#         - Implement all cuts and corrections (use Evan's)
# * __e16 exp-data__
#     * Need to make it consistent with __h10__ format used for __e1f-exp,e1f-sim__
# * __Simulation__
#     + Organizational
#         * Move to CentOS 6.2
#         * Reorganize configuration files and sim-dir
#     + Procedural
#         * Simulate entire Q2W area: __$\text{[2.0-5.0]}GeV^2 \text{by [1.3-3.0]GeV}$__ for
#             + e1f
#             + e16
#         * Make hY8D for sim-data on the farm
# * __Code Integrity__
#     * Need to think think of more appropriate terms for what I currently call __Top__ and __Varset__. Perhaps, __EvtSel__ and __Top__ respectively.
#     * Make the code adhere to Coding Convetions (noted in __zim__)
#     * Remove redundancy (for example, identify and remove redundancy in functions ProcYields::PlotPhi() and xsec-tools::plotphi()) 
#     * Encapsulate workign code into Functions.
# * __Code Infrastructure__
#     * Use TEntryList instead of making separate directories as per Q2W skim
#         + Currently this is only applicable for Experimental data. But ultimately, __if I Simulate a larger Q2W area__, this may be applicable for Simulation too and simplify the infrastructure
#         + Remove redundant environment variables from .tcshrc (for example, I now have $E1F_2PI_DATADIR and /e1f.2pi.datadir)

# <codecell>



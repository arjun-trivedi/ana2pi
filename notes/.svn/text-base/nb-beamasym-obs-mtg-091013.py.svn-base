# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# #Extracting beam helicity dependent observables from reaction: $p(\vec{e},e^{'}\pi^{+},\pi^{-})p$ 

# <markdowncell>

# ## At a particular photon virtuality ($Q^{2}, W$), the cross section can be written as
# ### $\frac{d^2\sigma}{dQ^{2}dW} = \Gamma\int\frac{d^5\sigma_{v}}{d\tau^5}d\tau^5$
# where:
# 
# * $\sigma_{v} = f(Q^{2}, W)$ is the cross section for the virtual photon and proton interaction
# * $\Gamma = f(Q^{2}, W)$ is the Virtual Photon Flux
# * $\tau$ is a symbol representing a point in the 5 dimensional kinematic phase space that uniquely identifies a particular configuration of the final state ($p\pi^{+}\pi^{-}$)
# * In this analysis, there are 3 sets of variables (varsets) that are used to represent this 5 dimensional space:
#     1. $ M_{p\pi^{+}}, M_{\pi^{+}\pi^{-}}, \theta_{\pi^{-}}, \phi_{\pi^{-}}, \alpha_{[p^{'}\pi^{+}][p\pi^{-}]}$
#     2. $ M_{p\pi^{+}}, M_{\pi^{+}\pi^{-}}, \theta_{p}, \phi_{p}, \alpha_{[\pi^{+}\pi^{-}][p^{'}p]}$
#     3. $ M_{p\pi^{+}}, M_{p\pi^{-}}, \theta_{\pi^{+}}, \phi_{\pi^{+}}, \alpha_{[p^{'}\pi^{-}][p\pi^{+}]}$

# <markdowncell>

# ## Cross sections $\longrightarrow$ Observables
# 
# ###1. Single-differential cross sections: $\frac{d\sigma}{dX^{ij}}$ $\longrightarrow$ $A_{\frac{1}{2}}, A_{\frac{3}{2}}, S_{\frac{1}{2}} $
# 
# where
# 
# * $i$ = index over the 3 varsets
# * $j$ = index over the 5 kinematic variables in each varset
# * Total = 9 (Ignoring 3 $\phi$ & 3 repeated invariant mass cross sections )
# 
# ###2. Helicity dependent $\phi$ distributions (direct cross sections, asymmetries) $\longrightarrow$ $\mathcal{R}_{LT^{'}}$ 
# 
# __When the beam helicity is known, from single pion electroproduction, we know the cross section to have the following form:__
# 
# (There are 2 d.o.f. for single pion production, namely, $\theta_{\pi}, \phi_{\pi}$):
# 
# #### $\Large \frac{d^2\sigma_{v}}{ {d\Omega} }^{h}_{\pi} = 
# \frac{p_{\pi}}{k_{\gamma^{*}}}
# ( \sigma_{unpol} + hP\sqrt{ 2\epsilon_{L}(1-\epsilon) }R_{{LT}^{'}}\sin\theta_{\pi}\sin\phi_{\pi} )$
# 
# where
# 
# * $p_{\pi^{-}} = f(Q^{2},W,\theta_{\pi^{-}})$
# * $k_{\gamma^{*}},\epsilon,\epsilon_{L} = f(Q^{2},W)$
# * $R_{L}, R_{T}, R_{LT}, R_{TT},R_{{LT}^{'}} = f(Q^{2},W,\theta)$; Define  $\bf R^{\theta}_{LT^{'}} \equiv R_{{LT}^{'}}(Q^{2},W,\theta)$
# * $h$ = helicity
# * $P$ = "degree" of Polarization
# * $\sigma_{unpol} = R_{T} + \epsilon_{L}R_{L} + \sqrt{ 2\epsilon_{L}(1+\epsilon) }R_{LT}\cos\phi_{\pi} + \epsilon R_{TT}\sin^2\theta_{\pi}\cos 2\phi_{\pi} $
# 
# __Motivated by the single pion channel, I state that for the double pion channel, we may expect the following equivalent form for the cross secions:__
# 
# For $X^{ij} \neq \phi^{i}$
# 
# #### $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}d\phi^{i}} }\right)^{h}  = A + B\cos\phi + C\cos2\phi + hPDR^{ij}_{LT^{'}}\sin\phi $$
# 
# where
# 
# * $A, B, C, D, R^{ij}_{LT^{'}}  = f(Q^{2},W, X^{ij})$
# 
# Define $\bf \mathcal{R}^{ij}_{LT^{'}} \equiv DR^{ij}_{LT^{'}}$ ; $ \text{If R is the classical Response Function, should I define} \mathcal{R} \text{ to be JM Model-Response Function ?} $
# 
# There, we now have
# 
# #### $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}d\phi^{i}} }\right)^{h}  = A + B\cos\phi + C\cos2\phi + hP\mathcal{R}^{ij}_{LT^{'}}\sin\phi $$

# <markdowncell>

# ## Extracting $\mathcal{R}^{ij}_{LT^{'}}$
# 
# For a particular bin of $X^{ij}$, project out histogram on $\phi$ axis:
# 
# ### $ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}d\phi^{i}} }\right)^{h} \longrightarrow \left(\frac{d^2\sigma_{v}}{ {dX^{ij}_{k}d\phi^{i}} }\right)^{h}$
# where $k$ = index over bins in $X^{ij}$
# 
# Now
# 
# $A, B, C, D, \mathcal{R}^{ij}_{LT^{'}}  \longrightarrow f(Q^{2},W, X^{ij}_{k})$; $k$ = index over bins in $X^{ij}$
# 
# Define
# 
# * $ A_{k}, B_{k}, C_{k} \equiv A, B, C(Q^{2},W, X^{ij}_{k}) $
# * $ \bf (\mathcal{R}^{ij}_{k})_{LT^{'}} \equiv \mathcal{R}^{ij}_{LT^{'}}(Q^{2},W, X^{ij}_{k}) $
# 
# Therefore, we can now write
# 
# #### $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}_{k}d\phi^{i}} }\right)^{h}  = A_{k} + B_{k}\cos\phi + C_{k}\cos2\phi + hP(\mathcal{R}^{ij}_{k})_{LT^{'}}\sin\phi $$
# 
# 
# 
# ### To get to $hP(\mathcal{R}^{ij}_{k})_{LT^{'}}$ 
# 
# ####1. Fit and extract $A_{k}, B_{k}, C_{k}, hP(\mathcal{R}^{ij}_{k})_{LT^{'}}$ as parameters of the fitting function
# 
# ####2. $\int f(\phi)\sin\phi d\phi \longrightarrow hP(\mathcal{R}^{ij}_{k})_{LT^{'}}\int\sin^2\phi d\phi = hP(\mathcal{R}^{ij}_{k})_{LT^{'}}\pi $
# 
# ####3. Calculate Asymmetry, $(A^{ij}_{k})_{LT^{'}}$
# ###$(A^{ij}_{k})_{LT^{'}} \equiv \left(\frac{d^2\sigma_{v}}{ {dX^{ij}_{k}d\phi^{i}} }\right)^{+} - \left(\frac{d^2\sigma_{v}}{ {dX^{ij}_{k}d\phi^{i}} }\right)^{-} = 2P(\mathcal{R}^{ij}_{k})_{LT^{'}}\sin\phi$
# 
# and either 
# 
# * fit $(A^{ij}_{k})_{LT^{'}}$
# * $\int (A^{ij}_{k})_{LT^{'}} \sin\phi d\phi$ 
# 
# to obtain $hPD_{k}(R^{ij}_{k})_{LT^{'}}$

# <markdowncell>

# #### Some advantages and disadvantages of the 3 different methods
# 
# * 1,2
#     * Hole and therefore, model dependent. 
# * 3
#     * Hole and therefore, model independent (contributions from Holes are "subtracted out")
#     * "More statistical errors": The process of subtracting $Y^{+} - Y^{-}$ "makes the statistical error bars larger"

# <markdowncell>

# ## Feedback from meeting
# 
# * Use Method 2. to begin with; Methods 1 & 3 for cross checks
# * Victor already has an indirect access to "${LT}^{'}$" observable from unpolarized cross section data. This is because in the JM Model, _all Reaction Amplitudes_ are incorporated to predict observed cross sections; Hence all observables have to be "tweaked". Therefore, an estimate on "${LT}^{'}$" is also obtained. 
#     * My analysis will be the first "direct" probe of this Observable

# <codecell>



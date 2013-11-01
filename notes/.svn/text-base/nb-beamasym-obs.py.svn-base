# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Definition of $A_{{LT}^{'}}$ (Beam Asymmetry Observable) for "my reaction"

# <markdowncell>

# [Source: Kijun's thesis (2006), page 164](file/KJPark_thesis.pdf)

# <markdowncell>

# ## Reaction
# * $\gamma^{*} p \rightarrow p \pi^{+} \pi^{-}$
#     * $\gamma^{*} = e - e^{'}$ 
#     * Polarization of $e$ is known 
# 
# ## Degrees of freedom
# * 5 $(M_{1},M_{2},\theta,\phi,\alpha)$
# * For example: $ M_{p\pi^{+}}, M_{\pi^{+}\pi^{-}}, \theta_{\pi^{-}}, \phi_{\pi^{-}}, \alpha_{[p^{'}\pi^{+}][p\pi^{-}]}$
# 
# ## At a particular photon virtuality i.e. particular $W,Q^{2}$ (in CMS, unless noted otherwise)
# * $\Large \frac{d\sigma}{ {d\Omega} }^{h}_{\pi^{-}} = 
# \frac{p_{\pi^{-}}}{k_{\gamma^{*}}}
# ( \sigma_{unpol} + hP\sqrt{ 2\epsilon_{L}(1-\epsilon) }R_{{LT}^{'}}\sin\theta_{\pi^{-}}\sin\phi_{\pi^{-}} )$
#     * h = helicity
#     * P = "degree" of Polarization
#     * $p_{\pi^{-}} = f(Q^{2},W,\theta_{\pi^{-}})$
#     * $\sigma_{unpol} = f(Q^{2},W,\theta_{\pi^{-}},\phi_{\pi^{-}})$
#     * $k_{\gamma^{*}},\epsilon,\epsilon_{L} = f(Q^{2},W)$
#     * $R_{{LT}^{'}} = f(Q^{2},W,\theta)$
#     
#     
# * $\Large \int{ \frac{d\sigma}{ {d\Omega} }^{h}_{\pi^{-}} d\theta_{\pi^{-}}} = \frac{d\sigma}{ {d\phi} }^{h}_{\pi^{-}}$
#     * $p_{\pi^{-}} \rightarrow \int{p_{\pi^{-}} d\theta_{\pi^{-}}} = f(Q^{2},W)$
#     * $\sigma_{unpol} \rightarrow \int{\sigma_{unpol} d\theta_{\pi^{-}}} = f(Q^{2},W,\phi_{\pi^{-}})$
#     * $R_{{LT}^{'}} \rightarrow \int{R_{{LT}^{'}} d\theta_{\pi^{-}}} = f(Q^{2},W)$
#     
#     
# * $\text{Assuming } P^{+} = P^{-}$
#     * $\frac{d\sigma}{ {d\phi} }^{+}_{\pi^{-}} - \frac{d\sigma}{ {d\phi} }^{-}_{\pi^{-}} = 2P\sqrt{ 2\epsilon_{L}(1-\epsilon) }R_{{LT}^{'}}\sin\phi_{\pi^{-}}$
#         
# 
# ## $\therefore d\sigma^{+} - d\sigma^{-} (\phi_{\pi^{-}} + d\phi_{\pi^{-}}) = \frac{p_{\pi^{-}}}{k_{\gamma^{*}}}P\int{(2\sqrt{ 2\epsilon_{L}(1-\epsilon) }R_{{LT}^{'}}\sin\phi_{\pi^{-}}) {d\phi_{\pi^{-}}}}$
# 
# ## $\text{Define }A_{{LT}^{'}}  = \frac{d\sigma^{+}-d\sigma^{-}}{P}$
# * $\text{Using } {d\sigma^{\pm}} = \frac{N^{\pm}}{L^{\pm}}$
# * $\text{Assuming } L^{+} = L^{-}$
# 
# ## $\therefore A_{{LT}^{'}} (\phi_{\pi^{-}} + d\phi_{\pi^{-}}) = \frac{N^{+}(\phi_{\pi^{-}} + d\phi_{\pi^{-}}) - N^{-}(\phi_{\pi^{-}} + d\phi_{\pi^{-}})}{P}$
# * Using this definition, we can extract $R_{{LT}^{'}} = f(Q^{2},W)$
# 
# * * *
# 
# ## Other definitions I have seen used, for example, by Kijun $(A_{{LT}^{'}} = \frac{(N^{+} - N^{-})}{P(N^{+} + N^{-})})$ and why it will not work for me:
# #1 
# >The above definition is useful when $A_{{LT}^{'}}$ is __directly__ calculated in the bin of interest $(\phi \text{ or }\Omega\text{ bin})$. In that case, the Acceptances do not need to be calculated, since they cancel out. In my case $N^{\pm}(\phi)$ is obtained by summing over bins in the remaining 4 d.o.f. In this case normalizing by $(N^{+} + N^{-})$ does not cancel Acceptances. 
# 
# * (__I have to carry out my analysis in each of the 5 d.o.f. bins to obtain Model independent results__ $\rightarrow$ verify my understanding of why 5D analysis in nb-proposal)
# * __Is my analysis the first attempt to extract $A_{{LT}^{'}}$ from double charged pion channel?__. If so, then these __new definitions__ may be important to note
# 
# > Therefore, in my case, acceptances have to be calculated in each bin, before $(N^{+} \pm N^{-})$ can be obtained
# 
# #2 
# > Even though after 1. I can use the noramlized definition of $A_{{LT}^{'}}$, the result will be Model dependent in that while the Holes drop out when calculating $(N^{+} - N^{-})$, they remain when calculating $(N^{+} + N^{-})$
# 
# ### To DO
# * Obtain P and verify $P^{+} = P^{-}$
# * BCA (= Charge Normalization ?)
# * RC
# * BCC
# 
# 
#    

# <markdowncell>

# 
# * $\Large k_{\gamma^{*}} = f(W;m_{p})$
# * $\Large \sigma_{unpol} = R_{T} + \epsilon_{L}R_{L} + \epsilon R_{TT} \sin^{2}\theta\cos2\phi + \sqrt{ 2\epsilon_{L}(1-\epsilon) }R_{LT}\cos\phi$ 
# * $\Large \epsilon = f(Q^{2},W)$
# * $\Large\epsilon_{L} = f(\epsilon, Q^2, W)$

# <markdowncell>

# from IPython.display import display, Math, Latex
# display(Math(r'F(k) = \int_{-\infty}^{\infty} f(x) e^{2\pi i k} dx'))

# <codecell>

from IPython.display import display, Math, Latex 
display(Math(r'F(k) = \int_{-\infty}^{\infty} f(x) e^{2\pi i k} dx'))
#display(\color{red}\text{test})
display(Math(r'\text{\color{red}test}'))

# <codecell>



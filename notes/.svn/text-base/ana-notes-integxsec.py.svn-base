# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Checklist to investigate why cross sections obtained by me are systematically lower than Evgeny's:
# ## Cuts:
# ### 1. <strike>'ProcQ2WSel' applied to exp-data and not sim-data. </strike> Moved to Odds &amp; Ends
# - 08-01-13: I do not expect this to have any significant effect, but have still taken this under consideration. 
# - 08-02-13:
#     - As I was setting up to investigate this effect, I realized my machinery was too heavy; I am carrying around all topologies (5: top0=(top1+top2+top3+top4), top1, top2, top3, top4) and each has 3 sets of variables (Varsets).I am not even using all of the Tops and Varsets therein currently. This excess baggage is really making ProcEvt2Pi very sluggish and high on resources.  
#     - Also after talking with Gothe, I feel like I need to start extracting polarization observables, even if the absolute cross section I obtain is still below Evgeny's. This is more in line with the spirit of the presentation I will make in DNP13.
# 
#         &rarr; I think I need to invest a couple of days to reduce the excess baggage so that I may reduce the time it takes to produce the minimum set of observables.
# - 08-03-13: 
#     - I have made the machinery much lighter; the performance statistics are noted in the Observations sections. 
#     - __The fact that ProcQ2WSel is not applied to sim-data makes absolutely no difference.__ However, for book keeping, the fact that ProcQ2WSel is not applied to sim-data may be problematic, at least while debugging. I could just apply it, since as I concluded that it makes no difference. However, ProcQ2WSel is not auto-configured to pick up the appropriate Q2W range and implementing it, at this stage, I do not find important. I am going to move this matter to the Odds &amp; Ends section 
#     
# ### 2. Fiducial cuts for Hadrons
# * 09-01-13
#     > I am testing out this Block Quote to see what it does. 
#     > Seems to work good!
# 
#     

# <markdowncell>

# ###Checklist for further checking (hopefully should have a lesser contribution)
# ####08-01-13
# * Simulate events in a Q2W region that is a little wider than the region in which the observables are extracted. 
#     * Currently for 'e1fs1' and 'e1fs2' the simulated Q2W region is not so 
#         * sim-q2w(e1fs1) = observable-q2w(e1fs1)
#         * sim-w(e1fs2)   = observable-w(e1fs2)
#     * Apply ProcQ2WSel to exp-data for this wider Q2W region.
# 
# * Apply ProcQ2WSel to exp-data after ProcMomCor

# <markdowncell>

# ### Odds &amp; Ends
# ####08-03-13
# * Apply ProcQ2WSel to sim-data

# <markdowncell>

# 
# 

# <codecell>



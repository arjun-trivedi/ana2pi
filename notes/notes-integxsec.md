notes-integxsec
===============
Here I keep track of how my extraction of the Integrated Cross-section is coming along.

#### [12-02-13] 
The story so far is that, although the shape of my Cross-sections look acceptable, they are factor of ~4 smaller than Evgeny's. The factors that could be affecting this, that I have so far identified are:

1. Fiducial Cuts for hadrons are not applied
2. Acceptance may currently be overestimated, since I am using `mctk` banks in the h10-Tree created after Reconstruction. At this stage, the `mctk' bank retains only those events that have been Reconstructed. Hence, I am not accounting for Acceptance using all the Thrown events.


[10-10-13]Analyzing MM cuts
============================
In the SIM data, the MM distributions seem to:
+ Be offset from the True value (EXP ~closer to True value)
+ Have a smaller width (EXP width is wider)
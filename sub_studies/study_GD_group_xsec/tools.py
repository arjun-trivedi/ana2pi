def get_bins(bincl,nround_digits):
    '''
    + Returns binl from bincl where bin=[le,ue] and binc=bin center
    + Based on algorithm where the midway point to the previous and next bin define
    the bin edges:
    + le_i= (binc_i) - ( (binc_i)-(binc_i-1) )/2
    + ue_i= (binc_i) + ( (binc_i+1)-(binc_i) )/2
    where 'i' is the bin index
    '''
    binl=[]
    for i in range(len(bincl)):
        #!print "******"
        #! debug
#         print "(binc_i-1),(binc_i),(binc_i+1):"
#         if i==0:
#             print bincl[i],bincl[i+1]
#         elif i==len(bincl)-1:
#             print bincl[i-1],bincl[i]
#         else:
#             print bincl[i-1],bincl[i],bincl[i+1]
        
        if i==len(bincl)-1: #! For last bin assume ue_dist=le_dist
            ue_dist=(bincl[i]-bincl[i-1])/2
        else:
            ue_dist=(bincl[i+1]-bincl[i])/2
        if i==0: #! For first bin assume le_dist = ue_dist
            le_dist=ue_dist
        else:
            le_dist=(bincl[i]-bincl[i-1])/2
        ue=bincl[i]+ue_dist
        le=bincl[i]-ue_dist
        #!print "ue,le=",round(le,nround_digits),round(ue,nround_digits)
        binl.append([round(le,nround_digits),round(ue,nround_digits)])
        #!print "******"
    return binl

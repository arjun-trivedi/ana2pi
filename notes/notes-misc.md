[10-16-13]`e1f` and `e16` datasets to use (as Evan found out that e16 pass2/v1 had "timing problems")

* 	`e1f`: pass2/v1
* 	`e16`: pass1/v2 (has Timing Correction, but not EC energy corrections) 
	* e16: pass2/v1 (has no Timing Correction, but EC energy loss corrections)

## Observations
* 08-03-13
    * Loading a root file (that contains THnSparses) that is 150MB on disk, takes up ~1GB in RAM. Why is this so?
    * SelH10 performance statistics on gothe14 (6 cores; 12GB MEM):
        * procorder = q2wsel:eid:efid:qskim:pid:top
        * TProof:4 wrkrs
        * Max MEM/wrkr = 25% (=100/4 workers)
        * Max CPU/wrkr = 17% (=100/6 cores) 
<TABLE border=1>
<CAPTION>Example of a simple data table
created using HTML markup.</CAPTION>
<TR>
<TD></TD>
<TH>e1fs1:100M:yield</TH>
<TH>e1fs1:100M:mcyield</TH>
<TH>e1fs2:400M:yield</TH>   
<TH>e1fs2:400M:mcyield</TH>
</TR>
<TR>
<TH>MEM/wrkr[%]</TH>
<TD>~3</TD>
<TD>~3</TD>
<TD>~10</TD>
    <TD>~<font color="red">20</font></TD>
</TR>
<TR>
<TH>&lt;Event rate&gt;[Kevt/s];[MB/s]</TH>
<TD>150;25</TD>
<TD>166;27</TD>
<TD>127;21</TD>
<TD>167;28</TD>
</TR>
<TR>
<TH>Disk R rate[MB/s]</TH>
<TD>~24</TD>
<TD>~24</TD>
<TD>~24</TD>
<TD>~30*</TD>
</TR>
<TR>
<TH>Total time[min]</TH>
<TD>10</TD>
<TD>10</TD>
<TD>50</TD>
<TD>38</TD>
</TR>
<TR>    
<TH>output file size[MB]</TH>
<TD>32</TD>
<TD>32</TD>
<TD>239</TD>
<TD>288</TD>
</TR>
</TABLE>
*Seems like Disk R rate increased at the expense of %CPU/wrkr. While in all other instances the job ran at max %CPU/wrkr ~ 16, here it reduced to 10. I think this is because making mcyield has just one Proc and that does not need as much processing time viz-a-viz the full "procorder" to make yields from cooked data. 

Total time to make all yields (e1fs1/e1fd1/e1fs2/e1fd2@08-05-13): 1216.136u 108.406s 8:42:26.76 4.2%      0+0k 4613600+11032408io 3594pf+0w

Total time to process all yields (e1fs1/e1fd1/e1fs2/e1fd2@08-05-13): 9710.398u 85.521s 2:44:10.30 99.4%      0+0k 5710272+3430456io 83pf+0w

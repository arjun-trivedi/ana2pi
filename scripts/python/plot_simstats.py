#!/usr/bin/python

# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
from __future__ import division
import os, sys, getopt, shutil, datetime, time
from pandas import *
import matplotlib.pyplot as plt

pandas.set_printoptions(max_columns=20)

#get options
global anadir
global outdir
global csvfname
print sys.argv[1:]
opts, args = getopt.getopt(sys.argv[1:],"h",["help", "e1fs1", "e1fs2", "hel="])

for opt, arg in opts:
    if opt in ('-h', '--help'):
        print 'plot_simstats.py --<e1fs1/e1fs2> --hel=<UNPOL/POS/NEG>'
        sys.exit()
    elif opt == "--e1fs1":
        anadir = os.environ['E1F_2PI_ANADIR1']
    elif opt == "--e1fs2":
        anadir = os.environ['E1F_2PI_ANADIR2']
    elif opt == "--hel":       
         if arg=='UNPOL': 
            csvfname = 'simstats_vm.csv'
            outdir = 'simstats'
         elif arg=='POS': 
            csvfname = 'simstats_vm_POS.csv'
            outdir = 'simstats_POS'
         elif arg=='NEG': 
            csvfname = 'simstats_vm_NEG.csv'
            outdir = 'simstats_NEG'
         else:
            print arg,' is not recognized'
            sys.exit()
            
print 'anadir = ', anadir
print 'outdir = ', outdir
print 'csvfilename = ', csvfname

#back up old outdir if already exists, else make it 
if os.path.isdir(os.path.join(anadir,outdir)):
    outdir_bk = outdir+'_bk_'+datetime.now().strftime('%Y-%m-%d_%H-%M')
    print 'going to backup', outdir, 'to', outdir_bk
    shutil.move(os.path.join(anadir,outdir), os.path.join(anadir,outdir_bk))
    os.makedirs(os.path.join(anadir,outdir))
else:
    os.makedirs(os.path.join(anadir,outdir))

csvf = os.path.join(anadir,csvfname)
dfss = pandas.read_csv(csvf, na_values=['\n'])
#print df

# <codecell>

dfss_grpd_q2wbin = dfss.groupby('Q2Wbin')
nq2wbins = len(dfss_grpd_q2wbin.groups)
ngridx = 4
ngridy = nq2wbins/ngridx
print ngridx, ngridy

# <codecell>

for q2wbin in dfss_grpd_q2wbin.groups:
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    #if (q2wbin==1):
    #    print df

# <codecell>

itops = {0,1,2,3,4} #NOTE, that these are indices of tops, not tops; vm_tops used for now
colors = ['b','r','g','m','k']
markers = ['s','o','^','v','<']
labels = ['Top1','Top2','Top3','Top4','Top5']
sels  = {}
gr_Fth = {}
gr_Fsr = {}
gr_Esr = {}
gr_EsrVsEac = {}
gr_Esr_frc = {} #Esr[top]/Fth
gr_Esr_rat = {} #Esr[top]/Esr[top=2]

for q2wbin in dfss_grpd_q2wbin.groups:
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    
    fg_Fth_name  = str.format('Fth_%02d'%q2wbin)
    fg_Fth_title = str.format('Fth[Q2W=%02d]'%q2wbin)
    fg_Fth = plt.figure(fg_Fth_name,figsize=(8,5))
    ax_Fth = fg_Fth.add_subplot(1,1,1,title=fg_Fth_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel=df.nFth.name)
    
    fg_Fsr_name  = str.format('Fsr_%02d'%q2wbin)
    fg_Fsr_title = str.format('Fsr[Q2W=%02d]'%q2wbin)
    fg_Fsr = plt.figure(fg_Fsr_name,figsize=(8,5))
    ax_Fsr = fg_Fsr.add_subplot(1,1,1,title=fg_Fsr_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel=df.nFsr.name)

    fg_Esr_name  = str.format('Esr_%02d'%q2wbin)
    fg_Esr_title = str.format('Esr[Q2W=%02d]'%q2wbin)
    fg_Esr = plt.figure(fg_Esr_name,figsize=(8,5))
    ax_Esr = fg_Esr.add_subplot(1,1,1,title=fg_Esr_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel=df.nEsr.name)
    
    fg_EsrVsEac_name  = str.format('EsrVsEac_%02d'%q2wbin)
    fg_EsrVsEac_title = str.format('EsrVsEac[Q2W=%02d]'%q2wbin)
    fg_EsrVsEac = plt.figure(fg_EsrVsEac_name,figsize=(8,5))
    ax_EsrVsEac = fg_EsrVsEac.add_subplot(1,1,1,title=fg_EsrVsEac_title,xlabel=df.nEac.name,ylabel=df.nEsr.name)

    fg_Esr_frc_name  = str.format('Esr_frc_%02d'%q2wbin)
    fg_Esr_frc_title = str.format('Esr_frc[Q2W=%02d]'%q2wbin)
    fg_Esr_frc = plt.figure(fg_Esr_frc_name,figsize=(8,5))
    ax_Esr_frc = fg_Esr_frc.add_subplot(1,1,1,title=fg_Esr_frc_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel='holes[top]/nFth')

    fg_Esr_rat_name  = str.format('Esr_rat_%02d'%q2wbin)
    fg_Esr_rat_title = str.format('Esr_rat[Q2W=%02d]'%q2wbin)
    fg_Esr_rat = plt.figure(fg_Esr_rat_name,figsize=(8,5))
    ax_Esr_rat = fg_Esr_rat.add_subplot(1,1,1,title=fg_Esr_rat_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel='holes[top]/holes[2]')


    for itop in itops:
        sels[itop] = (df['Top']==itop+1) & (df['Varset']==1)
        dfsel = df[sels[itop]]
        gr_Fth[itop] = ax_Fth.scatter(dfsel.Sim, dfsel.nFth,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])

        gr_Fsr[itop] = ax_Fsr.scatter(dfsel.Sim, dfsel.nFsr,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])

        gr_Esr[itop] = ax_Esr.scatter(dfsel.Sim, dfsel.nEsr,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])

        gr_EsrVsEac[itop] = ax_EsrVsEac.scatter(dfsel.nEac,dfsel.nEsr,
                                                s=50,c=colors[itop],marker=markers[itop],label=labels[itop])

        gr_Esr_frc[itop] = ax_Esr_frc.scatter(dfsel.Sim, dfsel.nEsr/dfsel.nFth,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])
        
        gr_Esr_rat[itop] = ax_Esr_rat.scatter(dfsel.Sim, dfsel.nEsr.values/df[(df['Top']==2) & (df['Varset']==1)].nEsr.values,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])
    
    #save plots
    plt.figure(fg_Fth_name)
    plt.legend(loc='lower right',prop={'size':6})
    fg_Fth.savefig("%s/%s/%s.jpg"%(anadir,outdir,fg_Fth_name))

    plt.figure(fg_Fsr_name)
    plt.legend(loc='upper left', prop={'size':6})
    fg_Fsr.savefig("%s/%s/%s.jpg"%(anadir,outdir,fg_Fsr_name))

    plt.figure(fg_Esr_name)
    plt.legend(loc='upper left', prop={'size':6})
    fg_Esr.savefig("%s/%s/%s.jpg"%(anadir,outdir,fg_Esr_name))

    plt.figure(fg_EsrVsEac_name)
    plt.legend(loc='upper right',prop={'size':6})
    fg_EsrVsEac.savefig("%s/%s/%s.jpg"%(anadir,outdir,fg_EsrVsEac_name))

    plt.figure(fg_Esr_frc_name)
    plt.legend(loc='center',prop={'size':6})
    fg_Esr_frc.savefig("%s/%s/%s.jpg"%(anadir,outdir,fg_Esr_frc_name))

    plt.figure(fg_Esr_rat_name)
    plt.legend(loc='center',prop={'size':6})
    fg_Esr_rat.savefig("%s/%s/%s.jpg"%(anadir,outdir,fg_Esr_rat_name))

print 'plot_simstats.py::Success!'

# <codecell>



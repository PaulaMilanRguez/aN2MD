#  File: Final_analysis.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "Final_analysis.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#

import re
import numpy as np
import pandas as pd
import math
import mdtraj as md
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec


# # File Analysis
#
# 1. Chose a type of analysis: structure (only .gro needed) or simulation (.gro + .xtc)
# 2. Convert .itp to dataframe (.itp needed)
# 3. Search itp NOEs in the structure/simulation
# 4. Calculate distances (option 0, structure/ option 1, simple mean/ option 2, mean 6)
# 5. Compare RMN distance with strucuture/simulation distances (attribute rang = error percentage in RMN distance)
#     30% by default
# 6. Extract violations
# 7. Create plot

def Search_violations(name_itp, name_gro, name_xtc= '', option = 0, rang = 30,
                      skip=200, strd=100, name_work='results'):
    
    df_noes = []
    
    if option == 0:
        df_noes = Ana_gro(name_itp, name_gro, rang)
    
    elif (option ==1 or option == 2 or option == 3):
        df_noes = Ana_xtc(name_itp, name_gro, name_xtc, option, rang, skip,
                          strd)

    
    pd.DataFrame.to_csv(df_noes, name_work)
    df_viol = Extract_viol(df_noes, name_work)
    Classify_viol(df_noes, name_work)
    return df_noes, df_viol


# In[32]:


def Read_itp(name_itp, name_gro):
    
    file_noe = open(name_itp, "r")
    
    comp ={'Res1': [],'At_name1':[],'Atom1':[], 'Res2':[], 'At_name2':[],
           'Atom2':[], 'Dist_RMN':[], 'Ind':[], 'Distance':[], 'Violation':[]}    
    gro = md.load(name_gro)
    top  = gro.topology
    cpt = 0
    for linenoe in file_noe:

        cpt += 1
        found1 = False
        found2 = False
        if cpt > 2:
            linen = linenoe.split()
            atom1 = int(linen[0]) -1
            atom2 = int(linen[1]) -1
            dist_rmn = float(linen[12][0:4])*10
            dist_rmn = round(dist_rmn, 2)
            index = int(linen[3])   

            comp['Res1'].append(str(top.atom(atom1)).split("-")[0])
            comp['At_name1'].append(str(top.atom(atom1)).split("-")[1])
            comp['Res2'].append(str(top.atom(atom2)).split("-")[0])
            comp['At_name2'].append(str(top.atom(atom2)).split("-")[1])
                
            comp['Atom1'].append(atom1)
            comp['Atom2'].append(atom2)
            comp['Dist_RMN'].append(dist_rmn)
            comp['Ind'].append(index)
            comp['Distance'].append(0)
            comp['Violation'].append('Viol')

    df_comp = pd.DataFrame.from_dict(comp)
    file_noe.close()
    return df_comp


def Ana_gro(name_itp, name_gro, rang):
    
    df_noes = Read_itp(name_itp, name_gro)
    gro = md.load(name_gro)
    top  = gro.topology
    
    for i in range(0, 302):
        couples = []

        for row in df_noes[df_noes.Ind == i].itertuples():
            couples.append([row.Atom1, row.Atom2])

        dist = md.compute_distances(gro, couples)
        moy = np.mean(dist)*10
        moy = round(moy, 2)   
        df_noes.loc[df_noes.Ind == i, 'Distance'] = moy
    
    for row in df_noes.itertuples():
        
        lower = row.Dist_RMN - (row.Dist_RMN*rang)/100
        
        upper = row.Dist_RMN + (row.Dist_RMN*rang)/100        
        
        if(row.Distance < upper):    
            df_noes.loc[df_noes.Ind == row.Ind, 'Violation'] = 'Res'
            
    
    return df_noes


def Ana_xtc(name_itp, name_gro, name_xtc, option, rang, skip, strd):
    
    df_noes = Read_itp(name_itp, name_gro)
    traj = md.load_xtc(name_xtc, name_gro, strd)
    top = traj.top

    for i in range(0, 302):
        couples = []

        for row in df_noes[df_noes.Ind == i].itertuples():
            couples.append([row.Atom1, row.Atom2])

        dist = md.compute_distances(traj, couples)
        
        if option == 1:
            moy = np.mean(dist)*10
        elif option == 2:
            moy = 0
            cpt = 0
            for d1 in dist:
                for d2 in d1:
                    moy = moy + (1/math.pow(d2, 6))
                    cpt += 1
            moy  = moy/cpt
            moy = (1/(moy**(1/6)))*10  
            
        elif option == 3:
            moy_r3 = []
            for couple in range(0, len(dist[0])):
                moy = 0
                for frame in dist[skip::]:
                    moy = moy + (1/math.pow(frame[couple], 3))
                moy = moy/len(dist[skip::])
                moy = (1/(moy**(1/3)))*10
                moy_r3.append(moy)
            moy = 0
            for moy3 in moy_r3:
                moy = moy + (1/math.pow(moy3, 6))
            moy = 1/(moy**(1/6))
        moy = round(moy, 2)    
        df_noes.loc[df_noes.Ind == i, 'Distance'] = moy
            
    for row in df_noes.itertuples():
        
        lower = row.Dist_RMN - (row.Dist_RMN*rang)/100
        upper = row.Dist_RMN + (row.Dist_RMN*rang)/100        
        
        if(row.Distance < upper):
            df_noes.loc[df_noes.Ind == row.Ind, 'Violation'] = 'Res'
    
    return df_noes


def Extract_viol(comp, title):

    v = comp[comp.Violation=='Viol']
    v = v[v.Distance > 5.5]
    pd.DataFrame.to_csv(v, title+'_violations')
    return v


def Classify_viol(comp, name_fig):

    comp.head(15)

    closedxV = []
    closedyV = []
    closedxR = []
    closedyR = []
    closedxO = []
    closedyO = []
    
    distantxV = []
    distantyV =[]
    distantxR = []
    distantyR = []
    distantxO = []
    distantyO = []
 

    ind = []

    for constr in comp.itertuples():
        if constr.Violation == 'Viol':
            res1 = int(constr.Res1[3:5])
            res2 = int(constr.Res2[3:5])
            if (res2 - res1) == 0 :
                if constr.Ind not in ind:
                    if constr.Distance <= 5.5:
                        closedxO.append(res1)
                        closedyO.append(res2)
                        ind.append(constr.Ind)
                    
                    else:
                        closedxV.append(res1)
                        closedyV.append(res2)
                        ind.append(constr.Ind)
            else:
                if constr.Ind not in ind:
                    if constr.Distance <= 5.5:
                        distantxO.append(res1)
                        distantyO.append(res2)
                        ind.append(constr.Ind)
                        
                    else:
                        distantxV.append(res1)
                        distantyV.append(res2)
                        ind.append(constr.Ind)
        else:
            res1 = int(constr.Res1[3:5])
            res2 = int(constr.Res2[3:5])
            if (res2 - res1) == 0 :
                if constr.Ind not in ind:
                    closedxR.append(res1)
                    closedyR.append(res2)
                    ind.append(constr.Ind)
            else:
                if constr.Ind not in ind:
                    distantxR.append(res1)
                    distantyR.append(res2)
                    ind.append(constr.Ind)

    closedV = closedxV + closedyV
    closedV = list(set(closedV))
    closedR = closedxR + closedyR
    closedR = list(set(closedR))
    closedO = closedxO + closedyO
    closedO = list(set(closedO))
    
    Violations_plot(distantxR, distantyR, distantxV, distantyV, distantxO, 
                    distantyO, closedR, closedV, closedO, name_fig)


def Violations_plot(x1, y1, x2, y2, xO, yO, cR, cV, cO, title):  
    
    #Figure creation
    
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 0.5])
    ax0 = plt.subplot(gs[0])
    title_text = plt.title('Comparison between the NMR data and '+title+'\n\n'
                           , loc='center')
    title_text = plt.title(title+'\n\n', loc='center')
    title_text.set_weight('heavy')
    
    ############### First plot ################

    
    #Plotting matches
    
    ax0.scatter(x1, y1, color='#498e1b')
  
    
    # Plotting oranges
    
    ax0.scatter(xO, yO, color='#ea950b')
   
    
    # Plotting mismatches

    ax0.scatter(x2, y2, color='#dd0000')
  
    
    ax0.legend(['Matches', 'False Mismatches', 'Mismatches'], bbox_to_anchor=
               (1.05, 0.75), loc=2, borderaxespad=0.)


    #NCP7 zones
    
    ax0.plot( [14.5, 14.5], [1, 55], c ='#bfbfbf')
    ax0.plot( [28.5, 28.5], [1, 55], c ='#bfbfbf')
    ax0.plot( [35.5, 35.5], [1, 55], c ='#bfbfbf')
    ax0.plot( [49.5, 49.5], [1, 55], c ='#bfbfbf')

    ax0.plot( [1, 55], [14.5, 14.5], c ='#bfbfbf')
    ax0.plot( [1, 55], [28.5, 28.5], c ='#bfbfbf')
    ax0.plot( [1, 55], [35.5, 35.5], c ='#bfbfbf')
    ax0.plot( [1, 55], [49.5, 49.5], c ='#bfbfbf')


    #Diagonal
    ax0.plot([1, 55],[1, 55], c='#5b5a5a')

    # Plot visualisation
    

    subt1 = ax0.set_title('\n\nInter residue', loc='left')
    subt1.set_style("italic")
    ax0.set_xlim(1,55)
    ax0.set_ylim(1,55)
    ax0.set_yticks([7, 21, 32, 42, 52])
    ax0.set_yticklabels(['NTerm', 'ZF1', 'Linker', 'ZF2', 'CTerm'])
    ax0.set_xticks([7, 21, 32, 42, 52])
    ax0.set_xticklabels(['NTerm', 'ZF1', 'Linker', 'ZF2', 'CTerm'])
    
    
    ############### Second plot ################
    
    ax1 = plt.subplot(gs[1])
    
    #Plotting matches
    
    vec1 = np.repeat(1, len(cR))
    ax1.scatter(cR, vec1, color='#498e1b')

    # Plotting oranges
    
    vecO = np.repeat(1, len(cO))
    ax1.scatter(cO, vecO, color='#ea950b')    
    
    # Plotting mismatches
    
    vec2 = np.repeat(1, len(cV))
    ax1.scatter(cV, vec2, color='#dd0000')
  

    #NCP7 zones
    
    ax1.plot( [14.5, 14.5], [0, 2], c ='#bfbfbf')
    ax1.plot( [28.5, 28.5], [0, 2], c ='#bfbfbf')
    ax1.plot( [35.5, 35.5], [0, 2], c ='#bfbfbf')
    ax1.plot( [49.5, 49.5], [0, 2], c ='#bfbfbf')

    # Plot visualisation
    
    subt2 = ax1.set_title('\nIntra residue', loc='left')
    subt2.set_style('italic')
    ax1.set_xlim(1,55)
    ax1.set_ylim((0,2))
    ax1.set_yticks([])
    ax1.set_xticks([7, 21, 32, 42, 52])
    ax1.set_xticklabels(['NTerm', 'ZF1', 'Linker', 'ZF2', 'CTerm'])

    # Save and quit :)
    
    plt.tight_layout()
    #plt.close()
    fig.savefig(title+".png", dpi=300, pad_inches = 0.1, 
                bbox_inches= 'tight' )


# df_noes, df_viol = Search_violations('noeres.itp', 'pi_1MFS.gro', 'pi10.2_1MFS.xtc'
#                                      ,option = 3, rang = 30, skip = 200,strd=100,
#                                      name_work='test_inutil')


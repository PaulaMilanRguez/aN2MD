
#  File: Final_analysis.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "searchviolations.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#

import re
import numpy as np
import pandas as pd
import math
import mdtraj as md


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
    nb_noes = list(set(df_noes.Ind))
    
    for i in range(0, len(nb_noes)):
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

    nb_noes = list(set(df_noes.Ind))
    
    for i in range(0, len(nb_noes)):
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



def Search_Violations(name_itp, name_gro, name_xtc= '', option = 0, rang = 30,
                      skip=200, strd=100, name_work='results'):

    df_noes = []
    
    if option == 0:
        df_noes = Ana_gro(name_itp, name_gro, rang)
    
    elif (option ==1 or option == 2 or option == 3):
        df_noes = Ana_xtc(name_itp, name_gro, name_xtc, option, rang, skip,
                          strd)

    
    pd.DataFrame.to_csv(df_noes, name_work)
    return df_noes





#  File: Final_analysis.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "violations_plot.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#

import numpy as np
import pandas as pd
import math
import mdtraj as md
from matplotlib import pyplot as plt
from matplotlib import gridspec



def Classify_viol(comp, title, ncp7, res_max):


    closedxV = []
    closedyV = []
    closedxR = []
    closedyR = []
    closedxO = []
    closedyO = []

    distantxV = []
    distantyV = []
    distantxR = []
    distantyR = []
    distantxO = []
    distantyO = []

    ind = []

    for constr in comp.itertuples():
        if constr.Violation == 'Viol':
            res1 = int(constr.Res1[3:5])
            res2 = int(constr.Res2[3:5])
            if (res2 - res1) == 0:
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
            if (res2 - res1) == 0:
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

    Plot_viol(distantxR, distantyR, distantxV, distantyV, distantxO,
                    distantyO, closedR, closedV, closedO, title, ncp7, res_max)


def Plot_viol(x1, y1, x2, y2, xO, yO, cR, cV, cO, title, ncp7, res_max):
    # Figure creation
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 0.5])
    ax0 = plt.subplot(gs[0])
    title_text = plt.title('Comparison between the NMR data and '+title+'\n\n',
                           loc='center')
    title_text = plt.title(title+'\n\n', loc='center')
    title_text.set_weight('heavy')

    if(ncp7):
        pass
    else:
        if(res_max==0):
            res_max = x1+ y1+ x2+ y2+ xO+ yO+ cR+ cV+ cO
            res_max = max(res_max)
        

    # -------------- First plot ----------------
    # Plotting matches

    ax0.scatter(x1, y1, color='#498e1b')

    # Plotting oranges

    ax0.scatter(xO, yO, color='#ea950b')

    # Plotting mismatches

    ax0.scatter(x2, y2, color='#dd0000')

    ax0.legend(['Matches', 'False Mismatches', 'Mismatches'],
               bbox_to_anchor=(1.05, 0.75), loc=2, borderaxespad=0.)

    # NCP7 zones

    if(ncp7):

        ax0.plot([14.5, 14.5], [1, 55], c='#bfbfbf')
        ax0.plot([28.5, 28.5], [1, 55], c='#bfbfbf')
        ax0.plot([35.5, 35.5], [1, 55], c='#bfbfbf')
        ax0.plot([49.5, 49.5], [1, 55], c='#bfbfbf')

        ax0.plot([1, 55], [14.5, 14.5], c='#bfbfbf')
        ax0.plot([1, 55], [28.5, 28.5], c='#bfbfbf')
        ax0.plot([1, 55], [35.5, 35.5], c='#bfbfbf')
        ax0.plot([1, 55], [49.5, 49.5], c='#bfbfbf')

    # Diagonal
 
        ax0.plot([1, 55], [1, 55], c='#5b5a5a')

    else:
        ax0.plot([1, res_max], [1, res_max], c='#5b5a5a')

    # Plot visualisation
    subt1 = ax0.set_title('\n\nInter residue', loc='left')
    subt1.set_style("italic")

    if(ncp7):
        ax0.set_xlim(1, 55)
        ax0.set_ylim(1, 55)
        ax0.set_yticks([7, 21, 32, 42, 52])
        ax0.set_yticklabels(['NTerm', 'ZF1', 'Linker', 'ZF2', 'CTerm'])
        ax0.set_xticks([7, 21, 32, 42, 52])
        ax0.set_xticklabels(['NTerm', 'ZF1', 'Linker', 'ZF2', 'CTerm'])
    else:
        ax0.set_xlim(1, res_max)
        ax0.set_ylim(1, res_max)
    

    # -------------- Second plot ----------------
    ax1 = plt.subplot(gs[1])

    # Plotting matches
    vec1 = np.repeat(1, len(cR))
    ax1.scatter(cR, vec1, color='#498e1b')

    # Plotting oranges
    vecO = np.repeat(1, len(cO))
    ax1.scatter(cO, vecO, color='#ea950b')

    # Plotting mismatches
    vec2 = np.repeat(1, len(cV))
    ax1.scatter(cV, vec2, color='#dd0000')


    # NCP7 zones
    if(ncp7):
        ax1.plot([14.5, 14.5], [0, 2], c='#bfbfbf')
        ax1.plot([28.5, 28.5], [0, 2], c='#bfbfbf')
        ax1.plot([35.5, 35.5], [0, 2], c='#bfbfbf')
        ax1.plot([49.5, 49.5], [0, 2], c='#bfbfbf')

    # Plot visualisation
    subt2 = ax1.set_title('\nIntra residue', loc='left')
    subt2.set_style('italic')
    ax1.set_ylim((0,2))
    ax1.set_yticks([])
    if(ncp7):
        ax1.set_xlim(1,55)
        ax1.set_xticks([7, 21, 32, 42, 52])
        ax1.set_xticklabels(['NTerm', 'ZF1', 'Linker', 'ZF2', 'CTerm'])
    else:
        ax1.set_xlim(1,res_max)        

    # Save and quit :)
    plt.tight_layout()
    fig.savefig(title, dpi=300, pad_inches=0.1,
                bbox_inches='tight')
    plt.close()


def Violations_plot(df_noes, title= 'Violations_plot.png', ncp7=False, res_max=0):

    Classify_viol(df_noes, title, ncp7, res_max)
    

    






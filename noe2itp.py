#  File: noe2itp.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "noe2itp.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#
import pandas as pd
import numpy as np


def Assing_index(df_noes):   
    
    ''' Assing an .itp index to each NOE distance (required by GROMACS)
    
    Input:
    
        - df_noes (pandas dataframe): NOE distances dataframe

    Output:
    
        - df_noes (pandas dataframe) : NOE distances modified dataframe
    
    '''

    # pre-assign dummy index
    df_noes.loc[:, "index"] = -1
    # keep track of current index
    cidx = 0

    for i in range(df_noes.shape[0]):
        row = df_noes.iloc[i, :]
        # 0) filter on residue
        I = df_noes.loc[:, "ResID1"] == row.ResID1
        I &= df_noes.loc[:, "ResID2"] == row.ResID2

        # 1) match atom1
        if len(row.Atom1) > 1:
            J = df_noes.loc[:, "Atom1"].str[:-1] == row.Atom1[:-1]  # match part
        else:
            J = df_noes.loc[:, "Atom1"] == row.Atom1  # match full
        J = J & I
        if J.sum() == 0:
            # Error matching Atom1 & Residues
            print(row)
            print(df_noes[I])
            raise ValueError("No match for atom1")

        # 2) match atom2
        if len(row.Atom2) > 1:
            K = df_noes.loc[:, "Atom2"].str[:-1] == row.Atom2[:-1]  # match part
        else:
            K = df_noes.loc[:, "Atom2"] == row.Atom2  # match full
        K = K & I
        if K.sum() == 0:  # no match!
            # Error matching Atom2 & Residues
            print(row)
            print(df_noes[I])
            raise ValueError("No match for atom2")

        # 3) must also have same distance !
        L = df_noes.loc[:, "Distance"] == row.Distance
        if L.sum() == 0:  # no match!
            # Error matching Distance
            print(row)
            raise ValueError("No match for distance")

        # 4) assing index to matching lines
        L = J & K & L
        if not row.name in df_noes[L].index:  # assert current row in match !
            print(row)
            print(df_noes[L])
            raise ValueError("Current row not in match !?")

        # 5) check indexes
        # indexes can be not naive (-1), in case row was matched earlier;
        # we could check index, but this way we pick up possible double matches.
        idxs = df_noes.loc[L, "index"].values
        if np.any(idxs >= 0):  # indexes are not naive
            if len(set(idxs)) > 1:  # multiple indexes: something's wrong
                print(row)
                print(df_noes[L])
                raise ValueError("Multiple indexes in the matching group !?")
            if idxs[0] != row["index"]:  # not current index ?
                print(row)
                print(df_noes[L])
                raise ValueError("Index mismatch !?")
        if row["index"] >= 0:
            continue
        df_noes.loc[L, "index"] = cidx
        cidx += 1

    # Write out
    df_noes.to_csv('noes_index.csv', sep='\t')
    return df_noes


def Change_format(df_noes, itp_file):
    
    ''' It creates an .itp file available for GROMACS
    
    Input:
        
        - df_noes (pandas dataframe): NOE distances dataframe
        - itp_file (str): name of the .itp file that will be created
    
    Output:
    
        NONE
    
    '''
    
    fic = open(itp_file, 'w')
    fic.write('[ distance_restraints ]\n; ai\taj\ttype\tindex\ttype’\tlow\tup1\tup2\tfac\n')
    
    for row in df_noes.itertuples():
        low = str(round((float(row.Distance) - float(row.Distance)/10)/10, 2))
        up1 = str(round((float(row.Distance) + float(row.Distance)/10)/10, 2))
        up2 = str(round((float(row.Distance) + 2*(float(row.Distance)/10))/10, 2))
        
        fic.write(str(row.AtomID1)+'\t')
        fic.write(str(row.AtomID2)+'\t')
        fic.write('1\t')
        fic.write(str(row.index)+'\t')
        fic.write('1\t')
        fic.write(low+'\t')
        fic.write(up1+'\t')
        fic.write(up2+'\t')
        fic.write('1.00\t')
        fic.write(";\t")
        fic.write(row.ResType1+str(row.ResID1)+"."+row.Atom1+"-"
                  +row.ResType2+str(row.ResID2)+"."+row.Atom2+'\t')
        fic.write(row.Distance+'\t\n')

    fic.close()


def noe2itp(df_noes, itp_file= None):
    
    ''' This function controles NOE2ITP scrip. It calls Assing_index() and 
    Change_format() if required.
    
    Input
    
        - df_noes (pandas dataframe): NOE distances dataframe
        - itp_file (str; default: None): name of the .itp file that will be created
    
    Output:
        
        - df_noes (pandas dataframe) : NOE distances modified dataframe        
    
    '''
    
    df_noes = Assing_index(df_noes)
    
    if itp_file != None:
        Change_format(df_noes, itp_file)
    return df_noes


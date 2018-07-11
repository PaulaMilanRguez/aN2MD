#  File: noecombine.py
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "noecombine.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#
import pandas as pd
import mdtraj as md
import numpy as np


def Combining(noes_df1, noes_df2):
    ''' Merge two NOE distances dataframe in a single one. A column "Origin" is
    assinged: value 0 for the first dataframe, 1 for the second one.

    Input:
        - noes_df1 (pandas dataframe): NOE distances dataframe
        - noes_df2 (pandas dataframe): NOE distances dataframe

    Output:
        - df_noes (pandas dataframe): dataframe containing the whole set of NOE
        distances

    '''

    origin_df1 = np.repeat(0, len(noes_df1))
    noes_df1["Origin"] = origin_df1

    origin_df2 = np.repeat(1, len(noes_df2))
    noes_df2["Origin"] = origin_df2

    df_noes = [noes_df1, noes_df2]
    df_noes = pd.concat(df_noes)
    return df_noes


def Extend_noes(df_noes, gro_file):
    ''' Ambiguous NOEs (M* and Q*) indicate multiple protons and are made
    explicit.

    Input:

        - df_noes (pandas dataframe): NOE distances dataframe
        - gro_file (str): topology file

    Output:

        - df_noes (pandas dataframe) : NOE distances modified dataframe
    '''

    top = md.load_topology(gro_file)

    # First Atom ##

    todel = []
    toadd = []

    for i, row in df_noes.iterrows():

        if row.Atom1.startswith("M"):
            selection = "resid "+str(int(row.ResID1)-1) + \
                        " and (name =~ 'H"+row.Atom1[1:]+".*')"

            sel = top.select(selection)
            assert len(sel) == 3, str(row)+str(sel)
            todel.append(i)
            for atom in [top.atom(a) for a in sel]:
                newrow = row.copy()
                newrow.Atom1 = atom.name
                newrow.AtomID1 = atom.index
                newrow.Origin += 10
                toadd.append(newrow)

        elif row.Atom1.startswith("Q"):
            selection = "resid "+str(int(row.ResID1)-1) + \
                        " and (name =~ 'H"+row.Atom1[1:]+".*')"

            sel = top.select(selection)
            assert len(sel) == 2, str(row)+str(sel)
            todel.append(i)
            for atom in [top.atom(a) for a in sel]:
                newrow = row.copy()
                newrow.Atom1 = atom.name
                newrow.AtomID1 = atom.index
                newrow.Origin += 10
                toadd.append(newrow)

    df_noes = df_noes.drop(todel)
    df_noes = df_noes.append(toadd, ignore_index=True)

    # Second Atom ##

    todel = []
    toadd = []

    for i, row in df_noes.iterrows():
        if row.Atom2.startswith("M"):
            selection = "resid "+str(int(row.ResID2)-1) + \
                        " and (name =~ 'H"+row.Atom2[1:]+".*')"

            sel = top.select(selection)
            assert len(sel) == 3, str(row)+str(sel)
            todel.append(i)
            for atom in [top.atom(a) for a in sel]:
                newrow = row.copy()
                newrow.Atom2 = atom.name
                newrow.AtomID2 = atom.index
                newrow.Origin += 10
                toadd.append(newrow)

        elif row.Atom2.startswith("Q"):
            selection = "resid "+str(int(row.ResID2)-1) + \
                        " and (name =~ 'H"+row.Atom2[1:]+".*')"

            sel = top.select(selection)
            assert len(sel) == 2, str(row)+str(sel)
            todel.append(i)
            for atom in [top.atom(a) for a in sel]:
                newrow = row.copy()
                newrow.Atom2 = atom.name
                newrow.AtomID2 = atom.index
                newrow.Origin += 10
                toadd.append(newrow)

    df_noes = df_noes.drop(todel)
    df_noes = df_noes.append(toadd, ignore_index=True)

    df_noes = df_noes.sort_values(['ResID1', 'ResID2'])

    return df_noes


def Search_atom_index(df_noes, gro_file):
    ''' Assing the topology index to each atom

    Input:

       - df_noes (pandas dataframe): NOE distances dataframe
       - gro_file (str): topology file

    Output:

        - df_noes (pandas dataframe) : NOE distances modified dataframe
    '''

    top = md.load_topology(gro_file)
    AtomID1 = []
    AtomID2 = []
    for noe in df_noes.itertuples():
        ai = top.select("resid "+str(int(noe[1])-1)+" and(name "+noe[3]+")")
        aj = top.select("resid "+str(int(noe[4])-1)+" and(name "+noe[6]+")")
        if len(ai) == 0:
            ai = -1
        else:
            ai = ai[0]
        if len(aj) == 0:
            aj = -1
        else:
            aj = aj[0]

        AtomID1.append(ai)
        AtomID2.append(aj)

    df_noes['AtomID1'] = AtomID1
    df_noes['AtomID2'] = AtomID2
    return df_noes


def Remove_repetitions(df_noes):
    ''' Deduplicate constraints. Only the fist to appear in the dataframe is
    saved.

    Input:

        - df_noes (pandas dataframe): NOE distances dataframe

    Output:

        - df_noes (pandas dataframe) : NOE distances modified dataframe

    '''

    todel = []
    for i, row1 in df_noes.iterrows():

        for j, row2 in df_noes.iterrows():
            if (row1.AtomID1 == row2.AtomID1 and
                row1.AtomID2 == row2.AtomID2 and
                i < j):
                todel.append(j)

    df_noes = df_noes.drop(todel)
    df_noes = df_noes.sort_values(['AtomID1', 'AtomID2'])
    return df_noes


def noecombine(noes_df1, noes_df2, gro_file, deduplicate=True):
    ''' This function controles NOECOMBINE scrip. It calls Combine(), 
    Search_atom_index(), Remove_repetitions() and Extend_noes() if required.

    Input:

        - noes_df1 (pandas dataframe): NOE distances dataframe
        - noes_df2 (pandas dataframe): NOE distances dataframe
        - gro_file (str): topology file
        - deduplicate (bool; default: True): deduplicate constraints

    Output:

        - df_noes (pandas dataframe) : NOE distances modified dataframe

    '''

    df_noes = Combining(noes_df1, noes_df2)

    if deduplicate:
        df_noes = Extend_noes(df_noes, gro_file)

    df_noes = Search_atom_index(df_noes, gro_file)
    df_noes = Remove_repetitions(df_noes)

    return df_noes

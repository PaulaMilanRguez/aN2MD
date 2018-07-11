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

import numpy as np

"""
"""


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
        if row.name not in df_noes[L].index:  # assert current row in match !
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


def write_itp(df_noes, itp_file):
    ''' It creates an .itp file available for GROMACS

    Input:

        - df_noes (pandas dataframe): NOE distances dataframe
        - itp_file (str): name of the .itp file that will be created

    Output:

        NONE
    '''

    ITP_FORMAT = """\
{:<5d} {:<5d} {:4d} {:7d}    {:4d}  {:7.2f} {:7.2f} {:7.2f} {:7.2f}"""
    COMMENT_FORMAT = """\
      ; {:3s}{}.{:<4s}-{:3s}{}.{:<4s}  {:7.2f}  {:>6d} {:>3d} """

    fo = open(itp_file, "w")
    fo.write("[ distance_restraints ]\n")
    fo.write("; ai   aj   type   index   type’      low     up1     up2     fac\n")

    # Settings
    sorters = ["index", "ResID1", "ResID2"]
    dfac = 0.1

    # write a comment with original names for future reference
    for i, row in df_noes.sort_values(sorters).iterrows():
        # 10     16      1       0       1      0.0     0.3     0.4     1.0
        dist = float(row.Distance) / 10
        fo.write(
            (ITP_FORMAT + COMMENT_FORMAT + "\n").format(
                row.AtomID1+1, row.AtomID2+1, 1, row["index"], 1,
                (1.-dfac)*dist, (1.+dfac)*dist, (1.+2*dfac)*dist, 1.0,
                row.ResType1, row.ResID1, row.Atom1,
                row.ResType2, row.ResID2, row.Atom2,
                dist, i, row.Origin
            ))
    fo.close()


def noe2itp(df_noes, itp_file=None):
    ''' This function controles NOE2ITP scrip. It calls Assing_index() and 
    Change_format() if required.

    Input

        - df_noes (pandas dataframe): NOE distances dataframe
        - itp_file (str; default: None): name of the .itp file that will be created

    Output:

        - df_noes (pandas dataframe) : NOE distances modified dataframe
    '''

    df_noes = Assing_index(df_noes)

    if itp_file is not None:
        write_itp(df_noes, itp_file)
    return df_noes

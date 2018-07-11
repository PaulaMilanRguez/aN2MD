#  File: str2noe.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "str2noe.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#
import pandas as pd
import re


def parse_str(str_file, skip_HB, eliminate_res):
    ''' This function parse de .str file and selects the lines containing NOE 
    distances information. Then, it creates a dataframe with the NOE distances
    set.

    Input:
        - str_file(str): name of the .str file
        - skip_HB (bool; default = True): skip Hidrogen Bonds
        - eliminate_res (list of int): ensemble of residues that doesn't have
        to be taken into account for the dataframe.

    Output:

        - df_noes(pandas dataframe): dataframe containing the NOE distances
        found in the .str file


    '''

    fic = open(str_file, "r")
    found = False
    lines = []
    eliminate = ["O", "N", "S", "SG"]
    res_start = None

    for line in fic:
        res_stop = re.match(r'.*stop_.*', line)
        if res_start is not None:
            found = True
        elif res_stop is not None:
            found = False
        if found:
            if line != "\n":
                info = line.split()

                if int(info[7]) > int(info[17]):
                    r1 = info[17]
                    t1 = info[18]
                    a1 = info[19]

                    info[17] = info[7]
                    info[18] = info[8]
                    info[19] = info[9]

                    info[7] = r1
                    info[8] = t1
                    info[9] = a1

                if skip_HB:
                    if (info[9] not in eliminate and
                        info[19] not in eliminate and
                        int(info[7]) not in eliminate_res and
                        int(info[17]) not in eliminate_res):
                        lines.append(info)
                else:
                    if (int(info[7]) not in eliminate_res and
                        int(info[17]) not in eliminate_res):
                        lines.append(info)

        res_start = re.match(r'.*_Gen_dist_constraint.Gen_dist_constraint_list_ID.*', line)
    fic.close()

    df_noes = pd.DataFrame(lines)
    df_noes = df_noes[[7, 8, 9, 17, 18, 19, 28]]
    df_noes = df_noes.rename(columns={7: "ResID1", 8: "ResType1", 9: "Atom1",
                                      17: "ResID2", 18: "ResType2", 19: "Atom2",
                                      28: "Distance"})

    return df_noes


def naming(df_noes, res_names, delta_resid, atom_names):
    ''' This function serves to change the naming of the chosen residues or
    atoms and eventually its numbering

    Input:
        - df_noes (pandas dataframe): dataframe containing NOE distances set
        - res_names (list of lists): each inner list must contain the old name
        of the residue and the new one. For example:

            res_names = [['CYS','CY3'], ['HIS', 'HD1']]

        - delta_resid (int) = variation in the numbering wanted for each residue

        - atom_name (list of list): each inner list must contain the old atom
        name and the new one. If the naming change of the atom should be done
        in just one residue, the name of the residue should be included. For ex:

            atom_name = [['H', 'HN'], ['MG', 'MG2', 'THR']]

    Output:

        -df_noes (pandas dataframe): modified dataframe
    '''

    for residue in res_names:
        df_noes.loc[df_noes.ResType1 == residue[0], 'ResType1'] = residue[1]
        df_noes.loc[df_noes.ResType2 == residue[0], 'ResType2'] = residue[1]

    for atom in atom_names:
        if len(atom) == 2:
            df_noes.loc[df_noes.Atom1 == atom[0], 'Atom1'] = atom[1]
            df_noes.loc[df_noes.Atom2 == atom[0], 'Atom2'] = atom[1]
        elif len(atom) == 3:
            df_noes.loc[(df_noes.Atom1 == atom[0]) & (df_noes.ResType1 == atom[2]),'Atom1'] = atom[1]
            df_noes.loc[(df_noes.Atom2 == atom[0]) & (df_noes.ResType2 == atom[2]),'Atom2'] = atom[1]

    for row in df_noes.itertuples():
        df_noes.loc[df_noes.index == row.Index, 'ResID1'] = int(row.ResID1) + delta_resid
        df_noes.loc[df_noes.index == row.Index, 'ResID2'] = int(row.ResID2) + delta_resid

    return df_noes


def str2noe(str_file, skip_HB=True, eliminate_res=[], change_naming=False,
            res_names=[], delta_resid=0, atom_name=[]):

    ''' Read a NOE data from a STAR file, optionally renaming entries.

    str2noe returns a pandas.Dataframe containing renamed NOE data. It calls
    functions Parse_str() to obtain the Dataframe and Change_naming() if
    required.

    Input:
        - str_file (str): name of the .str file
        - skip_HB (bool; default = True): skip Hidrogen Bonds
        - eliminate_res (list of int): ensemble of residues to skip
        - change_naming (bool; default = False): if True, perform renaming
        - res_names (list of lists): each list must contain the old name
          of the residue and the new one. For example:

            res_names = [['CYS','CY3'], ['HIS', 'HD1']]

        - delta_resid (int) = variation in the numbering wanted for each residue
        - atom_name (list of lists): each list must contain the old atom name
          and the new one. If the naming change of the atom should be done in
          just one residue, the name of the residue should be included. For ex:

            atom_name = [['H', 'HN'], ['MG', 'MG2', 'THR']]


    Output:

        - df_noes(pandas dataframe): dataframe containing the NOE distances
        found in the .str file

    '''
    df_noes = parse_str(str_file, skip_HB, eliminate_res)

    if (change_naming):
        df_noes = naming(df_noes, res_names, delta_resid, atom_name)

    return df_noes

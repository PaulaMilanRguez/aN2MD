#  File: __init__.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "__init__.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#
__author__ = 'Paula Milan Rodriguez, Marco Pasi'
__version__ = 'v0.1'
__license__ = 'GPLv2'


from str2noe import str2noe as s2n
from noecombine import noecombine as nc
from noe2itp import noe2itp as n2i
from searchviolations import Search_Violations as sv
from numresults import numresults as numr
from violations_plot import Violations_plot as viplot


def Create_itp(str_file1, str_file2, gro_file, skip_HB=True, eliminate_res1=[], eliminate_res2=[], change_naming=False,
            res_names=[], delta_resid1=0, delta_resid2 =0, atom_name=[], deduplicate=True, itp_file=None):

    ''' Creation of an .itp file available for GROMACS from two STAR files

    Input:
        - str_file1 (str): name of the first STAR file
        - str_file2 (str): name of the second STAR file
        - skip_HB (bool; default = True): skip Hidrogen Bonds
        - eliminate_res1 (list of int): ensemble of residues to skip in the first STAR file
        - eliminate_res2 (list of int): ensemble of residues to skip in the second STAR file
        - change_naming (bool; default = False): if True, perform renaming
        - res_names (list of lists): each list must contain the old name
          of the residue and the new one. For example:

            res_names = [['CYS','CY3'], ['HIS', 'HD1']]

        - delta_resid1 (int) = variation in the numbering wanted for each residue in the first STAR file
        - delta_resid2 (int) = variation in the numbering wanted for each residue in the second STAR file
        - atom_name (list of lists): each list must contain the old atom name
          and the new one. If the naming change of the atom should be done in
          just one residue, the name of the residue should be included. For ex:

            atom_name = [['H', 'HN'], ['MG', 'MG2', 'THR']]

        - gro_file (str): topology file
        - deduplicate (bool; default: True): deduplicate constraints
        - itp_file (str; default: None): name of the .itp file that will be created


    Output:

        - df_noes(pandas dataframe): dataframe containing the NOE distances
        found in the .str file

    '''

    
    df_noe1 = s2n.str2noe(str_file1, skip_HB, eliminate_res1, change_naming,
            res_names, delta_resid1, atom_name)

    df_noe2 = s2n.str2noe(str_file2, skip_HB, eliminate_res2, change_naming,
            res_names, delta_resid2, atom_name)

    df_noes = nc.noecombine(noes_df1, noes_df2, gro_file, deduplicate=True)

    df_noes = n2i.noe2itp(df_noes, itp_file=None)

    return df_noes



def Violation_analysis(name_itp, name_gro, name_xtc= '', option = 0, rang = 30,
                      skip=200, strd=100, limit_55= True, ncp7= False, res_max=0):

    ''' Analyze NOE distances violations in a simulation or structure

    Input:

       - name_itp (str): name of the .itp file contatining the NOE distances
       - name_gro (str): topology file
       - name_xtc (str): .xtc simulation file
       - option (int): mean type: 0 for analysing a structure; 1 for simple mean in a simulation; 2 for an 
        r6 mean in a simulation; 3 for a r3 mean and a r6 addition in a simulation
       - rang (int): error range considered for NMR results
       - skip (int): number of frames that should be removed in the begining of the simulation 
       - strd (int): Only read every stride-th frame
       - Limit_55 (bool; default = True): just consider violations with a distance value > 5.5 A
       - ncp7 (bool; default = False): working specifically with HIV-1 Nucleocapsid Protein 
       - res_max (int): protein length 
        

    Output:

       - df_viol(pandas dataframe): dataframe containing the NOE distances violations
        found in the structure/simulation 

    '''

    ana_noes = sv.Search_Violations(name_itp, name_gro, name_xtc, option, rang,
                      skip, strd, name_work)


    df_viol = numr.numresults(ana_noes, limit_55)


    viplot.Violations_plot(ana_noes)    

    return df_viol





	

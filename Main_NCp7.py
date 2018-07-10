
## Import

import str2noe as s2n
import noecombine as nc
import noe2itp as n2i
import searchviolations as sv

############################# MAIN ######################################

## This Main is an example of how to use aN2MD program. NCp7 protein is used in this case.

		
		#### HOW TO GET AN .ITP FILE ####


# Using STR2NOE module to obtain two dataframes from two .str: some residues are removed and 
#the naming and numbering modified


noes_1esk = s2n.str2noe('1esk_mr.str', eliminate_res = [1,13,15], change_naming = True,
 res_names=[['CYS', 'CY3'], ['HIS','HD1']], delta_resid = 11, atom_name=[['H','HN', 'CY3'],['MG', 'MG2', 'THR']])

noes_1mfs = s2n.str2noe('1mfs_mr.str', change_naming = True,
 res_names=[['CYS', 'CY3'], ['HIS','HD1']], atom_name=[['H','HN', 'CY3']])


# Then, we combine the two dataframes using NOECOMBINE module:  repetitions and ambiguities are 
#removed and the atom's topology index assigned 

noes = nc.noecombine(noes_1esk, noes_1mfs, gro_file = '1esk_md.gro')

# Finally, we obtain an .itp file available for GROMACS using NOE2ITP module

noes = n2i.noe2itp(noes, itp_file= 'ncp7_noes.itp')
print(noes.head(5))


		#### VIOLATIONS ANALYSIS ####

# First, we search for the distances that corresponds to the NOEs pairs in the simulation.
# We compare them with the NMR distance to determinate if the NOE is respected or violated (SEARCHVIOLATIONS MODULE)

noes = sv.Search_violations('ncp7_noes.itp', 'pi.gro', 'pi100.2.xtc',option = 3)

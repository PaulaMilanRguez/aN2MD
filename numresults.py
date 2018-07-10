#  File: Final_analysis.py 
#
#  Copyright (C) 2018 Paula Milan Rodriguez, Marco Pasi
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "numresults.py v0.1 (C) 2018 Paula Milan Rodriguez, Marco Pasi"
#


import numpy as np
import pandas as pd


def Extract_viol(df_noes, title, limit_55):

    df_viol = df_noes[df_noes.Violation == 'Viol']

    if(limit_55):
        df_viol = df_viol[v.Distance > 5.5]

    pd.DataFrame.to_csv(df_viol, title)
    
    return df_viol

def numresults(df_noes, title= 'violations.csv', limit_55=True):

    df_viol = Extract_viol(df_noes, title, limit_55)

    return df_viol

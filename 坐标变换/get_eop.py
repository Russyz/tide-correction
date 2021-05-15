import numpy as np
import os
import pandas as pd
path = os.getcwd()
def get_eops():
    """
   This function downloads the Earth Orientation Parameters (EOPs) from the IAU sources and returns them as a pandas
        dataframe; https://datacenter.iers.org/eop.php
    """
    url = 'ftp://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now'
    ds = np.DataSource(path)
    file = ds.open(url)
    array = np.genfromtxt(file, skip_header=14)
    headers = ['Year', 'Month', 'Day', 'MJD', 'x', 'y', 'UT1-UTC', 'LOD', 'dX',
               'dY', 'x Err', 'y Err', 'UT1-UTC Err', 'LOD Err', 'dX Err', 'dY Err']
    eop = pd.DataFrame(data=array, index=array[:, 3], columns=headers)
    return eop
print(get_eops())
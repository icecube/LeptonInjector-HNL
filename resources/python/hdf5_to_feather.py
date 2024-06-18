import tables

import pandas as pd

import logging
logging.basicConfig(level=logging.INFO)

from tqdm import tqdm

def hdf5_to_feather(
    infiles,
    outfilepath='/data/ana/BSM/HNL/MC/feather_files/test_files/test.feather',
    keys=None,
    logtofile=False,
):
    """
    Tool to convert a list of hdf5 files into a single feather file and return pandas dataframe.

    Parameters
    ----------
    infiles : list of strings
        Infile paths.
    outfilepath : string
        Path of the outfile.
    keys : list of strings
        Keys to extract, will take all keys in the first file if not given.

    Returns
    -------
    dataframe : pandas dataframe
        The dataframe with the combined data from all files.
    """

#     set up dict with correct keys
    datadict = dict()
    
    if keys:
        for key in keys:
            datadict[key] = list()
    else:
        first_file = tables.open_file(infiles[0])

        for child in first_file.root:
            try:datadict[child.name] = list()
            except tables.NoSuchNodeError:continue

        first_file.close()

    logging.info('Keys to be extracted: {}'.format(datadict.keys()))   
    
#     extract data
    for filepath in tqdm(infiles):
        f = tables.open_file(filepath, 'r')
        try:
            for key in datadict.keys():
                datadict[key].extend( f.get_node('/'+key).col('value') )
        except:
            print('Something went wrong and it happened with the file: {}'.format(filepath))
        f.close()

#     add number of infiles used
    datadict['n_files'] = len(infiles)

#     create dataframe and store it
    dataframe = pd.DataFrame(datadict)
    dataframe.to_feather(outfilepath)
    
    return dataframe

import h5py
from dimelo.utils.genome_regions import Region
import pandas as pd
import time
import os
import numpy as np

def load_region(
    region: Region,
    basemod: str,
    outDir: str,
    sampleName: str,
):
    merged_filepath = f'{outDir}{sampleName}_hdf5.h5'
    with h5py.File(merged_filepath,'r') as merged_file:
        chr_data = merged_file['chr'][:]
        chr_data_decoded = np.array([chrom.decode('utf-8') for chrom in chr_data])
        chr_indices = np.where(chr_data_decoded == region.chromosome)[0]
        del chr_data
        del chr_data_decoded
        pos_data = merged_file['pos'][:]
        end_data = merged_file['end'][:]
        pos_data_filtered = pos_data[chr_indices]
        end_data_filtered = end_data[chr_indices]
        range_indices = np.where(
            (end_data_filtered >= region.begin) & 
            (pos_data_filtered < region.end)
        )[0]
        final_indices = chr_indices[range_indices]
        reads_pos = pos_data[final_indices]
        reads_end = end_data[final_indices]
        mod_arrays = merged_file[f'mod{basemod}'][final_indices]
        valid_arrays = merged_file[f'val{basemod}'][final_indices]
    return reads_pos,reads_end,mod_arrays,valid_arrays
        
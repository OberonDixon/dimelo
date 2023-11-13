import pyBigWig
from dimelo.utils.genome_regions import Region
import numpy as np
import time

def load_region(
    region: Region,
    basemod: str,
    outDir: str,
    sampleName: str,
):
    modified_filepath = f'{outDir}{sampleName}_bigwig_mod{basemod}.bw'
    valid_filepath = f'{outDir}{sampleName}_bigwig_valid{basemod}.bw'
    
    with pyBigWig.open(modified_filepath) as mod_bw:
        mod_array = np.nan_to_num(mod_bw.values(region.chromosome,region.begin,region.end))
    with pyBigWig.open(valid_filepath) as valid_bw:
        valid_array = np.nan_to_num(valid_bw.values(region.chromosome,region.begin,region.end))
    
    return mod_array,valid_array
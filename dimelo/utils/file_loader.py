import os
import importlib.util
from dimelo.utils.genome_regions import Region

SINGLE_READ_DICT = {}
GENOMIC_TRACK_DICT = {}

def load_plugins(target_dict,directory):
    for file_name in os.listdir(directory):
        if file_name.endswith('.py'):
            module_name = file_name[:-3] # Strip suffix
            
            # Load the module
            spec = importlib.util.spec_from_file_location(module_name, os.path.join(directory, file_name))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Check if the module has both required functions
            if hasattr(module, 'load_region'):
                target_dict[module_name] = {
                    'load_region':module.load_region, 
                }
            else:
                print(f"Warning: module {directory}/{file_name} loader lacks load_region function.")
    
def check_for_duplicate_keys(dict1,dict2):
    # Get set of keys from both dicts
    single_read_keys = set(dict1.keys())
    genomic_track_keys = set(dict2.keys())

    # Find intersection (common elements)
    common_keys = single_read_keys & genomic_track_keys

    if common_keys:
        print(f"Warning: Single read file_saver plugins should not share a name with genomic track file_saver plugins: {common_keys}")
        return True
    else:
        return False   

def region_basemod_load_track(
    region: Region,
    basemod: str,
    outDir: str,
    sampleName: str,
    file_format: str,
):
    if file_format in GENOMIC_TRACK_DICT:
        return GENOMIC_TRACK_DICT[file_format]['load_region'](
            region,
            basemod,
            outDir,
            sampleName,
        )
    else:
        print(f'Warning: plugins folders contain no loader module {format_key}.py')
        return None
def region_basemod_load_reads(
    region: Region,
    basemod: str,
    outDir: str,
    sampleName: str,
    file_format: str,
):   
    if file_format in SINGLE_READ_DICT:
        return SINGLE_READ_DICT[file_format]['load_region'](
            region,
            basemod,
            outDir,
            sampleName,
        )
    else:
        print(f'Warning: plugins folders contain no loader module {file_format}.py')
        return None

single_reads_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'plugins/single_read_loaders')
genomic_tracks_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'plugins/genomic_track_loaders')
load_plugins(SINGLE_READ_DICT,single_reads_path)
load_plugins(GENOMIC_TRACK_DICT,genomic_tracks_path)
check_for_duplicate_keys(SINGLE_READ_DICT,GENOMIC_TRACK_DICT)
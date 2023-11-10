import os
import importlib.util
from dimelo.utils.genome_regions import SubregionTask, SubregionData, ProcesswiseTaskBuilder

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
            if hasattr(module, 'save_subregion_batch') and hasattr(module, 'merge_subregion_batches'):
                target_dict[module_name] = {
                    'save_subregion_batch':module.save_subregion_batch, 
                    'merge_subregion_batches':module.merge_subregion_batches,
                }
            else:
                print(f"Warning: module {file_name} lacks one or both of save_subregion_batch and merge_subregion_batches.")
    
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

def subregion_save_all(
    subregion_data: SubregionData,
    outDir_temp: str,
):
    for format_key in subregion_data.formats_list:
        if format_key in SINGLE_READ_DICT:
            SINGLE_READ_DICT[format_key]['save_subregion_batch'](
                subregion_data,
                outDir_temp,
                format_key,
            )
        elif format_key in GENOMIC_TRACK_DICT:
            GENOMIC_TRACK_DICT[format_key]['save_subregion_batch'](
                subregion_data,
                outDir_temp,
                format_key,
            )
        else:
            print(f'Warning: plugins folders contain no module {format_key}.py')
            
        print(format_key,'saved successfully',flush=True)
        
def merge_temp_files(
    processwise_tasks: ProcesswiseTaskBuilder,
    sampleName,
    outDir: str,
    outDir_temp: str,
    basemods,
):
    for format_key in processwise_tasks.formats_list:
        if format_key in SINGLE_READ_DICT:
            SINGLE_READ_DICT[format_key]['merge_subregion_batches'](
                processwise_tasks,
                sampleName,
                outDir,
                outDir_temp,
                format_key,
                basemods,
            )
        elif format_key in GENOMIC_TRACK_DICT:
            GENOMIC_TRACK_DICT[format_key]['merge_subregion_batches'](
                processwise_tasks,
                sampleName,
                outDir,
                outDir_temp,
                format_key,
                basemods,
            )  
        else:
            print(f'Warning: plugins folders contain no module {format_key}.py')
    
single_reads_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'plugins/single_read_savers')
genomic_tracks_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'plugins/genomic_track_savers')
load_plugins(SINGLE_READ_DICT,single_reads_path)
load_plugins(GENOMIC_TRACK_DICT,genomic_tracks_path)
check_for_duplicate_keys(SINGLE_READ_DICT,GENOMIC_TRACK_DICT)
import h5py
import pandas as pd
import time
import os
import numpy as np

def save_subregion_batch(
    subregion_data,
    outDir_temp,
    format_key,
):
#     print(subregion_data.read_dicts_list[0])
#     print(f'start writing batch {subregion_data.string} at {time.time()}')
    
    
    filepath = f"{outDir_temp}/{subregion_data.string}_temp_{format_key}.h5"
    if os.path.exists(filepath):
        os.remove(filepath)
    
    with h5py.File(filepath,'w') as hf:
        for key,value in subregion_data.read_lists_dict.items():
            if len(value)>0 and key in ['read_name','chr','pos','is_forward','ref_seq','valmA','modmA']:
                if isinstance(value[0],str):
                    dt = h5py.string_dtype(encoding='utf-8')
                    hf.create_dataset(key,data=value,dtype=dt)
                elif isinstance(value[0],np.ndarray):
                    if type(value[0][0]) in [float,np.float64]:
                        hf.create_dataset(key,
                                          data=value,
                                          dtype=h5py.vlen_dtype(np.float64))
                    elif type(value[0][0]) in [int,np.int64]:
                        hf.create_dataset(key,
                                          data=value,
                                          dtype=h5py.vlen_dtype(np.int64))                        
                else:
                    hf.create_dataset(key,data=value)
                    
            else:
                hf.create_dataset(key,data=(0,))
            
def merge_subregion_batches(
    processwise_tasks,
    sampleName,
    outDir,
    outDir_temp,
    format_key,
    basemods,
):
    merged_filepath = f'{outDir}/{sampleName}_{format_key}.h5'
    if os.path.exists(merged_filepath):
        os.remove(merged_filepath)
    with h5py.File(merged_filepath,'w') as target_file:
        for single_process_tasks in processwise_tasks.core_assignments.values():
            for subregion_tasks_list in single_process_tasks.tasks.values():
                for subregion_task in subregion_tasks_list:
                    source_filepath = f"{outDir_temp}/{subregion_task.string}_temp_{format_key}.h5"
                    with h5py.File(source_filepath,'r') as source_file:
                        for key in source_file.keys():
                            if key in target_file:
                                target_file[key].resize(
                                    (target_file[key].shape[0] + 
                                     source_file[key].shape[0],))
                                target_file[key][-source_file[key].shape[0]:] = source_file[key][:]
                            else:
                                target_file.create_dataset(
                                    key,
                                    data=source_file[key][:],
                                    maxshape=(None,)
                                )
        
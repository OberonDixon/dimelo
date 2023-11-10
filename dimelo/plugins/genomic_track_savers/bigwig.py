import pyBigWig
from dimelo.utils.file_saver import SubregionData,ProcesswiseTaskBuilder
import numpy as np
import time
from multiprocessing import Pool

def save_subregion_batch(
    subregion_data: SubregionData,
    outDir_temp: str,
    format_key: str,
):
    temp_file_prefix = f'{subregion_data.string}_temp_{format_key}'
    
    chr_sizes = {
        subregion_data.chromosome:subregion_data.end
    }
    
    rd_bw = pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_read_depth.bw","w")
    
    rd_bw.addHeader([(k, v) for k, v in chr_sizes.items()])
    
    rd_bw.addEntries(
        np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
        subregion_data.pile_coordinates,
        subregion_data.pile_coordinates+1,
        subregion_data.read_depth
    )
    
    rd_bw.close()
#     print('saved read depth',flush=True)
    
    for basemod in subregion_data.basemods:
        start_mod = time.time()
#         print(basemod,flush=True)
        with pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_mod{basemod}.bw","w") as mod_bw:        
            mod_bw.addHeader([(k, v) for k, v in chr_sizes.items()])        
            mod_bw.addEntries(
                np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
                subregion_data.pile_coordinates,
                subregion_data.pile_coordinates+1,   
                subregion_data.modified_pile_dict[basemod],
            )        
#         print('saved mod',flush=True)
        
        with pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_valid{basemod}.bw","w") as valid_bw:
            valid_bw.addHeader([(k, v) for k, v in chr_sizes.items()])       
            valid_bw.addEntries(
                np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
                subregion_data.pile_coordinates,
                subregion_data.pile_coordinates+1,   
                subregion_data.valid_pile_dict[basemod],
            )       
#         print('saved valid',flush=True)
#         print('no really though, I actually saved it')
        with pyBigWig.open(f"{outDir_temp}/{temp_file_prefix}_mod-over-valid{basemod}.bw","w") as mod_over_valid_bw:
#             print('created mod over valid bw',flush=True)
            mod_over_valid_bw.addHeader([(k, v) for k, v in chr_sizes.items()])
#             print('sized mod over valid bw',flush=True)
            mask = subregion_data.valid_pile_dict[basemod] != 0
            ratio = np.zeros_like(subregion_data.valid_pile_dict[basemod])
            ratio[mask] = subregion_data.modified_pile_dict[basemod][mask] / subregion_data.valid_pile_dict[basemod][mask]
#             print('ratio and mask lens',len(ratio),len(mask))
#             print('ratio and mask mean',np.mean(ratio),np.mean(mask))
#             print('pile coordinates min/max',np.min(subregion_data.pile_coordinates),np.max(subregion_data.pile_coordinates))
            mod_over_valid_bw.addEntries(
                np.repeat(subregion_data.chromosome,len(subregion_data.pile_coordinates)),
                subregion_data.pile_coordinates,
                subregion_data.pile_coordinates+1,   
                ratio,                      
            )        
#         print('saved mod over valid',flush=True)
           
def merge_subregion_batches(
    processwise_tasks: ProcesswiseTaskBuilder,
    sampleName: str,
    outDir: str,
    outDir_temp: str,
    format_key: str,
    basemods: list,
):
    merged_file_path = f'{outDir}{sampleName}_{format_key}_read_depth.bw'
#     print(merged_file_path,flush=True)
    with pyBigWig.open(merged_file_path,"w") as target_bw:
#         print([(k, v) for k, v in processwise_tasks.chr_ranges.items()])
        target_bw.addHeader([(k, v) for k, v in processwise_tasks.chr_ranges.items()])

        for single_process_tasks in processwise_tasks.core_assignments.values():
            for subregion_tasks_list in single_process_tasks.tasks.values():
                for subregion_task in subregion_tasks_list:
                    temp_file_path = f'{outDir_temp}/{subregion_task.string}_temp_{format_key}_read_depth.bw'
#                     print(temp_file_path,flush=True)
                    with pyBigWig.open(temp_file_path) as source_bw:
                        chrom = subregion_task.chromosome
                        
                        starts = np.arange(subregion_task.begin,subregion_task.end)
                        ends = starts+1
                        values = source_bw.values(chrom,subregion_task.begin,subregion_task.end) 

                        target_bw.addEntries(
                            np.repeat(chrom,len(starts)),
                            starts,
                            ends,
                            values,
                        )
        
    for basemod in basemods:
        for track_type in ['mod','valid','mod-over-valid']:
            merged_file_path = f'{outDir}{sampleName}_{format_key}_{track_type}{basemod}.bw'
#             print(merged_file_path,flush=True)
            with pyBigWig.open(merged_file_path,"w") as target_bw:
                target_bw.addHeader([(k, v) for k, v in processwise_tasks.chr_ranges.items()])
                for single_process_tasks in processwise_tasks.core_assignments.values():
                    for subregion_tasks_list in single_process_tasks.tasks.values():
                        for subregion_task in subregion_tasks_list:
                            temp_file_path = f'{outDir_temp}/{subregion_task.string}_temp_{format_key}_{track_type}{basemod}.bw'
#                             print(temp_file_path,flush=True)
                            with pyBigWig.open(temp_file_path) as source_bw:
                                chrom = subregion_task.chromosome

                                starts = np.arange(subregion_task.begin,subregion_task.end)
                                ends = starts+1
                                values = source_bw.values(chrom,subregion_task.begin,subregion_task.end)

                                target_bw.addEntries(
                                    np.repeat(chrom,len(starts)),
                                    starts,
                                    ends,
                                    values,
                                )

def create_bw(
    path:str,
    chr_sizes: dict,
    begin,
    end,
    values,
):
    pass
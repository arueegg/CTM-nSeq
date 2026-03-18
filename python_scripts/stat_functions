import os
import shutil
import pandas as pd
import numpy as np 
import subprocess

import helper_functions as hf
import main_functions as mf
import wrapper_functions as wf
    

def create_ROI_bed(input_bed,\
                   output_ROI_bed = True, \
                   output_bed="",\
                   return_unmodified = False,\
                   return_df = True):
    """ takes a .bed file containing the columns "chr", "start", "end", "name" and creates regions of interest (ROIs)
        for each "name" based on the minimum and maximum values found in "start" and "end", respectively.
    """
    bed_df = pd.read_csv(input_bed, sep = "\t", names = ["chr", "start", "end", "name"])
    if return_unmodified:
        return(bed_df)
    ROI_df = bed_df[bed_df['name'].str.contains("Cas")]
    minvalues = ROI_df.groupby("chr").min("start")[["start"]]
    maxvalues = ROI_df.groupby("chr").max("end")[["end"]]
    ROIs_bed = pd.merge(minvalues,maxvalues, on = "chr", how = "outer")
    ROIs_bed[["name"]] = "ROI"
    if output_ROI_bed:
        ROIs_bed.to_csv(output_bed, sep="\t", header = False, index = True)
    
    if return_df:
        ROIs_bed = ROIs_bed.reset_index()
        return(ROIs_bed)

def Sam_Stats(input_folder_name = "aligned", \
              searched_input_dir = "output", \
              extension = ".bam",\
              out_dir = "Sam_Stats",\
              extract_lengths = False,\
              generate_summary= True,\
              pre_calculated_stats = False):
    """

    """
    if extract_lengths:
        length_dir = os.path.join(out_dir, "Lengths")
        merged_length_dir = os.path.join(out_dir, "Merged_Lengths")
        hf.create_dir(length_dir)
        hf.create_dir(merged_length_dir)
        full_length_df = pd.DataFrame(columns=["size", "count",  "dataset", "barcode"])
    if generate_summary:
        summary_dir = os.path.join(out_dir, "Stat_Summaries")
        hf.create_dir(summary_dir)
        full_summary_df = pd.DataFrame(columns=["name","value", "comment","dataset", "barcode"])

    
    hf.create_dir(out_dir)
    hf.create_dir(os.path.join(out_dir, "SamStats"))
    all_bamfiles = hf.loop_output(searchterm = input_folder_name,\
                                  searched_input_dir = searched_input_dir, \
                                  file_extension = extension)
    
    for index, row in all_bamfiles.iterrows():
        name = row['name']
        path = row['path']
        dataset = "_".join([name.split("_")[0], name.split("_")[1]])
        stats_file = os.path.join(out_dir, "SamStats",name + "_stats.txt")
        function = "samtools stats"
        full_call = full_call = " ".join([function, path, ">", stats_file])
        if not pre_calculated_stats:
            subprocess.call(full_call, shell=True)
        if extract_lengths:
            try:
                bc="".join(["barcode", name.split("_barcode")[1][0:2]])
            except:
                bc="none"
            lengths_file = os.path.join(length_dir,name + "_lenghts.txt")
            lengths_call = " ".join(["grep ^RL", stats_file, "| cut -f 2-  >", lengths_file])
            subprocess.call(lengths_call, shell=True)
            length_df = pd.read_csv(lengths_file, sep = "\t", names = ["size", "count"])
            length_df['dataset']=dataset
            length_df['barcode'] = bc
            full_length_df= pd.concat([full_length_df, length_df], axis = 0, ignore_index=True)
        if extract_lengths:
            full_length_df.to_csv(os.path.join(merged_length_dir, dataset + "_full_df.txt"))
        if generate_summary:
            summary_file = os.path.join(summary_dir,name + "_summary.txt")
            summary_call = " ".join(["grep ^SN", stats_file, "| cut -f 2-  >", summary_file])
            subprocess.call(summary_call, shell=True)
            summary_df = pd.read_csv(summary_file, sep = "\t", names = ["name","value", "comment"],index_col=False)
            summary_df['dataset']=dataset
            summary_df['barcode'] = bc
            full_summary_df= pd.concat([full_summary_df, summary_df], axis = 0, ignore_index=True)
            full_summary_df.to_csv(os.path.join(summary_dir, dataset + "_merged_stat_sums.txt"))
        

def Depth_stats(ROI_bedfile,\
                input_folder_name = "aligned", \
                searched_input_dir = "output", \
                extension = ".bam",\
                out_dir = "Depth_Stats",
               combined_file_effix = ""):
    hf.create_dir(out_dir)
    hf.create_dir(os.path.join(out_dir, "TargetDepths"))
    hf.create_dir(os.path.join(out_dir, "TargetCoverage"))
    hf.create_dir(os.path.join(out_dir, "GenomeCoverage"))
    merged_depth_stats = os.path.join(out_dir, "Merged_Depths")
    hf.create_dir(merged_depth_stats)
    all_bamfiles = hf.loop_output(searchterm = input_folder_name,\
                                  searched_input_dir = searched_input_dir, \
                                  file_extension = extension)
    full_depth_df = pd.DataFrame(columns=["chr", "pos", "depth", "dataset", "barcode"])
    for index, row in all_bamfiles.iterrows():
        name = row['name']
        path = row['path']
        dataset = "_".join([name.split("_")[0], name.split("_")[1]])
        depth_file = os.path.join(out_dir, "TargetDepths", name + "_depths.txt")
        full_call = " ".join(["samtools depth -a -b", ROI_bedfile, "-o", depth_file, path])
        subprocess.call(full_call, shell=True)
        depth_df = pd.read_csv(depth_file, sep = "\t", names = ["chr", "pos", "depth"])
        depth_df['dataset']=dataset
        try:
            bc="".join(["barcode", name.split("_barcode")[1][0:2]])
        except:
            bc="none"
        depth_df['barcode'] = bc
        full_depth_df= pd.concat([full_depth_df, depth_df], axis = 0, ignore_index=True)
        cov_stats = os.path.join(out_dir, "TargetCoverage", name +  "_coverage.txt")
        ROIs_bed = pd.read_csv(ROI_bedfile, sep = "\t",  names = ["chr", "start", "end", "name"])
        for index, row in ROIs_bed.iterrows():
            region = row["chr"]+":"+str(row["start"])+"-"+str(row["end"])
            full_call = " ".join(["samtools coverage -r", region, path,  "| tee -a", cov_stats])
            subprocess.call(full_call, shell=True)
        genome_cov = os.path.join(out_dir, "GenomeCoverage", name + "_genome_cov.txt")
        full_call = " ".join(["samtools coverage", path,  ">", genome_cov])
        subprocess.call(full_call, shell=True)
    if combined_file_effix =="":
        from datetime import datetime
        combined_file_effix =datetime.now().strftime("%Y%m%d%H%M")+"_"
    full_depth_df.to_csv(os.path.join(merged_depth_stats, combined_file_effix+"_"+name+"_combined_depths.txt"), sep="\t", header = True, index = False)



    

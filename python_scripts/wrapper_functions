import os
import shutil
import pandas as pd
import numpy as np 
import subprocess
import time

import main_functions as mf
import helper_functions as hf
import stat_functions as sf

def core_analysis(input_df,\
                  name_col, \
                  ref_col, \
                  barcoding_col,\
                  input_col, \
                  out_dir_col, \
                  main_out_dir = "output", \
                  suffix="", \
                  overwrite =True,\
                  ###
                  basecall = True, \
                  basecall_dirname = "basecalled", \
                  quality = "hac", \
                  trim = "--no-trim", \
                  methylation = "", \
                  ###
                  filter_reads = True, \
                  filter_dirname = "filtered", \
                  seqlen =500, \
                  ###
                  demux = True, \
                  demux_dirname = "demuxed", \
                  dual_bc_only = False,\
                  demux_kit = 'HelperFiles/demux_files.txt', \
                  ###
                  align = True, \
                  align_dirname = "aligned", \
                  reference_locations = 'HelperFiles/references.txt',\
                  mmi_ref=""):
    """
    ***core_analysis*** 
    wrapper function to loop over input_df and feeding into all main functions
    """  
    error_filepath = os.path.join(main_out_dir, "_".join([time.strftime("%Y%m%d%H%M"), "run_log.txt"]))
    for index, row in input_df.iterrows():
        data = pd.DataFrame(row.to_frame().T)   
        reference = row[ref_col]
        if len(mmi_ref)==0:
            mmi_ref = reference
        demux_df = row[barcoding_col]
        name_base = row[name_col]
        #folder specification
        out_dir = os.path.join(main_out_dir, row[out_dir_col])
        input_folder = row[input_col]
        basecalled = os.path.join(out_dir + suffix, basecall_dirname)
        filtered_folder = os.path.join(out_dir + suffix,  filter_dirname)
        demuxed_folder = os.path.join(out_dir + suffix, demux_dirname)
        aligned_folder = os.path.join(out_dir + suffix, align_dirname)

        #generate an error file that collects potential errors
        
        
        next_input = basecalled
        
        if basecall:
            if os.path.exists(input_folder):
                try:
                    mf.basecall(name_base = name_base,\
                                input_folder = input_folder, \
                                quality =  quality, \
                                output_folder = basecalled, \
                                overwrite = overwrite, \
                                methylation = methylation)
                except subprocess.CalledProcessError as e:
                    print("Error occurred:", e)
                    hf.write_error_note(name_base, "basecalling",error_filepath, e)
                    continue
            else:
                hf.write_error_note(name_base, "basecalling",error_filepath, "- input not found")
                continue
        if filter_reads:
            if os.path.exists(next_input):
                try:
                    mf.filter_reads(name_base =name_base,\
                                    input_folder = next_input, \
                                    output_folder = filtered_folder, \
                                    overwrite =overwrite,\
                                    seqlen =seqlen)
                    next_input = filtered_folder
                except subprocess.CalledProcessError as e:
                    print("Error occurred:", e)
                    hf.write_error_note(name_base, "filtering",error_filepath, e)
                    continue
            else:
                hf.write_error_note(name_base, "filtering",error_filepath, "- input not found")
                continue
        if demux and demux_df:
            if os.path.exists(next_input):
                try:
                    mf.demux(name_base =name_base, \
                         input_folder = next_input, \
                         output_folder = demuxed_folder, \
                         overwrite =overwrite,\
                         dual_bc_only = dual_bc_only,\
                         demux_kit = demux_kit)
                    next_input = demuxed_folder
                except subprocess.CalledProcessError as e:
                    print("Error occurred:", e)
                    hf.write_error_note(name_base, "demultiplexing",error_filepath, e)
                    continue
            else:
                hf.write_error_note(name_base, "demultiplexing",error_filepath, "- input not found")
                continue
        if align:
            if os.path.exists(next_input):
                try:
                    mf.align(name_base =name_base,\
                         input_folder = next_input, \
                         output_folder = aligned_folder, \
                         overwrite =overwrite,\
                         reference_locations = reference_locations,\
                         mmi_ref=mmi_ref)
                except subprocess.CalledProcessError as e:
                    print("Error occurred:", e)
                    hf.write_error_note(name_base, "aligning",error_filepath, e)
                    continue
            else:
                hf.write_error_note(name_base, "aligning",error_filepath, "- input not found")
                continue
        out_message1 = " ".join([time.strftime("%Y%m%d%H%M"),":", name_base, "finished:"])
        bc ='basecalling' if basecall else 'no basecalling'
        filt = "".join(['filtering'," >", str(seqlen)," bp"]) if filter_reads else 'no filter'
        dem = 'demultiplexing' if demux else 'no demux '
        al = "".join(['aligned to', str(mmi_ref)]) if align else 'not aligned'
        out_message2 = ", ".join([bc, filt,dem,al])
        with open(error_filepath, "a") as f:
            f.write("".join([out_message1, " ", out_message2]))

#
def Stats_DFin(Input_DF,\
               Bedfile_All,\
               createROIbed, \
               out_dir_col,\
               input_folder_name, \
               main_out_dir, \
               extension,\
               out_dir_Depth,\
               out_dir_SamStats,\
               combined_file_effix = "",\
               extract_lengths = False,\
               generate_summary= True,\
               pre_calculated_stats = False):
               """
               Input_DF = dataframe with input file locations
               out_dir_col = column in Input_DF with output locations
               Bedfile_All = bed file with annotatioins and regions of interest (name column = "ROI")
               createROIbed = create a new bed file with only regions of interest based on Bedfile_All (true/false) 
               input_folder_name = name of folder with aligned files (subfolder of main_out_dir/dataset/)
               main_out_dir = name of main output folder containing aligned files
               extension = file extension to screen input_folder_name for
               out_dir_Depth = output folder for depth statistics
               out_dir_SamStats = output folder for genereal statistics
               combined_file_effix = effix for all files if left "", will be a timestamp
               extract_lengths = should length distribution statistics be collected for full dataset
               generate_summary= should summary files (for demuxed files) be generated
               pre_calculated_stats = use precalculated stats and only calculate summaries
               """
          
    if createROIbed:
        
        sf.create_ROI_bed(input_bed = Bedfile_All, output_bed="HelperFiles/all_ROI.bed")
        roibedfile = "HelperFiles/all_ROI.bed"
    else:
        roibedfile= Bedfile_All
    
    for index, row in Input_DF.iterrows():
        out_dir = os.path.join(main_out_dir, row[out_dir_col])
        print(input_folder_name, out_dir, extension, out_dir_SamStats)
        sf.Sam_Stats(input_folder_name, out_dir, extension, out_dir_SamStats , extract_lengths = True, generate_summary = True, pre_calculated_stats = False)
        sf.Depth_stats(ROI_bedfile = roibedfile,\
                input_folder_name = input_folder_name, \
                searched_input_dir = out_dir, \
                extension = extension,\
                out_dir = out_dir_Depth,\
               combined_file_effix = combined_file_effix)
        

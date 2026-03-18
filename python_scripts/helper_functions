import os
import shutil
import pandas as pd
import numpy as np 
import subprocess
import time

#import stat_functions as sf
import main_functions as mf
import wrapper_functions as wf
    

def create_dir(path):
    """creats a dir at path if it does not exist"""
    if not os.path.exists(path):
        os.makedirs(os.path.join(path))    

def empty_dir(folder):
    """removes all files and directories from folder"""
    for filename in os.listdir(folder):            
        filepath = os.path.join(folder, filename)
        if os.access(filepath, os.W_OK):
            try:
                shutil.rmtree(os.path.join(folder, filename), ignore_errors=True)
            except:
                print("cannot remove: ", os.path.join(folder, filename))

def bam_to_maindir(inputdir, name_base, outputdir, searchterm = "barcode", split = "_", input_alignment = ""):
    """
    ***bam_to_maindir*** restructures sub-folder structures generated through dorado tools.
    inputdir - directory with subfolders
    name_base - name to be used as prefix for all filenames
    outputdir - directory to output files found in subdirectories of inputdir
    searchterm - term to search for in subfolders and filenames to retain in all files
    split - separation characters found in original folders and to be added to new file names
    input_alignment - filepath, if set the filename found in filepath will be used to set the new file names
    """
    old_dirs = []
    new_bam_paths = []
    for dirpath, dirnames, filenames in os.walk(inputdir):
        old_dirs.append(dirpath)
        for dirname in dirnames:
            for bamfile in os.listdir(os.path.join(dirpath, dirname)):
                if bamfile.endswith('.bam') or bamfile.endswith(".bam.bai"):
                    original_path = os.path.join(dirpath, dirname, bamfile)
                    if len(input_alignment)>0:
                        if searchterm in input_alignment:
                            bc = split+searchterm+os.path.splitext(os.path.splitext(input_alignment.split(split+searchterm)[1])[0])[0]
                        else:
                            bc = "unclassified"
                    elif searchterm in dirname:
                        bc = split + dirname 
                    elif "unclassified" in dirname:
                        bc = split+ dirname 
                    else:
                        bc = ""
                    new_path = os.path.join(outputdir, name_base+ bc+".bam")
                    new_bam_paths.append(new_path)
                    if bamfile.endswith(".bam.bai"):
                        new_path = os.path.join(outputdir, name_base+ bc+".bam.bai")
                    os.rename(original_path, new_path)
    old_dirs = old_dirs[1:-1]
    for folder in old_dirs:
        shutil.rmtree(folder, ignore_errors=True)
    return(new_bam_paths)

def write_error_note(name_base, function ,error_filepath, error_message=""):
    errornote = " ".join([time.strftime("%Y%m%d%H%M%S"), ":", name_base, "had an error with", function, "\n", error_message])
    print(errornote)
    with open(error_filepath, "a") as f:
        f.write(errornote)

def loop_output(searchterm = "aligned", \
                searched_input_dir = "output", \
                file_extension = ".bam"):
    """
    output_stat - one of "ROI_depth", "genome_cov", "Sam_Stats"
    extract_readlength - can only be run in line with or if Sam_Stats was run previously
    ROI_bedfile file (bed format) with regions of interest. needed for "ROI_depth"
    """
    all_dirs = [x[0] for x in os.walk(searched_input_dir)]
    input_dirs = [x for x in all_dirs if searchterm in x]
    df_all = pd.DataFrame(data={'name': [], 'path': []})
    for folder in input_dirs:
        dataset = folder.split("/")[1]
        for bamfile in os.listdir(folder):
            if bamfile.endswith(file_extension):
                core_name = os.path.splitext(bamfile)[0]
                bam_path = os.path.join(folder,bamfile)
                df2 = pd.DataFrame(data={'name': [core_name], 'path': [bam_path]})
                df_all = pd.concat([df_all,df2])
    df_all = df_all.reset_index(drop = True)
    return(df_all)

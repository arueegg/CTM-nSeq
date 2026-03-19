import pandas as pd
import stat_functions as sf
import helper_functions as hf
import subprocess
import os

def methylation_analysis(dataset_name,\
                         aligned_folder,\
                         screen_for_foldername,\
                         screeened_outputfolder,\
                         barcode_list,\
                         reference_files,\
                         ref_name,\
                         main_output_folder,\
                         pile_up_folder,\
                         outfolder_full_tables,\
                         regions_bed,\
                         ROI_subset_filter,\
                         unmodified_ROI_bed = False,\
                         screened_extension = ".bam",\
                         combined_file_effix = ""
                         ):
                         """
                         dataset_name = name of dataset (bam files) to screen and select in aligned folder
                         aligned_folder = path of folder with aligned bam files
                         screen_for_foldername = base name of aligned_folder 
                         screeened_outputfolder = base name of  main output folder
                         
                         
                         barcode_list = list with barcodes used
                         reference_files =  a csv file containing the file locations for reference files, following columns:
                                 -name:   name of the reference - referenced through mmi_ref to find corresponding filelocations
                                 -fasta:   location of a fasta (.fa) file for the corresponding reference
                                 -fai:   location of a fasta index (.fa.fai) file for the corresponding reference
                                 -mmi:   location of a pre-computed minimap2 index file (.mmi) for the corresponding reference
                         ref_name = reference to select based on name in reference_files
                         main_output_folder = main folder for file output
                         pile_up_folder = subfolder to save pileup files
                         outfolder_full_tables = subfolder save full methylation tables
                         regions_bed = bed file with specified regions of interest
                         ROI_subset_filter = filter to apply (based on chr) if unmodified_ROI_bed == true
                         unmodified_ROI_bed = should regions_bed be screened or left unmodified (T/F)
                         screened_extension = file extension to screen for (.bam)
                         combined_file_effix = effix for all files if left "", will be a timestamp
                         """
    hf.create_dir(pile_up_folder)
    hf.create_dir(outfolder_full_tables)
    hf.create_dir(main_output_folder)
    if combined_file_effix =="":
        from datetime import datetime
        combined_file_effix =datetime.now().strftime("%Y%m%d_%H%M")+"_"
    ROIs = sf.create_ROI_bed(input_bed = regions_bed, output_bed="HelperFiles/all_ROI.bed", return_unmodified = unmodified_ROI_bed, return_df =True)
    all_bam = hf.loop_output(searchterm = screen_for_foldername, \
                searched_input_dir = screeened_outputfolder, \
                file_extension = screened_extension)
    indx=[]
    for i in range(len(all_bam)):
        for bc in barcode_list:
            if dataset_name in all_bam.name[i]:
                if bc in all_bam.name[i]:
                    indx.append(i)
    bam_test = all_bam.loc[indx]["path"]
    region_list = []
    names_list = []
    if ROI_subset_filter != "":
        region = "".join([ROIs.query("chr == '"+ROI_subset_filter+"'")["chr"].values[0], ":",\
                       str(ROIs.query("chr == '"+ROI_subset_filter+"'")["start"].values[0]), "-",\
                       str(ROIs.query("chr == '"+ROI_subset_filter+"'")["end"].values[0])])
        region_list.append(region)
        names_list = [ROI_subset_filter]
    else:
        for i in range(0, len(ROIs)):
            region = "".join([ROIs.iloc[i]["chr"], ":",\
                           str(ROIs.iloc[i]["start"]), "-",\
                           str(ROIs.iloc[i]["end"])])
            region_list.append(region)
            names_list.append(ROIs.iloc[i]["name"]+"_"+ROIs.iloc[i]["chr"])
    #print(region_list,  names_list)
    ref_df = pd.read_csv(reference_files, sep = ",")
    ref = ref_df.query("name == '"+ref_name+"'")["fasta"].values[0]
    threads=32
    
    for i in range(0, len(region_list)):
        region_modifier = "_"+names_list[i]+"_"
        file_name = main_output_folder+combined_file_effix+region_modifier+"CG_motifs.bed" 
        subprocess.call(" ".join(["modkit motif bed",ref, "CG 0 >",file_name]), shell=True)
        pileups = pd.DataFrame(columns=["pileup", "bc"])
        for file in bam_test:
            print(file)
            sample = file
            #create pileup tables
            pile_up = file.replace(".bam",region_modifier+"pileup.bed")
            pile_up = pile_up.replace(aligned_folder, pile_up_folder+region_modifier+combined_file_effix)
            function = "modkit pileup"
            call = " ".join([function, sample, pile_up, "--cpg --combine-strands","--region",region_list[i],"--ref", ref])
            subprocess.call(call, shell=True)
            if barcode_list ==["unclassified"]:
                bc = "unclassified"
            else:
                bc = "barcode"+pile_up.split("_barcode")[1].split("_")[0]    
            pileups.loc[len(pileups)] = [pile_up,bc]

            fullout = file.replace(".bam",region_modifier+"_full.tsv")
            fullout = fullout.replace(aligned_folder, outfolder_full_tables+region_modifier+combined_file_effix)
            function = "modkit extract full"
            full_call = " ".join([function, sample, fullout, "--region",region_list[i],"--ref", ref])
            subprocess.call(full_call, shell=True)

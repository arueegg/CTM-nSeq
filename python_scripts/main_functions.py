import os
import shutil
import pandas as pd
import numpy as np 
import subprocess

import helper_functions as hf

def basecall(name_base,\
             input_folder, \
             output_folder, \
             quality = "hac", \
             overwrite =True, \
             trim = "--no-trim", \
             #filtering = "--min-qscore 10", \
             methylation = "5mCG_5hmCG"):
    """
    ***basecall*** uses ont dorado for bascalling of raw pod5 files.
    name_base - name to be used as prefix for all filenames
    input_folder - folder with raw pod5 files
    quality - quality for dorado basecaller ("fast", "hac", "sup") - standard is "hac"
    output_folder - destination to write basecalled .bam file
    overwrite - True/False option to overwrite previous files in output_folder - standard is True
    trim - dorado option, use "--no-trim" to preserve barcodes
    """
    hf.create_dir(output_folder)
    if overwrite:
        hf.empty_dir(output_folder)
    temp_folder = os.path.join(output_folder, "temp")
    function = "dorado basecaller"
    if methylation != "":
        model = quality+","+ methylation #"--modified-bases-models "+ 
    else:
        model = quality
    full_call = " ".join([function, model, input_folder, trim, "--output-dir", temp_folder])
    print("********** " +"BASECALLING: "+  name_base + " **********" )
    subprocess.call(full_call, shell=True)
    hf.bam_to_maindir(temp_folder, name_base, output_folder)
    shutil.rmtree(temp_folder, ignore_errors=True)


def filter_reads(name_base,\
                  input_folder, \
                  output_folder, \
                  overwrite =True,\
                 seqlen =500,\
                 minqual = 90, \
                use_filtlong = False):
    """
    ***filter_reads*** uses samtools on all bam-files found in input folder to filter reads above a certain read length.
    name_base - name to be used as prefix for all filenames
    input_folder - folder with basecalled .bam
    output_folder - destination to write filtered .bam file
    overwrite - True/False option to overwrite previous files in output_folder - standard is True
    seqlen - filter for minimum read lengths
    """
    hf.create_dir(output_folder)
    temp_folder = os.path.join(output_folder, "temp")
    hf.create_dir(temp_folder)
    print(temp_folder)
    if overwrite:
        hf.empty_dir(output_folder)
    print("********** " +"FILTERING: " +  name_base+ " **********" )
    for bamfile in os.listdir(input_folder):
        if bamfile.endswith('.bam'):
            bam_path = os.path.join(input_folder, bamfile)
            base_bam = os.path.splitext(bamfile)[0]
            filt_bam = os.path.join(output_folder, base_bam+"_f.bam")
            if use_filtlong:
                fastq_temp = os.path.join(temp_folder, base_bam+".fastq")
                fastq_filt = os.path.join(temp_folder, base_bam+"filt.fastq")
                # convert to fastq
                subprocess.call(" ".join(["samtools fastq",bam_path, " > ", fastq_temp]), shell=True)
                #filter fastq
                function = "filtlong --min_length "+str(seqlen)+" --min_mean_q "+str(minqual)
                full_call = " ".join([function, fastq_temp, " > ", fastq_filt])
                subprocess.call(full_call, shell=True)
                print("********** " +"INPUT: "+bam_path+" OUTPUT: " +  filt_bam+ " **********" )
                subprocess.call(full_call, shell=True)
                #convert back to bam_file
                function = "samtools view -h -e 'length(seq)>="+str(seqlen)+"' -b"
                full_call = " ".join(["samtools import", fastq_filt, " -o ", filt_bam])
                subprocess.call(full_call, shell=True)
            else:
                function = "samtools view -h -e 'length(seq)>="+str(seqlen)+"' -b"
                full_call = " ".join([function, bam_path, " -o ", filt_bam])
                print("********** " +"INPUT: "+bam_path+" OUTPUT: " +  filt_bam+ " **********" )
                subprocess.call(full_call, shell=True)

def demux(name_base,\
          input_folder, \
          output_folder, \
          dual_bc_only = False,\
          overwrite =True,\
          demux_kit = 'HelperFiles/demux_files.txt'):
    """
    ***demux*** uses ont dorado to demultiplex all bamfiles found in input folder and outputs one file per barcode in the output_folder.
    name_base - name to be used as prefix for all filenames
    input_folder - folder with basecalled .bam
    output_folder - destination to write filtered .bam file
    dual_bc_only - only set to True, if only reads with barcodes on both ends should be counted.     
    overwrite - True/False option to overwrite previous files in output_folder - standard is True
    demux_kit = csv file containing kit name, arrangement (.toml), sequences (.fasta) for dorado demux
    """
    hf.create_dir(output_folder)
    if overwrite:
        hf.empty_dir(output_folder)
    #define custom parameters
    df = pd.read_csv(demux_kit, sep = ",", index_col=0, header = None)
    barcode_kit = df.loc["kit"][1]
    bc_arrangement = df.loc["arrangement"][1]
    bc_sequences = df.loc["sequences"][1]
    bc_opts = "--kit-name "+barcode_kit+" --barcode-arrangement "+bc_arrangement+" --barcode-sequences "+bc_sequences 
    temp_folder = os.path.join(output_folder, "temp")
    hf.create_dir(temp_folder)
    for bamfile in os.listdir(input_folder):
        if bamfile.endswith('.bam'):
            bam_path = os.path.join(input_folder, bamfile)   
            function = "dorado demux"
            if dual_bc_only:
                full_call = " ".join([function, bam_path, bc_opts, "--no-trim --sort-bam --output-dir", temp_folder,"--barcode-both-ends"])
            else:
                full_call = " ".join([function, bam_path, bc_opts, "--no-trim --sort-bam --output-dir", temp_folder])
            print("********** " +"DE-MULTIPLEXING: "+  name_base+ " **********" ) 
            subprocess.call(full_call, shell=True)
    hf.bam_to_maindir(temp_folder, name_base, output_folder)   
    shutil.rmtree(temp_folder, ignore_errors=True)

def align(name_base,\
          input_folder, \
          output_folder, \
          overwrite =True,\
          reference_locations = 'HelperFiles/references.txt',\
          mmi_ref="mm39",\
         create_new_mmi_index = False,\
         minmapqual =20,\
         use_only_dorado = True):
    """
    ***align*** uses ont dorado/minimap2 to align all bamfiles found in input folder to a reference file (specified through 
    reference_locations and mmi_ref)and outputs one file per barcode in the output_folder
    overwrite - True/False option to overwrite previous files in output_folder - standard is True
    reference_locations - a csv file containing the file locations for reference files, following columns:
        -name:   name of the reference - referenced through mmi_ref to find corresponding filelocations
        -fasta:   location of a fasta (.fa) file for the corresponding reference
        -fai:   location of a fasta index (.fa.fai) file for the corresponding reference
        -mmi:   location of a pre-computed minimap2 index file (.mmi) for the corresponding reference
    mmi_ref - name of the reference to be used (must be found as name in reference_locations)
    create_new_mmi_index - If set to true, a new minimap index file will be generated (temporarily)
    """
    hf.create_dir(output_folder)
    unifltered_dir = os.path.join("unfiltered_aligned", name_base)
    hf.create_dir(unifltered_dir)
    if overwrite:
        hf.empty_dir(output_folder)
        shutil.rmtree(unifltered_dir, ignore_errors=True)
        hf.create_dir(unifltered_dir)
    temp_folder = os.path.join(output_folder, "temp")
    hf.create_dir(temp_folder)
    ref_df = pd.read_csv(reference_locations, sep = ",")
    out_bam_paths = []
    if create_new_mmi_index:
        index_file = os.path.join(temp_folder,"temp.mmi")
        ref_file = ref_df[ref_df['name'].str.contains(mmi_ref)].iloc[0]['fasta']
        minimap_call = " ".join(["minimap2 -x map-ont -d ",index_file,ref_file])
        subprocess.call(minimap_call, shell=True)
        mmi = os.path.abspath(index_file)
    else:
        try:
            mmi = ref_df[ref_df['name'].str.contains(mmi_ref)].iloc[0]['mmi']
        except:
            mmi = os.path.abspath(mmi_ref)
    for bamfile in os.listdir(input_folder):
            if bamfile.endswith('.bam'):
                bam_path = os.path.join(input_folder, bamfile)
                print(bam_path)
                base_bam = os.path.splitext(bamfile)[0]
                ######### version with dorado #########
                
                if use_only_dorado:
                    ref_file = ref_df[ref_df['name'].str.contains(mmi_ref)].iloc[0]['fasta']
                    #mm2_opts = "--mm2-opts '-N 0'"
                    full_call = " ".join(["dorado aligner", ref_file, bam_path, "-o", temp_folder])
                    print("********** " +"ALIGNING: "+  name_base+ bamfile +" **********" )
                    subprocess.call(full_call, shell=True)
                else:
                    function = "dorado aligner"
                    mm_opts = "--mm2-opts '-k 15 -w 10 --secondary=no'"
                    full_call = " ".join([function, mmi, mm_opts, bam_path, "-o", temp_folder])
                    print("********** " +"ALIGNING: "+  name_base+ bamfile +" **********" )
                    subprocess.call(full_call, shell=True)
                    
                outpath = hf.bam_to_maindir(temp_folder, name_base, output_folder, input_alignment = bam_path)
                print(outpath)
                out_bam_paths += outpath
                print(out_bam_paths)

    ########### create bigwig #############
    for newbampath in list(set(out_bam_paths)):
        print("********** " +"FILTERING & GENERATING BIGWIGs"+" **********" )
        ########### filter #############
        function = "samtools view -b -h -F 0x900 -q" #  
        filt_bam_path = os.path.splitext(newbampath)[0]+"_ns"+".bam"
        filt_bigwig = os.path.splitext(filt_bam_path)[0]+".bw"
        full_call = " ".join([function, str(minmapqual), newbampath,"| samtools sort -o", filt_bam_path, "&& samtools index", filt_bam_path, "&& bamCoverage -b", filt_bam_path, "-o", filt_bigwig])
        subprocess.call(full_call, shell=True)
        #### move old bam and bai files
        if os.path.exists(newbampath):
            newbaipath = os.path.splitext(newbampath)[0]+".bam.bai"
            for file in [newbampath, newbaipath]:
                shutil.move(file, unifltered_dir, copy_function = shutil.copytree) 
    shutil.rmtree(temp_folder, ignore_errors=True)

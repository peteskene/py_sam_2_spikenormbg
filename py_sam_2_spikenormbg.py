def py_sam_2_spikenormbg(folder=None, spike_suffix='dm6', input_type = 'sam', ends=False, lengths_analysis=True, lengths_image=True, size_select_1=True, size_select_2=True, bedgraph=True, chrom_sizes='hg19',
                         big_wig=True, size_min_1 = 20, size_max_1=120, size_min_2=150, size_max_2=710, multiplying_factor=10000):
    
    """
    Script written by Pete Skene (peteskene@gmail.com). Free for academic use only.
    
    Script will take a folder of sam or bam files and generate bedgraphs. The sam/bam files need to be of properly matched pairs.
    Will always make bedgraphs using the entire inserts, can also make bedgraphs using just insert ends (see ends option)
    
    Expected file naming convention for data_files and spike_files, e.g. PS_HsDm_CTCF_1m.sam and PS_HsDm_CTCF_1m.sam.FLY
    ___________
    Parameters:
    -folder: path to folder containing sam or bam files. e.g. '/home/pskene/test_data/test_samfiles', 
            if 'None' then looks in current working directory
    -input_type: string describing either 'bam' or 'sam'. Default = 'sam'
    -ends: if True, will also make bed files and bedgraphs for just the ends of the inserts
    -lengths_analysis: calculate size distribution of all insert sizes in input sam or bam
    -lengths_image: if true, will generate and save a plot displaying length distribution for all data files 
    -size_select_1 or 2: Perform size selection (allows two size selections)
    -bedgraph: if True, will generate spike normalized bedgraphs from every generated bedfile,
            set to false, if you just want to make bedfiles
    -chrom_sizes: type genome build (e.g. hg19, mm9) to autogenerate chrom_sizes file, otherwise can provide file
    -big_wig: generate bigwig files from all the generated bedgraphs
    -multiplying_factor: stops bedgraph values being very small after normalizing to spikefile counts, defaults to 10000
    ___________
    Requirements (type the following before invoking jupyter):
    pip install pybedtools (only needs to be done once to install)
    pip install pysam (only needs to be done once to install)
    module load ucsc/2014_01_21 (needs to be every time, required to make big_wig, so can change to big_wig=False) 
    """
    
    import pybedtools 
    from pybedtools import BedTool
    import glob
    import os
    from subprocess import check_output
    import pandas as pd
    from pybedtools.contrib.bigwig import bedgraph_to_bigwig
    import matplotlib.pyplot as plt
    %matplotlib inline
    
    from datetime import datetime
    startTime = datetime.now()
    
    pybedtools.set_tempdir('/loc/scratch/' + user)
    
    #change directory as instructed
    if folder != None:
        os.chdir(folder)
        print 'will look for files in: ' + os.getcwd()    
    
    #generate list of names for data_files and spike_files from the Folder directory
    if input_type == 'sam':
        dataFiles = sorted(glob.glob('*.sam'))
    if input_type == 'bam':
        dataFiles = sorted(glob.glob('*.bam'))
      
    spikeFiles = sorted(glob.glob('*.' + spike_suffix))
    
    #print out the list of data_files and spike_files
    print 'WARNING: to make bigwig files the UCSC module must be loaded into rhino ' \
    'before starting jupyter kernel. Script will throw an error if big_wig=True and '\
    'module is not loaded!'
    print '\n'    
    print 'Data files imported as .sam or .bam :'
    print '\n'.join(dataFiles)
    print '\n'
    print 'Spike files imported as .spike_suffix'
    print '\n'.join(spikeFiles)
    print '\n'

    #check for equal numbers and matched data_files and spike_files. If incorrect, exit script
    if len(dataFiles)==0:
        return 'No datafiles imported. Exiting...'
    
    if len(dataFiles) != len(spikeFiles):
        return 'Unequal numbers of data files and spike files. Exiting...'
    
    for i in range(len(dataFiles)):
        if dataFiles[i] not in spikeFiles[i]:
            print 'Unmatched pairs and spike files.'
            print dataFiles[i] + ' not paired with ' + spikeFiles[i]            
            return 'Exiting...'
     
    print 'Properly matched data_files and spike_files. Continuing...'

    
    #if sam, then need to be converted to bam files
    if input_type == 'sam':
        print 'Files imported as "sam". Converting to bam format...'
        print '\n'
        
        #generate list of names for the bam files from the sam files
        bam_names = [f.replace('sam', 'bam') for f in dataFiles]
        
        #generate shell string for sam to bam conversion
        bam_string = []
        for i in range(len(bam_names)):
            bam_string.append('samtools view -b -S ' +  dataFiles[i] + ' > ' + bam_names[i])
   
        #calling samtools shell script to convert sam to bam for each file (will be saved to 'folder')
        for item in bam_string:
            check_output(item, shell = True)
        
        #replace dataFiles list with the bam_names list
        dataFiles = bam_names
 
    elif input_type == 'bam':
        print 'Files imported as "bam". Continuing...'
        print '\n'
        
    elif input_type != 'sam' or input_type != 'bam':
        return 'Unrecognized input type. Exiting...'
    

    
    ####################################################
    
    #generate bed file from the bam files (this assumes the bam files are just properly patched pairs)
    print 'Generating bed files representing whole insert from paired end reads in the data files'
    print '\n'
    
    if ends:
        #will make bed files plotting just the single base pair position of both ends of the inserts
        print 'Generating bed files plotting just the single base pair position of both ends of the inserts'
        print '\n'
        
    if size_select_1 or size_select_2:
        print 'Generating size selected bed files'
        print '\n'

    
    #generate bed file names from data file names (even if ends/size selection set to false)
    bed_names  = [f.replace('bam', 'bed') for f in dataFiles]
    
    bed_ends_names = [f.replace('bam', 'ends.bed') for f in dataFiles]
        
    size_selected_files_1 = [f.replace('bam', str(size_min_1) + '_' + str(size_max_1) + '.bed') for f in dataFiles]
        
    size_selected_files_1_ends = [f.replace('bam', str(size_min_1) + '_' + str(size_max_1) + '.ends.bed') for f in dataFiles]
  
    size_selected_files_2 = [f.replace('bam', str(size_min_2) + '_' + str(size_max_2) + '.bed') for f in dataFiles]
        
    size_selected_files_2_ends = [f.replace('bam', str(size_min_2) + '_' + str(size_max_2) + '.ends.bed') for f in dataFiles]
    
    all_beds = bed_names + bed_ends_names + size_selected_files_1 +size_selected_files_1_ends + size_selected_files_2 + size_selected_files_2_ends
    
   
    #generate filenames for length analysis as will perform on each datafile on fly rather than reloading
    lengths_names  = [f.replace('bam', 'lengths') for f in dataFiles]
    #create empty dataframe to be filled by each length analysis, used to plot lengths distribution 
    lengths_df = pd.DataFrame()

    #####################################################
    
    #generate bed files using bam_to_bed tool (makes bed12 format)
    for i in range(len(dataFiles)):
        temp_bed = BedTool(dataFiles[i]).bam_to_bed(bedpe=True).to_dataframe()
        
        #need to strip out the start and end position of the whole insert (bed12 is both sequenced reads)
        #note column names actually represent <chrom> <start of insert> <end of insert>
        temp_bed_stripped = temp_bed.iloc[:, [0,1,5]].sort_values(by = ['chrom', 'start', 'strand'])
        
        #calculate insert size and insert as column 4 and save file with bed_name
        #these bed files represent the entire insert
        temp_bed_stripped['length'] = temp_bed_stripped['strand'] - temp_bed_stripped['start']
        
        temp_bed_stripped.to_csv(bed_names[i], sep="\t", header=False, index=False)
        
        #perform analysis on the length of inserts sequenced
        if lengths_analysis:
            temp_lengths = temp_bed_stripped.groupby(by=['length'])['length'].count()
            
            temp_lengths.to_csv(lengths_names[i], sep='\t', header = [bed_names[i]], index=True, index_label='length')
            
            #add the lengths data from this datafile to the lengths_df dataframe, title each series with bed file name
            lengths_df = lengths_df.join(temp_lengths.rename(bed_names[i]), how='outer')
                      
        
        #generate size selected whole insert bed files
        if size_select_1:
            subset_1 = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=size_min_1)
                                         & (temp_bed_stripped.iloc[:,3]<=size_max_1)]
            
            subset_1.to_csv(size_selected_files_1[i], sep="\t", header=False, index=False)
        
        if size_select_2:
            subset_2 = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=size_min_2)
                                         & (temp_bed_stripped.iloc[:,3]<=size_max_2)]
            
            subset_2.to_csv(size_selected_files_2[i], sep="\t", header=False, index=False)
        
    
    
        if ends:
            #will make bed files plotting just the single base pair position of both ends of the inserts
            
            #rename columns for transparency and calculate base pair positions of insert ends
            #note left and right does not refer to sequenced strand
            temp_bed_stripped.columns = ['chrom', 'left_start', 'right_stop', 'length']
            
            temp_bed_stripped['left_stop'] = temp_bed_stripped['left_start']+1
            temp_bed_stripped['right_start'] = temp_bed_stripped['right_stop']-1
            
            #split into left and right ends 
            temp_bed_stripped_left = temp_bed_stripped[['chrom', 'left_start', 'left_stop', 'length']]
            temp_bed_stripped_right = temp_bed_stripped[['chrom', 'right_start', 'right_stop', 'length']]
            
            #rename columns to allow concatentation
            temp_bed_stripped_left.columns = ['chrom', 'start', 'stop', 'length']
            temp_bed_stripped_right.columns = ['chrom', 'start', 'stop', 'length']
            
            temp_ends = pd.concat([temp_bed_stripped_left, temp_bed_stripped_right],
                                  ignore_index=True).sort_values(by=['chrom', 'start', 'stop'])
            
            temp_ends.to_csv(bed_ends_names[i], sep="\t", header=False, index=False)
            
            #generate size selected ends bed files
            if size_select_1:
                subset_ends_1 = temp_ends[(temp_ends.iloc[:,3]>=size_min_1)
                                         & (temp_ends.iloc[:,3]<=size_max_1)]
            
                subset_ends_1.to_csv(size_selected_files_1_ends[i], sep="\t", header=False, index=False)
        
            if size_select_2:
                subset_ends_2 = temp_ends[(temp_ends.iloc[:,3]>=size_min_2)
                                         & (temp_ends.iloc[:,3]<=size_max_2)]
            
                subset_ends_2.to_csv(size_selected_files_2_ends[i], sep="\t", header=False, index=False)
    
    print 'finished generating bed files:'
    print '\n'
    print 'whole insert bed files:'+'\n'+'\n'.join(bed_names)
    print '\n'
    if ends:
        print 'insert ends bed files:'+'\n'+'\n'.join(bed_ends_names)
        print '\n'
    if size_select_1:
        print 'whole insert bed files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1)
        print '\n'
        if ends:
            print 'insert ends bed files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_ends)
            print '\n'
    if size_select_2:
        print 'whole insert bed files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2)
        print '\n'
        if ends:
            print 'insert ends bed files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_ends)
            print '\n'
    
  
    #####################################################
    #take all the bed files and generate spike normalized bedgraph files
        
    if bedgraph:
        print 'Generating spike normalized bedgraphs from all the bed files'
        print '\n'
    
        #generate file names for the bedgraphs
        bg_names  = [f.replace('bed', 'bg') for f in bed_names]
    
        bg_ends_names = [f.replace('bed', 'bg') for f in bed_ends_names]
        
        size_selected_files_1_bg = [f.replace('bed', 'bg') for f in size_selected_files_1]
         
        size_selected_files_1_ends_bg = [f.replace('bed', 'bg') for f in size_selected_files_1_ends]
  
        size_selected_files_2_bg = [f.replace('bed', 'bg') for f in size_selected_files_2]
        
        size_selected_files_2_ends_bg = [f.replace('bed', 'bg') for f in size_selected_files_2_ends]
        
        all_bg = [f.replace('bed', 'bg') for f in all_beds]
            
        #calculate the number of reads in each of the spike sam files
        #spikeCount is a lit of reads populated by the 'samtools view -c' shell command for each spikeFile
        spikeCount = []

        spike_string = []

        for item in spikeFiles:
            spike_string.append('samtools view -c -S ' + item)

        for item in spike_string:
            spikeCount.append(int(check_output(item, shell =True))/2)
            
        print 'spike counts:' 
        print spikeCount
        
        #calculating list of scaling factors
        scaling_factor = []
    
        for item in spikeCount:
            scaling_factor.append(float(10000)/item)
        
        #run bedtools genomecov to generate bedgraph files
        for i in range(len(bg_names)):
            BedTool(bed_names[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(bg_names[i])
    
        if ends:
            for i in range(len(bg_ends_names)):
                BedTool(bed_ends_names[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(bg_ends_names[i])
                
        if size_select_1:
            for i in range(len(size_selected_files_1_bg)):
                BedTool(size_selected_files_1[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(size_selected_files_1_bg[i])
                
            if ends:
                for i in range(len(size_selected_files_1_ends_bg)):
                    BedTool(size_selected_files_1_ends[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(size_selected_files_1_ends_bg[i])
                    
        if size_select_2:
            for i in range(len(size_selected_files_2_bg)):
                BedTool(size_selected_files_2[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(size_selected_files_2_bg[i])
                
            if ends:
                for i in range(len(size_selected_files_2_ends_bg)):
                    BedTool(size_selected_files_2_ends[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(size_selected_files_2_ends_bg[i])
             
        print 'finished generating bedgraph files:'
        print '\n'
        print 'whole insert bedgraph files:'+'\n'+'\n'.join(bg_names)
        print '\n'
        if ends:
            print 'insert ends bedgraph files:'+'\n'+'\n'.join(bg_ends_names)
            print '\n'
        if size_select_1:
            print 'whole insert bedgraph files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_bg)
            print '\n'
            if ends:
                print 'insert ends bedgraph files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_ends_bg)
                print '\n'
        if size_select_2:
            print 'whole insert bedgraph files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_bg)
            print '\n'
            if ends:
                print 'insert ends bedgraph files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_ends_bg)
                print '\n'
    
    #####################################################
    #make bigwig files from all the bedgraphs generated
        
    if big_wig:
        print 'Generating big_wig files from each of the bedgraphs'
        
        if big_wig==True and bedgraph==False:
            return 'WARNING: no bedgraphs to make into big_wig files'
            
        #generate file names for the bigwigs
        bw_names  = [f.replace('bg', 'bw') for f in bg_names]
    
        bw_ends_names = [f.replace('bg', 'bw') for f in bg_ends_names]
        
        size_selected_files_1_bw = [f.replace('bg', 'bw') for f in size_selected_files_1_bg]
         
        size_selected_files_1_ends_bw = [f.replace('bg', 'bw') for f in size_selected_files_1_ends_bg]
  
        size_selected_files_2_bw = [f.replace('bg', 'bw') for f in size_selected_files_2_bg]
        
        size_selected_files_2_ends_bw = [f.replace('bg', 'bw') for f in size_selected_files_2_ends_bg]
        
        all_bw = [f.replace('bg', 'bw') for f in all_beds]
        
                
        #run bedgraph_to_bigwig tool
        for i in range(len(bg_names)):
            bedgraph_to_bigwig(BedTool(bg_names[i]), chrom_sizes, bw_names[i])
    
        if ends:
            for i in range(len(bg_ends_names)):
                bedgraph_to_bigwig(BedTool(bg_ends_names[i]), chrom_sizes, bw_ends_names[i])
                
        if size_select_1:
            for i in range(len(size_selected_files_1_bg)):
                bedgraph_to_bigwig(BedTool(size_selected_files_1_bg[i]), chrom_sizes, size_selected_files_1_bw[i])
                
            if ends:
                for i in range(len(size_selected_files_1_ends_bg)):
                    bedgraph_to_bigwig(BedTool(size_selected_files_1_ends_bg[i]), chrom_sizes, size_selected_files_1_ends_bw[i])
                    
        if size_select_2:
            for i in range(len(size_selected_files_2_bg)):
                bedgraph_to_bigwig(BedTool(size_selected_files_2_bg[i]), chrom_sizes, size_selected_files_2_bw[i])
                
            if ends:
                for i in range(len(size_selected_files_2_ends_bg)):
                    bedgraph_to_bigwig(BedTool(size_selected_files_2_ends_bg[i]), chrom_sizes, size_selected_files_2_ends_bw[i])
       
        print 'finished generating bigwig files:'
        print '\n'
        print 'whole insert bigwig files:'+'\n'+'\n'.join(bw_names)
        print '\n'
        if ends:
            print 'insert ends bigwig files:'+'\n'+'\n'.join(bw_ends_names)
            print '\n'
        if size_select_1:
            print 'whole insert bigwig files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_bw)
            print '\n'
            if ends:
                print 'insert ends bigwig files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_ends_bw)
                print '\n'
        if size_select_2:
            print 'whole insert bigwig files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_bw)
            print '\n'
            if ends:
                print 'insert ends bigwig files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_ends_bw)
                print '\n'
    
    if lengths_image:
        if lengths_analysis==False:
            print 'lengths analysis set to false, so no image to display'
            
        else:
            print 'saving combined lengths distribution to file: ' + os.getcwd().split('/')[-1] + str('.lengths')
            
            temp_name = os.getcwd().split('/')[-1] + str('.lengths')
            lengths_df.to_csv(temp_name, sep='\t', header=True, index=True, index_label='bp')
            
            print 'generating image of lengths distribution'
            
            temp_plot_name = temp_name + str('_plot.png')
            
            fig = plt.figure(figsize=(12,6))
            ax  = fig.add_subplot(111)
            ax.set_position([0.1,0.1,0.5,0.8])
            ax.plot(lengths_df)
            leg = ax.legend(lengths_df.columns.values.tolist(), loc = 'center left', bbox_to_anchor = (1.0, 0.5))
            plt.title('Length Distribution')
            plt.xlabel('Insert Lengths (bp)')
            plt.ylabel('Count')
            fig.savefig(temp_plot_name)

            
                
    print 'Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
    return 'Finished'
    
 

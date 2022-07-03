#!/usr/bin/python
# -*- coding: utf-8 -*-
# CONDA ENV SOFTWARE: python (3.7), trimmomatic (v0.39), picard (v2.21.1), bwa (v0.7.17), samtools (v1.9), gatk4 (v4.1.4.0), vcftools (v0.1.16)

import os, sys, operator, time, argparse, subprocess, re, collections, statistics, itertools # required modules
from core import CAPTURE, SBATCH, REVIEW, ezSub, FindSupplementaryFile, JoinNeatly

pipe_script = os.path.realpath(__file__)
buddy_script = f'{os.path.dirname(pipe_script)}/sub_buddy.py'

parser = argparse.ArgumentParser(description='Next Generation Sequencing Pipe', prog=os.path.basename(__file__), usage='%(prog)s [options]', epilog='see the readme file for further details')

modes = parser.add_mutually_exclusive_group(required=True) # run modes
modes.add_argument('-setup', action='store_true', help='setup initial directories')
modes.add_argument('-index', action='store_true', help='index reference genome')
modes.add_argument('-stage', metavar='<number>', dest='stage', type=int, choices=range(1,15), help='specify processing stage')

additional = parser.add_argument_group(description='additional modes:') # additional modes
additional.add_argument('-pipe', action='store_true', help='submit pipe as task')
additional.add_argument('-review', action='store_true', help='review task accounting data')
additional.add_argument('-test', action='store_true', help='test script locally')

inputs = parser.add_argument_group(description='user inputs:') # user inputs
inputs.add_argument('-u', metavar='<name>', type=str, help='specify user name')
inputs.add_argument('-p', metavar='</path/to/>', type=str, help='specify path')
inputs.add_argument('-d', metavar='<name>', type=str, default='ngs_pipe', help='specify working directory')
inputs.add_argument('-m', metavar='<name>', type=str, help='specify user email address')
inputs.add_argument('-e', metavar='<name>', type=str, default='ngs_env', help='specify conda environment')
inputs.add_argument('-l', metavar='<number>', type=int, default=100, help='specify parallel task limit')
inputs.add_argument('-x', metavar='<number>', type=int, help='specify ploidy for calling haplotypes')

hpcc = parser.add_argument_group(description='hpcc settings:') # overide hpcc settings
hpcc.add_argument('-pa', metavar='<name>', type=str, choices=['defq','shortq','devq','mmemq','hmemq','voltaq','visq'], help='specify partition')
hpcc.add_argument('-no', metavar='<number>', type=int, help='specify nodes')
hpcc.add_argument('-nt', metavar='<number>', type=int, help='specify ntasks')
hpcc.add_argument('-me', metavar='<number[units]>', type=str.lower, help='specify memory [k|m|g|t]')
hpcc.add_argument('-wt', metavar='<HH:MM:SS>', type=str, help='specify walltime')
hpcc.add_argument('-sb', action='store_true', help='assistance from sub buddy')

setup, indexing, stage, submitting_self, reviewing, testing, user, path, wrk_dir, address, environment, limit, ploidy, *hpcc_overides, subbuddy = vars(parser.parse_args()).values() # define user inputs

*_, memory, walltime = hpcc_overides # define hpcc overides
if memory and (sum(memory.count(unit) for unit in ['k','m','g','t']) != 1 or not memory.rstrip('kmgt').isdigit()): parser.error(f"argument -me: invalid str format: '{memory}'") # check memory format
if walltime and ( not 8 <= len(walltime) <= 9 or walltime[-6:-2:3].count(':') != 2 or not walltime.replace(':','').isdigit()): parser.error(f"argument -wt: invalid str format: '{walltime}'") # check walltime format

if stage == 5 and not ploidy: parser.error(f"the argument -x is required for stage {stage}") # check required input
if stage != 5 and ploidy: parser.error(f"the argument -x is not required for stage {stage}") # check required input

(pipe_flag, *_), *_ = [ info.option_strings for info in additional._group_actions ] # extract pipe flag



def sort_via_L(name):
    # extract non-overlapping lane references
    Lx = re.findall('_L[0-9]+_', name)
    # substitute no matches with zero
    if not Lx:
        L = '0'
        print(f'Sequencing lane reference not found in {name} - files may not be concatenated in correct order')
    # extract last lane reference match
    else:
        print(Lx)
        *_, L = Lx

    # remove non-integer text integers from lane reference
    L = L.strip('L_')
    
    return int(L)


def SplitRegions(region_list):
    if len(region_list) < 50: # when less than 50 scaffolds process each seperately
        grouped = [ [info] for info in region_list ]
    else: # when more than 50 scaffolds process groups in blocks
        *_, block = grouped = [[]]
        for info in region_list: 
            if not block or sum(int(length) for scaffold,length in block) < 50000000: 
                block.append(info) # group scaffolds into blocks no more than 50 M BP
            else: 
                grouped.append([]) # create new block
                *_, block = grouped # reset block variable
                block.append(info) # stored region info in new block
    return [ list(scaffolds) for scaffolds,*_ in [zip(*regions) for regions in grouped] ] # extract scaffold names

directory_name, in_suffix, out_suffix = category_labels = 'directory_name','in_suffix','out_suffix' # specify subdictionary labels
stage_formats = { # stage-specific formatting
  'i': [ 'ref.gen',            ('.fasta','.fa','.faa','.fas')                                                                   ],
    1: [ '00.fastqs',          ('.fq.gz', '.fastq.gz'),         ('fastq.gz')                                                    ],
    2: [ '01.merged',          ('.fastq.gz',),                  ('fastq.gz')                                                    ],
    3: [ '02.trimmed',         ('.fastq.gz',),                  ('aln-pe.sam','srtd.bam','flagstat.txt')                        ],
    4: [ '03.mapped',          ('.srtd.bam',),                  ('mrkd.bam','metrics.txt','mrkd.nmd.bam')                       ],
    5: [ '04.marked',          ('.nmd.bam',),                   ('vcf.gz', 'bam')                                               ],
    6: [ '05.haplotyped',      ('.vcf.gz',),                    ('combined.vcf.gz')                                             ],
    7: [ '06.combined',        ('.vcf.gz',),                    ('genotyped.vcf.gz')                                            ],
    8: [ '07.genotyped',       ('.vcf.gz',),                    ('recombined.vcf.gz')                                           ],
    9: [ '08.recombined',      ('.vcf.gz',),                    ('filtered.biallelic.vcf.gz','filtered.best.vcf.gz','filtered') ],
   10: [ '09.filtered.best',   ('best.vcf.gz','.ldepth'),       ('depth.mask.vcf.gz','filtered.F3.vcf.gz','filtered.F4.vcf.gz') ],
   11: [ '10.filtered.depth'  ]}
stage_formats = { key:  {label:info for label,info in zip(category_labels,value)} for key,value in stage_formats.items() } # creatre subdictionaries

resource_labels = 'partition','nodes','ntasks-per-node', 'memory', 'walltime' # specify subdictionary labels
hpcc_settings = { # default hpcc settings 
  'i': ['defq','1','1','8g', '24:00:00'],   1: ['defq','1','1','4g' ,'24:00:00'],   2: ['defq','1','1','4g', '24:00:00'],   3: ['defq','1','1','32g','48:00:00'],  4: ['defq','1','1','32g','24:00:00'],
    5: ['defq','1','1','64g','96:00:00'],   6: ['defq','1','1','64g','96:00:00'],   7: ['defq','1','1','64g','96:00:00'],   8: ['defq','1','1','32g','48:00:00'],  9: ['defq','1','1','32g','24:00:00'],  
   10: ['defq','1','1','32g','24:00:00'], '+': ['defq','1','1','4g','168:00:00'] }
hpcc_settings = { key:  {label:info for label,info in zip(resource_labels,value)} for key,value in hpcc_settings.items() } # creatre subdictionaries

print(f'\n{"SETUP" if setup else f"REVIEWING " if reviewing else "SUBMITTING " if stage else ""}{"INDEXING REFERENCE GENOME" if indexing else f"STAGE {stage} TASKS" if stage else ""}\n')

user, path = [ os.getenv(bash) if not argument else argument for argument,bash in [ (user,'USER'),(path,'HOME')] ] # find user & path & directory
wrk_dir = f'/{path.strip("/")}/{wrk_dir.strip("/")}'  # ensure correct path format
system_script, *arguments = sys.argv # extract script location & arguments

if setup: 
    initial_dirs = [ f'{wrk_dir}/{subdict[directory_name]}' for key,subdict in stage_formats.items() if str(key).isalpha() or key <= 1] # extract initial directory names
    [ os.makedirs(path, exist_ok=True) for path in initial_dirs ]; print('SETUP COMPLETE\n'); exit() # create initial directories

else: # proceed with alternative pipeline mode    
    assert(os.path.exists(wrk_dir)), f'Problem finding the working directory "{wrk_dir}".'# check path exists
    assert(environment in CAPTURE('conda info --env')), f'Problem finding the conda environment "{environment}".' # check conda environment exists

    relevant = 'i' if indexing else stage # specify stage key

    stage_info, next_stage_info = operator.itemgetter(relevant, relevant if indexing else relevant+1)(stage_formats) # specify stage formats
    batch_labels = [ f'x.slurm/{"pipe" if submitting_self else "index" if indexing else f"stage.{stage:02}"}/{label}' for label in ['scripts','out.err','ids'] ] # specify batch formats        
    
    *_, sh_dir, oe_dir, id_dir = in_dir, out_dir, *stage_dirs = [ f'{wrk_dir}/{label}' for label in [stage_info[directory_name], next_stage_info[directory_name], *batch_labels] ] # specify stage directories
    id_file = f'{id_dir}/task.ids' # specify slurm id file

    if reviewing: # review task accounting data
        if not os.path.exists(id_file): print(error_no_files.format(f'\"{id_file}\"')) # check id file exists
        else: REVIEW(id_file) # review task accounting data via job ids

    else: # proceed with indexing or stage
        [ os.makedirs(stage_dir, exist_ok=True) for stage_dir in stage_dirs ] # make stage-specific directories as required
        
        hpcc_settings[relevant].update({ label:setting for label,setting in zip(resource_labels, hpcc_overides) if setting }) # update hpcc settings as required
        partition, nodes, ntasks, memory, walltime = resources = [  hpcc_settings['+' if submitting_self else relevant][label] for label in resource_labels ] # extract hpcc settings

        setting_info = zip(['user','working directory','conda environment','ploidy',*resource_labels],[user,wrk_dir,environment,ploidy,*resources])
        [ print(f'{label}: {info}') for label,info in setting_info ]


        if indexing or stage >= 2: # find supplementary files
            print('\nFINDING SUPPLEMENTARY FILES...')
            ref_dir = f'{wrk_dir}/{stage_formats["i"][directory_name]}' 
            fasta_file = FindSupplementaryFile(search_dir=ref_dir, search_suffix=stage_formats['i'][in_suffix] ) # find reference genome
            if indexing:
                if not fasta_file.endswith('.fasta'): # check suffix appropriate for GATK
                    print('\tRENAMING REFERENCE GENOME')
                    fasta_prefix, *_ = fasta_file.rsplit('.', 1) # extract reference prefix
                    new_fasta_file = f'{fasta_prefix}.fasta' # specify new file name
                    os.rename(fasta_file,new_fasta_file) # rename original file
                    fasta_file = new_fasta_file # update reference genome variable
                chr_names = [ name.lstrip('>') for name in CAPTURE(f'grep \> {fasta_file}').split('\n') ] # extract chromosome names
                dodgy_names = any('|' in name for name in chr_names) # check if any problematic characters in chromosome names
                proceeding = False if dodgy_names else True
                while not proceeding: # allow user to manually change scaffold names
                    response = input('\nSome scaffold names contain pipe characters that can cause problems with the software used here - would you like to quit & remove these? (y/n): ') # instruct task resubmission
                    proceeding = response in ['y','n']
                    if response == 'n': print('\nOk but you\'ve been warned.\n')
                    if response == 'y': print('\nExiting.\n'); sys.exit(0)
            if stage == 7: 
                fai_file = FindSupplementaryFile(search_dir=ref_dir, search_suffix=('.fai',)) # find regions file
            print('\tSUPPLEMENTARY FILES FOUND!')


        if indexing: 
            to_process = ['bwa index','samtools faidx','gatk CreateSequenceDictionary -R'] # define index stages
        if stage: # find input files
            print('\nFINDING TASK FILES...')
            assert(os.path.isdir(in_dir)), f'Problem finding the input directory \"{in_dir}\".'
            identified_inputs = [ (root,files) for root,*_,files in os.walk(in_dir, topdown=False) if any(name.endswith(stage_info[in_suffix]) for name in files) ] # find relevant inputs   
            assert(identified_inputs), f'No files found in \"{in_dir}\" ending with {JoinNeatly(stage_info[in_suffix])}.'
            if 'vcf' in stage_info[in_suffix]: 
                accompanied_by_indexes = [ files for *_,files in identified_inputs if not any(name.endswith('.tbi') for name in files) ] # find vcf files without accompanying tbi fil
                assert(accompanied_by_indexes), 'Not all vcf files are accompanied by ".tbi" indexes.' # check index file (TBI) exists (i.e. previous stage completed)

            file_to_process, *_ = files_to_process = sorted([ (root, [ f'{root}/{name}' for name in files if name.endswith(stage_info[in_suffix]) ]) for root,files in identified_inputs ]) # filter relevant files 

            if stage == 7: # genomic region to process
                regions_info = [ (scaffold,length) for scaffold,length,*_ in [line.split('\t') for line in open(fai_file, 'r').readlines()] ] # extract scaffolds
                regions_to_process = SplitRegions(region_list=regions_info) # group scaffolds into manageable parts (i.e. single scaffolds or groups of less than 50,000,000 bp)      

            to_process = files_to_process if stage < 6 else regions_to_process if stage == 7 else [1] # exceptions from stage 6 onwards
            print(f'\t{len(to_process)} TASK{"S" if len(to_process) > 1 else ""} FOUND!')

        scripts = []
        paired, single = categories = ['paired','single']
        processing_mode = { category:[] for category in categories }

        single_output_stages = [6,8,9,10]

        for i,task in enumerate(to_process, 1):
            
            if indexing: tool, *_ = task.split(' ', 1) # extract software name from command

            elif stage != 6 and stage != 8: # not required for stages 6 & 8 nor indexing
                *_, (*_, seq_file) = in_subdir, task_files = file_to_process if stage > 6 else task # specify input subdirectories & files

                if stage < 9:
                    sample = f'part{i:02}' if stage == 7 else os.path.basename(in_subdir) # specify sample
                    out_subdir = f'{out_dir}/{sample}' # specify output sub directory
                    #os.makedirs(out_subdir, exist_ok=True) # make output sub directories
            
            task_id = tool if indexing else 'merging' if 'combine' in out_dir else 'filtering' if 'filter' in out_dir else sample # specify task id
            
            if ploidy: task_id += f'-{ploidy}x' # specify task id ploidy info

            sh_file = f'{sh_dir}/{task_id}.sh' # specify sh file
            
            scripts.append(sh_file) # record sh file

            with open(sh_file, 'w') as sh:
                        
                hpcc_directives = SBATCH(job_id=task_id, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, out_err=oe_dir, conda_env=environment) # specify hpcc directives (slurm)
                sh.write(hpcc_directives)                

                sh.write('echo TASK STARTED `date`\n') # log start time

                # INDEX     ref.gen
                
                if indexing: sh.write(f'{task} {fasta_file}\n') # index FASTA reference sequence
                else:

                    task_dir = out_dir if stage in single_output_stages else out_subdir

                    sh.write(
                        f'rm -rf {task_dir}\n'+
                        f'mkdir -p {task_dir}\n'
                        )

                # CHECK READ TYPE

                    if stage <= 3:

                        raw_suffixes = '|'.join([ re.escape(suffix) for suffix in stage_info[in_suffix] ]) # specify regex search pattern
                        *_, P, U = S,  *PU  = ['S','P','U'] # specify read identifiers
                        files_by_member = { member: [ name for name in task_files if re.search(f'_[R]?({member})[{P}]?({raw_suffixes})$', name) ] for member in [1,2] } # categorise by read pair member if present
                        
                        all_reads_assigned_member = set(itertools.chain(*files_by_member.values())).issuperset(task_files) # establish if all reads belong to members
                        read_pairs_found = all(files_by_member.values()) and all_reads_assigned_member if stage == 1 else all(files_by_member.values()) # check that both members of read pair present & all reads allocated a member

                        read_type = paired if read_pairs_found else single
                        processing_mode[read_type].append(sample)
                        
                        relevant_reads = files_by_member if read_pairs_found else {S: task_files}

                        sh.write(f'echo PROCESSING AS {read_type.upper()} END\n')


                        # STAGE 1   00.fastqs --> 01.merged

                        if stage == 1:

                            for label, reads in relevant_reads.items():
                                
                                # added to sort paired reads by lane
                                if read_pairs_found:
                                    
                                    reads = sorted( reads, key=lambda name: sort_via_L(name) )

                                in_reads = ' '.join(reads)
                                out_reads = f'{out_subdir}/{sample}_{label}.{stage_info[out_suffix]}'
                                sh.write(f'cat {in_reads} > {out_reads}\n') # merging seperate reads 
                        
                        else: # process reads as appropriate
                            relevant_members = [1,2] if read_pairs_found else [S] # specify relevant read members
                            in_reads = ' '.join([ relevant_reads[member].pop()  for member in relevant_members ]) # extract & organise read files


                            # STAGE 2   01.merged --> 02.trimmed
            
                            if stage == 2:
                                
                                identifiers = [ f'{member}{identifier}' for member in relevant_members for identifier in PU ] if read_pairs_found else relevant_members # specify relevant read labels
                                out_reads = ' '.join([ f'{out_subdir}/{sample}_{identifier}.{stage_info[out_suffix]}' for identifier in identifiers ]) # specify output files
                                mode = 'PE' if read_pairs_found else 'SE'# specify relevant software mode                            
                                sh.write(f'trimmomatic {mode} -phred33 {in_reads} {out_reads} ' # use phred+33 quality scores
                                'LEADING:10 '+ # cut bases at start if quality below 10 
                                'TRAILING:10 '+ # cut bases at end if quality below 10 
                                'SLIDINGWINDOW:4:15 '+ # cut within 4-bp sliding window if average quality below 15
                                'MINLEN:50\n') # drop reads shorter than 50 bp


                            # STAGE 3   02.trimmed --> 03.mapped

                            if stage == 3:

                                aligned, organised, stats = [ f'{out_subdir}/{sample}.{suffix}' for suffix in stage_info[out_suffix] ]
                                sh.write(f'bwa mem -t10 {fasta_file} {in_reads} > {aligned}\n'+ # align 70-bp to 1-Mbp query sequences via BWA-MEM algorithm with 10 threads in paired-end mode
                                f'samtools view -bS {aligned} | samtools sort -T {sample} > {organised}\n'+ # sort alignments by chromosomal coordinates with @HD-SO tag; writes temporary files to PREFIX.nnnn.bam
                                f'samtools index {organised}\n'+ # index sorted BAM for fast random access
                                f'samtools flagstat {organised} > {stats}\n'+ # calculate BAM QC statistics
                                f'rm {aligned}\n') # delete intermediary sam


                    # STAGE 4   03.mapped --> 04.marked

                    if stage == 4:

                        marked, metrics, named = [ f'{out_subdir}/{sample}.{suffix}' for suffix in stage_info[out_suffix] ]
                        sh.write('picard -Xmx30g MarkDuplicates \\\n'+ # locate & tag duplicate reads
                        f'I={seq_file} \\\n'+ # input (BAM; sorted)
                        f'O={marked} \\\n'+ # output (BAM; marked)
                        f'M={metrics} \\\n'+ # output (duplication metrics)
                        'VALIDATION_STRINGENCY=SILENT \\\n'+ # SAM validation stringency; silenced for BAM to improve performance as variable-length data not decoded
                        'ASSUME_SORTED=true \\\n'+ # assume sorted by chromosomal coordinates
                        'REMOVE_DUPLICATES=true \n\n'+ # duplicates not written to output
                        'picard AddOrReplaceReadGroups \\\n'+ # assign reads a single read-groups
                        f'I={marked} \\\n'+ # input (BAM; marked)
                        f'O={named} \\\n'+ # output (BAM; named)
                        'SORT_ORDER=coordinate \\\n'+ # output sort order
                        'CREATE_INDEX=True \\\n'+ # index sorted BAM for fast random access (BAI)
                        'VALIDATION_STRINGENCY=LENIENT \\\n'+ # SAM validation stringency
                        f'RGLB={sample} \\\n'+ # read-group library
                        'RGPL=illumina \\\n'+ # read-group platform
                        'RGPU=AAAAAA \\\n'+ # read-group platform unit (eg. run barcode)
                        f'RGSM={sample}\n') # read-group sample name
                    

                    # STAGE 5   04.marked --> 05.haplotyped
                    
                    if stage == 5 and ploidy:

                        vcf, bam = [ f'{out_subdir}/{sample}.{ploidy}x.{suffix}' for suffix in stage_info[out_suffix] ]
                        sh.write('gatk HaplotypeCaller \\\n'+ # identify SNPs/indels relative to reference genome via local de-novo re-assembly of haplotypes in regions of variation
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-I {seq_file} \\\n'+ # input (BAM; marked)
                        f'-O {vcf} \\\n'+ # output (gVCF; raw SNPs/indel calls)
                        #f'-bamout {bam} \\\n'+ # output (bam; raw SNPs/indel calls)
                        '--emit-ref-confidence BP_RESOLUTION \\\n'+ # emit reference confidence scores site by site
                        f'--sample-ploidy {ploidy} \\\n'+ # sample ploidy (or pooled sample number)
                        '--min-base-quality-score 25 \\\n'+ # minimum base quality to be considered for calling
                        '--minimum-mapping-quality 25\n') # minimum read mapping quality to be considerd for calling (MappingQualityReadFilter)
                        # N.B. HaplotypeCaller applies NotDuplicateReadFilter & MappingQualityReadFilter by default to filter out reads marked as duplicate & with specific mapping qualities


                    # STAGE 6   05.haplotyped --> 06.combined
                    
                    if stage == 6:

                        out_file = f'{out_dir}/{stage_info[out_suffix]}'
                        sh.write('gatk CombineGVCFs \\\n'+ # combine per-sample gVCF into multi-sample gVCF
                        f'-R {fasta_file} \\\n') # reference genome (fasta)
                        for *_, (part,*_) in files_to_process: sh.write(f'-V {part} \\\n') # input (gVCF; raw SNPs/indel calls)
                        sh.write(f'-O {out_file} \n') # output (gVCF; combined)
                    
                    
                    # STAGE 7   06.combined --> 07.genotyped
                    
                    if stage == 7:

                        out_file = f'{out_subdir}/{sample}.{stage_info[out_suffix]}'
                        sh.write('gatk GenotypeGVCFs \\\n'+ # perform joint genotyping
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-V {seq_file} \\\n') # input (gVCF; combined)
                        for region in task: sh.write(f'-L {region} \\\n') # indvidual regions within group
                        sh.write(f'-O {out_file} \\\n'+ # output (gVCF; genotyped)
                        '-G StandardAnnotation \\\n'+ # annotation classes/groups applied to variant calls
                        '--include-non-variant-sites true \n') # include non-variant loci
                    

                    # STAGE 8  07.genotyped --> 08.recombined

                    if stage == 8:

                        recombined = f'{out_dir}/{stage_info[out_suffix]}'
                        sh.write('gatk GatherVcfs \\\n') # gather per-region gVCF into single gVCF
                        for *_,(part,*_) in files_to_process: sh.write(f'-I {part} \\\n') # input (gVCF; genotyped)
                        sh.write(f'-O {recombined}\n'+ # output (gVCF; gathered)
                        'gatk IndexFeatureFile \\\n'+ # index gathered gVCF for fast random access (TBI)
                        f'-F {recombined}\n') # input (gVCF; gathered)


                    # STAGE 9   08.recombined --> 09.filtered.best
                    
                    if stage == 9:
                        F1, F2, depth_per_site = [ f'{out_dir}/{suffix}' for suffix in stage_info[out_suffix] ]
                        sh.write('gatk SelectVariants \\\n'+ # select variant subset
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-V {seq_file} \\\n'+ # input (gVCF; genotyped)
                        f'-O {F1} \\\n'+ # output (gVCF; filtered F1 biallelic)
                        '--select-type-to-exclude INDEL \\\n'+ # do not select insertion/deletion variants
                        '--select-type-to-exclude MIXED \\\n'+ # do not select mixed (combination of SNPs & insertion/deletion at single position) variant
                        '--restrict-alleles-to BIALLELIC \n'+ # select only biallelic variants
                        'gatk VariantFiltration \\\n'+ # filter variants based on info/format annotations
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-V {F1} \\\n'+ # input (gVCF; filtered F1 biallelic)
                        f'-O {F2} \\\n'+ # output (gVCF; filtered F2 best practice)
                        '--filter-name "QD" \\\n'+
                        '--filter-expression "QD < 2.0" \\\n'+ # filter with QUAL score normalized by unfiltered allele depth (QUAL/AD) < 2
                        '--filter-name "FS" \\\n'+ # fishers exact test to detect forward/reverse strand allele bias i.e. artifact variant due to sequencing error
                        '--filter-expression "FS > 60.0" \\\n'+ # filter when strand bias as per phred-scaled p-value
                        '--filter-name  "MQ" \\\n'+
                        '--filter-expression "MQ < 40.0" \\\n'+ # filter with root mean square read mapping quality < 40
                        '--filter-name "MQRS" \\\n'+ # rank sum test to detect read mapping quality allele bias i.e. artifact variant due to sequencing/mapping error
                        '--filter-expression "MQRankSum < -12.5" \\\n'+ # filter when ALT mapping quality lower than REF
                        '--filter-name "RPRS" \\\n'+ # rank sum test to detect relative read position allele bias i.e. artifact variant due to sequencing error
                        '--filter-expression "ReadPosRankSum < -8.0" \\\n'+ # filter when ALT more often at read end than REF 
                        '--filter-name "HS" \\\n'+ # site consistency when strictly two segregating haplotypes i.e. artifact variant due to alignment error
                        '--filter-expression "HaplotypeScore < 13.0" \n'+ # filter when multiple segregating haplotypes (REDUNDANT?)
                        f'vcftools --gzvcf {F2} '+ # input (gVCF; filtered F2 best practice)
                        f'--out {depth_per_site} '+ # output (depth statistics)
                        '--site-depth \n') # calculate depth per site summed across all individuals
                    

                    # STAGE 10  09.filtered.best --> 10.filtered.depth (depth mask; F3: masked; F4: depth)
                    
                    if stage == 10:

                        F2, depth_per_site = sorted(task_files)
                        mask, F3, F4 = [ f'{out_dir}/{suffix}' for suffix in stage_info[out_suffix] ]
                        depth_input = open(depth_per_site, 'r').readlines()[1:]
                        depths = [ int(line.split('\t')[2]) for line in depth_input ]
                        mean_depth = statistics.mean(depths) # calculate mean depth
                        cut_off = int(round(mean_depth * 1.6)) # # calculate upper limit N.B. JEXL expression value type must match VCF annotation field (i.e. DP = int)
                        sh.write('gatk SelectVariants \\\n'+ # select variant subset
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-V {F2} \\\n'+ # input (gVCF; filtered F2 best practice)
                        f'-O {mask} \\\n'+ # output (gVCF; depth mask)
                        f'--selectExpressions "DP < {cut_off}" \n'+ # select with depth of coverage < mean x 1.6
                        'gatk VariantFiltration \\\n'+ # filter variants based on info/format annotations
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-V {F2} \\\n'+ # input (gVCF; filtered F2 best practice)
                        f'-O {F3} \\\n'+ # output (gVCF; filtered F3)
                        f'--mask {mask} \\\n'+ # input (gVCF; depth mask)
                        '--filter-not-in-mask \n'+ # filter variants not in mask
                        'gatk SelectVariants \\\n'+ # select variant subset
                        f'-R {fasta_file} \\\n'+ # reference genome (fasta)
                        f'-V {F3} \\\n'+ # input (gVCF; filtered F3)
                        f'-O {F4} \\\n'+ # input (gVCF; filtered F4)
                        '--exclude-filtered true \n') # exclude sites marked as filtered

                sh.write('echo TASK COMPLETED `date` \n') # log end time 

        if stage and stage <= 3:
            for category, samples in processing_mode.items():
                print(f'\nProcessing as {category}-end:')
                print(*[f' - {entry}' for entry in samples] if samples else [f'   NA'], sep='\n')

        if subbuddy: # companion script to manage slurm when it misbehaves; requires user & id file generated by pipeline
            buddy_name = f'buddy-{stage}'
            sh_file = f'{sh_dir}/{buddy_name}.sh'
            buddy_found = int(CAPTURE(f'squeue -u {user} | grep -c "{buddy_name}"')) > 0

            if buddy_found: print('\nSUB BUDDY ALREADY RUNNING FOR CURRENT STAGE')
            else: # submit sub buddy
                with open(f'{sh_file}', 'w') as sh: 
                    partition, nodes, ntasks, memory, walltime = resources = [  hpcc_settings['+'][label] for label in resource_labels ] # extract hpcc settings
                    hpcc_directives = SBATCH(job_id=buddy_name, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, out_err=oe_dir, email=address) # specify hpcc directives (slurm)
                    sh.write(f'{hpcc_directives}\n'+
                    'echo BUDDY CALLED `date`\n'+
                    f'python {buddy_script} -u {user} -i {id_file}\n')
                scripts.append(sh_file)

        if testing:
            print('\nIndividual Task Scripts:\n')
            for sh_file in scripts: print(open(sh_file, 'r').read()) 
            print('\nTESTING COMPLETE\n')
        else: # submit tasks to hpcc

            if submitting_self: # submit pipeline to trickle tasks from hpcc
                pipe_description = f'pipe-{relevant}'
                sh_script = f'{sh_dir}/{pipe_description}.sh'
                with open(sh_script, 'w') as sh:  
                    hpcc_directives = SBATCH(job_id=pipe_description, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, out_err=oe_dir, conda_env=environment) # specify hpcc directives (slurm)
                    sh.write(hpcc_directives)                
                    arguments.remove(pipe_flag) # remove self submission flag
                    print('python', pipe_script, *arguments, sep=' ', file=sh) # submit self as task
                pipe_id = CAPTURE(f'sbatch {sh_script}')
                print(f'\nPIPELINE SUBMITTED ({pipe_id})\n')
            
            else: # submit tasks

                *_, final_stage, final_stage_output = all_stages = list(stage_formats) # extract final stage keys
                next_stage = all_stages[all_stages.index(relevant)+1] # extract next stage key
                final_stage = stage == final_stage # check if last stage

                print('\nSUBMITTING TASKS...')                    
                with open(id_file, 'a+') as id_output:
                    for i,script in enumerate(scripts,1): 
                        ezSub(i=i, check=600, user=user, limit=limit) # maintain tasks below parellel task limit
                        sub_id = CAPTURE(f'sbatch -d singleton {script}') # submit task
                        if os.path.basename(script).startswith('buddy'): continue
                        else: print(sub_id, script, sep='\t', file=id_output, flush=True) # record task job id & shell script
                print('\tALL TASKS SUBMITTED\n')

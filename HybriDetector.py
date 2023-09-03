#execution script to initialize snakemake pipeline of HybriDetector
import argparse
import os
import sys
import pandas as pd
from snakemake.shell import shell


def main(sysargs = sys.argv[1:]):
    # Initialize parser
    parser = argparse.ArgumentParser(description='HybriDetector.py: run HybriDetector snakemake workflows', usage='python HybriDetector.py -i <input_sample1,input_sampleX> -l <read_length> -u <is_umi> -g <map_perc_single_genomic> -s <map_perc_softclip> -c <cores> -r <ram> ') 
 
 
    # Adding argments
    parser.add_argument('-i', '--input_sample', help='Fastq file names here written without the ".fastq.gz" suffix separated by comma.')
    parser.add_argument('-l', '--read_length', default=75, type=int, help='Maximal length of sequenced read within input fastq file (integer).')
    parser.add_argument('-u', '--is_umi', default=False, help='Specify whether your library contains extracted UMIs in the read header (boolean).')
    parser.add_argument('-g', '--map_perc_single_genomic', default=0.85, type=float, help='Specify the percentage when the alignment of single genomic non-chimeric reads will be reported. Only if the ratio of "aligned length / read length" is higher than or equal to this value, the alignemnt will be output.')
    parser.add_argument('-s', '--map_perc_softclip', default=0.75, type=float, help='Specify the percentage when the alignment of genomic part of the chimeric reads will be reported. Only if the ratio of "aligned length / read length" is higher than or equal to this value, the alignemnt will be output.')
    parser.add_argument('-c', '--cores', default=6, type=int, help='Number of provided CPUs.')
    parser.add_argument('-r', '--ram', default=24, type=int, help='Number of provided RAM in GB.')

    # Read arguments from command line
    args = parser.parse_args(sysargs)
    print('--------')
    print('details!')
    print('\tInput samples: {}'.format(args.input_sample))
    print('\tRead length: {}'.format(args.read_length))
    print('\tIs UMI: {}'.format(args.is_umi))
    print('\tMapped percentage for single genomic reads: {}'.format(args.map_perc_single_genomic))
    print('\tMapped percentage for the softclipped parts of hybrid reads: {}'.format(args.map_perc_softclip))
    print('\tNumber of required cores: {}'.format(args.cores))
    print('\tNumber of required RAM: {}'.format(args.ram))
    print('--------')
    config_sheet = prepare_sheet(args)
    config_file = write_config_file(config_sheet)
    run_snakemake(config_sheet)

def prepare_sheet(args):
    config_table = pd.DataFrame()
    sample_query = args.input_sample.split(",")
    #print(sample_query)
    for sample in sample_query:
        single_sample_dict = {}
        single_sample_dict['Sample'] = sample
        single_sample_dict['map_perc_single_genomic'] = args.map_perc_single_genomic
        single_sample_dict['map_perc_softclip'] = args.map_perc_softclip
        single_sample_dict['is_umi'] = str(args.is_umi).upper()
        single_sample_dict['read_length'] = args.read_length 
        single_sample_dict['cores'] = args.cores
        single_sample_dict['ram'] = args.ram
        if len(config_table) == 0:
            config_table = pd.DataFrame(single_sample_dict, index=[0])
        else:
            config_table.loc[len(config_table.index)] = single_sample_dict
    return(config_table)

def write_config_file(config_sheet):
    work_dir = os.getcwd()
    f = open(f'{work_dir}/config.json', "w")
    f.write(f'{{')
    for i, (key,value) in enumerate(config_sheet.items()):
        if i == len(config_sheet.keys())-1:
            f.write(f""""{key}":[{','.join(f'"{str(i)}"' for i in value.tolist())}]\n""")
        else:    
            f.write(f""""{key}":[{','.join(f'"{str(i)}"' for i in value.tolist())}],\n """)
    f.write(f'}}')
    f.close()

def run_snakemake(config_sheet):
    work_dir = os.getcwd()
    snakefile = f'{work_dir}/HybriDetector.smk'
    directory = work_dir
    configfile = f'{work_dir}/config.json'
    ram = config_sheet["ram"].min()
    cpu = config_sheet["cores"].min()

    command = "snakemake --snakefile " + snakefile + " --directory " + directory + " --configfile " + configfile + " --use-conda --conda-frontend mamba -p --res mem=" + str(ram) + " --jobs " + str(cpu)
    print("Executed command:\n")
    print(command)
    print("\nTask is being processed!")

    shell(command)

    #snakemake --snakefile HybriDetector.smk --directory /home/vasek/Documents/CLASH_snakemake_github/HybriDetector --configfile configuration_file.json  --use-conda --conda-frontend mamba -p --res mem=100 -n -r
    #status = snakemake.snakemake(snakefile, configfile=paramsfile,
    #                             targets=[target], printshellcmds=True,
    #                             dryrun=args.dry_run, forceall=args.force,
    #                             config=config)
    #
    #if status: # translate "success" into shell exit code of 0
    #   return 0
    #return 1
if __name__ == '__main__':
    main()  


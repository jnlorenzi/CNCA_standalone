import os
import sys
import argparse
from dotenv import dotenv_values
from pathlib import Path
import re
import subprocess
import datetime
import pandas


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='command line example: python3 ./pipeline_generator.py -i ./config_server.env')

    # Positional mandatory arguments
    parser.add_argument('-i', '--config_file_path', help="absolute or relative path to the config file", type=str)
    parser.add_argument('-t', '--temp_dir_path', help="absolute or relative path to the temp dir", type=str)

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.00 - 02/08/2023')

    # Parse arguments
    args = parser.parse_args()

    return args


def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]


def edit_path_for_windows(old_path):
    handler = old_path.split('\\')
    new_path = '/mnt/c/' + '/'.join(handler[1::])
    return new_path


def gbk_parsing_call_maker(config, temp_dir):
    for key in config:
        if key == 'genome_path':
            genome_path = config['genome_path']
    try:
        cgi_bin_path = config['cgi_bin_path']
        gbk_parsing_call_string = f'python3 {cgi_bin_path}/gbk_parsing.py ' \
                                  f'-i {temp_dir}/{genome_path} ' \
                                  f'-o {temp_dir}/{genome_path} ' \
                                  f'-f {temp_dir}/{genome_path}'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing')

    gbk_parsing_call_string = re.sub(r'-[a-z,A-Z]* {2,}', '', gbk_parsing_call_string)
    gbk_parsing_call_string = re.sub(r'/{2,}', '/', gbk_parsing_call_string)

    return gbk_parsing_call_string


def aligner_prep_call_maker(config, temp_dir, data_type):
    try:
        cgi_bin_path = config['cgi_bin_path']
        aligner_prep_call_string = f'python3 {cgi_bin_path}/aligner_preparator.py'
        for key in config:
            if key == 'genome_path':
                genome_path = config['genome_path']
                aligner_prep_call_string = f'{aligner_prep_call_string} -d {temp_dir}/{genome_path}'
            if key == 'infile_align':
                infile_align = config['infile_align']
                aligner_prep_call_string = f'{aligner_prep_call_string} -i {temp_dir}/{infile_align}/{data_type}'
            if key == 'outfile_align':
                outfile_align = config['outfile_align']
                aligner_prep_call_string = f'{aligner_prep_call_string} -o {temp_dir}/{outfile_align}/{data_type}'
            if key == 'sequence_to_use_path':
                sequence_to_use_path = config['sequence_to_use_path']
                if sequence_to_use_path:
                    aligner_prep_call_string = f'{aligner_prep_call_string} -s {temp_dir}/{sequence_to_use_path}'
            if key == 'aligner_method':
                aligner_method = config['aligner_method']
                aligner_prep_call_string = f'{aligner_prep_call_string} -m {aligner_method}'
            if key == 'tail_head':
                tail_head = config['tail_head']
                if tail_head:
                    aligner_prep_call_string = f'{aligner_prep_call_string} -th {tail_head}'
        aligner_prep_call_string = f'{aligner_prep_call_string} -t {data_type}'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing')

    aligner_prep_call_string = re.sub(r'-[a-z,A-Z]* {2,}', '', aligner_prep_call_string)
    aligner_prep_call_string = re.sub(r'/{2,}', '/', aligner_prep_call_string)

    return aligner_prep_call_string


def aligner_command_maker(config, temp_dir, parameters, data_type):
    try:
        if config['aligner_method'] == 'mafft':
            # retrieve mafft parameters
            strategy = {
                'auto': 'mafft --auto',
                'L-INS-i': 'mafft --localpair --maxiterate 1000',
                'G-INS-i': 'mafft --globalpair --maxiterate 1000',
                'E-INS-i': 'mafft --ep 0 --genafpair --maxiterate 1000',
                'FFT-NS-i': 'mafft --retree 2 --maxiterate 1000',
                'FFT-NS-1': 'mafft --retree 1 --maxiterate 0',
                'FFT-NS-2': 'mafft --retree 2 --maxiterate 2'
            }

            dist_matrix = {
                'aa': {
                    'bl 30': '--amino --bl 30',
                    'bl 45': '--amino --bl 45',
                    'bl 62': '--amino --bl 62',
                    'bl 80': '--amino --bl 80',
                    'jtt 100': '--amino --jtt 100',
                    'jtt 200': '--amino --jtt 200'
                },
                'nuc': {
                    'kimura 1': '--nuc --jtt 1',
                    'kimura 20': '--nuc --jtt 20',
                    'kimura 200': '--nuc --jtt 200'
                }
            }

            if data_type == "nuc":
                score_matrix_type = "scorematrix_nt"
            else:
                score_matrix_type = "scorematrix_aa"

            # aligner_prep_call_string = f'mafft --quiet --auto --clustalout --thread -1'
            aligner_prep_call_string = f'{strategy[parameters["strategy"]]} ' \
                                       f'--op {parameters["gapopeningpenalty"]} ' \
                                       f'--ep {parameters["offsetvalue"]} ' \
                                       f'{dist_matrix[data_type][parameters[score_matrix_type]]} ' \
                                       f'--clustalout --thread -1'
            infile_mafft = config['infile_align'] + '/' + data_type + '/'
            outfile_mafft = config['outfile_align'] + '/' + data_type + '/'
            if config['windows'] == 'True':
                aligner_prep_call_string = f'wsl.exe {aligner_prep_call_string} '
                aligner_prep_call_string = f'{aligner_prep_call_string} {edit_path_for_windows(temp_dir)}/' \
                                           f'{infile_mafft}/multiple.fasta > {temp_dir}/{outfile_mafft}/multiple.aln'
            else:
                aligner_prep_call_string = f'{aligner_prep_call_string} {temp_dir}/{infile_mafft}/multiple.fasta > ' \
                                           f'{temp_dir}/{outfile_mafft}/multiple.aln'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing')

    return aligner_prep_call_string


def crossing_call_maker(config, temp_dir):
    try:
        cgi_bin_path = config['cgi_bin_path']
        crossing_call_string = f'python3 {cgi_bin_path}/crossing_align_annot.py '
        for key in config:
            if key == 'genome_path':
                genome_path = config['genome_path']
                crossing_call_string = f'{crossing_call_string} -c {temp_dir}/{genome_path}'
                crossing_call_string = f'{crossing_call_string} -o {temp_dir}/{genome_path}'
            if key == 'outfile_align':
                outfile_align = config['outfile_align'] + '/nuc/'
                crossing_call_string = f'{crossing_call_string} -i {temp_dir}/{outfile_align}'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing in crossing_call_maker()')

    return crossing_call_string


def alignment_patching_call_maker(config, temp_dir):
    try:
        cgi_bin_path = config['cgi_bin_path']
        if config['windows'] == 'True':
            alignment_patching_call_string = f'Rscript {cgi_bin_path}/nuc_aa_match_table.R'
        else:
            alignment_patching_call_string = f'Rscript {cgi_bin_path}/nuc_aa_match_table.R'
        align_path = '/'.join(config['infile_align'].split('/')[0:-2])
        json_path = config['genome_path']
        method = config['aligner_method']

        alignment_patching_call_string = f'{alignment_patching_call_string}' \
                                         f' -i {temp_dir}/{align_path}' \
                                         f' -n nuc' \
                                         f' -a aa' \
                                         f' -j {temp_dir}/{json_path}' \
                                         f' -m {method}'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing in alignment_patching_call_maker()')

    return alignment_patching_call_string


def patching_call_maker(config, temp_dir):
    try:
        cgi_bin_path = config['cgi_bin_path']
        if config['windows'] == 'True':
            patching_call_string = f'Rscript {cgi_bin_path}/patching.R'
        else:
            patching_call_string = f'Rscript {cgi_bin_path}/patching.R'
        align_path = '/'.join(config['infile_align'].split('/')[0:-2])

        patching_call_string = f'{patching_call_string}' \
                               f' -i {temp_dir}/{align_path}' \
                               f' -n nuc' \
                               f' -a aa'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing in patching_call_maker()')

    return patching_call_string


def check_call_maker(config, temp_dir):
    try:
        cgi_bin_path = config['cgi_bin_path']
        if config['windows'] == 'True':
            check_call_string = f'Rscript {cgi_bin_path}/check_patching.R'
        else:
            check_call_string = f'Rscript {cgi_bin_path}/check_patching.R'
        align_path = '/'.join(config['infile_align'].split('/')[0:-2])
        method = config['aligner_method']

        check_call_string = f'{check_call_string}' \
                            f' -i {temp_dir}/{align_path}' \
                            f' -n nuc' \
                            f' -a aa' \
                            f' -m {method}'
    except ValueError:
        sys.exit('Pb during the config_server.env file parsing in check_call_maker()')

    return check_call_string


def main():
    args = parseArguments()

    config_file_path = args.config_file_path
    temp_dir = args.temp_dir_path

    dotenv_path = Path(config_file_path)
    config = dotenv_values(dotenv_path)

    parameters_path = Path(f'{temp_dir}/parameters.csv')
    with open(parameters_path, 'r') as infile_parameters:
        parameters = {line.split('=')[0]: line.split('=')[1].rstrip() for line in infile_parameters.readlines()}

    handler_call_string = []

    # parse the user data
    gbk_parsing_call_string = gbk_parsing_call_maker(config, temp_dir)
    handler_call_string.append(gbk_parsing_call_string)

    aligner_prep_call_string_nuc = aligner_prep_call_maker(config, temp_dir, data_type='nuc')
    handler_call_string.append(aligner_prep_call_string_nuc)

    aligner_call_string_nuc = aligner_command_maker(config, temp_dir, parameters, data_type='nuc')
    handler_call_string.append(aligner_call_string_nuc)

    crossing_call_string = crossing_call_maker(config, temp_dir)
    handler_call_string.append(crossing_call_string)

    aligner_prep_call_string_aa = aligner_prep_call_maker(config, temp_dir, data_type='aa')
    handler_call_string.append(aligner_prep_call_string_aa)

    aligner_call_string_aa = aligner_command_maker(config, temp_dir, parameters, data_type='aa')
    handler_call_string.append(aligner_call_string_aa)

    alignment_patching_call_string = alignment_patching_call_maker(config, temp_dir)
    handler_call_string.append(alignment_patching_call_string)

    patching_call_string = patching_call_maker(config, temp_dir)
    handler_call_string.append(patching_call_string)

    check_call_string = check_call_maker(config, temp_dir)
    handler_call_string.append(check_call_string)

    with open(f'{temp_dir}/run.sh', 'a') as outfile:
        outfile.write(' && '.join(handler_call_string))

    progress_file = f'{temp_dir}/progress.txt'

    progress = {
        'data loading': 'pending',
        'nucleotide alignment': 'pending',
        'protein alignment': 'pending',
        'patching': 'pending'
    }
    with open(f'{temp_dir}/log.txt', 'w') as log_file:
        # pipeline calling
        log_file.write(f'### START --- {datetime.datetime.now()}  --- gbk_parsing.py\n')
        progress['data loading'] = 'ongoing'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)
        subprocess.run(gbk_parsing_call_string, shell=True, stdout=log_file)
        log_file.write(f'### END --- {datetime.datetime.now()}  --- gbk_parsing.py\n\n')
        if os.path.exists(f'{temp_dir}/genome/gbk_parsing_sumup.csv'):
            with open(f'{temp_dir}/genome/gbk_parsing_sumup.csv', 'r') as infile:
                if infile.readline().rstrip() == 'multiple':
                    progress['data loading'] = 'error'
                    progress['nucleotide alignment'] = 'halted'
                    progress['protein alignment'] = 'halted'
                    progress['patching'] = 'halted'
                    pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)
                    sys.exit()
        else:
            sys.exit()
        progress['data loading'] = 'done'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)

        log_file.write(f'### START --- {datetime.datetime.now()}  --- prepare nuc alignment\n')
        progress['nucleotide alignment'] = 'ongoing'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)
        subprocess.run(aligner_prep_call_string_nuc, shell=True, stdout=log_file)
        log_file.write(f'### END --- {datetime.datetime.now()}  --- prepare nuc alignment\n\n')

        log_file.write(f'### START --- {datetime.datetime.now()}  --- run nuc alignment\n')
        subprocess.run(aligner_call_string_nuc, shell=True, stdout=log_file)
        log_file.write(f'### END --- {datetime.datetime.now()}  --- run nuc alignment\n\n')
        progress['nucleotide alignment'] = 'done'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)

        log_file.write(f'### START --- {datetime.datetime.now()}  --- prepare aa alignment\n')
        progress['protein alignment'] = 'ongoing'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)
        subprocess.run(aligner_prep_call_string_aa, shell=True, stdout=log_file)
        log_file.write(f'### END --- {datetime.datetime.now()}  --- prepare aa alignment\n\n')

        log_file.write(f'### START --- {datetime.datetime.now()}  --- run aa alignment\n')
        subprocess.run(aligner_call_string_aa, shell=True, stdout=log_file)
        log_file.write(f'### END --- {datetime.datetime.now()}  --- run aa alignment\n\n')
        progress['protein alignment'] = 'done'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)

        log_file.write(f'### START --- {datetime.datetime.now()}  --- run patching process\n')
        progress['patching'] = 'ongoing'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)
        subprocess.run(crossing_call_string, shell=True, stdout=log_file)
        subprocess.run(alignment_patching_call_string, shell=True, stdout=log_file)
        subprocess.run(patching_call_string, shell=True, stdout=log_file)
        subprocess.run(check_call_string, shell=True, stdout=log_file)
        log_file.write(f'### END --- {datetime.datetime.now()}  --- run patching process\n\n')
        progress['patching'] = 'done'
        pandas.DataFrame(progress, index=[0]).to_csv(progress_file, index=False)
        log_file.write(f'### ALL DONE ###\n')


if __name__ == "__main__":
    main()

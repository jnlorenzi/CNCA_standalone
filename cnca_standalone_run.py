import os
import sys
import argparse
import shutil
import subprocess


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Execute the CNCA pipeline with the default parameters, see https://cnca.ijm.fr/ for other options')

    # Positional mandatory arguments
    parser.add_argument('-i', '--gbk_path', help="absolute or relative path to the genbank files", type=str,
                        required=True)
    parser.add_argument('-w', '--working_dir', help="absolute or relative path to the working directory", type=str,
                        required=True)

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.00')
    # Parse arguments
    args = parser.parse_args()

    return args


def main():
    args = parseArguments()
    gbk_path = args.gbk_path
    working_dir = args.working_dir

    # create working directory, stop process if already exist (dot not overwrite previous run)
    if os.path.exists(working_dir):
        sys.exit(f'Error: the working directory "{working_dir}" already exist')
    else:
        os.makedirs(working_dir)

    count_gbk_file = 0
    for gbk_file in os.listdir(gbk_path):
        # copy file in working directory
        genome_dir = f'{working_dir}/genome/'
        os.makedirs(genome_dir, exist_ok=True)
        if os.path.splitext(gbk_file)[1] == ".gbk" and not gbk_file.startswith('._'):
            shutil.copy(f'{gbk_path}/{gbk_file}', genome_dir)
            count_gbk_file += 1
    if count_gbk_file == 0:
        sys.exit(f'No genbank file (.gbk) in the given directory: "{gbk_path}"')

    base_path = os.path.dirname(os.path.abspath(__file__))
    # copy mafft parameter file to the working directory
    shutil.copy(f'{base_path}/parameters.csv', working_dir)

    call_pipeline_string = f'python3 {base_path}/pipeline_generator.py ' \
                           f'-i {base_path}/config.env ' \
                           f'-t {working_dir}'

    subprocess.run(call_pipeline_string, shell=True)

    # copy final file to root of the working directory
    for res_file in os.listdir(f'{working_dir}/mafft/'):
        if os.path.splitext(res_file)[1] == ".nexus" or not os.path.splitext(res_file)[1] == '.fasta':
            shutil.copy(f'{working_dir}/mafft/{res_file}', working_dir)


if __name__ == "__main__":
    main()


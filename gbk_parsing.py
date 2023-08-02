import argparse, os, sys
import shutil
import tempfile
import Bio
from Bio import SeqIO
import json
import pandas


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='command line example: python gbk_parsing.py ~/Documents/SARS-cov2/table/genome/SARS-cov2-vs-nr_1_gl3_236-38/ ~/Documents/SARS-cov2/table/cds/SARS-cov2-vs-nr_1_gl3_236-38/')

    # Positional mandatory arguments
    parser.add_argument('-i', '--gbk_path',
                        help='absolute or relative path to the gbk files (generated by retrieve_genome.py)',
                        type=str)
    parser.add_argument('-o', '--path_cds_outfile', help='absolute or relative path to the cds json file (created)',
                        type=str)
    parser.add_argument('-f', '--path_fasta_outfile', help='absolute or relative path to the fasta file (created)',
                        type=str)

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.5 - 30.08.2022')

    # Parse arguments
    args = parser.parse_args()

    return args


def get_num(x):
    return int(''.join(ele for ele in x if ele.isdigit()))


def format_position(position):
    start = get_num(position.split(':')[0].split('[')[-1])
    end = get_num(position.split(':')[-1].split(']')[0])
    return start, end


def adjust_conding_sequence(sequence):
    while len(sequence) % 3 != 0:
        sequence = sequence + 'N'
    return sequence


def get_all_cds(seq_record):
    """
  Function to extract all the CDS and their features
  """
    all_cds = {'cds': {}}
    c = 2
    multiple_annotation = False
    annotated_range = {}
    for feature in seq_record.features:
        rna_seq = str(seq_record.seq)
        all_cds['complete_sequence'] = rna_seq
        # Select on CDS feature
        if feature.type == 'CDS':
            # skip pseudogene
            if 'pseudo' in list(feature.qualifiers.keys()):
                continue
            position = str(feature.location)
            if not 'product' in list(feature.qualifiers):
                if 'gene' in list(feature.qualifiers):
                    cds_id = feature.qualifiers['gene'][0]
                elif 'protein_id' in list(feature.qualifiers):
                    cds_id = feature.qualifiers['protein_id'][0]
                else:
                    sys.exit('Impossible to set a cds id')
            else:
                cds_id = feature.qualifiers['product'][0]

            # If there is a join position, skip it on CNCA servor
            if position.startswith('join'):
                if len(position.split(',')) == 2:
                    pos_a, pos_b = position.split(',')
                    pos_a_start, pos_a_end = format_position(pos_a)
                    pos_b_start, pos_b_end = format_position(pos_b)
                    if abs(pos_a_end - pos_b_start) > 1:
                        sys.exit('Join situation different from ORF1ab, abort')
                    a_seq = rna_seq[pos_a_start: pos_a_end].lower().replace('u', 't')
                    b_seq = rna_seq[pos_b_start: pos_b_end].lower().replace('u', 't')
                    prot_a_start = 0
                    prot_a_end = int(((pos_a_end + 1) / 3) - ((pos_a_start + 1) / 3))
                    if not 'translation' in list(feature.qualifiers):
                        a_aa_seq = Bio.Seq.translate(a_seq, to_stop='False')
                        b_aa_seq = Bio.Seq.translate(b_seq, to_stop='False')
                    else:
                        aa_seq = feature.qualifiers['translation'][0]
                        a_aa_seq = aa_seq[prot_a_start:prot_a_end]
                        aa_b_start = int(((pos_b_start + 1) / 3) - ((pos_a_start + 1) / 3)) + 1
                        b_aa_seq = aa_seq[aa_b_start::]
                    all_cds['cds'][cds_id + '_1'] = {}
                    all_cds['cds'][cds_id + '_1']['start'] = pos_a_start
                    all_cds['cds'][cds_id + '_1']['end'] = pos_a_end
                    all_cds['cds'][cds_id + '_1']['seq'] = a_seq
                    all_cds['cds'][cds_id + '_1']['aa_seq'] = a_aa_seq
                    all_cds['cds'][cds_id + '_2'] = {}
                    all_cds['cds'][cds_id + '_2']['start'] = pos_b_start
                    all_cds['cds'][cds_id + '_2']['end'] = pos_b_end
                    all_cds['cds'][cds_id + '_2']['seq'] = b_seq
                    all_cds['cds'][cds_id + '_2']['aa_seq'] = b_aa_seq

                    for x in range(pos_a_start, pos_b_end + 1):
                        if x in annotated_range:
                            multiple_annotation = True
                        else:
                            annotated_range.update({x: True})
                else:
                    print(cds_id)
                    sys.exit('join involving more than 2 join CDS, this situation is not managed, abort')
            # No join position (normal case)
            else:
                if cds_id in all_cds:
                    if cds_id.lower() == 'hypothetical protein':
                        cds_id = 'hypothetical protein ' + str(c)
                        c += 1
                    else:
                        sys.exit(cds_id + " already exist, unanticipated case")
                pos_start, pos_end = format_position(position)
                cds_seq = adjust_conding_sequence(rna_seq[pos_start: pos_end]).lower().replace('u', 't')
                all_cds['cds'][cds_id] = {}
                all_cds['cds'][cds_id]['start'] = pos_start
                all_cds['cds'][cds_id]['end'] = pos_end
                all_cds['cds'][cds_id]['seq'] = cds_seq
                if not 'translation' in list(feature.qualifiers):
                    all_cds['cds'][cds_id]['aa_seq'] = Bio.Seq.translate(cds_seq, to_stop='False')
                else:
                    all_cds['cds'][cds_id]['aa_seq'] = feature.qualifiers['translation'][0]
                for x in range(pos_start, pos_end + 1):
                    if x in annotated_range:
                        # print(x, cds_id, pos_start)
                        multiple_annotation = True
                    else:
                        annotated_range.update({x: True})
    return all_cds, multiple_annotation


def insert_line_front(insert_filename, to_insert):
    with open(insert_filename) as src, tempfile.NamedTemporaryFile(
            'w', dir=os.path.dirname(insert_filename), delete=False) as dst:
        # Discard first line
        src.readline()

        # Save the new first line
        dst.write(to_insert + '\n')

        # Copy the rest of the file
        shutil.copyfileobj(src, dst)

    # remove old version
    os.unlink(insert_filename)

    # rename new version
    os.rename(dst.name, insert_filename)

    return ()


def main():
    args = parseArguments()
    gbk_path = args.gbk_path
    path_cds_outfile = args.path_cds_outfile
    path_fasta_outfile = args.path_fasta_outfile

    os.makedirs(os.path.split(path_cds_outfile)[0], exist_ok=True)

    csv_data = {'Organism': [], 'Nb CDS': [], 'CDS name (ordered)': []}

    for gbk_file in os.listdir(gbk_path):
        if os.path.splitext(gbk_file)[1] == ".gbk" and not gbk_file.startswith('._'):
            genome_record = SeqIO.read(gbk_path + '/' + gbk_file, "genbank")
            all_cds, multiple_annotation = get_all_cds(genome_record)
            if multiple_annotation:
                break
            organism = os.path.splitext(os.path.basename(gbk_file))[0]
            with open(path_cds_outfile + '/' + organism + '.txt',
                      'w') as outfile:
                json.dump(all_cds, outfile)
            # create fasta file
            fasta_file = f'{path_fasta_outfile}/{organism}.fa'
            SeqIO.convert(f'{gbk_path}/{gbk_file}', 'genbank', fasta_file, 'fasta')
            # change first line
            insert_line_front(fasta_file, '>' + os.path.splitext(os.path.basename(gbk_file))[0])
            # update csv_data
            csv_data['Organism'].append(organism)
            csv_data['Nb CDS'].append(len(all_cds['cds']))
            temp_cds_name = [(x, all_cds['cds'][x]['start']) for x in all_cds['cds']]
            temp_cds_name.sort(key=lambda x: x[1])
            csv_data['CDS name (ordered)'].append(' - '.join([x[0] for x in temp_cds_name]))
    sumup_file = f'{path_fasta_outfile}/gbk_parsing_sumup.csv'
    if multiple_annotation:
        with open(sumup_file, 'w') as infile:
            infile.write('multiple')
    else:
        df = pandas.DataFrame(csv_data)
        df.to_csv(sumup_file, index=False)


if __name__ == "__main__":
    main()
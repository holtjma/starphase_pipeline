'''
Utility script that will auto-populate a TSV file based on the input cohorts. 
This is later parsed in our collation script to aggregate the data.
'''

import argparse as ap
import csv

from PipelineConfig import PIPELINE_FOLDER
from SharedParsing import load_batches

if __name__ == '__main__':
    description = 'Creates a TSV file for data aggregation using expected pipeline locations'
    p = ap.ArgumentParser(description=description, formatter_class=ap.RawTextHelpFormatter)

    # this is really all we need, all inputs are from the PipelineConfig
    p.add_argument('-o', '--output-tsv', dest='output_tsv', required=True, help='the output file path (TSV)')
    
    args = p.parse_args()

    # load the batch information
    batch_data = load_batches()

    # open the TSV file, which just needs the prefixes
    fp = open(args.output_tsv, 'w+')
    tsv_writer = csv.DictWriter(fp, delimiter='\t', fieldnames=[
        'unique_id', 'starphase_prefix', 'peddy_prefix'
    ])
    tsv_writer.writeheader()

    # add a row for each entry
    for sample_id in batch_data:
        row = {
            'unique_id' : sample_id,
            'starphase_prefix' : f'{PIPELINE_FOLDER}/starphase/{sample_id}',
            'peddy_prefix' : f'{PIPELINE_FOLDER}/peddy/{sample_id}'
        }
        tsv_writer.writerow(row)

    fp.close()
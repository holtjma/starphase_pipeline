
import csv

from PipelineConfig import *

# if we need to ignore any datasets, this is a simple work-around while keeping our input files the same
IGNORE_SET = set([
    # 'unique_id' column would go here
])

def load_batches():
    '''
    This will load the batch information
    '''
    # return is a key -> row dictionary
    ret = {}
    for bf in BATCH_FILES:
        print(f"Parsing {bf}...")
        fp = open(bf, 'r')
        tsv_reader = csv.DictReader(fp, delimiter='\t')
        for d in tsv_reader:
            unique_id = d['unique_id']
            if unique_id in IGNORE_SET:
                continue
            if unique_id in ret:
                if ret[unique_id] == d:
                    # it's exactly the same, so it's okay
                    pass
                else:
                    raise Exception(f'Unique ID collision: {unique_id}')
            ret[unique_id] = d
        fp.close()
    return ret

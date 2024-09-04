'''
Loads all the data and aggregates it, reporting only ancestry level summaries.
'''

import argparse as ap
import csv
import json

#####################################################
# Config
#####################################################

# controls whether we collect the HLA delta metrics, which loads a semi-large JSON
ENABLE_HLA_DELTAS = True
MIN_COVERAGE = 15.0

#####################################################
# Data loading
#####################################################

def loadInputTsv(fn):
    '''
    Loads the collation file so we can prepare to parse everything
    '''
    fp = open(fn, 'r')
    tsv_reader = csv.DictReader(fp, delimiter='\t')
    ret = {}
    for row in tsv_reader:
        sample_id = row['unique_id']
        if sample_id in ret:
            raise Exception(f'Duplicate unique_id detected: {sample_id}')
        ret[sample_id] = row
    fp.close()
    return ret

def loadCohortData(cohort_datasets):
    '''
    This will load all of the sample data for the cohort
    '''
    ret = {}
    for k in cohort_datasets:
        print(f'\tLoading data for {k}...')

        # first, load the coverage to determine if this dataset is getting filtered
        mosdepth_prefix = cohort_datasets[k]['mosdepth_prefix']
        mosdepth_summary_fn = f'{mosdepth_prefix}.mosdepth.summary.txt'
        mean_coverage = loadMeanCoverage(mosdepth_summary_fn)
        if mean_coverage < MIN_COVERAGE:
            print(f'\t\tSkipping {k} => mean_coverage = {mean_coverage}')
            continue
        else:
            print(f'\t\tmean_coverage = {mean_coverage}')

        # load the sex for the sample
        peddy_prefix = cohort_datasets[k]['peddy_prefix']
        peddy_sex_fn = f'{peddy_prefix}.sex_check.csv'
        sex = loadPeddySex(peddy_sex_fn)
        print(f'\t\tPredicted sex = {sex}')

        # load the ancestry for the sample
        peddy_anc_fn = f'{peddy_prefix}.het_check.csv'
        (ancestry, anc_prob) = loadPeddyAncestry(peddy_anc_fn)
        print(f'\t\tPredicted ancestry = {ancestry} ({anc_prob})')

        # load the StarPhase results
        starphase_prefix = cohort_datasets[k]['starphase_prefix']
        starphase_fn = f'{starphase_prefix}.diplotypes.json'
        diplotypes = loadStarphase(starphase_fn)
        if ENABLE_HLA_DELTAS:
            starphase_debug_fn = f'{starphase_prefix}_debug/hla_debug.json'
            hla_deltas = loadHlaDeltas(starphase_debug_fn, diplotypes)
        else:
            hla_deltas = {}

        # add it to the result
        ret[k] = {
            'sex' : sex,
            'ancestry' : ancestry,
            'ancestry_prob' : anc_prob,
            'diplotypes' : diplotypes,
            'hla_deltas' : hla_deltas
        }
    
    return ret

def loadMeanCoverage(fn):
    '''
    File format:
    chrom	length	bases	mean	min	max
    chr1	248956422	8709049892	34.98	0	21554
    ...
    total	3099813762	98832345199	31.88	0	35443
    '''
    ret = None
    fp = open(fn, 'r')
    tsv_reader = csv.DictReader(fp, delimiter='\t')
    for row in tsv_reader:
        if row['chrom'] == 'total':
            ret = float(row['mean'])
            break

    fp.close()
    return ret

def loadPeddySex(fn):
    '''
    File format:
    sample_id,ped_sex,hom_ref_count,het_count,hom_alt_count,het_ratio,predicted_sex,error
    {sample},unknown,0,7107,2856,2.488,female,True
    '''
    fp = open(fn, 'r')
    csv_reader = csv.DictReader(fp)
    sex = None
    for row in csv_reader:
        assert(sex == None) # verifies we only iterate on one row
        sex = row['predicted_sex']
    fp.close()
    return sex

def loadPeddyAncestry(fn):
    '''
    File format:
    sample_id,depth_outlier,het_count,het_ratio,idr_baf,mean_depth,median_depth,p10,p90,sampled_sites,call_rate,ancestry-prediction,ancestry-prob,PC1,PC2,PC3,PC4
    {sample},False,2941,0.5762,0.2546,27.71,27,0.375,0.6296,5088,1,EAS,0.9919,1.478,1.868,-0.134,-0.05689
    '''
    fp = open(fn, 'r')
    csv_reader = csv.DictReader(fp)
    ancestry = None
    anc_prob = 0.0
    for row in csv_reader:
        assert(ancestry == None) # verifies we only iterate on one row
        ancestry = row['ancestry-prediction']
        anc_prob = float(row['ancestry-prob'])
    fp.close()
    return (ancestry, anc_prob)

def loadStarphase(data_fn):
    '''
    Loads all the starphase diplotype results
    '''
    # load the dict
    fp = open(data_fn, 'r')
    j = json.load(fp)
    fp.close()

    # parse the genes, mostly copy paste
    gene_results = {}
    for gene in j['gene_details']:
        diplotypes = j['gene_details'][gene]['diplotypes']
        gene_results[gene] = {'diplotypes' : diplotypes}
    
    # tie it to the sample
    return gene_results

def loadHlaDeltas(data_fn, diplotypes):
    '''
    This will load the match information for each haplotype so we can generate some metrics on identity to DB
    '''
    # load the dict
    fp = open(data_fn, 'r')
    j = json.load(fp)
    fp.close()

    # parse the genes, mostly copy paste
    genes = ['HLA-A', 'HLA-B']
    ret = {}
    for gene in genes:
        # get the haplotypes and strip the '*' at the start
        dips = diplotypes[gene]['diplotypes']
        assert(len(dips) == 1)
        hap1 = dips[0]['hap1'][1:]
        hap2 = dips[0]['hap2'][1:]

        # make sure they match what is in the debug file
        con1_dict = j['read_mapping_stats'][gene]['consensus1']
        if hap1 == hap2 and 'consensus2' not in j['read_mapping_stats'][gene]:
            # homozygous
            con2_dict = con1_dict
        else:
            con2_dict = j['read_mapping_stats'][gene]['consensus2']
        
        if hap1 == hap2 and con1_dict['best_match_star'] != hap1:
            # con1 failed the CDF or MAF test, copy con2 into con1
            con1_dict = con2_dict
        
        if hap1 == hap2 and con2_dict['best_match_star'] != hap2:
            # con2 failed the CDF or MAF test, copy con1 into con2
            con2_dict = con1_dict

        # verify these match after all the checks above
        assert(con1_dict['best_match_star'] == hap1)
        assert(con2_dict['best_match_star'] == hap2)

        # get the mapping differences, preferring DNA over cDNA
        delta_dict = {}
        iter_set = [
            ('hap1_dna_delta', con1_dict['mapping_stats'][hap1]['dna_mapping']),
            ('hap1_cdna_delta', con1_dict['mapping_stats'][hap1]['cdna_mapping']),
            ('hap2_dna_delta', con2_dict['mapping_stats'][hap2]['dna_mapping']),
            ('hap2_cdna_delta', con2_dict['mapping_stats'][hap2]['cdna_mapping']),
        ]
        for (k, d) in iter_set:
            if d == None:
                delta_dict[k] = None
            else:
                delta_dict[k] = d['nm'] + d['unmapped']

        # TODO: length is stored in 'query_len' if we need ratios
        ret[gene] = delta_dict
    
    # tie it to the sample
    return ret

#####################################################
# Data aggregation
#####################################################

def collectByAncestry(all_data):
    '''
    Go through and reorganize by (gene, ancestry, sex, diplotype) -> count
    '''
    # (gene, ancestry, sex, sorted_diplotype/value) -> count
    ancestry_counts = {}

    for sample_id in all_data:
        # pull out the relevant info we might filter on
        sex = all_data[sample_id]['sex']
        ancestry = all_data[sample_id]['ancestry']
        anc_prob = all_data[sample_id]['ancestry_prob']
        diplotypes = all_data[sample_id]['diplotypes']
        hla_deltas = all_data[sample_id]['hla_deltas']

        # TODO: any ancestry filtering? Seems like there is a built-in for Unknown, so we will ignore for now

        for gene in diplotypes:
            dip_calls = diplotypes[gene]['diplotypes']
            num_diplotypes = len(dip_calls)

            if num_diplotypes > 1:
                print(f'Skipping {sample_id} + {gene} with {num_diplotypes} diplotypes')
                continue
            
            # get the haplotypes
            hap1 = dip_calls[0]['hap1']
            hap2 = dip_calls[0]['hap2']

            # sort it just to get a consistent order
            ordered_diplotype = '/'.join(sorted([hap1, hap2]))

            # add one for this observation
            data_key = (gene, ancestry, sex, ordered_diplotype)
            ancestry_counts[data_key] = ancestry_counts.get(data_key, 0) + 1

            # check if this is an HLA gene with delta values
            if gene in hla_deltas:
                # first do cDNA
                # can be an int or None, so convert to str for sorting
                cdna_deltas = sorted([
                    str(hla_deltas[gene]['hap1_cdna_delta']),
                    str(hla_deltas[gene]['hap2_cdna_delta'])
                ])
                c_delta_dip = '/'.join(cdna_deltas)
                data_key = (f'{gene}_cdna_delta', ancestry, sex, c_delta_dip)
                ancestry_counts[data_key] = ancestry_counts.get(data_key, 0) + 1

                # now do DNA
                dna_deltas = sorted([
                    str(hla_deltas[gene]['hap1_dna_delta']),
                    str(hla_deltas[gene]['hap2_dna_delta'])
                ])
                delta_dip = '/'.join(dna_deltas)
                data_key = (f'{gene}_dna_delta', ancestry, sex, delta_dip)
                ancestry_counts[data_key] = ancestry_counts.get(data_key, 0) + 1
    
    return ancestry_counts

def writeAggregateTsv(aggregate_data, output_fn):
    '''
    Converts the aggregate dict into a simple TSV file
    '''
    # open file
    fp = open(output_fn, 'w+')
    tsv_writer = csv.DictWriter(fp, delimiter='\t', fieldnames=[
        'gene', 'ancestry', 'sex', 'starphase_diplotype', 'count'
    ])
    tsv_writer.writeheader()

    # write each aggregate value
    for k in sorted(aggregate_data.keys()):
        (gene, anc, sex, dip) = k
        count = aggregate_data[k]

        row = {
            'gene' : gene,
            'ancestry' : anc,
            'sex' : sex,
            'starphase_diplotype' : dip,
            'count' : count
        }
        tsv_writer.writerow(row)

    fp.close()

if __name__ == '__main__':
    description = 'Loads all indicated data and aggregates it at the ancestry level'
    p = ap.ArgumentParser(description=description, formatter_class=ap.RawTextHelpFormatter)

    p.add_argument('-i', '--input-tsv', dest='input_tsv', required=True, help='the input collation file (TSV)')
    p.add_argument('-o', '--output-tsv', dest='output_tsv', required=True, help='the output aggregate file (TSV)')
    
    args = p.parse_args()

    # load the data files we need to parse
    print(f'Loading input collation file {args.input_tsv}...')
    data_collection = loadInputTsv(args.input_tsv)

    # now load the actual data
    print('Loading all data files...')
    sample_data = loadCohortData(data_collection)

    # aggregate everything, we do NOT want to reduce alleles at this point; that will happy on all data
    # output here is a dict with (gene, ancestry, sex, diplotype) -> count; no sample info exists anymore
    aggregate_data = collectByAncestry(sample_data)

    # finally, save everything for the greater aggregation
    print(f'Saving aggregate data to {args.output_tsv}...')
    writeAggregateTsv(aggregate_data, args.output_tsv)
    print('Done.')

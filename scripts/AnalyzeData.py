'''
Final analysis file that loads all aggregated data files and generates all output data / figures
'''

import csv
import glob
import numpy as np
import os

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from PipelineConfig import *

#####################################################
# Config
#####################################################

# Constants for the HLA delta fields
ENABLE_HLA_DELTAS = True # if True, this will run the HLA delta figures as well
HLA_MAX_DELTA = 5 # this is a cut-off between minor/major ED
HLA_MISSING = "Missing" # no entry in DB, should only happen to DNA
HLA_EQUAL = "Equal" # exact sequence match in overlap
HLA_OFF_BY_ONE = "Off-by-one (ED=1)" # 1bp delta, could be homo-polymer error
HLA_MINOR_DELTA = f"Minor delta (ED<={HLA_MAX_DELTA})" # small delta, but greater than 1
HLA_MAJOR_DELTA = f"Major delta (ED>{HLA_MAX_DELTA})" # big delta

# controls the CYP2D6 stuff
IMPACT_MODE = 'status' # 'status' or 'score'
D6_POOR = 'Poor'
D6_INTERMEDIATE = 'Intermediate'
D6_NORMAL = 'Normal'
D6_ULTRA = 'Ultrarapid'
D6_INDETERMINATE = 'Indeterminate'

def loadD6Impact():
    '''
    Loads the various D6 allele impact scores
    '''
    fn = f'{DATA_FOLDER}/CYP2D6_allele_functionality_reference.csv'
    fp = open(fn, 'r')
    fp.readline() # skip a line
    csv_reader = csv.DictReader(fp)
    ret = {}
    for row in csv_reader:
        allele_id = row['Allele/cDNA/rsID']
        score = row['Activity Value (Optional)']
        status = row['Allele Clinical Functional Status (Required)']
        ret[allele_id] = {
            'score' : score,
            'status' : status
        }

    ret['*68x2'] = {
        'score' : '0.0',
        'status' : 'No function'
    }
    
    fp.close()
    return ret

D6_IMPACT_VALUES = loadD6Impact()

#####################################################
# Data loading
#####################################################

def loadAllAggregateFiles(filenames):
    '''
    Goes through all files and loads them as a single cohort, original cohort ID is masked
    '''
    ret = {}
    print('Loading all aggregate files...')
    for fn in filenames:
        print(f'\tLoading aggregate file: {fn}')
        fp = open(fn, 'r')
        tsv_reader = csv.DictReader(fp, delimiter='\t')
        for row in tsv_reader:
            # pull out all these values as a key
            gene = row['gene']
            ancestry = row['ancestry']
            sex = row['sex']
            diplotype = row['starphase_diplotype']
            key = (gene, ancestry, sex, diplotype)

            # observed totals from this cohort
            count = int(row['count'])

            # add the count
            ret[key] = ret.get(key, 0) + count

        fp.close()
    print('Done loading aggregate files.')
    return ret

#####################################################
# Data aggregation and translation
#####################################################

def collectByAncestry(all_data, reduce_to_core):
    '''
    Go through and reorganize by gene, then ancestry
    '''
    # [gene][ancestry][haplotype] -> count
    ancestry_counts = {}

    for data_key in all_data:
        # parse the data key and get the count
        (gene, ancestry, sex, diplotype) = data_key
        count = all_data[data_key]

        # diplotype is always a single entry in our aggregates
        
        if gene not in ancestry_counts:
            ancestry_counts[gene] = {}
            if gene == 'CYP2D6':
                # extra values parsed for D6
                ancestry_counts['CYP2D6_cn'] = {}
                ancestry_counts['CYP2D6_impact'] = {}
                ancestry_counts['CYP2D6_dip_func'] = {}

        if ancestry not in ancestry_counts[gene]:
            ancestry_counts[gene][ancestry] = {}
            if gene == 'CYP2D6':
                # extra values parsed for D6
                ancestry_counts['CYP2D6_cn'][ancestry] = {}
                ancestry_counts['CYP2D6_impact'][ancestry] = {}
                ancestry_counts['CYP2D6_dip_func'][ancestry] = {}
        
        # get the haplotypes
        hap1, hap2 = diplotype.split('/')

        # wrapper for later
        hap_calls = {
            'hap1' : hap1,
            'hap2' : hap2
        }

        if gene == 'MT-RNR1':
            # special handling for mitochrondria and heteroplasmy; there may be a handful of these
            if hap1 == hap2:
                final_hap = hap1
            else:
                final_hap = 'heteroplasmy'
            ancestry_counts[gene][ancestry][final_hap] = ancestry_counts[gene][ancestry].get(final_hap, 0) + count
        
        elif gene == 'G6PD' and sex == 'male':
            # males should only count once here, need to verify hom.
            if hap1 != hap2:
                raise Exception('handle het male')
            
            final_hap = hap1
            ancestry_counts[gene][ancestry][final_hap] = ancestry_counts[gene][ancestry].get(final_hap, 0) + count
        
        elif gene.endswith('dna_delta'):
            # special HLA DNA/cDNA delta handling, we need to translate the delta value
            for hap in ['hap1', 'hap2']:
                hap_val = hap_calls[hap]            
                category = categoryHlaDelta(hap_val)
                ancestry_counts[gene][ancestry][category] = ancestry_counts[gene][ancestry].get(category, 0) + count

        else:
            # handle hap reductions
            if reduce_to_core and gene in ['HLA-A', 'HLA-B']:
                # modify hap1 and hap2 by stripping down to 2-field
                hap1 = reduceHlaHap(hap1)
                hap2 = reduceHlaHap(hap2)
            
            elif reduce_to_core and gene == 'CYP2D6':
                # remove sub-alleles; e.g *4.001 -> *4
                hap1 = reduceCyp2d6Hap(hap1)
                hap2 = reduceCyp2d6Hap(hap2)

            # all others are "normal" diplotypes, each gets count per hap
            ancestry_counts[gene][ancestry][hap1] = ancestry_counts[gene][ancestry].get(hap1, 0) + count
            ancestry_counts[gene][ancestry][hap2] = ancestry_counts[gene][ancestry].get(hap2, 0) + count

            if gene == 'CYP2D6':
                # get the haplotype copy number
                hap_scores = []
                for hap in ['hap1', 'hap2']:
                    hap_cn = reduceCyp2d6CN(hap_calls[hap])
                    ancestry_counts['CYP2D6_cn'][ancestry][hap_cn] = ancestry_counts['CYP2D6_cn'][ancestry].get(hap_cn, 0) + count

                    # get the haplotype impact
                    hap_impact = scoreCyp2d6(hap_calls[hap], IMPACT_MODE)
                    ancestry_counts['CYP2D6_impact'][ancestry][hap_impact] = ancestry_counts['CYP2D6_impact'][ancestry].get(hap_impact, 0) + count

                    # also get the score, which is always just a number or n/a
                    hap_score = scoreCyp2d6(hap_calls[hap], 'score')
                    hap_scores.append(hap_score)
                
                # lastly, get the diplotype function
                dip_function = scoreDipCyp2d6(hap_scores[0], hap_scores[1])
                ancestry_counts['CYP2D6_dip_func'][ancestry][dip_function] = ancestry_counts['CYP2D6_dip_func'][ancestry].get(dip_function, 0) + count
    
    return ancestry_counts

def reduceHlaHap(hap):
    '''
    Reduces an N-field HLA allele into a 2-field allele
    '''
    hap_frags = hap.split(':')
    assert(len(hap_frags) >= 2)
    return ':'.join(hap_frags[:2])

def reduceCyp2d6Hap(hap):
    '''
    Removes fractions from the allele defs
    '''
    d6_frags = hap.split(' + ')
    reduc = []
    for frag in d6_frags:
        copy_frags = frag.split('x')
        assert(len(copy_frags) <= 2)

        allele = copy_frags[0][1:] # remove *
        float_val = float(allele)
        int_val = int(np.floor(float_val))

        copy_frags[0] = str(int_val)
        reduc.append('*' + 'x'.join(copy_frags))
    return ' + '.join(reduc)

def reduceCyp2d6CN(hap):
    '''
    Parses a haplotype into full allele and hybrid counts
    '''
    d6_frags = hap.split(' + ')
    cn_count = 0
    hybrid_count = 0
    for frag in d6_frags:
        copy_frags = frag.split('x')
        assert(len(copy_frags) <= 2)

        allele = copy_frags[0][1:] # remove *
        float_val = float(allele)
        int_val = int(np.floor(float_val))

        if len(copy_frags) == 1:
            cn = 1
        else:
            cn = int(copy_frags[1])

        if int_val == 5:
            assert(cn == 1)
            assert(len(d6_frags) == 1)
            cn_count = 0

        # hybrids from page 9 of https://a.storyblok.com/f/70677/x/169927da7a/cyp2d6_structural-variation_v3-4.pdf
        elif int_val in [13, 36, 61, 63, 68, 83] or allele == '4.013':
            hybrid_count += cn

        else:
            cn_count += cn

    if hybrid_count > 0:
        ret = f'CYP2D6x{cn_count}, Hybridx{hybrid_count}'
    else:
        if cn_count == 0:
            ret = f'CYP2D6*5 (x0)'
        else:
            ret = f'CYP2D6x{cn_count}'
    
    # debug
    # print(hap, '->', ret)
    return ret

def scoreCyp2d6(hap, impact_mode):
    '''
    Calculates the impact score or category for the provided haplotype
    '''
    d6_frags = hap.split(' + ')
    impacts = []
    for frag in d6_frags:
        copy_frags = frag.split('x')
        assert(len(copy_frags) <= 2)

        allele = copy_frags[0][1:] # remove *
        float_val = float(allele)
        int_val = int(np.floor(float_val))

        if len(copy_frags) == 2:
            copy_count = int(copy_frags[1])
        else:
            copy_count = 1

        if copy_count >= 3:
            # special matching here
            frag_lookup = f'*{int_val}x≥3'
        else:
            # normal system with x2
            copy_frags[0] = str(int_val)
            frag_lookup = '*' + 'x'.join(copy_frags)

        impacts.append(D6_IMPACT_VALUES[frag_lookup])
    
    if len(impacts) > 1:
        # see if we can filter down the impacts to something easy by removing all defunct ones
        filtered_impacts = []
        for imp in impacts:
            if imp['status'] == 'No function':
                # can be ignored
                pass
            else:
                filtered_impacts.append(imp)
        
        if len(filtered_impacts) == 0:
            # everything is defunct, so populate with No function
            impacts = [{'score' : '0.0', 'status' : 'No function'}]
        else:
            # we had something left over
            impacts = filtered_impacts

    if len(impacts) > 1:
        # raise Exception(f'Handle {impacts}')
        print(f"ERROR: Handle {hap} -> {impacts}")
        return 'complex'
    else:
        if impact_mode == 'score':
            return impacts[0]['score']
        else:
            return impacts[0]['status']

def scoreDipCyp2d6(score1, score2):
    '''
    Evaluates the diplotype function based on total score:
    'n/a' + anything = indeterminate
    0.0 => poor
    [0.25, 1.0] => intermediate
    [1.25, 2.25] => normal
    >2.25 => Ultrarapid
    '''
    if score1 == 'n/a' or score2 == 'n/a':
        return D6_INDETERMINATE
    
    if score1.startswith('≥'):
        v1 = float(score1[1:])
    else:
        v1 = float(score1)
    if score2.startswith('≥'):
        v2 = float(score2[1:])
    else:
        v2 = float(score2)
    
    total = v1+v2
    if total == 0:
        return D6_POOR
    elif total >= 0.25 and total <= 1.0:
        return D6_INTERMEDIATE
    elif total >= 1.25 and total <= 2.25:
        return D6_NORMAL
    elif total > 2.25:
        return D6_ULTRA
    else:
        raise Exception(f'Handle total: {score1}+{score2} = {total}')

def categoryHlaDelta(delta_value):
    '''
    Categorizes an HLA delta score based on the value
    '''
    # first, check for missing
    if delta_value == 'None':
        return HLA_MISSING
    
    # everything else is an 'int' type
    delta_value = int(delta_value)
    if delta_value == 0:
        return HLA_EQUAL
    elif delta_value == 1:
        return HLA_OFF_BY_ONE
    elif delta_value <= HLA_MAX_DELTA:
        return HLA_MINOR_DELTA
    else:
        return HLA_MAJOR_DELTA

#####################################################
# Data analysis / visualization
#####################################################

def generateAncestryPlots(ancestry_data):
    '''
    This is where we can create a bunch of plots
    @param ancestry_data - [gene][ancestry][hap] -> count
    '''
    image_folder = f'{RESULTS_FOLDER}/ancestry_images'
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)
    
    print(f'Generating images at {image_folder}...')

    for gene in sorted(ancestry_data.keys()):
        print(f'\tCreating image for {gene}...')

        # contains the total counts across all ancestries
        total_counts = {}
        ancestry_totals = {}

        for ancestry in ancestry_data[gene]:
            for hap in ancestry_data[gene][ancestry]:
                # accumulate totals
                total_counts[hap] = total_counts.get(hap, 0) + ancestry_data[gene][ancestry][hap]

                # count the denominator for each ancestry
                ancestry_totals[ancestry] = ancestry_totals.get(ancestry, 0) + ancestry_data[gene][ancestry][hap]
        
        if gene.endswith('dna_delta'):
            # specify the order for these
            ordered_keys = [(k, total_counts.get(k, 0)) for k in [
                HLA_EQUAL, HLA_OFF_BY_ONE, HLA_MINOR_DELTA, HLA_MAJOR_DELTA, HLA_MISSING
            ]]
        elif gene == 'CYP2D6_dip_func':
            ordered_keys = [(k, total_counts.get(k, 0)) for k in [
                D6_POOR, D6_INTERMEDIATE, D6_NORMAL, D6_ULTRA, D6_INDETERMINATE
            ]]
        else:
            # order from most count to least count
            ordered_keys = [(k, total_counts[k]) for k in total_counts]
            ordered_keys.sort(key=lambda t: (t[1], t[0]), reverse=True)

        # get the total denom
        total_hap_count = sum([t[1] for t in ordered_keys])

        # figure out if we need to restrict our display
        max_colors = 20 # max number of categories
        cycle_length = 10 # if greater than this, we add a hatch
        if len(ordered_keys) > max_colors:
            # we need to truncate to 9, and then have an "everything else"
            low_values = ordered_keys[max_colors-1:]
            high_values = ordered_keys[:max_colors-1]

            other_set = set(t[0] for t in low_values)
            tail_count = sum(t[1] for t in low_values)
            ordered_keys = high_values+[('Other', tail_count)]
        else:
            other_set = set([])

        # create the main figure        
        fig = plt.figure()
        bottoms = np.full((len(ancestry_totals)+1, ), 100.0)
        label_order = ['Combined'] + sorted(ancestry_totals.keys())
        for (i, (hap, count)) in enumerate(ordered_keys):
            values = []
            axes_labels = []
            for ancestry in label_order:
                if ancestry == 'Combined':
                    v = count
                    denom = total_hap_count
                else:
                    if hap == 'Other':
                        v = 0
                        for h in ancestry_data[gene][ancestry]:
                            if h in other_set:
                                v += ancestry_data[gene][ancestry][h]
                    else:
                        v = ancestry_data[gene][ancestry].get(hap, 0)
                    denom = ancestry_totals[ancestry]
                # we watch a percentage
                values.append(100.0 * v / denom)
                axes_labels.append(f'{ancestry}\n(N={denom})')
            
            bottoms -= np.array(values)
            max_len = 30
            if len(hap) > max_len:
                hap = hap[:max_len-3]+'...'
            
            if i >= 2*cycle_length:
                raise Exception('need additional hatches')
            elif i >= cycle_length:
                hatch = '//'
            else:
                hatch = None
            plt.bar(axes_labels, values, width=0.6, label=hap, bottom=bottoms, hatch=hatch)
        
        # customize title for some of the "weird" plots
        if gene == 'CYP2D6_cn':
            title = 'Copy number distribution for CYP2D6 by ancestry'
            legend_title = 'Top copy numbers'
        elif gene == 'CYP2D6_impact':
            title = 'CYP2D6 haplotype function by ancestry'
            legend_title = 'Functional categories'
        elif gene == 'CYP2D6_dip_func':
            title = 'CYP2D6 diplotype predicted metabolizer phenotype'
            legend_title = 'Metabolizer categories'
        elif gene.endswith('dna_delta'):
            g, tag, d = gene.split('_')
            if tag == 'dna':
                tag = 'DNA'
            else:
                tag = 'cDNA'
            title = f'Differences between DB and consensus for {g} {tag}'
            legend_title = 'Difference category'
        else:
            title = f'Haplotype distribution for {gene} by ancestry'
            legend_title = 'Top haplotypes'

        plt.grid(axis='y')
        # plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left') # anchors to right side
        plt.legend(title=legend_title, bbox_to_anchor=(1.0, 0.5), loc='center left')

        plt.ylim(0, 100)
        # plt.xticks(rotation=60) #no rotation really needed here
        plt.ylabel('Percentage')
        plt.xlabel('Ancestry (peddy)')

        plt.title(title)

        plt.savefig(f'{image_folder}/{gene}_distribution.png', bbox_inches='tight')
        plt.close()
    
    print('Image generation complete.')

if __name__ == '__main__':
    # first, load all the data into a single dict
    aggregate_files = glob.glob(f'{DATA_FOLDER}/cohort_aggregate_files/*.tsv')
    all_loaded_data = loadAllAggregateFiles(aggregate_files)
    
    # consolidate everything by ancestry
    ancestry_collection = collectByAncestry(all_loaded_data, True)

    # now generate ancestry images
    generateAncestryPlots(ancestry_collection)

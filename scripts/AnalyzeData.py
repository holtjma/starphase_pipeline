'''
Final analysis file that loads all aggregated data files and generates all output data / figures
'''

import csv
import glob
import gzip
import json
import numpy as np
import os

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from AggregateData import writeAggregateTsv
from PipelineConfig import *

#####################################################
# Config
#####################################################

# controls a bunch of fonts
USE_LATEX = True # if True, requires a local install of latex to generate fonts
if USE_LATEX:
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "Bookman" # best so far: "Helvetica" or "Bookman"

# Constants for the HLA delta fields
ENABLE_HLA_DELTAS = True # if True, this will run the HLA delta figures as well
HLA_MAX_DELTA = 5 # this is a cut-off between minor/major ED
HLA_MISSING = "Missing" # no entry in DB, should only happen to DNA
HLA_EQUAL = "Exact match (ED=0)" # exact sequence match in overlap
HLA_OFF_BY_ONE = "Off-by-one (ED=1)" # 1bp delta, could be homo-polymer error
HLA_MINOR_DELTA = f"Minor delta (ED<={HLA_MAX_DELTA})" # small delta, but greater than 1
HLA_MAJOR_DELTA = f"Major delta (ED>{HLA_MAX_DELTA})" # big delta
GENERATE_JOINT_FIGURES = True # enables a single joint image for some of the paper figures
USE_CUSTOM_COLORS = True # enables a non-default color scheme
MAKE_TITLES = False # if True, include titles
CUSTOM_LEGEND = False # if True, shift the legend location below and make it 2-column

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

    # while not specified in the file, additional copies should also have no function
    ret['*68x2'] = {
        'score' : '0.0',
        'status' : 'No function'
    }
    ret['*68x≥3'] = ret['*68x2']

    # sentinel for our undefined hybrids
    ret['*None'] = {
        'score' : 'n/a',
        'status' : 'Uncertain function'
    }
    
    fp.close()
    return ret

D6_IMPACT_VALUES = loadD6Impact()

def loadPopFreq(pop_freq_folder):
    '''
    Loads the population frequency information from the AllOfUs v7 populations.
    Indeterminate haplotypes are ignored and removed from the frequency count.
    Return values is a dict[gene][ancestry/"ALL"][allele/"TOTAL"] => count.
    @param pop_freq_folder - the "AllOfUs_Frequencies_v7/allele" folder path
    '''
    ret = {}

    filenames = sorted(glob.glob(f'{pop_freq_folder}/*_allele.tsv'))
    for fn in filenames:
        gene = fn.split('/')[-1][:-11]
        gene_vals = {}
        all_anc = {
            'TOTAL' : 0
        }

        fp = open(fn, 'r')
        tsv_reader = csv.DictReader(fp, delimiter='\t')
        for row in tsv_reader:
            anc = row['biogeographic_group'].upper()
            allele = row['allele']
            if allele == 'Indeterminate':
                # ignore these
                continue
            hap_count = int(row['n_haplotype'])

            if anc not in gene_vals:
                gene_vals[anc] = {
                    'TOTAL' : 0
                }
            assert(allele not in gene_vals[anc])

            #add the counts for the ancestry and the totals
            gene_vals[anc][allele] = hap_count
            gene_vals[anc]['TOTAL'] += hap_count
            all_anc[allele] = all_anc.get(allele, 0) + hap_count
            all_anc['TOTAL'] += hap_count
        
        gene_vals['ALL'] = all_anc
        ret[gene] = gene_vals
        fp.close()

    return ret

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
        gene_count = {}
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
            gene_count[gene] = gene_count.get(gene, 0) + count

        fp.close()

        observed_counts = set([])
        for g in sorted(gene_count.keys()):
            c = gene_count[g]
            if c not in observed_counts:
                print(f'\t\t{g} -> {c}')
                observed_counts.add(c)
        
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

    # summary stats
    het_dict = {}
    hom_dict = {}
    hemi_dict = {}
    no_match_dict = {}
    mt_hom_dict = {}
    mt_het_dict = {}

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
                mt_hom_dict[gene] = mt_hom_dict.get(gene, 0) + count
            else:
                final_hap = 'heteroplasmy'
                mt_het_dict[gene] = mt_het_dict.get(gene, 0) + count
            ancestry_counts[gene][ancestry][final_hap] = ancestry_counts[gene][ancestry].get(final_hap, 0) + count
        
        elif gene == 'G6PD' and sex == 'male':
            # males should only count once here, need to verify hom.
            if hap1 != hap2:
                raise Exception('handle het male')
            
            # hemi-zygous
            hemi_dict[gene] = hemi_dict.get(gene, 0) + count

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

            # these are the only ones with het/hom info
            if hap1 == 'NO_MATCH' or hap2 == 'NO_MATCH':
                # no counting it
                no_match_dict[gene] = no_match_dict.get(gene, 0) + count
            elif hap1 == hap2:
                hom_dict[gene] = hom_dict.get(gene, 0) + count
            else:
                het_dict[gene] = het_dict.get(gene, 0) + count

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
    
    # add it up
    total_hets = sum(het_dict.values())
    total_homs = sum(hom_dict.values())
    total_hemi = sum(hemi_dict.values())
    total_no_match = sum(no_match_dict.values())
    total_mt_hom = sum(mt_hom_dict.values())
    total_mt_het = sum(mt_het_dict.values())
    combined = total_hets + total_homs + total_hemi + total_no_match + total_mt_hom + total_mt_het

    combined_counts = {
        'hets' : het_dict,
        'homs' : hom_dict,
        'hemis' : hemi_dict,
        'no_matches' : no_match_dict,
        'mt_homs' : mt_hom_dict,
        'mt_hets' : mt_het_dict,
    }

    # print it
    print(f'\tTotal hets: {total_hets}')
    print(f'\tTotal homs: {total_homs}')
    print(f'\tHet/hom ratio (includes REF alleles): {total_hets / total_homs}')
    print(f'\tTotal hemi: {total_hemi}')
    print(f'\tTotal NO_MATCH: {total_no_match} ({total_no_match / combined})')
    print(f'\tTotal MT homs: {total_mt_hom}')
    print(f'\tTotal MT heteroplasmy: {total_mt_het}')

    return ancestry_counts, combined_counts

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
        if allele.startswith('CYP2D6::CYP2D7::'):
            # this is an unknown hybrid, we do not need to convert anything here
            copy_frags[0] = allele
        elif allele == '4.013':
            # do not collapse this one, it's a special case that is different from *4
            copy_frags[0] = allele
        else:
            float_val = float(allele)
            int_val = int(np.floor(float_val))
            copy_frags[0] = str(int_val)
        
        reduc.append('*' + 'x'.join(copy_frags))
    return ' + '.join(reduc)

def decomposeCyp2d6Hap(hap):
    '''
    This will return a set of all observed allele in a D6 haplotype
    '''
    d6_frags = hap.split(' + ')
    reduc = []
    for frag in d6_frags:
        copy_frags = frag.split('x')
        assert(len(copy_frags) <= 2)
        reduc.append(copy_frags[0]) # if it is *4x2, we just want the *4 component
    return set(reduc)

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
        is_undefined_hybrid = False
        if allele.startswith('CYP2D6::CYP2D7::'):
            # this is an unknown hybrid, so count as a hybrid later
            is_undefined_hybrid = True
            int_val = None
        elif allele == '4.013':
            # do not collapse this one, it's a special case that is different from *4
            int_val = None
        else:
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
        elif int_val in [13, 36, 61, 63, 68, 83] or allele == '4.013' or is_undefined_hybrid:
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

    # first interate over the haplotypes and combine those with the same base integer label; e.g. *1.001 + *1.002 = *1x2
    hap_copy = {}
    for frag in d6_frags:
        copy_frags = frag.split('x')
        assert(len(copy_frags) <= 2)

        allele = copy_frags[0][1:] # remove *
        if allele.startswith('CYP2D6::CYP2D7::'):
            # this is an unknown hybrid, we will not have an impact score here
            int_val = None
        #elif allele == '4.013': # we do not currently need a special case for this one here
        else:
            float_val = float(allele)
            int_val = int(np.floor(float_val))

        if len(copy_frags) == 2:
            copy_count = int(copy_frags[1])
        else:
            copy_count = 1
        
        hap_copy[int_val] = hap_copy.get(int_val, 0) + copy_count
    
    # now go through the counts and add the impact values
    for int_val in hap_copy:
        copy_count = hap_copy[int_val]

        if copy_count == 1:
            frag_lookup = f'*{int_val}'
        elif copy_count == 2:
            frag_lookup = f'*{int_val}x2'
        elif copy_count >= 3:
            # special matching here
            frag_lookup = f'*{int_val}x≥3'
        else:
            raise Exception(f'unhandled copy_count={copy_count}')

        impacts.append(D6_IMPACT_VALUES[frag_lookup])
    
    # if we have more than one impact, make sure that all the others are defunct
    if len(impacts) > 1:
        # see if we can filter down the impacts to something easy by removing all defunct ones
        filtered_impacts = []
        uncertain_detected = False
        for imp in impacts:
            if imp['status'] == 'No function':
                # can be ignored
                pass
            elif imp['status'] == 'Uncertain function':
                uncertain_detected = True
                filtered_impacts.append(imp)
            else:
                filtered_impacts.append(imp)
        
        if len(filtered_impacts) == 0:
            # everything is defunct, so populate with No function
            impacts = [{'score' : '0.0', 'status' : 'No function'}]
        elif uncertain_detected:
            # we detected at least one that is uncertain, so all will be uncertain
            impacts = [{'score' : 'n/a', 'status' : 'Uncertain function'}]
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
# Basic writers
#####################################################

def writeAggregateHaps(aggregate_data, output_fn):
    '''
    Converts the aggregate haplotype dict into a simple TSV file
    '''
    # open file
    fp = open(output_fn, 'w+')
    tsv_writer = csv.DictWriter(fp, delimiter='\t', fieldnames=[
        'gene', 'ancestry', 'starphase_haplotype', 'count'
    ])
    tsv_writer.writeheader()

    # write each aggregate value
    for gene in sorted(aggregate_data.keys()):
        for anc in sorted(aggregate_data[gene].keys()):
            for hap in sorted(aggregate_data[gene][anc].keys()):
                count = aggregate_data[gene][anc][hap]
                row = {
                    'gene' : gene,
                    'ancestry' : anc,
                    'starphase_haplotype' : hap,
                    'count' : count
                }
                tsv_writer.writerow(row)

    fp.close()

#####################################################
# Data analysis / visualization
#####################################################

def generateAncestryPlots(ancestry_data, popfreqs):
    '''
    This is where we can create a bunch of plots
    @param ancestry_data - dict[gene][ancestry][hap] -> count
    @param popfreqs - dict[gene][ancestry/"ALL"][allele/"TOTAL"] => count
    '''
    image_folder = f'{RESULTS_FOLDER}/ancestry_images'
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)
    
    print(f'Generating ancestry images at {image_folder}...')

    if USE_CUSTOM_COLORS:
        # colors can be any length, but more than 10 is usually bad for visuals
        colors = [
            '#1383C6', #bright blue,
            # '#E16A2C', #bright orange
            '#F99D41', #light orange
            '#009D4E', #bright green; default palette has red next, which is bad for RG colorblind
            '#5F249F', #purple
            '#FF66CC', #magenta
            'lightgrey',
            #'#6ABF6A', #light green
            #'red',
            #'green',
            #'blue'
        ]
        
    else:
        # get the default color cycle
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']

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
        # cycle_length = 10 # if greater than this, we add a hatch
        # max_colors = 20 # max number of categories
        # assert(cycle_length == len(colors))
        
        cycle_length = len(colors)
        max_colors = 2 * cycle_length

        if len(ordered_keys) > max_colors:
            # we need to truncate to 9, and then have an "everything else"
            low_values = ordered_keys[max_colors-1:]
            high_values = ordered_keys[:max_colors-1]

            other_set = set(t[0] for t in low_values)
            tail_count = sum(t[1] for t in low_values)
            ordered_keys = high_values+[('Other', tail_count)]
            other_plotted = True
        else:
            other_set = set([])
            other_plotted = False

        # create the main figure        
        fig = plt.figure()
        width = 0.3
        ind = np.arange(len(ancestry_totals)+1)

        bottoms = np.full((len(ancestry_totals)+1, ), 100.0)
        bottoms_exp = np.full((len(ancestry_totals)+1, ), 100.0)

        gene_in_pop = (popfreqs != None and (gene in popfreqs))

        label_order = ['Combined'] + sorted(ancestry_totals.keys())
        for (i, (hap, count)) in enumerate(ordered_keys):
            values = []
            exp_values = []
            axes_labels = []
            for (a_ind, ancestry) in enumerate(label_order):
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

                if hap == 'NO_MATCH' and ancestry == 'Combined':
                    print('', '', ancestry, hap, 100 * v / denom, sep='\t')

                # we watch a percentage
                values.append(100.0 * v / denom)
                axes_labels.append(f'{ancestry}\n(N={denom})')

                if gene_in_pop and hap != 'NO_MATCH':
                    if ancestry == 'Combined':
                        anc_translate = "ALL"
                    elif ancestry == 'UNKNOWN':
                        anc_translate = 'OTH'
                    else:
                        anc_translate = ancestry
                    
                    if hap == 'Other':
                        # we use everything left
                        exp_v = bottoms_exp[a_ind]
                        exp_denom = 100.0 # it's already in fraction form since we're looking at the bottom
                    else:
                        exp_v = popfreqs[gene][anc_translate].get(hap, 0)
                        exp_denom = popfreqs[gene][anc_translate]["TOTAL"]
                    exp_values.append(100.0 * exp_v / exp_denom)
                else:
                    # if the gene is not in the population OR we have a NO_MATCH, just use 0.0 placehold
                    exp_values.append(0.0)
            
            bottoms -= np.array(values)
            bottoms_exp -= np.array(exp_values)

            max_len = 30
            if len(hap) > max_len:
                hap = hap[:max_len-3]+'...'
            
            if USE_LATEX:
                hap = hap.replace('>', '$>$').replace('<=', '$\\le$')
            
            if i >= 2*cycle_length:
                raise Exception('need additional hatches')
            elif i >= cycle_length:
                hatch = '//'
            else:
                hatch = None

            if gene_in_pop:
                plt.bar(ind + width / 2.0, values, width=width, label=hap, bottom=bottoms, hatch=hatch, color=colors[i % cycle_length], edgecolor='black', linewidth=1)
                plt.bar(ind - width / 2.0, exp_values, width=width, bottom=bottoms_exp, hatch=hatch, color=colors[i % cycle_length], alpha=0.5, edgecolor='black', linewidth=1)
            else:
                # plt.bar(axes_labels, values, width=0.6, label=hap, bottom=bottoms, hatch=hatch)
                plt.bar(ind, values, width=2.0*width, label=hap, bottom=bottoms, hatch=hatch, color=colors[i % cycle_length], edgecolor='black', linewidth=1)
        
        if gene_in_pop and (not other_plotted) and np.any(bottoms_exp > 0.0001):
            # only plot this special one IF no "Other" was already plotted AND popfreqs is enabled AND we have pop unaccounted for
            plt.bar([0], [0], width=width, label='Other', hatch='//', color=colors[-1], edgecolor='black', linewidth=1)
            plt.bar(ind - width / 2.0, bottoms_exp, width=width, hatch='//', color=colors[-1], alpha=0.5, edgecolor='black', linewidth=1)

        plt.xticks(ind, axes_labels)
            
        # customize title for some of the "weird" plots
        if gene == 'CYP2D6_cn':
            title = 'Copy number distribution for CYP2D6 by ancestry'
            if USE_LATEX:
                title = title.replace("CYP2D6", "\\textit{CYP2D6}")
            legend_title = 'Top copy numbers'
        elif gene == 'CYP2D6_impact':
            title = 'Predicted CYP2D6 haplotype function by ancestry'
            if USE_LATEX:
                title = title.replace("CYP2D6", "\\textit{CYP2D6}")
            legend_title = 'Functional categories'
        elif gene == 'CYP2D6_dip_func':
            title = 'Predicted CYP2D6 diplotype metabolizer phenotype by ancestry'
            if USE_LATEX:
                title = title.replace("CYP2D6", "\\textit{CYP2D6}")
            legend_title = 'Metabolizer categories'
        elif gene.endswith('dna_delta'):
            g, tag, d = gene.split('_')
            if USE_LATEX:
                g = f'\\textit{{{g}}}'
            if tag == 'dna':
                tag = 'DNA'
            else:
                tag = 'cDNA'
            title = f'Differences between DB and consensus for {g} {tag}'
            legend_title = 'Difference category'
        else:
            if USE_LATEX:
                title = f'Haplotype distribution for \\textit{{{gene}}} by ancestry'
            else:
                title = f'Haplotype distribution for {gene} by ancestry'
            legend_title = 'Top haplotypes'

        ax = plt.gca()
        ax.set_axisbelow(True)
        plt.grid(axis='y')
        # plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left') # anchors to right side
        if CUSTOM_LEGEND:
            leg_handles, leg_labels = ax.get_legend_handles_labels()
            # order = [0, 2, 4, 6, 8, 1, 3, 5, 7, 9]
            order = range(0, len(leg_handles))
            #plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
            plt.legend(
                [leg_handles[idx] for idx in order],
                [leg_labels[idx] for idx in order], 
                title=legend_title,
                #bbox_to_anchor=(1.0, 1.01, 1., .102), ncol=3, loc='lower right', borderaxespad=0.)
                bbox_to_anchor=(0.5, -0.2), ncol=3, loc='upper center', borderaxespad=0.)
        else:
            plt.legend(title=legend_title, bbox_to_anchor=(1.0, 0.5), loc='center left')

        plt.ylim(0, 100)
        # plt.xticks(rotation=60) #no rotation really needed here
        plt.ylabel('Percentage')
        plt.xlabel('Ancestry (peddy)')

        if MAKE_TITLES:
            plt.title(title)

        plt.savefig(f'{image_folder}/{gene}_distribution.png', bbox_inches='tight')
        plt.close()
    
    if GENERATE_JOINT_FIGURES:
        categories = [
            'HLA-A_cdna_delta', 'HLA-B_cdna_delta',
            'HLA-A_dna_delta', 'HLA-B_dna_delta'
        ]
        # basic
        # fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)

        # more complicated, but we can control spacing here
        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(2, 2, hspace=0.1, wspace=0.05)
        axes = gs.subplots(sharex=True, sharey=True)

        for (j, gene) in enumerate(categories):
            # set the figure
            ax = axes[j // 2, j % 2]

            # contains the total counts across all ancestries
            total_counts = {}
            ancestry_totals = {}

            for ancestry in ancestry_data[gene]:
                for hap in ancestry_data[gene][ancestry]:
                    # accumulate totals
                    total_counts[hap] = total_counts.get(hap, 0) + ancestry_data[gene][ancestry][hap]

                    # count the denominator for each ancestry
                    ancestry_totals[ancestry] = ancestry_totals.get(ancestry, 0) + ancestry_data[gene][ancestry][hap]
            
            # specify the order for these
            ordered_keys = [(k, total_counts.get(k, 0)) for k in [
                HLA_EQUAL, HLA_OFF_BY_ONE, HLA_MINOR_DELTA, HLA_MAJOR_DELTA, HLA_MISSING
            ]]
            
            # get the total denom
            total_hap_count = sum([t[1] for t in ordered_keys])

            # figure out if we need to restrict our display
            # max_colors = 20 # max number of categories
            # cycle_length = 10 # if greater than this, we add a hatch
            # assert(cycle_length == len(colors))
            cycle_length = len(colors)
            max_colors = 2 * cycle_length

            other_set = set([])
            other_plotted = False

            # create the main figure        
            width = 0.3
            ind = np.arange(len(ancestry_totals)+1)
            # bottoms = np.full((len(ancestry_totals)+1, ), 100.0)
            bottoms = np.zeros(shape=(len(ancestry_totals)+1, ))
            
            label_order = ['Combined'] + sorted(ancestry_totals.keys())
            for (i, (hap, count)) in enumerate(ordered_keys):
                values = []
                exp_values = []
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
                    # axes_labels.append(ancestry)
                
                # bottoms -= np.array(values)
                if USE_LATEX:
                    hap = hap.replace('>', '$>$').replace('<=', '$\\le$')
                ax.bar(ind, values, width=2.0*width, label=hap, bottom=bottoms, color=colors[i % cycle_length], edgecolor='black', linewidth=1)
                bottoms += np.array(values)
            
            ax.set_xticks(ind, axes_labels)
                
            # customize title for some of the "weird" plots
            if gene.endswith('dna_delta'):
                g, tag, d = gene.split('_')
                if USE_LATEX:
                    g = f'\\textit{{{g}}}'
                if tag == 'dna':
                    tag = 'DNA'
                else:
                    tag = 'cDNA'
                # title = f'Differences between DB and consensus for {g} {tag}'
                title = f'{g} {tag}'
                legend_title = 'Difference category'

            ax.set_axisbelow(True)
            ax.grid(axis='y')
            # plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left') # anchors to right side
            # plt.legend(title=legend_title, bbox_to_anchor=(1.0, 0.5), loc='center left')

            ax.set_ylim(84, 100)
            # plt.xticks(rotation=60) #no rotation really needed here
            ax.set_title(title, fontsize=16)
            # ax.set_ylabel('Percentage')
            # ax.set_xlabel('Ancestry (peddy)')
            ax.label_outer()

        fig.supxlabel('Ancestry (peddy)', fontsize=16, y=0.02)
        fig.supylabel('Percentage', fontsize=16, x=0.06)
        if MAKE_TITLES:
            fig.suptitle('Difference between database sequence and StarPhase consensus', fontsize=18, y=0.95)
        # plt.legend(title=legend_title, bbox_to_anchor=(1.0, 1.1), loc='center left', fontsize=16, title_fontsize=16)
        plt.legend(title=legend_title, bbox_to_anchor=(0.20, 1.4), loc='center left', fontsize=16, title_fontsize=16, framealpha=1.0)
        
        plt.savefig(f'{image_folder}/HLA_joint.png', bbox_inches='tight')
        plt.close()
    
        # now do the other figures a bit more systematically
        figure_dict = {
            'cpic_stack' : ('VKORC1', 'SLCO1B1'),
            'CYP2D6_stack' : ('CYP2D6', 'CYP2D6_dip_func')
        }
        for figure_name in sorted(figure_dict.keys()):
            # basic
            # fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)

            #is_shared_labels = figure_name != 'CYP2D6_stack'
            is_shared_labels = False
            
            # more complicated, but we can control spacing here
            if is_shared_labels:
                h_space = 0.1
            else:
                h_space = 0.25
            fig = plt.figure(figsize=(7, 10))
            gs = fig.add_gridspec(2, 1, hspace=h_space, wspace=0.05)

            axes = gs.subplots(sharex=is_shared_labels)

            for (j, gene) in enumerate(figure_dict[figure_name]):
                # set the figure
                #ax = axes[i // 2, i % 2]
                ax = axes[j]

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
                # max_colors = 20 # max number of categories
                # cycle_length = 10 # if greater than this, we add a hatch
                # assert(cycle_length == len(colors))
                cycle_length = len(colors)
                max_colors = 2 * cycle_length

                if len(ordered_keys) > max_colors:
                    # we need to truncate to 9, and then have an "everything else"
                    low_values = ordered_keys[max_colors-1:]
                    high_values = ordered_keys[:max_colors-1]

                    other_set = set(t[0] for t in low_values)
                    tail_count = sum(t[1] for t in low_values)
                    ordered_keys = high_values+[('Other', tail_count)]
                    other_plotted = True
                else:
                    other_set = set([])
                    other_plotted = False

                # create the main figure        
                # fig = plt.figure()
                width = 0.3
                ind = np.arange(len(ancestry_totals)+1)

                bottoms = np.full((len(ancestry_totals)+1, ), 100.0)
                bottoms_exp = np.full((len(ancestry_totals)+1, ), 100.0)

                gene_in_pop = (popfreqs != None and (gene in popfreqs))

                label_order = ['Combined'] + sorted(ancestry_totals.keys())
                for (i, (hap, count)) in enumerate(ordered_keys):
                    values = []
                    exp_values = []
                    axes_labels = []
                    for (a_ind, ancestry) in enumerate(label_order):
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

                        if gene_in_pop and hap != 'NO_MATCH':
                            if ancestry == 'Combined':
                                anc_translate = "ALL"
                            elif ancestry == 'UNKNOWN':
                                anc_translate = 'OTH'
                            else:
                                anc_translate = ancestry
                            
                            if hap == 'Other':
                                # we use everything left
                                exp_v = bottoms_exp[a_ind]
                                exp_denom = 100.0
                            else:
                                exp_v = popfreqs[gene][anc_translate].get(hap, 0)
                                exp_denom = popfreqs[gene][anc_translate]["TOTAL"]
                            exp_values.append(100.0 * exp_v / exp_denom)
                        else:
                            # if the gene is not in the population OR we have a NO_MATCH, just use 0.0 placehold
                            exp_values.append(0.0)
                    
                    bottoms -= np.array(values)
                    bottoms_exp -= np.array(exp_values)

                    max_len = 30
                    if len(hap) > max_len:
                        hap = hap[:max_len-3]+'...'

                    if USE_LATEX:
                        hap = hap.replace('>', '$>$').replace('<=', '$\\le$')
                    
                    if i >= 2*cycle_length:
                        raise Exception('need additional hatches')
                    elif i >= cycle_length:
                        hatch = '//'
                    else:
                        hatch = None

                    if gene_in_pop:
                        ax.bar(ind + width / 2.0, values, width=width, label=hap, bottom=bottoms, hatch=hatch, color=colors[i % cycle_length], edgecolor='black', linewidth=1)
                        ax.bar(ind - width / 2.0, exp_values, width=width, bottom=bottoms_exp, hatch=hatch, color=colors[i % cycle_length], alpha=0.5, edgecolor='black', linewidth=1)
                    else:
                        # plt.bar(axes_labels, values, width=0.6, label=hap, bottom=bottoms, hatch=hatch)
                        ax.bar(ind, values, width=2.0*width, label=hap, bottom=bottoms, hatch=hatch, color=colors[i % cycle_length], edgecolor='black', linewidth=1)
                
                if gene_in_pop and (not other_plotted) and np.any(bottoms_exp > 0.0001):
                    # only plot this special one IF no "Other" was already plotted AND popfreqs is enabled AND we have pop unaccounted for
                    ax.bar([0], [0], width=width, label='Other', hatch='//', color=colors[-1], edgecolor='black', linewidth=1)
                    ax.bar(ind - width / 2.0, bottoms_exp, width=width, hatch='//', color=colors[-1], alpha=0.5, edgecolor='black', linewidth=1)

                ax.set_xticks(ind, axes_labels)
                    
                # customize title for some of the "weird" plots
                if gene == 'CYP2D6_cn':
                    title = 'Copy number distribution for CYP2D6 by ancestry'
                    legend_title = 'Top copy numbers'
                elif gene == 'CYP2D6_impact':
                    title = 'Predicted CYP2D6 haplotype function by ancestry'
                    legend_title = 'Functional categories'
                elif gene == 'CYP2D6_dip_func':
                    if USE_LATEX:
                        title = 'Predicted \\textit{CYP2D6} diplotype metabolizer phenotype by ancestry'
                    else:
                        title = 'Predicted CYP2D6 diplotype metabolizer phenotype by ancestry'
                    legend_title = 'Metabolizer categories'
                elif gene.endswith('dna_delta'):
                    g, tag, d = gene.split('_')
                    if USE_LATEX:
                        g = f'\\textit{{{g}}}'
                    if tag == 'dna':
                        tag = 'DNA'
                    else:
                        tag = 'cDNA'
                    title = f'Differences between DB and consensus for {g} {tag}'
                    legend_title = 'Difference category'
                else:
                    if USE_LATEX:
                        title = f'Haplotype distribution for \\textit{{{gene}}} by ancestry'
                    else:
                        title = f'Haplotype distribution for {gene} by ancestry'
                    legend_title = 'Top haplotypes'

                # ax = plt.gca()
                ax.set_axisbelow(True)
                ax.grid(axis='y')
                # plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left') # anchors to right side
                ax.legend(title=legend_title, bbox_to_anchor=(1.0, 0.5), loc='center left')

                ax.set_ylim(0, 100)
                # plt.xticks(rotation=60) #no rotation really needed here
                # ax.set_ylabel('Percentage')
                # ax.set_xlabel('Ancestry (peddy)')

                ax.set_title(title)
            
            fig.supxlabel('Ancestry (peddy)', fontsize=16, y=0.04)
            fig.supylabel('Percentage', fontsize=16, x=0.04)
            # fig.suptitle('Difference between database sequence and StarPhase consensus', fontsize=18, y=0.95)    

            plt.savefig(f'{image_folder}/{figure_name}.png', bbox_inches='tight')
            plt.close()

    print('Ancestry image generation complete.')

def loadDatabaseHaplotypes(fn, reduce_to_core):
    '''
    Loads our database haplotype set for each gene.
    '''
    fp = gzip.open(fn, 'r')
    j = json.load(fp)
    fp.close()

    ret = {}
    for gene in j['gene_entries']:
        hap_dict = j['gene_entries'][gene]['defined_haplotypes']
        hap_set = set(hap_dict.keys())
        ret[gene] = hap_set

    # now do HLA
    ret['HLA-A'] = set([])
    ret['HLA-B'] = set([])
    for hap_id in j['hla_sequences']:
        gene = j['hla_sequences'][hap_id]['gene_name']
        star_allele = '*'+':'.join(j['hla_sequences'][hap_id]['star_allele'])
        if reduce_to_core:
            star_allele = reduceHlaHap(star_allele)
        ret[gene].add(star_allele)

    # lastly, CYP2D6; we manually add the ones that are not built into the DB, hybrids & deletions
    ret['CYP2D6'] = set([
        '*5', # deletion
        '*13', '*61', '*63', '*68' # hybrids annotated in StarPhase
    ])
    for raw_allele in j['cyp2d6_gene_def'].keys():
        # strip the prefix
        assert(raw_allele.startswith('CYP2D6*'))
        allele = raw_allele[6:]

        # reduce if needed and add
        if reduce_to_core:
            allele = reduceCyp2d6Hap(allele)
        ret['CYP2D6'].add(allele)
    
    return ret

def generateDbRep(db_haps, ancestry_data):
    '''
    Compares the DB size to our observations to see how many alleles we have observed
    @param db_haps - [gene] -> set of haplotypes in DB
    @param ancestry_data - [gene][ancestry][hap] -> count
    '''
    debug = False
    PLOT_UNOBSERVED = False
    PLOT_POP_LIMIT = False
    COLOR_BY_FRACTION = True

    gene_order = sorted(db_haps.keys())
    gene_labels = []
    observed_fractions = []
    observed_counts = []
    missing_counts = []
    total_observations = []
    for gene in gene_order:
        # build the sets
        db_gene = db_haps[gene]
        obs_gene = set([])
        total_obs = 0
        for anc in ancestry_data[gene]:
            if gene == 'CYP2D6':
                for haplotype in ancestry_data[gene][anc].keys():
                    # split the haplotype into individual alleles
                    allele_composition = decomposeCyp2d6Hap(haplotype)
                    obs_gene |= allele_composition
            else:
                obs_gene |= set(ancestry_data[gene][anc].keys())
            
            for h in ancestry_data[gene][anc]:
                total_obs += ancestry_data[gene][anc][h]

        # calculate overlaps
        shared = db_gene & obs_gene
        missing = db_gene - obs_gene
        unknown = obs_gene - db_gene
        frac_shared = 100.0 * len(shared) / len(db_gene)
        
        # output and verify everything we observe is known
        if debug:
            print(gene, total_obs, f'{frac_shared:.2f}', len(shared), len(missing), len(unknown), sep='\t')

        if len(unknown) > 0 and unknown != set(['NO_MATCH']) and unknown != set(['NO_MATCH', 'heteroplasmy']):
            # if this happens, we may need to encode more D6 haplotypes
            # OR if someone ran it with a wrong DB, we could get weird mismatches
            print(f'Unknown haps encountered in {gene}:', unknown)
            all_fine = True
            for hap in unknown:
                if hap.startswith('*CYP2D6::CYP2D7::'):
                    # this is just an unknown hap, so it's fine
                    pass
                else:
                    all_fine = False

            if not all_fine:
                raise Exception('Unknown haps encountered')
        
        # save the data for plotting
        observed_fractions.append(frac_shared)
        observed_counts.append(len(shared))
        missing_counts.append(len(missing))
        total_observations.append(total_obs)
        
    print(f'\tUnique alleles: {observed_counts}')
    print(f'\tSum = {np.sum(observed_counts)}')
    
    # we just have observed and unobserved
    plt.figure(figsize=(8,5)) # adding legend made dimensions weird, this gets back to something semi-normal looking
    width = 0.3
    ind = np.arange(len(gene_order))

    if PLOT_UNOBSERVED:
        first_bar = plt.bar(ind - width / 2.0, observed_counts, width=width, label='Observed', edgecolor='black', linewidth=1)
        last_bar = plt.bar(ind + width / 2.0, missing_counts, width=width, label='Unobserved', edgecolor='black', linewidth=1) #, bottom=observed_counts)
    else:
        if COLOR_BY_FRACTION:
            # cmap = plt.get_cmap('RdYlGn') # this one is hard for color blind
            # winter - too bright
            # PRGn - a little dark, but pretty good
            # PiYG - is this good for color blind? this seems to say yes: https://www.aptech.com/releases/gauss18/graphics-updates/color-brewer-palettes/
            cmap = plt.get_cmap('PiYG')
            bar_colors = [cmap(obs / 100.0) for obs in observed_fractions]
            first_bar = plt.bar(ind, observed_counts, width=2*width, label='Observed', color=bar_colors, edgecolor='black', linewidth=1)
            last_bar = [None]*len(first_bar)

            # custom color bar placement
            norm = matplotlib.colors.Normalize(vmin=0.0, vmax=100.0, clip=False)
            cb = plt.colorbar(
                matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
                format="%.0f\\%%", # %% for format; \\ for latex,
                # location='left',
                # orientation='horizontal',
                shrink=0.2, aspect=5,
                anchor=(-0.7, 0.96),
                ticks=[0, 50, 100]
            )
            cb.ax.yaxis.set_ticks_position('left')
        else:
            first_bar = plt.bar(ind, observed_counts, width=2*width, label='Observed', edgecolor='black', linewidth=1)
            last_bar = [None]*len(first_bar)

    # LaTeX doesn't like to give us a buffer space, so instead we're going to just disable for this one and then re-enable
    with plt.rc_context({"text.usetex": False, "font.family" : "DejaVu Sans"}):
        for (i, (fb, lb, obs_frac)) in enumerate(zip(first_bar, last_bar, observed_fractions)):
            # X-offset, can be derived (less precise) or fixed by enumerate
            # offset = fb.get_x() + fb.get_width() + 0.1 # shift it right a little
            offset = i + width + 0.1

            # Y-offset, a few options here
            # placed at the top - gets cut off
            if lb == None:
                height = fb.get_height()
            else:
                height = max(fb.get_height(), lb.get_height())

            text = f' {obs_frac:.1f}%'
            plt.annotate(text, (offset, height), ha='right', va='bottom', rotation='vertical', annotation_clip=False)

    if USE_LATEX:
        gene_ticks = [f'\\textit{{{g}}}' for g in gene_order]
    else:
        gene_ticks = gene_order
    plt.xticks(ind, gene_ticks, rotation=90)
    ax = plt.gca()
    ax.set_axisbelow(True)
    plt.grid(axis='y')
    plt.yscale('symlog')
    if PLOT_UNOBSERVED:
        plt.ylim([0, 25000]) # hard-coded just to make the numbers fit in the axis
    else:
        plt.ylim([0, 1000])
    plt.ylabel('Allele count')
    if MAKE_TITLES:
        plt.title('Count/percentage of database alleles observed in cohort')

    if PLOT_POP_LIMIT:
        max_obs = max(total_observations)
        plt.axhline(y=max_obs, linestyle='--', label='Cohort limit')

    if PLOT_POP_LIMIT or PLOT_UNOBSERVED:
        # we only need the legend if we're plotting more than one thing
        plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left')

    out_fn = f'{RESULTS_FOLDER}/observed_haps.png'
    print(f'Saving figure to {out_fn}...')
    plt.savefig(out_fn, bbox_inches='tight')
    plt.close()

def generateZygFig(combined_zyg_counts):
    '''
    Compares the DB size to our observations to see how many alleles we have observed
    @param combined_zyg_counts - [zygosity] -> dict of counts per gene
    '''
    debug = False
    colors = [
        '#1383C6', #bright blue,
        # '#E16A2C', #bright orange
        '#F99D41', #light orange
        '#009D4E', #bright green; default palette has red next, which is bad for RG colorblind
        '#5F249F', #purple
        '#FF66CC', #magenta
        'lightgrey',
        #'#6ABF6A', #light green
        #'red',
        #'green',
        #'blue'
    ]
    
    label_translate = {
        'hets' : 'Heterozygous',
        'homs' : 'Homozygous',
        'hemis' : 'Hemizygous',
        'mt_hets' : 'Heteroplasmic',
        'mt_homs' : 'Homoplasmic',
        'no_matches' : 'No match diplotype'
    }
    stack_order = ['hets', 'homs', 'hemis', 'mt_hets', 'mt_homs', 'no_matches']
    all_genes = set([])
    for s in stack_order:
        all_genes |= combined_zyg_counts[s].keys()
    gene_order = sorted(all_genes)
    gene_labels = []
    
    # we just have observed and unobserved
    plt.figure(figsize=(8,5)) # adding legend made dimensions weird, this gets back to something semi-normal looking
    width = 0.6
    ind = np.arange(len(gene_order))

    if debug:
        # this plot just print the hom / total values
        for g in gene_order:
            het_count = combined_zyg_counts['hets'].get(g, 0)
            hom_count = combined_zyg_counts['homs'].get(g, 0)
            total = het_count+hom_count
            if total > 0:
                print(g, het_count, hom_count, hom_count / total)

    bottoms = np.zeros(shape=(len(gene_order), ))
    for (i, stack_key) in enumerate(stack_order):
        counts = []
        for gene in gene_order:
            value = combined_zyg_counts[stack_key].get(gene, 0)
            counts.append(value)
        
        # plot it and shift the bottoms
        plt.bar(ind, counts, width=width, label=label_translate[stack_key], edgecolor='black', linewidth=1, bottom=bottoms, color=colors[i])
        bottoms += np.array(counts)

    if USE_LATEX:
        gene_ticks = [f'\\textit{{{g}}}' for g in gene_order]
    else:
        gene_ticks = gene_order
    plt.xticks(ind, gene_ticks, rotation=90)
    ax = plt.gca()
    ax.set_axisbelow(True)
    plt.grid(axis='y')
    # plt.yscale('symlog')
    plt.ylabel('Zygosity count')
    if MAKE_TITLES:
        plt.title('Count of zygosity in cohort')
    
    plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left') # anchors to right side

    out_fn = f'{RESULTS_FOLDER}/observed_zygs.png'
    print(f'Saving figure to {out_fn}...')
    plt.savefig(out_fn, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    # first, load all the data into a single dict
    aggregate_files = glob.glob(f'{DATA_FOLDER}/cohort_aggregate_files/*.tsv')
    all_loaded_data = loadAllAggregateFiles(aggregate_files)

    # save the combined aggregate to a file
    aggregate_fn = f'{RESULTS_FOLDER}/combined_aggregate.tsv'
    print(f'Saving combined aggregate to {aggregate_fn}...')
    writeAggregateTsv(all_loaded_data, aggregate_fn)

    display_expected = True # need to have the annotations for this to work
    if display_expected:
        pop_freq_folder = f'{DATA_FOLDER}/AllOfUs_Frequencies_v7/allele'
        popfreqs = loadPopFreq(pop_freq_folder)
        print(f'Loaded population frequencies for {len(popfreqs)} genes.')
    else:
        # empty if we are not loading it
        popfreqs = None
    
    # consolidate everything by ancestry
    reduce_to_core = True
    ancestry_collection, combined_zyg_counts = collectByAncestry(all_loaded_data, reduce_to_core)

    generateZygFig(combined_zyg_counts)
    
    hap_aggregate_fn = f'{RESULTS_FOLDER}/haplotype_aggregate.tsv'
    print(f'Saving combined haplotype aggregate to {hap_aggregate_fn}...')
    writeAggregateHaps(ancestry_collection, hap_aggregate_fn)

    # now generate ancestry images
    generateAncestryPlots(ancestry_collection, popfreqs)

    # now lets look at database representation
    database_fn = f'{DATA_FOLDER}/starphase_db/v0.14.1/pbstarphase_20240826.json.gz'
    db_haps = loadDatabaseHaplotypes(database_fn, reduce_to_core)

    generateDbRep(db_haps, ancestry_collection)

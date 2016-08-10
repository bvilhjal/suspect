"""
Various prediction code... for now.
"""

import sys
import os
import ldpred
from ldpred import ld
from ldpred import LDpred
from ldpred import validate
import scipy as sp
from scipy import stats
import pandas as pd
import cPickle
from plinkio import plinkfile
import h5py
import gzip


def predict_chrons_disease(genot_file = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/genoabc_maf0.05', 
                           phenot_file ='/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc/phenadj.pheno',
                           res_dir = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc',
                           ld_radius=50, N=4100):
    """
    Prediction and validation for Doug Speed.
    
    Since this is a 10-fold cross validation, I will use the full data as LD reference..  which is sort of cheating!!!!!!
    """
    
    data_dir = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns'
    
    #Coordinate data??? 
    #Generate LD reference file
    
    
    train_filter_files = [data_dir+'/predabc/train%d'%i for i in range(1,11)]
    test_filter_files = [data_dir+'/predabc/test%d'%i for i in range(1,11)]
    #Iterate over train-test datasets to perform LDpred
    for i in range(10):
        # - Run LDpred on summary stats.. 
        #     1. Generate a coordinated data set.
        coord_filename = res_dir+'/coord_data%i.hdf5'%i
        if not os.path.isfile(coord_filename):
            print "Performing PLINK GWAS"
            # - Calculate summary statistics...
            run_plink_log_reg_assoc(genot_file=genot_file, indiv_filter_file=train_filter_files[i], phenot_file=phenot_file, 
                                    out_file = res_dir+'/assoc_test%d'%(i+1))
            res_file = res_dir+'/assoc_test%d.qassoc'%(i+1)  
            generate_coordinated_data(res_file, phenot_file, genot_file, coord_filename)

        else:
            print 'Coordinated file found at: %s'%coord_filename
        
            #     2. Run LDpred
        ld_file_prefix = res_dir+'/ld_tab%d'%(i+1)
        out_file_prefix = res_dir + '/ldpred_res%d'%(i+1) 
#         
        print "Running LDpred on the %d'th train/test data set."%(i+1)
        run_ldpred(coord_filename, ld_file_prefix, ld_radius, out_file_prefix, N, )
    
        #     3. Validate predictions
        snp_weights_prefix = out_file_prefix
        out_file_prefix = res_dir + '/ldpred_prs%d'%(i+1) 
        validate_pred(snp_weights_prefix, phenot_file, genot_file, test_filter_files[i], out_file_prefix)
        
    #Summarize results
    

def _get_chrom_dict_(loci, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%d'%chrom
        chr_dict[chr_str] = {'sids':[],'snp_indices':[],'positions':[], 'nts':[]}
     
    for i, l in enumerate(loci):
        chrom = l.chromosome
        pos = l.bp_position
        chr_str = 'chrom_%d'%chrom
        chr_dict[chr_str]['sids'].append(l.name)
#         chr_dict[chr_str]['sids'].append('%d_%d'%(chrom,pos))
        chr_dict[chr_str]['snp_indices'].append(i)
        chr_dict[chr_str]['positions'].append(pos)
        chr_dict[chr_str]['nts'].append([l.allele1,l.allele2])
     
    print 'Genotype dictionary filled'
    return chr_dict


def _parse_plink_snps_(genotype_file, snp_indices):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    num_snps = len(snp_indices)
    raw_snps = sp.empty((num_snps,num_individs),dtype='int8')
    #If these indices are not in order then we place them in the right place while parsing SNPs.
    snp_order = sp.argsort(snp_indices)
    ordered_snp_indices = list(snp_indices[snp_order])
    ordered_snp_indices.reverse()
    print 'Iterating over file to load SNPs'
    snp_i = 0
    next_i = ordered_snp_indices.pop()
    line_i = 0
    max_i = ordered_snp_indices[0]
    while line_i <= max_i:
        if line_i < next_i:
            plinkf.next()
        elif line_i==next_i:
            line = plinkf.next()
            snp = sp.array(line, dtype='int8')
            bin_counts = line.allele_counts()
            if bin_counts[-1]>0:
                mode_v = sp.argmax(bin_counts[:2])
                snp[snp==3] = mode_v
            s_i = snp_order[snp_i]
            raw_snps[s_i]=snp
            if line_i < max_i:
                next_i = ordered_snp_indices.pop()
            snp_i+=1
        line_i +=1
    plinkf.close()
    assert snp_i==len(raw_snps), 'Failed to parse SNPs?'
    num_indivs = len(raw_snps[0])
    freqs = sp.sum(raw_snps,1, dtype='float32')/(2*float(num_indivs))
    return raw_snps, freqs


def generate_coordinated_data(res_file,
                              phen_file,
                              genotype_file,
                              hdf5_filename,
                              indiv_filter_file=None,
                              freq_file=None):
    """
    Assumes plink BED files.  Imputes missing genotypes.
    """

    #Parse indiv filter
    #indiv_tab = pd.read_table(indiv_filter_file,delim_whitespace=True, header=None)

    #Parse phenotypes
    phen_tab = pd.read_table(phen_file,delim_whitespace=True, header=None)
        
    #iids = sp.array(indiv_tab[0])
    #phen_tab = phen_tab[phen_tab[0].isin(iids)]
    
    #Parse frequency file
#     freq_tab = pd.read_table(freq_file,delim_whitespace=True)
    
    #Parse summary stats
    res_tab = pd.read_table(res_file,delim_whitespace=True)
    assert sp.all(res_tab.columns == pd.Index([u'CHR', u'SNP', u'BP', u'NMISS', u'BETA', u'SE', u'R2', u'T', u'P'], dtype='object')), 'Plink results are not in the expected format.'

    
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    Y = sp.array(phen_tab[2])
    fids = [s.fid for s in samples]
    iids = [s.iid for s in samples]
    unique_phens = sp.unique(Y)
    if len(unique_phens)==1:
        print 'Unable to find phenotype values.'
        has_phenotype=False
    elif len(unique_phens)==2:
        cc_bins = sp.bincount(Y)
        assert len(cc_bins)==2, 'Problems with loading phenotype'
        print 'Loaded %d controls and %d cases'%(cc_bins[0], cc_bins[1])
        has_phenotype=True
    else:
        print 'Found quantitative phenotype values'
        has_phenotype=True
    risk_scores = sp.zeros(num_individs)
    rb_risk_scores = sp.zeros(num_individs)
    num_common_snps = 0


    h5f = h5py.File(hdf5_filename)
    h5f.create_dataset('y', data=Y)
    h5f.create_dataset('fids', data=fids)
    h5f.create_dataset('iids', data=iids)
    cdg = h5f.create_group('cord_data')

    #Figure out chromosomes and positions by looking at SNPs.  
    loci = plinkf.get_loci()
    plinkf.close()
    gf_chromosomes = [l.chromosome for l in loci] 

    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()
    chr_dict = _get_chrom_dict_(loci, chromosomes)
    
    tot_num_non_matching_nts = 0
    for chrom in chromosomes:
        chrom_res_tab = res_tab.loc[res_tab['CHR']==chrom]
#         chrom_freq_tab = freq_tab.loc[freq_tab['CHR']==chrom]
        chr_str = 'chrom_%d'%chrom
        print 'Working on chromosome: %s'%chr_str
        
        chrom_d = chr_dict[chr_str]

        g_sids = chrom_d['sids']
        g_sid_set = set(g_sids)
        assert len(g_sid_set) == len(g_sids), 'Some duplicates?'

        ss_sids = sp.array(chrom_res_tab['SNP'])
        ss_sid_set = set(ss_sids)
        assert len(ss_sid_set) == len(ss_sids), 'Some duplicates?'
        
        log_odds = sp.array(chrom_res_tab['BETA'])
        ps = sp.array(chrom_res_tab['P'])
        ns = sp.array(chrom_res_tab['NMISS'])
        betas = sp.sign(log_odds) * stats.norm.ppf(ps/2.0)/sp.sqrt(ns)
        positions = sp.array(chrom_d['positions'])
        g_nts = sp.array(chrom_d['nts'])

        assert sp.all(ss_sids==g_sids), 'Uncoordinated data'


        order = sp.argsort(positions)
        betas = betas[order]
        log_odds = log_odds[order]
        ps = ps[order]
        sids = ss_sids[order]                 
        positions = positions[order]
        nts = g_nts[order]
        
        #Parse SNPs
        snp_indices = sp.array(chrom_d['snp_indices'])
        snp_indices = snp_indices[order] #Pinpoint where the SNPs are in the file.
        raw_snps, freqs = _parse_plink_snps_(genotype_file, snp_indices)
#         mafs = sp.minimum(freqs,1-freqs)
        print 'raw_snps.shape=', raw_snps.shape

        snp_stds = sp.sqrt(2*freqs*(1-freqs)) #sp.std(raw_snps, 1) 
        snp_means = freqs*2 #sp.mean(raw_snps, 1)
        print snp_stds, snp_means
        
            
        
        rb_prs = sp.dot(sp.transpose(raw_snps), log_odds)
        if has_phenotype:
            print 'Normalizing SNPs'
            snp_means.shape = (len(raw_snps),1)
            snp_stds.shape = (len(raw_snps),1)
            snps = (raw_snps - snp_means) / snp_stds
            assert snps.shape==raw_snps.shape, 'Aha!'
            snp_stds = snp_stds.flatten()
            snp_means = snp_means.flatten()
            prs = sp.dot(sp.transpose(snps), betas)
            corr = sp.corrcoef(Y, prs)[0, 1]
            print 'PRS correlation for chromosome %d was %0.4f' % (chrom, corr)
            rb_corr = sp.corrcoef(Y, rb_prs)[0, 1]
            print 'Raw effect sizes PRS correlation for chromosome %d was %0.4f' % (chrom, rb_corr)
        
        
        print 'Now storing coordinated data to HDF5 file.'
        ofg = cdg.create_group('chrom_%d' % chrom)
        ofg.create_dataset('raw_snps_ref', data=raw_snps, compression='lzf')
        ofg.create_dataset('snp_stds_ref', data=snp_stds)
        ofg.create_dataset('snp_means_ref', data=snp_means)
        ofg.create_dataset('freqs_ref', data=freqs)
        ofg.create_dataset('ps', data=ps)
        ofg.create_dataset('positions', data=positions)
        ofg.create_dataset('nts', data=nts)
        ofg.create_dataset('sids', data=sp.array(sids.tolist()))
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        ofg.create_dataset('log_odds_prs', data=rb_prs)
        if has_phenotype:
            risk_scores += prs
        rb_risk_scores += rb_prs
        num_common_snps += len(betas)

    if has_phenotype:
        # Now calculate the prediction r^2
        corr = sp.corrcoef(Y, risk_scores)[0, 1]
        rb_corr = sp.corrcoef(Y, rb_risk_scores)[0, 1]
        print 'PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (corr ** 2,corr)
        print 'Log-odds (effects) PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (rb_corr ** 2, rb_corr)
    print 'There were %d SNPs in common' % num_common_snps
    print 'In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts
    print 'Done coordinating genotypes and summary statistics datasets.'



    

    
    
def run_plink_log_reg_assoc(genot_file = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/genoabc_maf0.05', 
                            indiv_filter_file = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc/train1',
                            phenot_file ='/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc/phenadj.pheno',
                            out_file = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc/assoc_test',
                            plink_path = '/Users/bjv/Dropbox/Cloud_folder/programs/plink1.9/plink'):
    """
    Perform GWAS on the set of individuals 
    """
    from subprocess import call
    
    
    plink_args = ['--bfile',genot_file,
                  '--keep',indiv_filter_file,
                  '--freq',
                  '--pheno',phenot_file,
                  '--allow-no-sex',
                  '--assoc',
                  '--out', out_file,
                  ]

    call([plink_path]+plink_args)
    

        
def run_ldpred(coord_file, ld_file_prefix, ld_radius, out_file_prefix, N, 
               ps=[1,0.3,0.1,0.03,0.01,0.003,0.001], num_iter=50, h2=None, 
               verbose=False):
    local_ld_dict_file = '%s_ldradius%d.pickled.gz'%(ld_file_prefix, ld_radius)
    
    print """
Note: For maximal accuracy all SNPs with LDpred weights should be included in the validation data set.
If they are a subset of the validation data set, then we suggest recalculate LDpred for the overlapping SNPs. 
"""
    if not os.path.isfile(local_ld_dict_file):
        df = h5py.File(coord_file)
                 
        chrom_ld_scores_dict = {}
        chrom_ld_dict = {}
        chrom_ref_ld_mats = {}
        ld_score_sum = 0
        num_snps = 0
        print 'Calculating LD information w. radius %d'% ld_radius

        cord_data_g = df['cord_data']

        for chrom_str in cord_data_g.keys():
            print 'Working on %s'%chrom_str
            g = cord_data_g[chrom_str]
            if 'raw_snps_ref' in g.keys():
                raw_snps = g['raw_snps_ref'][...]
                snp_stds = g['snp_stds_ref'][...]
                snp_means = g['snp_means_ref'][...]
            
            
            #Filter monomorphic SNPs
            ok_snps_filter = snp_stds>0
            ok_snps_filter = ok_snps_filter.flatten()
            raw_snps = raw_snps[ok_snps_filter]
            snp_means = snp_means[ok_snps_filter]
            snp_stds = snp_stds[ok_snps_filter]

            n_snps = len(raw_snps)
            snp_means.shape = (n_snps,1)   
            snp_stds.shape = (n_snps,1)   
            
            
            # Normalize SNPs..
            snps = sp.array((raw_snps - snp_means)/snp_stds,dtype='float32')
            assert snps.shape==raw_snps.shape, 'Array Shape mismatch'
            ret_dict = ld.get_LDpred_ld_tables(snps, ld_radius=ld_radius, ld_window_size=2*ld_radius)
            chrom_ld_dict[chrom_str] = ret_dict['ld_dict']
            chrom_ref_ld_mats[chrom_str] = ret_dict['ref_ld_matrices']
            ld_scores = ret_dict['ld_scores']
            chrom_ld_scores_dict[chrom_str] = {'ld_scores':ld_scores, 'avg_ld_score':sp.mean(ld_scores)}
            ld_score_sum += sp.sum(ld_scores)
            num_snps += n_snps
        avg_gw_ld_score = ld_score_sum / float(num_snps)
        ld_scores_dict = {'avg_gw_ld_score': avg_gw_ld_score, 'chrom_dict':chrom_ld_scores_dict}    
        
        print 'Done calculating the LD table and LD score, writing to file:', local_ld_dict_file
        print 'Genome-wide average LD score was:', ld_scores_dict['avg_gw_ld_score']
        ld_dict = {'ld_scores_dict':ld_scores_dict, 'chrom_ld_dict':chrom_ld_dict, 'chrom_ref_ld_mats':chrom_ref_ld_mats}
        f = gzip.open(local_ld_dict_file, 'wb')
        cPickle.dump(ld_dict, f, protocol=2)
        f.close()
        print 'LD information is now pickled.'
    else:
        print 'Loading LD information from file: %s'%local_ld_dict_file
        f = gzip.open(local_ld_dict_file, 'r')
        ld_dict = cPickle.load(f)
        f.close()
    LDpred.ldpred_genomewide(data_file=coord_file, out_file_prefix=out_file_prefix, ps=ps, ld_radius=ld_radius, 
                      ld_dict = ld_dict, n=N, num_iter=num_iter, h2=h2, verbose=verbose)
            


def validate_pred(res_file_prefix, phenotype_file, val_gt_file, indiv_filter_file, out_file_prefix, 
                  ps=[1,0.3,0.1,0.03,0.01,0.003,0.001], ts=[1,0.3,0.1,0.03,0.01,0.003,0.001], cov_file=None, res_format='LDPRED'):
    """
    Validation...
    """
    #Parse indiv filter
    indiv_tab = pd.read_table(indiv_filter_file,delim_whitespace=True, header=None)
    
    

    #Parse phenotypes FIXME!!
    phen_tab = pd.read_table(phenotype_file,delim_whitespace=True, header=None)
    iids = sp.array(indiv_tab[0])
    phen_tab = phen_tab[phen_tab[0].isin(iids)]
    
    phen_map = {}
    for i, r in phen_tab.iterrows():
        phen_map[r[0]] = {'phen':float(r[2])}

    iids = set(phen_map.keys())            
    
                        
#     non_zero_chromosomes = set()

    if cov_file is not None:
        print 'Parsing additional covariates'
        #FIXME parse using pandas!
#         with open(cov_file,'r') as f:
#             num_missing = 0
#             for line in f:
#                 l = line.split()
#                 iid = l[0]
#                 if iid in phen_map:
#                     covariates = map(float,l[1:])
#                     phen_map[iid]['covariates']=covariates
#                 else:
#                     num_missing +=1
#             if num_missing>0:
#                 print 'Unable to find %d iids in phen file!'%num_missing
                    
                    
    
    num_individs = len(phen_map)
    assert num_individs>0, 'No phenotypes were found!' 
    
    if res_format=='LDPRED':
        weights_file = '%s_LDpred-inf.txt'%(res_file_prefix)
        if os.path.isfile(weights_file):
            print ''
            print 'Calculating LDpred-inf risk scores'
            rs_id_map = validate.parse_ldpred_res(weights_file)       
            out_file = '%s_LDpred-inf.txt'%(out_file_prefix)
            validate.calc_risk_scores(val_gt_file, rs_id_map, phen_map, out_file=out_file, split_by_chrom=False, 
                                      adjust_for_sex=False, adjust_for_covariates=False, 
                                      adjust_for_pcs=False)        
        
        for p in ps:
            weights_file = '%s_LDpred_p%0.4e.txt'%(res_file_prefix, p)
            if os.path.isfile(weights_file):
                print ''
                print 'Calculating LDpred risk scores using p=%0.3e'%p
                rs_id_map = validate.parse_ldpred_res(weights_file)       
                out_file = '%s_LDpred_p%0.4e.txt'%(out_file_prefix, p)
                validate.calc_risk_scores(val_gt_file, rs_id_map, phen_map, out_file=out_file, split_by_chrom=False, 
                                      adjust_for_sex=False, adjust_for_covariates=False, 
                                      adjust_for_pcs=False)        
            
        #Plot results?

    elif res_format=='P+T':
        weights_file = '%s_all_snps.txt'%(res_file_prefix)
        if os.path.isfile(weights_file):
            print ''
            print 'Calculating risk scores using all SNPs'
            rs_id_map = validate.parse_ldpred_res(weights_file)       
            out_file = '%s_all_snps.txt'%(out_file_prefix)
            validate.calc_risk_scores(val_gt_file, rs_id_map, phen_map, out_file=out_file, split_by_chrom=False, 
                                      adjust_for_sex=False, adjust_for_covariates=False, 
                                      adjust_for_pcs=False)        
        
        for p_thres in ts:
            weights_file = '%s_P+T_p%0.4e.txt'%(res_file_prefix, p_thres)
            print weights_file
            if os.path.isfile(weights_file):
                print ''
                print 'Calculating P+T risk scores using p-value threshold of %0.3e'%p_thres
                rs_id_map = validate.parse_pt_res(weights_file)       
                out_file = '%s_P+T_p%0.4e.txt'%(out_file_prefix, p_thres)
                validate.calc_risk_scores(val_gt_file, rs_id_map, phen_map, out_file=out_file, split_by_chrom=False, 
                                          adjust_for_sex=False, adjust_for_covariates=False, 
                                          adjust_for_pcs=False)            
        #Plot results?
    else:
        raise NotImplementedError('Results file format missing or unknown: %s'%res_format)
    
    

    
    
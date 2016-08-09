"""
Various prediction code... for now.
"""

import sys
import ldpred
import scipy as sp
from scipy import stats
import pandas as pd
from plinkio import plinkfile
import h5py


def predict_chrons_disease(genot_file = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/genoabc_maf0.05', 
                           phenot_file ='/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc/phenadj.pheno',
                           res_dir = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc'):
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
        # - Calculate summary statistics...
        print "Running LDpred on the %d'th train/test data set."%(i+1)
        print "Performing PLINK GWAS"
        run_plink_log_reg_assoc(genot_file=genot_file, indiv_filter_file=train_filter_files[i], phenot_file=phenot_file, 
                                out_file = res_dir+'/assoc_test%d'%(i+1))
        res_file = res_dir+'/assoc_test%d.qassoc'%(i+1)  
        hdf5_filename = res_dir+'/coord_data%i.hdf5'%i
       
        
        # - Run LDpred on summary stats.. 
        #     1. Generate a coordinated data set.
        generate_coordinated_data(res_file, phenot_file, genot_file, hdf5_filename)
        
        #     2. Run LDpred
        
    
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
        betas = sp.sign(log_odds) * stats.norm.ppf(ps/2.0)


        if not sp.all(ss_sids==g_sids):
            print 'Coordinating SNPs'
            g_filter = sp.in1d(g_sids,ss_sids)
            ss_filter = sp.in1d(ss_sids,g_sids)

            #Order by SNP IDs
            g_order = sp.argsort(g_sids)
            ss_order = sp.argsort(ss_sids)

            g_indices = []
            for g_i in g_order:
                if g_filter[g_i]:
                    g_indices.append(g_i)
    
            ss_indices = []
            for ss_i in ss_order:
                if ss_filter[ss_i]:
                    ss_indices.append(ss_i)

            g_nts = chrom_d['nts']
            snp_indices = chrom_d['snp_indices']
            betas = sp.array(chrom_res_tab['BETA'])
            assert not sp.any(sp.isnan(betas)) and not sp.any(sp.isinf(betas)), 'Some betas are inf or nan.'
            
            #Convert p-values to betas...FIXME

            print 'Found %d SNPs present in both datasets'%(len(g_indices))

            ok_nts = g_nts[g_indices]


        #Resorting by position
        positions = sp.array(chrom_d['positions'])[g_indices]
        order = sp.argsort(positions)
        g_indices = list(sp.array(g_indices)[order])
        ss_indices = list(sp.array(ss_indices)[order])
        positions = positions[order]
        nts = ok_nts[order]
        
        #Parse SNPs
        snp_indices = sp.array(chrom_d['snp_indices'])
        snp_indices = snp_indices[g_indices] #Pinpoint where the SNPs are in the file.
        raw_snps, freqs = _parse_plink_snps_(genotype_file, snp_indices)
#         mafs = sp.minimum(freqs,1-freqs)
        print 'raw_snps.shape=', raw_snps.shape

        snp_stds = sp.sqrt(2*freqs*(1-freqs)) #sp.std(raw_snps, 1) 
        snp_means = freqs*2 #sp.mean(raw_snps, 1)

        betas = betas[ss_indices]
        log_odds = log_odds[ss_indices]
        ps = ps[ss_indices]
        sids = ss_sids[ss_indices]                 
        
            
        
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
        ofg.create_dataset('sids', data=sids)
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
                            out_file = '/Users/bjv/Dropbox/Cloud_folder/Data/bjarni_crohns/predabc/assoc_test'):
    """
    Perform GWAS on the set of individuals 
    """
    from subprocess import call
    plink_path = './../../programs/plink1.9/plink' 
    
    plink_args = ['--bfile',genot_file,
                  '--keep',indiv_filter_file,
                  '--freq',
                  '--pheno',phenot_file,
                  '--allow-no-sex',
                  '--assoc',
                  '--out', out_file,
                  ]

    call([plink_path]+plink_args)
    

        
            
            
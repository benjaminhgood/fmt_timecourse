import parse_midas_data
import config
import numpy
import gene_diversity_utils

treatment_subjects = set(['FAT_006','FAT_015','FAT_008','FAT_020','FAT_012'])
control_subjects = set(['FAT_010','FAT_014','FAT_017','FAT_023','FAT_024'])

pretty_name_map = { 'FAT_006': 'FMT1',
                    'FAT_015': 'FMT2',
                    'FAT_008': 'FMT3',
                    'FAT_020': 'FMT4',
                    'FAT_012': 'FMT5',
                    'FAT_010': 'SELF1',
                    'FAT_014': 'SELF2',
                    'FAT_017': 'SELF3',
                    'FAT_023': 'SELF4',
                    'FAT_024': 'SELF5' }

pretty_time_map = {0:'Host',
                   1:'Donor',
                   2:'Day 2',
                   14:'Day 14',
                   42:'Day 42',
                   84:'Day 84'}

subject_donor_map = { 'FAT_006': 'FAT_DON_11',
                      'FAT_015': 'FAT_DON_11',
                      'FAT_008': 'FAT_DON_11',
                      'FAT_020': 'FAT_DON_19',
                      'FAT_012': 'FAT_DON_8',
                      'FAT_010': 'FAT_010',
                      'FAT_014': 'FAT_014',
                      'FAT_017': 'FAT_017',
                      'FAT_023': 'FAT_023',
                      'FAT_024': 'FAT_024' }
 
time_map = {0:0, 1:1, 2:2, 14: 3, 42: 4, 84:5} 
                      
def parse_sample_map():
    
    sample_map = {'control': {}, 'treatment': {}, 'donor': {}}
    sample_strs = set([])
    
    file = open(parse_midas_data.scripts_directory+"li_etal_metadata.txt","r")
    file.readline()
    for line in file:
        items = line.split("\t")
        sample_str = items[2].strip()
        sample_strs.add(sample_str)
        
    file.close()
    
    
    for sample_str in sample_strs:
        
        sample_items = sample_str.split("-")
        
        subject = sample_items[0].strip()
    
        # first item is FAT
        if "DON" in subject:
            # donor strain!
            sample_map['donor'][subject] = sample_str
            
        else:
            # FMT strain        
            
            
            time = long(sample_items[2])
            
        
            if subject in control_subjects:
                # control!
                treatment = 'control'
            else:
                treatment = 'treatment'
                
            if subject not in sample_map[treatment]:
                sample_map[treatment][subject] = {}
            
            sample_map[treatment][subject][time] = sample_str
    
    
    # add in donor
    for subject in sample_map['control'].keys():
        sample_map['control'][subject][1] = sample_map['control'][subject][0]
        
    for subject in sample_map['treatment'].keys():    
        sample_map['treatment'][subject][1] = sample_map['donor'][subject_donor_map[subject]]
            
    return sample_map
            
    
    
###############################################################################
#
# Loads a subset of "core" genes using copynum information in the genes/ folder 
#
###############################################################################   
def load_core_timecourse_genes(desired_species_name, min_copynum=0.3, min_prevalence=0.9, min_marker_coverage=20):

    # Load subject and sample metadata
    subject_sample_map = parse_subject_sample_map()
    sample_time_map = parse_sample_time_map()
    
    desired_samples = set(subject_sample_map[focal_patient].keys())
    
    # Load reference genes
    reference_genes = parse_midas_data.load_reference_genes(desired_species_name)
    
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(desired_species_name)
    
    gene_names = numpy.array(gene_names)
   
    reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])

    sample_idxs = numpy.array([sample_name in desired_samples for sample_name in gene_samples])*(marker_coverages>=min_marker_coverage)
    
    if sample_idxs.sum()>0:   
        prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,sample_idxs], marker_coverages[sample_idxs], min_copynum)
        core_gene_idxs = reference_gene_idxs*(prevalences>=min_prevalence)  
    else:
        sys.stderr.write("Not enough samples for core genome!\n")
        reference_gene_idxs

    return set(gene_names[core_gene_idxs])

def calculate_timecourse_idxs(sample_map, desired_subject, samples):
    
    sample_list = list(samples)
    

    # TODO
    idxs = []
    times = []
    
    if desired_subject in sample_map['control']:
        treatment = 'control'
    else:
        treatment = 'treatment'
        
    for time in sorted(sample_map[treatment][desired_subject].keys()):
        sample = sample_map[treatment][desired_subject][time]
        if sample in sample_list:
            idx = sample_list.index(sample)
            idxs.append(idx)
            times.append(time)
    
    new_times = [time_map[t] for t in times]
    
    return numpy.array(new_times), numpy.array(idxs)
    
                
if __name__=='__main__':
    sample_map = parse_sample_map()
    print sample_map 
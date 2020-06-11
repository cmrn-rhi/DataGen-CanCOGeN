
"""
Generate COVID19 CanCOGen Metadata

Version: 1.6
COVID-19 Vocab Version #39
Date: 2020/06/10

Python v 3.8.3

@author: Rhiannon Cameron
         orcid.org/0000-0002-9578-0788
"""

import library as lib
import pandas as pd
import os
from sys import argv

# generate dictionary with CanCOGeN labels as keys.
covid_vars = lib.generate_covid19_vars()
# generate list of variable names.
covid_var_names = []
# associate variable names values with CanCOGeN label keys.
for key in covid_vars:
  covid_var_names.append(covid_vars[key])


def make_row(invalid_data=[]):
  """
  Generate Row of CanCOGeN Metadata
  
  return: a list of one rows of fake metadata.
  """

  ### Database Identifiers ###
  
  # specimen_collector_sample_id moved to "Dependent IDs"
  phac_sample_id = lib.random_phac_id()
  # irida_sample_name moved to "Dependent IDs"
  umbrella_bioproject_accession = lib.umbrella_bioproject_accession()
  bioproject_accession = lib.random_bioproject_accession()
  biosample_accession = lib.random_biosample_accession()
  sra_accession = lib.random_sra_accession()
  genbank_accession = lib.random_genbank_accession()
  gisaid_accession = lib.random_gisaid_accession()
  
  ### Sample Collection and Processing ###
  
  sample_collected_by = lib.random_agency()
  sample_collector_contact_email = lib.random_email()
  sample_collector_contact_address = lib.random_address()
  sequence_submitter_contact_email = lib.random_email()
  sequence_submitter_contact_address = lib.random_address()
  sample_collection_date = lib.random_date()
  sample_received_date = lib.random_date()
  geo_loc_name_country = lib.random_country()
  geo_loc_name_province_territory = lib.random_province_territory()
  geo_loc_name_city = lib.random_city()
  organism = lib.random_organism()
  # isolate moved to "Dependent IDs"
  purpose_of_sampling = lib.random_purpose_of_sampling()
  anatomical_material = lib.random_anatomical_material()
  anatomical_part = lib.random_anatomical_part()
  body_product = lib.random_body_product()
  environmental_material = lib.random_environmental_material()
  environmental_site = lib.random_environmental_site()
  collection_device = lib.random_collection_device()
  collection_method = lib.random_collection_method()
  collection_protocol = lib.fake_protocol()
  specimen_processing = lib.random_specimen_processing()
  lab_host = lib.random_lab_host()
  passage_number = lib.random_passage_number()
  passage_method = lib.passage_method_text()
  biomaterial_extracted = lib.random_biomaterial_extracted()
  
  ### Host Information ###
  
  host_common_name = lib.random_host_common_name()
  host_scientific_name = lib.random_host_scientific_name()
  host_health_state = lib.random_host_health_state()
  host_health_status_details = lib.random_host_health_status_details()
  host_disease = lib.random_host_disease()
  host_age, host_age_error = lib.random_host_age(invalid_data)
  host_gender = lib.random_host_gender()
  host_origin_geo_loc_name_country = lib.random_country()
  host_subject_id = lib.random_host_subject_id()
  symptom_onset_date = lib.random_date()
  signs_and_symptoms = lib.random_signs_symptoms()
  
  ### Host Exposure Information ###
  
  location_of_exposure_geo_loc_name_country = lib.random_country()
  travel_history = lib.random_travel_history()
  exposure_event = lib.random_exposure_event()
  
  ### Sequencing ###
  
  minion_barcode = lib.random_minIon_barcode()
  sequencing_instrument = lib.random_seq_instrument()
  sequencing_protocol_name = lib.fake_protocol()
  sequencing_protocol_source = lib.random_seq_protocol_source()
  sequencing_kit_number = lib.random_alphanumeric()
  amplicon_pcr_primers_filename = lib.random_txt_filename()
  
  ### Bioinformatics and QC metrics ###
  
  raw_sequence_data_processing = lib.random_seq_process()
  sequencing_depth_average = lib.random_seq_depth()
  # assembly_name moved to "Dependent IDs"
  assembly_method = lib.random_assembly_software()
  assembly_coverage_breadth = lib.random_assembly_coverage_breadth()
  assembly_coverage_depth = lib.random_seq_depth()
  r1_fastq_filename = lib.random_fastq_filename()
  r2_fastq_filename = lib.random_fastq_filename()
  r1_fastq_filepath = lib.random_filepath(r1_fastq_filename)
  r2_fastq_filepath = lib.random_filepath(r2_fastq_filename)
  fast5_filename = lib.random_fast5_filename()
  fast5_filepath = lib.random_filepath(fast5_filename)
  fasta_filename = lib.random_fasta_filename()
  fasta_filepath = lib.random_filepath(fasta_filename)
  number_base_pairs = lib.random_bp_num()
  consensus_genome_length = lib.random_genome_length()
  mean_contig_length = lib.random_contig_length()
  n50 = lib.random_n50()
  ns_per_100_kbp = lib.random_ns_100kbp()
  reference_genome_accession = lib.random_ref_genome()
  consensus_sequence_id = lib.random_consensus_seq_id()
  consensus_sequence_method = lib.random_consensus_seq_method()
  consensus_sequence_filename = lib.random_fasta_filename()
  consensus_sequence_filepath = lib.random_filepath(consensus_sequence_filename)
  annotation_feature_table_filename = lib.random_feature_table_filename()
  bioinformatics_protocol = lib.bioinformatics_protocol()
  
  ### Pathogen Diagnostic Testing ###
  
  gene_name_1 = lib.random_gene()
  diagnostic_pcr_protocol_1 = lib.fake_protocol()
  diagnostic_pcr_ct_value_1 = lib.random_pcr_ct_val()
  gene_name_2 = lib.random_gene()
  diagnostic_pcr_protocol_2 = lib.fake_protocol()
  diagnostic_pcr_ct_value_2 = lib.random_pcr_ct_val()
  
  ### Contributor Acknowledgement ###
  
  authors = lib.authors()
  
  ### Dependent IDs ###
  
  specimen_collector_sample_id = lib.random_specimen_collector_sample_id(
                                 geo_loc_name_country, 
                                 geo_loc_name_province_territory,
                                 geo_loc_name_city)
  irida_sample_name = specimen_collector_sample_id
  sequence_submitted_by = sample_collected_by
  library_id = lib.random_library_id(specimen_collector_sample_id)
  isolate = specimen_collector_sample_id
  assembly_name = lib.random_assembly_name(specimen_collector_sample_id)

  cols = [specimen_collector_sample_id, 
          phac_sample_id, 
          irida_sample_name, 
          umbrella_bioproject_accession, 
          bioproject_accession, 
          biosample_accession, 
          sra_accession, 
          genbank_accession, 
          gisaid_accession, 
          sample_collected_by, 
          sample_collector_contact_email, 
          sample_collector_contact_address, 
          sequence_submitted_by, 
          sequence_submitter_contact_email, 
          sequence_submitter_contact_address, 
          sample_collection_date, 
          sample_received_date, 
          geo_loc_name_country, 
          geo_loc_name_province_territory, 
          geo_loc_name_city, 
          organism, 
          isolate, 
          purpose_of_sampling, 
          anatomical_material, 
          anatomical_part, 
          body_product, 
          environmental_material, 
          environmental_site, 
          collection_device, 
          collection_method, 
          collection_protocol, 
          specimen_processing, 
          lab_host, 
          passage_number, 
          passage_method, 
          biomaterial_extracted, 
          host_common_name, 
          host_scientific_name, 
          host_health_state, 
          host_health_status_details, 
          host_disease, 
          host_age, 
          host_gender, 
          host_origin_geo_loc_name_country, 
          host_subject_id, 
          symptom_onset_date, 
          signs_and_symptoms, 
          location_of_exposure_geo_loc_name_country, 
          travel_history, 
          exposure_event, 
          library_id, 
          minion_barcode, 
          sequencing_instrument, 
          sequencing_protocol_name, 
          sequencing_protocol_source, 
          sequencing_kit_number, 
          amplicon_pcr_primers_filename, 
          raw_sequence_data_processing, 
          sequencing_depth_average, 
          assembly_name, 
          assembly_method, 
          assembly_coverage_breadth, 
          assembly_coverage_depth, 
          r1_fastq_filename, 
          r2_fastq_filename, 
          r1_fastq_filepath, 
          r2_fastq_filepath, 
          fast5_filename, 
          fast5_filepath, 
          fasta_filename, 
          fasta_filepath, 
          number_base_pairs, 
          consensus_genome_length, 
          mean_contig_length, 
          n50, 
          ns_per_100_kbp, 
          reference_genome_accession, 
          consensus_sequence_id, 
          consensus_sequence_method, 
          consensus_sequence_filename, 
          consensus_sequence_filepath, 
          annotation_feature_table_filename, 
          bioinformatics_protocol, 
          gene_name_1, 
          diagnostic_pcr_protocol_1, 
          diagnostic_pcr_ct_value_1, 
          gene_name_2, 
          diagnostic_pcr_protocol_2, 
          diagnostic_pcr_ct_value_2, 
          authors]

  error_grid = [specimen_collector_sample_id, 
          phac_sample_id, 
          irida_sample_name, 
          umbrella_bioproject_accession, 
          bioproject_accession, 
          biosample_accession, 
          sra_accession, 
          genbank_accession, 
          gisaid_accession, 
          sample_collected_by, 
          sample_collector_contact_email, 
          sample_collector_contact_address, 
          sequence_submitted_by, 
          sequence_submitter_contact_email, 
          sequence_submitter_contact_address, 
          sample_collection_date, 
          sample_received_date, 
          geo_loc_name_country, 
          geo_loc_name_province_territory, 
          geo_loc_name_city, 
          organism, 
          isolate, 
          purpose_of_sampling, 
          anatomical_material, 
          anatomical_part, 
          body_product, 
          environmental_material, 
          environmental_site, 
          collection_device, 
          collection_method, 
          collection_protocol, 
          specimen_processing, 
          lab_host, 
          passage_number, 
          passage_method, 
          biomaterial_extracted, 
          host_common_name, 
          host_scientific_name, 
          host_health_state, 
          host_health_status_details, 
          host_disease, 
          host_age_error, 
          host_gender, 
          host_origin_geo_loc_name_country, 
          host_subject_id, 
          symptom_onset_date, 
          signs_and_symptoms, 
          location_of_exposure_geo_loc_name_country, 
          travel_history, 
          exposure_event, 
          library_id, 
          minion_barcode, 
          sequencing_instrument, 
          sequencing_protocol_name, 
          sequencing_protocol_source, 
          sequencing_kit_number, 
          amplicon_pcr_primers_filename, 
          raw_sequence_data_processing, 
          sequencing_depth_average, 
          assembly_name, 
          assembly_method, 
          assembly_coverage_breadth, 
          assembly_coverage_depth, 
          r1_fastq_filename, 
          r2_fastq_filename, 
          r1_fastq_filepath, 
          r2_fastq_filepath, 
          fast5_filename, 
          fast5_filepath, 
          fasta_filename, 
          fasta_filepath, 
          number_base_pairs, 
          consensus_genome_length, 
          mean_contig_length, 
          n50, 
          ns_per_100_kbp, 
          reference_genome_accession, 
          consensus_sequence_id, 
          consensus_sequence_method, 
          consensus_sequence_filename, 
          consensus_sequence_filepath, 
          annotation_feature_table_filename, 
          bioinformatics_protocol, 
          gene_name_1, 
          diagnostic_pcr_protocol_1, 
          diagnostic_pcr_ct_value_1, 
          gene_name_2, 
          diagnostic_pcr_protocol_2, 
          diagnostic_pcr_ct_value_2, 
          authors]


  return cols, error_grid
  #return cols, None


def check_file_name(file_name):
  """
  Check For Existing Filename
  credit: lastro
  
  given: string containing file name.
  return: string with additional number if file already exists with given 
          filename; otherwise returns original string.
  """
  if os.path.isfile(file_name):
      expand = 1
      
      while True:
          expand += 1      
          
          if ".tsv" in file_name:
            new_file_name = (file_name.split(".tsv")[0] + "-" + str(expand) + 
                             ".tsv")
            if os.path.isfile(new_file_name):
                continue
            else:
                return new_file_name

          elif ".csv" in file_name:
            new_file_name = (file_name.split(".csv")[0] + "-" + str(expand) + 
                             ".csv")
            if os.path.isfile(new_file_name):
                continue
            else:
                return new_file_name
              
  return file_name
              

def generate_data_files(file_name, rows, delimiter, invalid_data=[]):
  """
  Generate Sample Metadata File
  
  given: string containing file name, number of rows to generate, and delimiter.
  return: csv or tsv delimited data file.
  """
  # check that filename will not overwrite existing file, append number to the 
  # end if filename already exists.
  if delimiter == 'tab':
    # save as tsv if using tab delimiter.
    file_name = check_file_name(file_name + '.tsv')
    # create filename for error grid.
    error_file_name = (file_name.split(".tsv")[0] + "_error-grid.tsv")
  else:
    # save as csv if using comma delimiter.
    file_name = check_file_name(file_name + '.csv')
    # create filename for error grid.
    error_file_name = (file_name.split(".csv")[0] + "_error-grid.csv")
  
  # CanCOGeN category headers.
  category_headers = ['Database Identifiers','','','','','','','','',
                      'Sample collection and processing','','','','','','','',
                      '','','','','','','','','','','','','','','','','',
                      '','','Host Information','','','','','','','','','',
                      '','Host exposure information','','','Sequencing','',
                      '','','','','','Bioinformatics and QC metrics','','','',
                      '','','','','','','','','','','','','','','','','','','',
                      '','','','Pathogen diagnostic testing','','','','','',
                      'Contributor acknowledgement']
  # Error Grid header row.
  error_grid_header = [file_name.capitalize() + ' Error Grid','','','','',
                       '','','','','','','','','','','','','','','','','','',
                       '','','','','','','','','','','','','','','','','','',
                       '','','','','','','','','','','','','','','','','','',
                       '','','','','','','','','','','','','','','','','','',
                       '','','','','','','','','','','','','']
  
  # make category headers first row for data file.
  df = pd.DataFrame(columns=category_headers)
  # make error_grid_header first row for error grid file.
  error_df = pd.DataFrame(columns=error_grid_header)
  
  # make column headers second row for both files.
  df.loc[0] = covid_var_names  
  error_df.loc[0] = covid_var_names  
    
  # generate rows of data, add one to avoid overwriting column headers row.
  for i in range(rows):
    df.loc[i+1], error_df.loc[i+1] = make_row(invalid_data)
  
  if delimiter == 'tab':
    return (df.to_csv(file_name, sep='\t', index=False), 
            error_df.to_csv(error_file_name, sep='\t', index=False))
  else:
    return (df.to_csv(file_name, sep=',', index=False),
            error_df.to_csv(error_file_name, sep=',', index=False))

# for running script from command line.
#generate_data_files(argv[1], int(argv[2]), argv[3])
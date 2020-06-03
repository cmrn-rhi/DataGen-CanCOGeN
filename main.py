# -*- coding: utf-8 -*-
"""
Generate COVID19 CanCOGen Metadata

Version: 1.0
COVID-19 Template Version: 0.95
Date: 2020/06/01

@author: Rhiannon Cameron
         orcid.org/0000-0002-9578-0788
"""

import library as lib
import pandas as pd
import random
import os

# generate dictionary with CanCOGeN labels as keys.
covid_vars = lib.generate_covid19_vars()
# generate list of variable names.
covid_var_names = []
# associate variable names values with CanCOGeN label keys.
for key in covid_vars:
  covid_var_names.append(covid_vars[key])


def make_row():
  """
  Generate Row of Metadata
  
  return: a list of one rows of fake CanCOGeN metadata.
  """
  # Database Identifiers
  
  # specimen_collector_sample_id moved to "Dependent IDs"
  phac_sample_id = lib.random_phac_id()
  # irida_sample_name moved to "Dependent IDs"
  umbrella_bioproject_accession = lib.umbrella_bioproject_accession()
  bioproject_accession = lib.random_bioproject_accession()
  biosample_accession = lib.random_biosample_accession()
  sra_accession = lib.random_sra_accession()
  genbank_accession = lib.random_genbank_accession()
  gisaid_accession = lib.random_gisaid_accession()
  
  # Sample Collection and Processing
  
  sample_collected_by = lib.random_name()
  sample_collector_contact_email = lib.random_email()
  sample_collector_contact_address = lib.random_address()
  sequence_submitted_by = lib.random_name()
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
  
  # Host Information
  
  host_common_name = lib.random_host_common_name()
  host_scientific_name = lib.random_host_scientific_name()
  host_health_state = lib.random_host_health_state()
  host_health_status_details = lib.random_host_health_status_details()
  host_disease = lib.random_host_disease()
  host_age = lib.random_host_age()
  host_gender = lib.random_host_gender()
  host_origin_geo_loc_name_country = lib.random_country()
  host_subject_id = lib.random_host_subject_id()
  symptom_onset_date = lib.random_date()
  signs_and_symptoms = lib.random_signs_symptoms()
  
  # Host Exposure Information
  
  location_of_exposure_geo_loc_name_country = lib.random_country()
  travel_history = lib.random_travel_history()
  exposure_event = lib.random_exposure_event()
  
  # Sequencing
  
  minion_barcode = lib.random_minIon_barcode()
  sequencing_instrument = lib.random_seq_instrument()
  sequencing_protocol_name = lib.fake_protocol()
  sequencing_protocol_source = lib.random_seq_protocol_source()
  sequencing_kit_number = lib.random_alphanumeric()
  amplicon_pcr_primers_filename = lib.random_txt_filename()
  
  # Bioinformatics and QC metrics
  
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
  
  # Pathogen Diagnostic Testing
  
  gene_name_1 = lib.random_gene()
  diagnostic_pcr_protocol_1 = lib.fake_protocol()
  diagnostic_pcr_ct_value_1 = lib.random_pcr_ct_val()
  gene_name_2 = lib.random_gene()
  diagnostic_pcr_protocol_2 = lib.fake_protocol()
  diagnostic_pcr_ct_value_2 = lib.random_pcr_ct_val()
  
  # Contributor Acknowledgement
  
  authors = lib.authors()
  
  # Dependent IDs
  
  specimen_collector_sample_id = lib.random_specimen_collector_sample_id(
                                 geo_loc_name_country, 
                                 geo_loc_name_province_territory,
                                 geo_loc_name_city)
  irida_sample_name = specimen_collector_sample_id
  library_id = lib.random_library_id(specimen_collector_sample_id)
  isolate = specimen_collector_sample_id
  assembly_name = lib.random_assembly_name(specimen_collector_sample_id)

  return [specimen_collector_sample_id, 
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
          new_file_name = (file_name.split(".csv")[0] + "-" + str(expand) + 
                           ".csv")
          if os.path.isfile(new_file_name):
              continue
          else:
              return new_file_name
              
  return file_name


def generate_data_file(file_name, rows, delimiter):
  """
  Generate Sample Metadata File
  
  given: string containing file name, number of rows to generate, and delimiter..
  return: csv delimited data file.
  """
  
  file_name = check_file_name(file_name)
  
  df = pd.DataFrame(columns=covid_var_names)
  
  for i in range(rows):
    df.loc[i] = make_row()
  
  if delimiter == 'tab':
    return df.to_csv(file_name, sep='\t', index=False)
  else:
    return df.to_csv(file_name, sep=',', index=False)
  

#generate_data_file('COVID19_Metadata_v0.95_Sample_Complete_v3.0.csv', 
  #                 50, 'comma')
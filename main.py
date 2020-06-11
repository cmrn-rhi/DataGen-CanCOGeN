
"""
Generate COVID19 CanCOGen Metadata

Version: 1.7
COVID-19 Vocab Version #39
Date: 2020/06/11

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
  
  # phac_sample_id, valid/invalid error-grid information.
  phac_sample_id, phac_sample_id_error = lib.random_phac_id(invalid_data)
  # umbrella_bioproject_accession, valid/invalid error-grid information.
  ub_accession, ub_accession_error = lib.umbrella_bioproject_accession(invalid_data)
  # bioproject_accession, valid/invalid error-grid information.
  bp_accession, bp_accession_error = lib.random_bioproject_accession(invalid_data)
  # biosample_accession, valid/invalid error-grid information.
  bs_accession, bs_accession_error = lib.random_biosample_accession(invalid_data)
  # sra_accession, valid/invalid error-grid information.
  sra_accession, sra_accession_error = lib.random_sra_accession(invalid_data)
  # genbank_accession, valid/invalid error-grid information.
  gb_accession, gb_accession_error = lib.random_genbank_accession(invalid_data)
  # gisaid_accession, valid/invalid error-grid information.
  gisaid_accession, gisaid_accession_error = lib.random_gisaid_accession(invalid_data)
  
  ### Sample Collection and Processing ###
  
  # sample_collected_by, valid/invalid error-grid information.
  samp_col_by, samp_col_by_error = lib.random_agency(invalid_data)
  # sample_collector_contact_email, valid/invalid error-grid information.
  samp_col_email, samp_col_email_error = lib.random_email(invalid_data)
  # sample_collector_contact_address, valid/invalid error-grid information.
  samp_col_address, samp_col_address_error = lib.random_address(invalid_data)
  # sequence_submitter_contact_email, valid/invalid error-grid information.
  seq_sub_email, seq_sub_email_error = lib.random_email(invalid_data)
  # sequence_submitter_contact_address, valid/invalid error-grid information.
  seq_sub_address, seq_sub_address_error = lib.random_address(invalid_data)
  # sample_collection_date, valid/invalid error-grid information.
  samp_col_date, samp_col_date_error = lib.random_date(invalid_data)
  # sample_received_date, valid/invalid error-grid information.
  samp_rec_date, samp_rec_date_error = lib.random_date(invalid_data)
  # geo_loc_name_country, valid/invalid error-grid information.
  geo_loc_country, geo_loc_country_error = lib.random_country(invalid_data)
  # geo_loc_name_province_territory, valid/invalid error-grid information.
  geo_loc_prov_ter, geo_loc_prov_ter_error = lib.random_province_territory(invalid_data)
  # geo_loc_name_city, valid/invalid error-grid information.
  geo_loc_city, geo_loc_city_error = lib.random_city(invalid_data)
  # organism, valid/invalid error-grid information.
  organism, organism_error = lib.random_organism(invalid_data)
  # purpose_of_sampling, valid/invalid error-grid information.
  p_o_sampling, p_o_sampling_error = lib.random_purpose_of_sampling(invalid_data)
  # anatomical_material, valid/invalid error-grid information.
  anat_material, anat_material_error = lib.random_anatomical_material(invalid_data)
  # anatomical_part, valid/invalid error-grid information.
  anat_part, anat_part_error = lib.random_anatomical_part(invalid_data)
  # body_product, valid/invalid error-grid information.
  body_product, body_product_error = lib.random_body_product(invalid_data)
  # environmental_material, valid/invalid error-grid information.
  envi_material, envi_material_error = lib.random_environmental_material(invalid_data)
  # environmental_site, valid/invalid error-grid information.
  envi_site, envi_site_error = lib.random_environmental_site(invalid_data)
  # collection_device, valid/invalid error-grid information.
  col_device, col_device_error = lib.random_collection_device(invalid_data)
  # collection_method, valid/invalid error-grid information.
  col_method, col_method_error = lib.random_collection_method(invalid_data)
  # collection_protocol, valid/invalid error-grid information.
  col_protocol, col_protocol_error = lib.fake_protocol(invalid_data)
  # specimen_processing, valid/invalid error-grid information.
  spec_process, spec_process_error = lib.random_specimen_processing(invalid_data)
  # lab_host, valid/invalid error-grid information.
  lab_host, lab_host_error = lib.random_lab_host(invalid_data)
  # passage_number , valid/invalid error-grid information.
  passage_num, passage_num_error = lib.random_passage_number(invalid_data)
  # passage_method, valid/invalid error-grid information.
  passage_method, passage_method_error = lib.passage_method_text(invalid_data)
  # biomaterial_extracted, valid/invalid error-grid information.
  biom_extract, biom_extract_error = lib.random_biomaterial_extracted(invalid_data)
  
  ### Host Information ###
  
  # host_common_name, valid/invalid error-grid information.
  host_com_name, host_com_name_error = lib.random_host_common_name(invalid_data)
  # host_scientific_name, valid/invalid error-grid information.
  host_sci_name, host_sci_name_error = lib.random_host_scientific_name(invalid_data)
  # host_health_state, valid/invalid error-grid information.
  host_health_state, host_health_state_error = lib.random_host_health_state(invalid_data)
  # host_health_status_details, valid/invalid error-grid information.
  host_health_status, host_health_status_error = lib.random_host_health_status_details(invalid_data)
  # host_disease, valid/invalid error-grid information.
  host_disease, host_disease_error = lib.random_host_disease(invalid_data)
  # host_age, valid/invalid error-grid information.
  host_age, host_age_error = lib.random_host_age(invalid_data)
  # host_gender, valid/invalid error-grid information.
  host_gender, host_gender_error = lib.random_host_gender(invalid_data)
  # host_origin_geo_loc_name_country, valid/invalid error-grid information.
  host_loc_country, host_loc_country_error = lib.random_country(invalid_data)
  # host_subject_id, valid/invalid error-grid information.
  host_sub_id, host_sub_id_error = lib.random_host_subject_id(invalid_data)
  # symptom_onset_date, valid/invalid error-grid information.
  symp_onset_date, symp_onset_date_error = lib.random_date(invalid_data)
  # signs_and_symptoms, valid/invalid error-grid information.
  signs_symptoms, signs_symptoms_error = lib.random_signs_symptoms(invalid_data)
  
  ### Host Exposure Information ###
  
  # location_of_exposure_geo_loc_name_country, valid/invalid error-grid information.
  loc_exp_country, loc_exp_country_error = lib.random_country(invalid_data)
  # travel_history, valid/invalid error-grid information.  
  trav_history, trav_history_error = lib.random_travel_history(invalid_data)
  # exposure_event, valid/invalid error-grid information.
  exp_event, exp_event_error = lib.random_exposure_event(invalid_data)
  
  ### Sequencing ###
  
  # minion_barcode, valid/invalid error-grid information.
  minion_barcode, minion_barcode_error = lib.random_minIon_barcode(invalid_data)
  # sequencing_instrument, valid/invalid error-grid information.
  seq_instrument, seq_instrument_error = lib.random_seq_instrument(invalid_data)
  # sequencing_protocol_name, valid/invalid error-grid information.
  seq_prot_name, seq_prot_name_error = lib.fake_protocol(invalid_data)
  # sequencing_protocol_source, valid/invalid error-grid information.
  seq_prot_source, seq_prot_source_error = lib.random_seq_protocol_source(invalid_data)
  # sequencing_kit_number, valid/invalid error-grid information.
  seq_kit_num, seq_kit_num_error = lib.random_seq_kit_num(invalid_data)
  # amplicon_pcr_primers_filename, valid/invalid error-grid information.
  amp_pcr_filename, amp_pcr_filename_error = lib.random_txt_filename(invalid_data)
  
  ### Bioinformatics and QC metrics ###
  
  # raw_sequence_data_processing, valid/invalid error-grid information.
  raw_seq_process, raw_seq_process_error = lib.random_seq_process(invalid_data)
  # sequencing_depth_average, valid/invalid error-grid information.
  seq_depth_avg, seq_depth_avg_error = lib.random_seq_depth(invalid_data)
  # assembly_method, valid/invalid error-grid information.
  assemb_method, assemb_method_error = lib.random_assembly_software(invalid_data)
  # assembly_coverage_breadth, valid/invalid error-grid information.
  assemb_cov_breadth, assemb_cov_breadth_error = lib.random_assembly_coverage_breadth(invalid_data)
  # assembly_coverage_depth, valid/invalid error-grid information.
  assemb_cov_depth, assemb_cov_depth_error = lib.random_seq_depth(invalid_data)
  # r1_fastq_filename, valid/invalid error-grid information.
  r1_filename, r1_filename_error = lib.random_fastq_filename(invalid_data)
  # r2_fastq_filename, valid/invalid error-grid information.
  r2_filename, r2_filename_error = lib.random_fastq_filename(invalid_data)
  # r1_fastq_filepath, valid/invalid error-grid information.
  r1_filepath, r1_filepath_error = lib.random_filepath(r1_filename, invalid_data)
  # r2_fastq_filepath, valid/invalid error-grid information.
  r2_filepath, r2_filepath_error = lib.random_filepath(r2_filename, invalid_data)
  # fast5_filename, valid/invalid error-grid information.
  fast5_filename, fast5_filename_error = lib.random_fast5_filename(invalid_data)
  # fast5_filepath, valid/invalid error-grid information.
  fast5_filepath, fast5_filepath_error = lib.random_filepath(fast5_filename, invalid_data)
  # fasta_filename, valid/invalid error-grid information.
  fasta_filename, fasta_filename_error = lib.random_fasta_filename(invalid_data)
  # fasta_filepath, valid/invalid error-grid information.
  fasta_filepath, fasta_filepath_error = lib.random_filepath(fasta_filename, invalid_data)
  # number_base_pairs, valid/invalid error-grid information.
  num_bp, num_bp_error = lib.random_bp_num(invalid_data)
  # consensus_genome_length, valid/invalid error-grid information.
  cons_genome_len, cons_genome_len_error = lib.random_genome_length(invalid_data)
  # mean_contig_length, valid/invalid error-grid information.
  mean_contig_len, mean_contig_len_error = lib.random_contig_length(invalid_data)
  # n50, valid/invalid error-grid information.
  n50, n50_error = lib.random_n50(invalid_data)
  # ns_per_100_kbp, valid/invalid error-grid information.
  ns_100kbp, ns_100kbp_error = lib.random_ns_100kbp(invalid_data)
  # reference_genome_accession, valid/invalid error-grid information.
  ref_genome_accession, ref_genome_accession_error = lib.random_ref_genome(invalid_data)
  # consensus_sequence_id, valid/invalid error-grid information.
  cons_seq_id, cons_seq_id_error = lib.random_consensus_seq_id(invalid_data)
  # consensus_sequence_method, valid/invalid error-grid information.
  cons_seq_method, cons_seq_method_error = lib.random_consensus_seq_method(invalid_data)
  # consensus_sequence_filename, valid/invalid error-grid information.
  cons_seq_filename, cons_seq_filename_error = lib.random_fasta_filename(invalid_data)
  # consensus_sequence_filepath, valid/invalid error-grid information.
  cons_seq_filepath, cons_seq_filepath_error = lib.random_filepath(cons_seq_filename, invalid_data)
  # annotation_feature_table_filename, valid/invalid error-grid information.
  annot_table_filename, annot_table_filename_error = lib.random_feature_table_filename(invalid_data)
  # bioinformatics_protocol, valid/invalid error-grid information.
  biof_protocol, biof_protocol_error = lib.bioinformatics_protocol(invalid_data)
  
  ### Pathogen Diagnostic Testing ###
  
  # gene_name_1, valid/invalid error-grid information.
  gene_1, gene_1_error = lib.random_gene(invalid_data)
  # diagnostic_pcr_protocol_1, valid/invalid error-grid information.
  pcr_protocol_1, pcr_protocol_1_error = lib.fake_protocol(invalid_data)
  # diagnostic_pcr_ct_value_1, valid/invalid error-grid information.
  pcr_ct_1, pcr_ct_1_error = lib.random_pcr_ct_val(invalid_data)
  # gene_name_2, valid/invalid error-grid information.
  gene_2, gene_2_error = lib.random_gene(invalid_data)
  # diagnostic_pcr_protocol_2, valid/invalid error-grid information.
  pcr_protocol_2, pcr_protocol_2_error = lib.fake_protocol(invalid_data)
  # diagnostic_pcr_ct_value_2, valid/invalid error-grid information.
  pcr_ct_2, pcr_ct_2_error = lib.random_pcr_ct_val(invalid_data)
  
  ### Contributor Acknowledgement ###
  
  # authors, valid/invalid error-grid information.
  authors, authors_error = lib.authors(invalid_data)
  
  ### Dependent IDs ###
  
  # specimen_collector_sample_id.
  spec_col_sample_id, spec_col_sample_id_error = lib.random_specimen_collector_sample_id(
                                                            geo_loc_country, 
                                                            geo_loc_prov_ter,
                                                            geo_loc_city,
                                                            invalid_data)
  # irida_sample_name.
  irida_sample_id = spec_col_sample_id
  # sequence_submitted_by.
  seq_sub_by = samp_col_by
  # library_id, valid/invalid error-grid information.
  library_id, library_id_error = lib.random_library_id(spec_col_sample_id, invalid_data)
  # isolate.
  isolate = spec_col_sample_id
  # assembly_name, valid/invalid error-grid information.
  assemb_name, assemb_name_error = lib.random_assembly_name(spec_col_sample_id, invalid_data)

  # Row of generated data organised by column.
  cols = [spec_col_sample_id,
          phac_sample_id,
          irida_sample_id,
          ub_accession,
          bp_accession,
          bs_accession,
          sra_accession,
          gb_accession,
          gisaid_accession,
          samp_col_by,
          samp_col_email,
          samp_col_address,
          seq_sub_by,
          seq_sub_email,
          seq_sub_address,
          samp_col_date,
          samp_rec_date,
          geo_loc_country,
          geo_loc_prov_ter,
          geo_loc_city,
          organism,
          isolate,
          p_o_sampling,
          anat_material,
          anat_part,
          body_product,
          envi_material,
          envi_site,
          col_device,
          col_method,
          col_protocol,
          spec_process,
          lab_host,
          passage_num,
          passage_method,
          biom_extract,
          host_com_name,
          host_sci_name,
          host_health_state,
          host_health_status,
          host_disease,
          host_age,
          host_gender,
          host_loc_country,
          host_sub_id,
          symp_onset_date,
          signs_symptoms,
          loc_exp_country,
          trav_history,
          exp_event,
          library_id,
          minion_barcode,
          seq_instrument,
          seq_prot_name,
          seq_prot_source,
          seq_kit_num,
          amp_pcr_filename,
          raw_seq_process,
          seq_depth_avg,
          assemb_name,
          assemb_method,
          assemb_cov_breadth,
          assemb_cov_depth,
          r1_filename,
          r2_filename,
          r1_filepath,
          r2_filepath,
          fast5_filename,
          fast5_filepath,
          fasta_filename,
          fasta_filepath,
          num_bp,
          cons_genome_len,
          mean_contig_len,
          n50,
          ns_100kbp,
          ref_genome_accession,
          cons_seq_id,
          cons_seq_method,
          cons_seq_filename,
          cons_seq_filepath,
          annot_table_filename,
          biof_protocol,
          gene_1,
          pcr_protocol_1,
          pcr_ct_1,
          gene_2,
          pcr_protocol_2,
          pcr_ct_2,
          authors]

  # Error grid valid/invalid (error specific) information.
  grid = ['-', # spec_col_sample_id_error placeholder
          phac_sample_id_error,
          '-', # irida_sample_id_error placeholder
          ub_accession_error,
          bp_accession_error,
          bs_accession_error,
          sra_accession_error,
          gb_accession_error,
          gisaid_accession_error,
          samp_col_by_error,
          samp_col_email_error,
          samp_col_address_error,
          '-', # seq_sub_by_error placeholder
          seq_sub_email_error,
          seq_sub_address_error,
          samp_col_date_error,
          samp_rec_date_error,
          geo_loc_country_error,
          geo_loc_prov_ter_error,
          geo_loc_city_error,
          organism_error,
          '-', # isolate_error placeholder
          p_o_sampling_error,
          anat_material_error,
          anat_part_error,
          body_product_error,
          envi_material_error,
          envi_site_error,
          col_device_error,
          col_method_error,
          col_protocol_error,
          spec_process_error,
          lab_host_error,
          passage_num_error,
          passage_method_error,
          biom_extract_error,
          host_com_name_error,
          host_sci_name_error,
          host_health_state_error,
          host_health_status_error,
          host_disease_error,
          host_age_error,
          host_gender_error,
          host_loc_country_error,
          host_sub_id_error,
          symp_onset_date_error,
          signs_symptoms_error,
          loc_exp_country_error,
          trav_history_error,
          exp_event_error,
          library_id_error,
          minion_barcode_error,
          seq_instrument_error,
          seq_prot_name_error,
          seq_prot_source_error,
          seq_kit_num_error,
          amp_pcr_filename_error,
          raw_seq_process_error,
          seq_depth_avg_error,
          assemb_name_error,
          assemb_method_error,
          assemb_cov_breadth_error,
          assemb_cov_depth_error,
          r1_filename_error,
          r2_filename_error,
          r1_filepath_error,
          r2_filepath_error,
          fast5_filename_error,
          fast5_filepath_error,
          fasta_filename_error,
          fasta_filepath_error,
          num_bp_error,
          cons_genome_len_error,
          mean_contig_len_error,
          n50_error,
          ns_100kbp_error,
          ref_genome_accession_error,
          cons_seq_id_error,
          cons_seq_method_error,
          cons_seq_filename_error,
          cons_seq_filepath_error,
          annot_table_filename_error,
          biof_protocol_error,
          gene_1_error,
          pcr_protocol_1_error,
          pcr_ct_1_error,
          gene_2_error,
          pcr_protocol_2_error,
          pcr_ct_2_error,
          authors_error]

  # return row of data for data file and row of data validity for error grid.
  return cols, grid


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
              

def generate_data_file(file_name, rows, delimiter, invalid_data=[]):
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
generate_data_file(argv[1], int(argv[2]), argv[3])
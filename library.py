
"""
Generate COVID19 CanCOGen Metadata Function Library

Version: 1.7
COVID-19 Vocab Version #39
Date: 2020/06/11

Python v 3.8.3

@author: Rhiannon Cameron
         orcid.org/0000-0002-9578-0788
"""

import pandas as pd
import random
from faker import Faker
import datetime
import string
import uuid
from num2words import num2words


# Error grid output indicating valid data
global_valid_data = '-'


def generate_covid19_dict(covid19_vocab):
  """
  Generate a dictionary of the CanCOGeN COVID19 vocabulary

  given: csv file where the first row contains column headers, and each column
         lists items from the CanCOGeN pick-lists or theoretical answers for
         a CanCOGeN metadata field.
  return: data dictionary where the keys are the column headers and the values
          are series of the data from said column. 
  """
  data = pd.read_csv(covid19_vocab, header=0, encoding='utf-8')
  covid19_vocab_dict = {}

  for col in data.columns:
    # clean up header names, 
    sanitized = (col.lower().replace(' ','_').replace('(','').replace(')','')
                  .replace('/','_').replace('ï»¿',''))
    # make headers the dict keys, columns series dict values, and drop all 
    # empty/NaN cells.
    covid19_vocab_dict[sanitized] = data.loc[:,col].dropna()

  return covid19_vocab_dict

# load and generate metadata dictionary
covid19_vocab_dict = generate_covid19_dict(open('vocab-lists.csv','r'))


def generate_covid19_vars():
  """
  Generate a dictionary of CanCOGeN COVID19 variable names
  
  given: csv file containing all current CanCOGeN metadata labels, comma
         seperated.
  return: dictionary where keys are CanCOGeN metadata labels and values are the
          labels formatted as variable names.
  """
  data = pd.read_csv('variables.csv', header=0)
  covid19_vars = {}
  
  for col in data.columns:
    # clean up names
    sanitized = col.lower().replace(' ','_').replace('(','').replace(')','')
    covid19_vars[sanitized] = str(col)
    
  return covid19_vars


#------------------------------------------------------------------------------
#                                    MODULES
#------------------------------------------------------------------------------


# PHAC sample ID
def random_phac_id(invalid_data):
  """
  Generate Random PHAC ID
  
  return: string formatted to imitate a PHAC sample ID.
  """
  return ("PHAC_" + str(random.randint(99, 2001))), global_valid_data


# umbrella bioproject accession
def umbrella_bioproject_accession(invalid_data):
  """
  Umbrella BioProject Accession
  
  return: string containing CanCOGeN umbrella bioproject accession number. 
          Different provinces will have their own BioProjects, however these 
          BioProjects will be linked under one umbrella BioProject.
  """
  return ("PRJNA623807"), global_valid_data


# bioproject accession
def random_bioproject_accession(invalid_data):
  """
  Generate Random BioProject Accession
  
  return: string formatted to imitate a bioproject accession number. A valid 
          NCBI BioProject accession has prefix PRJNA, followed by six numerals.
  """
  return ("PRJNA" + str(random.randint(400000,1000000))), global_valid_data


# biosample accession
def random_biosample_accession(invalid_data):
  """
  Generate Random BioSample Accession
  
  return: string formatted to imitate a biosample accession number. NCBI 
          BioSamples will have the prefix SAMN.
  """
  return ("SAMN" + str(random.randint(10000000,20000001))), global_valid_data


# SRA accession
def random_sra_accession(invalid_data):
  """
  Generate Random SRA Accession
  
  return: string formatted to imitate a SRA accession. NCBI-SRA accessions 
          start with SRR.
  """
  return ("SRR" + str(random.randint(10000000,90000001))), global_valid_data


# GenBank accession
def random_genbank_accession(invalid_data):
  """
  Generate Random GenBank Accession
  
  return: string formatted to imitate a GenBank submission accession for viral
          genome assemlby. 
  """
  # prefixs for GenBank direct submissions, excluding single letter 'U' option.
  prefix = ['AF', 'AY', 'DQ', 'EF', 'EU', 'FJ', 'GQ', 'GU', 'HM', 'HQ', 'JF', 
            'JN', 'JQ', 'JX', 'KC', 'KF', 'KJ', 'KM', 'KP', 'KR', 'KT', 'KM', 
            'KP', 'KR', 'KT', 'KU', 'KX', 'KY', 'MF', 'MG', 'MH', 'MK', 'MN', 
            'MT']
  # suffix is 6 numerals.
  return ((random.choice(prefix) + str(random.randint(100000,1000000))), 
          global_valid_data)


# GISAID accession
def random_gisaid_accession(invalid_data):
  """
  Generate Random GISAID Accession
  
  return: string formatted to imitate the accession returned from the GISAID
          submission. GISAID currently uses the prefix "EPI_ISL_" for 
          SARS-CoV-2 submissions.
  """
  return ("EPI_ISL_" + str(random.randint(400000,600000))), global_valid_data


# reference genome accession
def random_ref_genome(invalid_data):
  """
  Generate Random Genome Accession 
  
  return: string containing a real SARS-CoV-2 reference genome accession or a
          randomly generated one.
  """
  # choose from list of commonly used SARS-CoV-2 reference genomes.
  real = random.choice(covid19_vocab_dict.get('reference_genome_accession'))
  # choose from random accession formats defined in library, 
  # omitting data validity information.
  fake = random.choice([random_genbank_accession(invalid_data)[0], 
                        random_sra_accession(invalid_data)[0], 
                        random_biosample_accession(invalid_data)[0], 
                        random_gisaid_accession(invalid_data)[0]])
  
  return random.choice([real, fake]), global_valid_data


# sample collected by
# sequence submitter by
def random_agency(invalid_data):
  """
  Generate Random Sample Collection Agency
  
  return: string referencing a Canadian Public Health agency.
  """
  agency = random.choice(covid19_vocab_dict.get('agencies'))
  
  return agency, global_valid_data


Faker('en_GB').name()
# sample collector contact email
# sequence submitter contact email
def random_email(invalid_data):
  """
  Generate Random E-mail
  
  return: string containing imitation email address from random domains.
  """
  
  return (random.choice([Faker().email(), Faker().company_email()]), 
          global_valid_data)


# sample collector contact address
# sequence submitter contact address
def random_address(invalid_data):
  """
  Generate Random Address
  
  return: string containing imitation postal address.
  """
  fake = Faker(['en_CA']) # localized to Canada
  
  return fake.address().replace('\n',', '), global_valid_data


# sample collection date
# sample received date
# symptom onset date
def random_date(invalid_data):
#def random_date(invalid_data):
  """
  Generate Random Date
  
  error: boolean; True will generate invalid data, False will generate valid 
         data.
  return: string containing valid random date, in ISO 8601 standard 
          "YYYY-MM-DD", from 2019-12 to the present year; or a string of an 
          invalid date format and/or a non-existant day (e.g. Jun 31).
          
  """
  # start date around time of initial COVID19 cases.
  start_date = datetime.date(year=2019, month=12, day=1)
  # generate date between start date and today.
  date = Faker().date_between(start_date, end_date='today')
  
  # output invalid date format.
  if 'date_error' in invalid_data:
   
    while True:
                 
      joiner = random.choice(['-','/',' ',''])   
      year = random.choice([str(date)[0:4],str(date)[2:4]])
      month_names = ['January','February','March','April','May','June','July',
                     'August','September','October','November','December']
      month = random.choice([str(date)[5:7], random.choice(month_names), 
                             random.choice(month_names)[0:3], 
                             random.choice(month_names)[0]])      
      day_names = ['Monday','Tuesday','Wednesday','Thursday','Friday',
                   'Saturday','Sunday']
      day = random.choice([str(date)[8:10], random.choice(day_names), 
                           random.choice(day_names)[0:3]])      
      options = [year, month, day]      
      # randomize order of items in options.
      random.shuffle(options)
      # join date components together.
      invalid_date = joiner.join(options)
      # check that date isn't in valid YYYY-MM-DDf ormat.
      valid_dates = [str(date)]
      
      if invalid_date not in valid_dates:
        return invalid_date
  
  # output valid ISO 8601 standard date format.
  else:
    # return date in YYYY-MM-DD.
    return random.choice([date]), global_valid_data


# geo_loc_name_country
# host origin geo_loc name (country)
# location of exposure geo_loc name (country)
def random_country(invalid_data):
  """
  Generate Random Country
  
  return: string containing country name from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('geo_loc_name_country')), 
          global_valid_data)


# geo_loc_name_province/territory
def random_province_territory(invalid_data):
  """
  Generate Random Canadian Province/Territory
  
  return: string containing province/territory name from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('geo_loc_name_province_territory')), 
          global_valid_data)


# geo_loc_name_city
def random_city(invalid_data):
  """
  Generate Random Canadian City
  
  return: string containing city name.
  """
  fake = Faker(['en_CA']) # localized to Canada
  
  return fake.city(), global_valid_data


# organism
def random_organism(invalid_data):
  """
  Generate Random Organism
  
  return: string containing "organism" name from CanCOGeN vocabulary.
  """
  return random.choice(covid19_vocab_dict.get('organism')), global_valid_data


# purpose of sampling
def random_purpose_of_sampling(invalid_data):
  """
  Generate Random Purpose of Sampling
  
  return: string containing "purpose of sampling" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('purpose_of_sampling')), 
          global_valid_data)


# anatomical material
def random_anatomical_material(invalid_data):
  """
  Generate Random Anatomical Material
  
  return: string containing "anatomical material" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('anatomical_material')), 
          global_valid_data)


# anatomical part
def random_anatomical_part(invalid_data):
  """
  Generate Random Anatomical Part
  
  return: string containing "anatomical part" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('anatomical_part')), 
          global_valid_data)


# body product
def random_body_product(invalid_data):
  """
  Generate Random Body Product
  
  return: string containing "body part" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('body_product')), 
          global_valid_data)


# environmental material
def random_environmental_material(invalid_data):
  """
  Generate Random Environmental Material
  
  return: string containing "environmental material" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('environmental_material')), 
          global_valid_data)


# environmental site
def random_environmental_site(invalid_data):
  """
  Generate Random Environmental Site
  
  return: string containing "environmental site" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('environmental_site')), 
          global_valid_data)


# collection device
def random_collection_device(invalid_data):
  """
  Generate Random Collection Device
  
  return: string containing "collection device" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('collection_device')), 
          global_valid_data)


# collection method 
def random_collection_method(invalid_data):
  """
  Generate Random Collection Device
  
  return: string containing "collection device" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('collection_method')), 
          global_valid_data)


# collection protocol
# sequencing_protocol_name
# diagnostic_pcr_protocol 1
# diagnostic_pcr_protocol 2
def fake_protocol(invalid_data):
  """
  Generate Fake Protocol Name
  
  return: string containing a fake protocol name.
  """
  prefix = Faker().text(20).replace('.','').replace(' ','')
  prefix = random.choice([prefix, prefix.upper(), prefix.lower()])
  affix = (random.choice(['_','-','/','|',]) + 
           random.choice(['v','V','v.','V.']) + 
           random.choice(['',' ','_']))
  suffix = str(Faker().random_int(1,11))

  return (prefix + affix + suffix), global_valid_data

# specimen processing
def random_specimen_processing(invalid_data):
  """
  Generate Random Specimen Processing
  
  return: string containing "specimen processing" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('specimen_processing')), 
          global_valid_data)


# lab_host
def random_lab_host(invalid_data):
  """
  Generate Random Lab Host Cell Line
  
  return: string containing "lab host" cell line from CanCOGeN vocabulary.
  """
  return random.choice(covid19_vocab_dict.get('lab_host')), global_valid_data


# passage number
def random_passage_number(invalid_data):
  """
  Generate Random Viral Passage Number
  
  return: int or string if int value is zero.
  """
  passage = Faker().random_int(0,25) # maximum can be increased.
  if passage == 0:
    passage = 'not applicable'
    
  return passage, global_valid_data


# passage method
def passage_method_text(invalid_data):
  """
  Generate Viral Passage Method Name
  
  return: string containing sentence of random words.
  """
  words = random.randint(2,10) # less than 10 words.
  text = (' '.join(Faker().words(words)) + '.').capitalize()
  
  return random.choice([text, 'not applicable']), global_valid_data


# biomaterial extracted
def random_biomaterial_extracted(invalid_data):
  """
  Generate Random Extracted Biomaterial
  
  return: string containing "biomaterial extracted" from CanCOGeN vocabulary.
          vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('biomaterial_extracted')),
          global_valid_data)


# host common name
def random_host_common_name(invalid_data):
  """
  Generate Random Host - Common Name
  
  return: string containing "host (common name)" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('host_common_name')),
          global_valid_data)


# host scientific name
def random_host_scientific_name(invalid_data):
  """
  Generate Random Host - Scientific Name
  
  return: string containing "host (scientific name)" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('host_scientific_name')),
          global_valid_data)


# host health state
def random_host_health_state(invalid_data):
  """
  Generate Random Host Health State
  
  return: string containing "host health state" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('host_health_state')),
          global_valid_data)


# host health status details
def random_host_health_status_details(invalid_data):
  """
  Generate Random Host Health Status Details
  
  return: string containing "host health status details" from CanCOGeN 
          vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('host_health_status_details')), 
          global_valid_data)


# host disease
def random_host_disease(invalid_data):
  """
  Generate Random Host Disease
  
  return: string containing "host disease" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('host_disease')), 
          global_valid_data)


# host age
def random_host_age(invalid_data):
  """
  Generate Random Host Age
  
  error: boolean; True will generate invalid data, False will generate valid 
         data.
  return: string containing a valid numerical age; or an invalid 'range of ages'
          or age as text.
          
  """
  age = random.randint(0,105)
  
  # output invalid date format.
  if 'age_error' in invalid_data:
    age_range = random.choice(covid19_vocab_dict.get('host_age'))
    age_text = num2words(age)
    return random.choice([age_range, age_text]), 'age_error'
  
  # output valid date format.
  return str(age), global_valid_data


# host gender
def random_host_gender(invalid_data):
  """
  Generate Random Host Gender
  
  return: string containing "host gender" from CanCOGeN vocabulary.
  """
  return random.choice(covid19_vocab_dict.get('host_gender')), global_valid_data


# host subject ID
def random_host_subject_id(invalid_data):
  """
  Generate Random Host Subject ID
  
  return: string containing imitation user-defined "host subject ID".
  """
  prefix = Faker().pystr(min_chars=1, max_chars=6)
  affix = ''
  suffix = str(Faker().random_int(1,11))
  with_affix = random.choice([True,False])

  if with_affix == True:
    affix = random.choice(['_','-','/','|',':',';','+'])

  return (prefix + affix + suffix), global_valid_data


# signs symptoms
def random_signs_symptoms(invalid_data):
  """
  Generate Random Signs/Symptoms
  
  return: string containing "signs and symptoms" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('signs_symptoms')), 
          global_valid_data)


# travel history
def random_travel_history(invalid_data):
  """
  Generate Random Travel History
  
  return: string containing at least one hypothetical country and city, and up 
          to four additional countries (cities optional).
  """
  locations = (str(random_country(invalid_data)) + ', ' 
               + str(random_city(invalid_data)))

  for i in range(random.randint(0,5)):

    with_city = random.choice([True,False])

    if with_city == True:
      # City known; currently only outputting Canadian cities, regardless of country.
      locations = locations + str('; ' + str(random_country(invalid_data)) 
                                  + ', ' + str(random_city(invalid_data)))
    else:
      # City unknown.
      locations = locations + str('; ' + str(random_country(invalid_data)))

  return locations, global_valid_data


# exposure event
def random_exposure_event(invalid_data):
  """
  Generate Random Exposure Event
  
  return: string containing "exposure event" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('exposure_event')), 
          global_valid_data)


# library_ID
def random_library_id(specimen_collector_sample_id, invalid_data):
  """
  Generate Random Sequencing Library ID
  
  return: string containing imitation user-defined "library ID".
  """
  return ((specimen_collector_sample_id
          + random.choice(['_','-','/','|',':',';',' ']) 
          + str(random.randint(0,9))), 
          global_valid_data)


# MinIon barcode
def random_minIon_barcode(invalid_data):
  """
  Generate Random minIon Barcode
  
  return: int representing imitation barcode of a MinIon sequencing unit.
  """
  return random.randint(1000000000,10000000000), global_valid_data


# sequencing instrument
def random_seq_instrument(invalid_data):
  """
  Generate Random Sequencing Instrument
  
  return: string containing "sequencing instrument" from CanCOGeN vocabulary.
  """
  choices = covid19_vocab_dict.get('sequencing_instrument')
  
  return random.choice(choices), global_valid_data


# sequencing_protocol_source
def random_seq_protocol_source(invalid_data):
  """
  Generate Random Sequencing Protocol Source Organization/Authors
  
  return: string containing random sequencing protocol source.
  """
  source = random.choice(covid19_vocab_dict.get('sequencing_protocol_source'))
  decimal = random.choice([random.randint(0,2),None])
  version = random.choice([' v',' V',' v.',' V.',' v ',' V ',' v. ',' V. '])
  version_num = round(random.uniform(1,15), decimal)

  return source + version + str(version_num), global_valid_data


def random_alphanumeric():
  """
  Generate Random Alphanumeric
  
  return: string of random numbers and letters.
  """
  alpha_num = (''.join(random.choice(string.ascii_uppercase + 
                                     string.ascii_lowercase + 
                                     string.digits) for _ in 
                       range(random.randint(6,20))))
    
  return (random.choice([alpha_num,alpha_num.upper(),alpha_num.lower()]),
          global_valid_data)


# seq_kit_number
def random_seq_kit_num(invalid_data):

  return random_alphanumeric()

def random_filename(invalid_data):
  """
  Generate Random Filename
  
  return: string containing imitation file name.
  """
  joiner = random.choice([' ','_','-'])
  return joiner.join(Faker().words(random.randint(1,4))), global_valid_data


# amplicon_pcr_primers_filename
def random_txt_filename(invalid_data):
  """
  Generate Random TXT Filename
  
  return: string containing imitation filename in ".txt" format.
  """
  # call random_filename, ignoring invalid_data information.
  return (random_filename(invalid_data)[0] + '.txt'), global_valid_data


# r1 fastq filename
# r2 fastq filename
def random_fastq_filename(invalid_data):
  """
  Generate Random FASTQ Filename
  
  return: string containing imitation filename in ".fastq.gz" format.
  """
  # call random_filename, ignoring invalid_data information.
  return (random_filename(invalid_data)[0] + '.fastq.gz'), global_valid_data


# fast5 filename
def random_fast5_filename(invalid_data):
  """
  Generate Random FAST5 Filename
  
  return: string containing imitation filename in ".fast5" format.
  """
  # call random_filename, ignoring invalid_data information.
  return (random_filename(invalid_data)[0] + '.fast5'), global_valid_data
  
  
# fasta filename
# consensus sequence filename
def random_fasta_filename(invalid_data):
  """
  Generate Random FASTA Filename
  
  return: string containing imitation filename in ".fasta" format.
  """
  # call random_filename, ignoring invalid_data information.
  return (random_filename(invalid_data)[0] + '.fasta'), global_valid_data


# annotation feature table filename
def random_feature_table_filename(invalid_data):
  """
  Generate Random Feature Table Filename
  
  return: string containing imitation filename in ".tbl" format.
  """
  # call random_filename, ignoring invalid_data information.
  return (random_filename(invalid_data)[0] + '.tbl'), global_valid_data


# raw sequence data processing
def random_seq_process(invalid_data):
  """
  Generate Random Raw Sequence Data Process Name.
  
  return: string containing imitation processing software name.
  """
  real_options = random.choice(covid19_vocab_dict.get(
      'raw_sequence_data_processing'))
  fake_options = ''.join(Faker().words(random.randint(1,4))).capitalize()
  prefix = random.choice([real_options, fake_options])

  decimal = random.choice([random.randint(0,2),None])
  version = random.choice([' v',' V',' v.',' V.',' v ',' V ',' v. ',' V. '])
  version_num = round(random.uniform(1,15), decimal) 
  suffix = version + str(random.randint(0,10)) + '.' + str(version_num)

  return (prefix + suffix), global_valid_data


# sequence depth (average)
# assembly coverage depth
def random_seq_depth(invalid_data):
  """
  Generate Random Average Sequencing Depth
  
  return: string containing random sequencing depth number.
  """
  return (str(random.randint(1,150)) + 'x'), global_valid_data


# assembly name
def random_assembly_name(specimen_collector_sample_id, invalid_data):
  """
  Generate Random Sequence Assembly Name
  
  return: string containing imitation assembly name/version.
  """
  return (specimen_collector_sample_id.lower() + str(random.randint(0,100)) + 
          'assembly.fasta'), global_valid_data


# assembly software/method
def random_assembly_software(invalid_data):
  """
  Generate Random Assembly Software Name
  
  return: string containing "assembly software" from CanCOGeN vocabulary.
  """
  return (random.choice(covid19_vocab_dict.get('assembly_software')), 
          global_valid_data)


# assembly coverage breadth
def random_assembly_coverage_breadth(invalid_data):
  """
  Generate Random Assembly Coverage Breadth
  
  return: string containing random percentage.
  """
  decimal = random.choice([random.randint(1,2),None])
  value = random.choice([random.randint(1,101), round(random.uniform(1,101), decimal)])
  
  return (str(value) + '%'), global_valid_data


# r1 fastq filepath
# r2 fastq filepath
# fast5 filepath
# fasta filepath
# consensus sequence filepath
def random_filepath(filename, invalid_data):
  """
  Generate Random Filepath
  
  return: string containing imitation filepath.
  """
  fake = Faker()
  os = random.choice(['Windows','MacOS','Linux'])
  os = 'Linux'

  username = random.choice([fake.name(),fake.first_name(),fake.last_name()])
  folder_1 = random.choice(['Documents','Downloads','Google Drive','Desktop'])
  folder_2 = random.choice([fake.word(),fake.word().capitalize(),''])

  # User may not include the filename in the filepath.
  filename = random.choice([filename,''])

  # User is using Windows OS.
  if os == 'Windows':
    # Option of commonly used 'C:', an alternative volume, a volume name or path
    # that exceeds Windows 255 char limit, or a random drive mounted on C:.
    prefix = random.choice(['C:', random.choice(string.ascii_letters).upper() + 
                            ':', "\\\\?\\Volume{" + str(uuid.uuid4()) + "}", 
                            'C:\Mount' + 
                            random.choice(string.ascii_letters).upper()])
    
    # Windows puts double quotes around folders with spaces on the command line.
    # Make it random whether the use includes them.
    if any(" " in s for s in username):
      user_quotes = random.choice([True,False])
      if user_quotes == True:
        username = '"' + username + '"'

    filepath_list = [prefix, 'Users', username, folder_1, folder_2, filename]

    # Windows filepaths defaults to '\' but do work with '/'.
    return '\\'.join(filepath_list).replace('\\\\','\\'), global_valid_data

  # User is using Mac OS.
  elif os == 'MacOS':
    filepath_list = ['/Users', username, folder_1, folder_2, filename]
    return '/'.join(filepath_list).replace('//','/'), global_valid_data

  # User is using Linux OS.
  else:
    # some linux filesystems allow the use of unusual special characters, have
    # included a sampling to test and added to the end since some are not allowed    !!!!!!!!!!!!!!!!!!! /:
    # at the beginning. 
    folder_2 = folder_2 + random.choice(['\\',';','&','"','.','#','$','-','*'])
    filepath_list = ['/Users', username, folder_1, folder_2, filename]
    return '/'.join(filepath_list).replace('//','/'), global_valid_data


# number base pairs
def random_bp_num(invalid_data):
  """
  Generate Random Base Pair Number
  
  return: int.
  """
  return random.randint(30000,400000), global_valid_data


# consensus genome length
def random_genome_length(invalid_data):
  """
  Generate Random Genome Length
  
  return: int.
  """
  return random.randint(2500,35000), global_valid_data


# mean contig length
def random_contig_length(invalid_data):
  """
  Generate Random Contig Length
  
  return: int.
  """
  return random.randint(10,50000), global_valid_data


# N50
def random_n50(invalid_data):
  """
  Generate Random N50
  
  return: int.
  """
  return round(random.randint(10,5000), 2), global_valid_data

# Ns per 100 kbp
def random_ns_100kbp(invalid_data):
  """
  Generate Random Ns per 100 kbp Value
  
  return: float, rounded to two decimal places.
  """
  return round(random.uniform(0,1000), 2), global_valid_data


# consensus sequence ID
def random_consensus_seq_id(invalid_data):
  """
  Generate Random Consensus Sequence ID
  
  return: string containing imitation sequence ID.
  """
  # Could not find template format, those generated are completely fake.
  prefix = random.choice(['SCV2','Prov','ConsensusSeq','SARS2','SARS-CoV-2'])
  alpha_num = (''.join(random.choice(string.ascii_uppercase + 
                                     string.ascii_lowercase + 
                                     string.digits) for _ in range(random.randint(3,3))))
  return ((prefix + '_' + alpha_num.upper() + str(random.randint(100,1000))), 
          global_valid_data)


# consensus sequence method
def random_consensus_seq_method(invalid_data):
  """
  Generate Consensus Sequence Method Name
  
  return: string containing existing or imitation names and version numbers of 
          sequence method software.
  """
  source = random.choice(covid19_vocab_dict.get('consensus_sequence_method'))
  decimal = random.choice([random.randint(0,2),None])
  version = random.choice([' v',' V',' v.',' V.',' v ',' V ',' v. ',' V. '])
  version_num = round(random.uniform(1,15), decimal)

  return ((random.choice([source, Faker().word()]) + version + str(version_num)), 
          global_valid_data)


# bioinformatics protocol
def bioinformatics_protocol(invalid_data):
  """
  Generate Random Bioinformatics Protcol Name
  
  return: string containing a fake bioinformatics protocol name.
  """
  word_list = ['SARS-CoV-2', 'Enrichment', 'Sequencing', 'by', 'Spiked', 'RNA',
               'Primer', 'MSSPR', 'method', 'testing', 'automated', 'kit-free',
               'Extraction', 'for', 'Genome', 'Long', 'Pooled', 'Amplicons',
               'a', 'Platoforms', 'using', 'reads', 'nanopore', 'protocol',
               'sampling', 'viral', 'concentration', 'Detection', 'assay', 
               'PCR']
  return (random.choice([fake_protocol(invalid_data)[0], # omit invalid_data error info.
                         Faker().sentence(ext_word_list=word_list).capitalize()]),
          global_valid_data)


# gene name 1
# gene name 2
def random_gene(invalid_data):
  """
  Generate Random Gene Name
  
  return: string containing a random SARS-CoV-2 gene name.
  """
  orf = random.choice(['orf','ORF'])
  x = random.choice(['1a','1b','1ab','3','3a','6','7a','7b','8','10'])
  orfs = orf + random.choice([x.upper(), x.lower()])

  gene_dict = {'E':['envelope protein', 'CoV envelope protein', 
                    'envelope small membrane protein', 'sM protein'], 
              'M':['membrane protein', 'M protein', 'E1 glycoprotein',
                    'Matrix glycoprotein','membrane glycoprotein', 'orf5', 
                    'ORF5'], 
              'N':['nulceocapsid phosphoprotein', 'Nucleoprotein', 'Protein N',
                    'NC', 'orf9a', 'ORF9a', 'ORF9A'],
              'RdRP':['RNA dependent RNA Polymerase', 
                      'RNA-dependent RNA polymerase'],
              'S':['Spike glycoprotein', 'S glycoprotein', 'E2', 
                    'Peplomer protein', 'orf2', 'ORF2']}

  gene_letter = random.choice(list(gene_dict.keys()))
  gene_names = [item for sublist in [*gene_dict.values()] for item in sublist]

  prefix = random.choice(['SARS-CoV-2 ','2019-nCoV ', 'COVID-19-', ''])
  affix_1 = random.choice([gene_letter,'(' + gene_letter + ')'])
  affix_2 = random.choice(gene_names)
  suffix = random.choice([' gene', '-gene', ''])

  return ((random.choice([orfs, prefix + orfs, prefix + 
                          random.choice([affix_1, affix_2])]) + suffix),
          global_valid_data)


# diagnostic_pcr_Ct_value 1
# diagnostic_pcr_Ct_value 2
def random_pcr_ct_val(invalid_data):
  """
  Generate Random PCR Ct Value
  
  return: string.
  """
  return str(random.randint(1,40)), global_valid_data


def random_name():
  """
  Generate Random Name
  
  return: string containing random first and last name, occasionally with title.
  """
  # localized Faker provider packages without accents.
  local_ascii_no_accents = ['ar_EG','bs_BA','en_AU','en_CA','en_GB',
                            'en_IN','en_NZ','en_US']
  
  return Faker(local_ascii_no_accents).name()


# authors
def authors(invalid_data):
  """
  Generate Random Authors
  
  return: string listing names of 1-10 fake authors. 
  """
  list_authors = []
  num_authors = random.randint(1,11)
  for i in range(num_authors):
    list_authors.append(random_name())
    
  return (', '.join(list_authors)), global_valid_data


# specimen collector sample ID
# IRIDA sample name
# isolate
def random_specimen_collector_sample_id(country, province_territory, city, 
                                        invalid_data):
  """
  Generate Sample/Isolate ID
  
  return: string containing and immitation user-define sample ID.
  """
  opt1 = str(random.choice([country, province_territory, city, 
                          random_environmental_site(invalid_data)[0]])) # omit invalid_data error info.
  opt2 = str(random.choice([random_environmental_material(invalid_data)[0], # omit invalid_data error info.
                         random_organism(invalid_data)[0], # omit invalid_data error info.
                         random_collection_device(invalid_data)[0]])) # omit invalid_data error info.
  opt3 = str(random.choice([random.randint(1,1000), 
                            random_alphanumeric()[0:4]])) # omit invalid_data error info.

  joiner = random.choice(['_','-',''])

  return joiner.join([opt1, opt2, str(opt3)]).replace(' ',''), global_valid_data

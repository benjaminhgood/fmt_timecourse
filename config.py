###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 

data_directory = os.path.expanduser("~/fmt_data/")
analysis_directory = os.path.expanduser("~/fmt_analysis/")
scripts_directory = os.path.expanduser("~/fmt_scripts/")
patric_directory = os.path.expanduser("~/patric_db/")
midas_directory = os.path.expanduser("~/midas_db/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 10
good_species_min_prevalence = 5
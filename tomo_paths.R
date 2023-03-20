# ----------------------------------------------------------------------------------------
# PATHS
# ----------------------------------------------------------------------------------------

# Number of cores used in MUMPS depends on the machine
# Hence given here, not in the "paramaters" file
ncores <- 24
results_directory <- "results_remote_9"
data_directory    <- "data"

SAVE_RESULTS <- TRUE

# Additional_label is added to all filenames
label <- ""

LEO_directory       <- NULL#paste(data_directory, "/beacon_data/", sep = "")
GNSS_directory      <- paste(data_directory, "/GNSS_data/",  sep = "")
IONOSONDE_directory <- paste(data_directory, "/ionosonde_data/",  sep = "")

RO_directory        <- paste(data_directory, "/RO_data/",  sep = "")
RO_pattern          <- "podTec"

ISR_directory       <- NULL #paste(data_directory, "/EISCAT_data/",  sep = "")
ISR_paths           <- NULL #c(
#  paste(ISR_directory, "EISCAT_2018_plasmacalibrated/2018-11-09_bella_120_tomoscand@uhf_plcalib",  sep = ""),
#  paste(ISR_directory, "EISCAT_2018_reanalysed/2018-11-09-tomoscand/2018-11-09_beata_120_tomoscand@vhf",  sep = ""),
#  paste(ISR_directory, "EISCAT_2018_reanalysed/2018-11-09-tomoscand/2018-11-09_tau7_240_tomoscand@32m",  sep = ""),
#  paste(ISR_directory, "EISCAT_2018_reanalysed/2018-11-09-tomoscand/2018-11-09_tau7_240_tomoscand@42m",  sep = ""))

#NeQuick_directory   <- "/Users/norberg/Dropbox/Projects/NeQuick/R-REC-P/NeQuick2_P531-12/"
NeQuick_directory   <- "/home/ubuntu/R-REC-P/NeQuick2_P531-12/"

tomoscand_bg_path <- "results_cl10_100_090/tomo/"
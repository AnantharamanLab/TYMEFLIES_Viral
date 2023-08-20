# Read rds file
temp <- readRDS("10-water_temp_profiles_with_thermocline_stats.rds")

# Export table
write.table(temp, "temp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Read rds file
secchi <- readRDS("5-secchi_all_sources_combined-only_1995-2020_has_manual_annotation.rds")

# Export table
write.table(secchi, "secchi.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Read rds file
DO <- readRDS("14-DO_profiles_with_more_oxycline_stats-colnames_updated.rds")

# Export table
write.table(DO, "DO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

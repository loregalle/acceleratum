bench <- read.csv("data-raw/26147-cal-2025-acc.csv")
bench$timestamp <- as.POSIXct(bench$timestamp, "UTC")
# sanitise dimnames
colnames(bench) <- gsub("\\.", "_", colnames(bench))

usethis::use_data(bench, overwrite = TRUE)

bench <- read.csv("data-raw/bench.csv",
                  colClasses = c("POSIXct", "character", "character"))

usethis::use_data(bench, overwrite = TRUE)

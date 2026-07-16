muskox <- read.csv("data-raw/muskox.csv",
                  colClasses = c("POSIXct", "character", "character"))

usethis::use_data(muskox, overwrite = TRUE)

# the original dataset is not included as it is over 10 GB in size and
# is more than 1 year of data. For more information, contact
# Prof. Niels Martin Schmidt.
# original data file name: 26147-34-2019-acc.csv
muskox <- read.csv("./data-raw/26147-34-2019-acc_3daysubset.csv")

muskox$timestamp <- as.POSIXct(muskox$timestamp, tz = "UTC")
filtermask <- muskox$timestamp >= as.POSIXct("2019-09-28 00:00:00", tz = "UTC") &
  muskox$timestamp < as.POSIXct("2019-10-01 00:00:00", tz = "UTC")
muskox <- muskox[filtermask,]
# sanitise dimnames
rownames(muskox) <- NULL
colnames(muskox) <- gsub("\\.", "_", colnames(muskox))

usethis::use_data(muskox, overwrite = TRUE)

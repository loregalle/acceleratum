library("tidyverse")
library("data.table")
load_all()

dat <- fread("dev/scratch/calib/26147-cal-2025-acc.csv")

brst <- new_burst(dat, "accelerations-raw", "timestamp")
validate_burst(brst)

brst <- burst(dat, "accelerations-raw", "xyz", "timestamp")
brst
brst[1:3,]

write_burst(brst, "dev/scratch/calib/test_brst_save.csv", row.names = F)

xx <- accelerometry(matrix(rnorm(90), ncol = 3), start_time = 0, sampling_rate = 4)
xx
aa <- xx[1:3,]

xx[4:6,]

a2b(xx, 1)

brst
a2b(b2a(brst), burst_size = 40)

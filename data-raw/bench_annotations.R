bench_annotations <- read.csv("data-raw/bench_annotations.csv")

# Note: annotations were originally collected manually and interactively from an
# aclrtm_accelerometery-class representation of the bench dataset and the
# range_select() function.
# Here provided as .csv output for repeatability and example purposes.
#
# Original code:
#
# data(bench)
# bench <- burst(bench, "accelerations_raw", "xyz", "timestamp")
# bench <- burst_to_accelerometery(bench)
# bench_annotations <- segment_select(bench)

usethis::use_data(bench_annotations, overwrite = TRUE)

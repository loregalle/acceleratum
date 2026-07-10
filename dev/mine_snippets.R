library("tidyverse")
library("data.table")
load_all()

path <- "dev/scratch_data/calib/26147-cal-2025-acc.csv"
dat <- fread(path) |> burst("accelerations-raw", "xyz", "timestamp")
dat |> colnames()
acc <- dat |> burst_to_accelerometry()
head(acc)


mined_list <- mine_reservoir(path, "accelerations-raw", "timestamp", "xyz",
                             window_sec = 3,
                             vedba_thresh = 0.02,
                             sr = 8,
                             scratch_path = "dev/tmp")
mined_list

selected_list <- select_fps(mined_list, 6)
selected_list

reconstructed <- reconstruct_selected(mined_list, selected_list$selected_idx, TRUE)
reconstructed

reconstructed[[6]]$data
norm <- reconstructed[[6]]$data |> colMeans() |> `^`(x=_, y = 2) |> sum() |> sqrt()
reconstructed[[6]]$data |> colMeans() |> `/`(x=_, y = norm)
selected_list$orientations[6,]

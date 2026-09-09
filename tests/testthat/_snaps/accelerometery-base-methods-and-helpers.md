# print.aclrtm_accelerometry displays header and matrix

    Code
      print(x)
    Output
      <aclrtm_accelerometry>  [3 samples × 3 axes (xyz)]
        sampling_rate : 5 Hz
        start_time : 2024-01-01
           x y z
      [1,] 1 4 7
      [2,] 2 5 8
      [3,] 3 6 9

# print.aclrtm_accelerometry handles missing attributes gracefully

    Code
      print(x)
    Output
      <aclrtm_accelerometry>  [3 samples × 2 axes (xy)]
           x y
      [1,] 1 4
      [2,] 2 5
      [3,] 3 6


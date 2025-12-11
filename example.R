library(spdep)
library(ggplot2)


### Samples are first processed here to generate original and permutated CLQs for each cell to cell pair in a given sample.

### Input file name example: “TAFs1_cell_type_assignment.csv”


files <- (Sys.glob("*cell_type_assignment.csv"))

for (f in files) {
  print(f)
  filename_c <- f

  count_file <- get_counts(filename = filename_c)

  ### CLQ_permutated_matrix_gen function is using iteration number, filename and the count_file generated in the previous step. This function is dependent on multiple functions [get_CLQ(), KNN_neighbors(), find_cell_type_neighbors()]
  ### iternum is the iteration number for permutation analysis, which is set to 500.


  CLQ_permutated <- CLQ_permutated_matrix_gen(
    iternum = 500,
    filename = filename_c,
    df_c = count_file
  )
}


###
### Then, significance of CLQs are calculated based on their percentile.
### The original CLQs and permutated CLQs are retrieved for each sample through the functions [CLQ_matrix_gen(), CLQ_permutated_matrix_gen_caller(),get_counts_caller() ]

files <- (Sys.glob("*cell_type_assignment.csv"))

for (f in files) {
  print(f)

  filename_c <- f

  CLQ_matrix <- CLQ_matrix_gen(filename = filename_c)

  CLQ_permutated <- CLQ_permutated_matrix_gen_caller(filename = filename_c)

  count_file <- get_counts_caller(filename = filename_c)

  ###  “list_of_matrices” is the 500 different CLQ sets for each iteration.

  significance_matrices <- significance_matrix_gen(
    iternum = 500,
    filename = filename_c,
    list_of_matrices = CLQ_permutated,
    CLQ_matrix_original = CLQ_matrix,
    df_c = count_file
  )

  ### plot_gen generates plot for each CLQ for a pair of cell type A to cell type B. It retrieves the permutated CLQ values from CLQ_permutated, and the original CLQ values from CLQ_matrix.

  plot_gen(
    iternum = 500,
    filename = filename_c,
    list_of_matrices = CLQ_permutated,
    CLQ_matrix_original = CLQ_matrix
  )
}

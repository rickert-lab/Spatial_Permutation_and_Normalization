## These packages need to be installed in R (OS-independent):
# install.packages(c('dplyr', 'ggplot2', 'spdep', 'tidyr'))
## package dependencies might require installation of system libraries (OS-dependent):
# brew install pkg-config udunits gdal

# message from 'spdep': "To access larger datasets in this package, install the spDataLarge package with:"
# install.packages('spDataLarge', repos='https://nowosad.github.io/drat/', type='source')

# load custom library with imports
source("./Spatial_Perm_and_Normalization.R")


### Samples are first processed here to generate original and permutated CLQs for each cell to cell pair in a given sample.

### Input file name example: “TAFs1_cell_type_assignment.csv”

if (FALSE) {
  files <- (Sys.glob("./input/*cell_type_assignment.csv"))
  if (length(files) == 0) {
    stop("No files found at path. Check file path and pattern!")
  }

  for (file in files) {
    print(file)

    count_file <- write_counts(sample_path = file, out_dir = "./output")


    ### CLQ_permutated_matrix_gen function is using iteration number, filename and the count_file generated in the previous step.
    ### This function is dependent on multiple functions [get_CLQ(), KNN_neighbors(), find_cell_type_neighbors()]
    ### iternum is the iteration number for permutation analysis, which is set to 500.


    CLQ_permutated <- CLQ_permutated_matrix_gen1(
      iternum = 500,
      sample_path = file,
      counts_path = count_file,
      out_dir = "./output"
    )
  }
}

if (TRUE) {
  ###
  ### Then, significance of CLQs are calculated based on their percentile.
  ### The original CLQs and permutated CLQs are retrieved for each sample through the functions [CLQ_matrix_gen(),
  ### CLQ_permutated_matrix_gen_caller(),get_counts_caller() ]

  files <- (Sys.glob("./input/*cell_type_assignment.csv"))
  if (length(files) == 0) {
    stop("No files found at path. Check file path and pattern!")
  }

  for (file in files) {
    print(file)

    CLQ_matrix <- CLQ_matrix_gen(sample_path = file)

    CLQ_permutated <- CLQ_permutated_matrix_gen2(sample_path = file)

    stop("STOP!")

    count_file <- get_counts_caller(filename = filename_c)

    ###  “list_of_matrices” is the 500 different CLQ sets for each iteration.

    significance_matrices <- significance_matrix_gen(
      iternum = 500,
      filename = filename_c,
      list_of_matrices = CLQ_permutated,
      CLQ_matrix_original = CLQ_matrix,
      df_c = count_file
    )

    ### plot_gen generates plot for each CLQ for a pair of cell type A to cell type B. It retrieves the permutated CLQ values
    ### from CLQ_permutated, and the original CLQ values from CLQ_matrix.

    plot_gen(
      iternum = 500,
      filename = filename_c,
      list_of_matrices = CLQ_permutated,
      CLQ_matrix_original = CLQ_matrix
    )
  }
}

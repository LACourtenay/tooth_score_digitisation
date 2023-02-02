

directory = list.files("C:/", # insert directory where all the profile images are saved
                       pattern = ".jpg",
                       full.names = TRUE)

options(locatorBell = FALSE) # optional: remove windows sound when each landmark is placed

# This function is for the digitisation of multiple images in a directory, this function can be used to place both fixed landmarks
# and semilandmarks along a curve.

digitise_images_2d(directory,
                   n_landmarks = 7, # number of fixed landmarks
                   n_curves = 0, # number of curves
                   n_semilandmarks = 0, # number of semilandmarks on the curve
                   scale = 2.5, # indicate the scale of the scale bar
                   verbose = FALSE, # if verbose is TRUE, then the software will ask if you want to save the last curve or point placed
                   output = "txt" # save the results in a .txt format
) 

# This function is for the digitisation of score profiles using a single semilandmark curve, while experimenting with the number of
# semilandmarks as well. The user defines the directory where images are obtained and then selects the different models of semilandmarks
# they wish to use. The user then has to trace out the entire profile per image and the algorithm will automatically save the different
# configurations with varying numbers of landmarks

extract_bsm_profile_smlm(directory,
                         n_semilandmarks = c(10, 25, 50), # insert the number of semilandmarks you want to digitise on each curve
                         scale = 0.5, # indicate the scale of the scale bar
                         verbose = FALSE # if verbose is TRUE then with every step you will be asked if you want to save your progress
)

# This function is for the digitisation of score profiles using a two semilandmark curves, the user defines how many semilandmarks they
# wish to place on each side of the score, then they trace using two different curves each side. The algorithm will then automatically
# calculate the position of the semilandmarks and save the results

new_lmsmlm_method(directory,
                  n_semilandmarks = c(7, 12, 17), # insert the number of semilandmarks you want to digitise
                  scale = 0.5, # indicate the scale of the scale bar
                  external_window = TRUE # visualise in an external window
)

# This function is to save all of the results into a morphologika file to be imported into R

generate_morphologika_file(
  "C:/", # define the directory where information is contained
  7 # state the number of landmarks in the model
)


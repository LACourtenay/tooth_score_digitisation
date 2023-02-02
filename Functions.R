
`%!in%` = Negate(`%in%`)

generate_morphologika_file <- function(folder_name, n_landmarks) {
  
  directory_contents <- list.dirs(folder_name, recursive = FALSE)
  
  landmarks <- array(numeric(),
                     dim = c(n_landmarks, 2, 0))
  labels <- c()
  
  for (folder in directory_contents)  {
    
    folder_contents <- list.files(folder, recursive = FALSE,
                                  full.names = TRUE)
    
    for (file in folder_contents) {
      
      single_file_contents <- read.table(file, sep = ",", head = FALSE)
      
      if (nrow(single_file_contents) != n_landmarks) {
        stop(
          paste0("File ", file, " has an incorrect number of landmarks")
        )
      }
      
      landmarks <- abind::abind(
        landmarks,
        single_file_contents,
        along = 3
      )
      
    }
    
    labels <- c(labels, rep(
      substr(folder, nchar(folder_name) + 2, nchar(folder)),
      length(list.files(folder))
    ))
  }
  
  GraphGMM::write_morphologika_file(paste0("BSM_morphologika_file_", n_landmarks, "_smlms"),
                                    landmarks, as.factor(labels))
  
}

scale_value <- function(scale_bar, scale) {
  
  x <- scale_bar$x
  y <- scale_bar$y
  
  pixel_distance <- sqrt(
    (x[2] - x[1])^2 + (y[2] - y[1])^2
  )
  
  scale_factor <- scale / pixel_distance
  
  return(scale_factor)
  
}

convert_list_to_coordinate_matrix <- function(coordinate_list) {
  
  x <- coordinate_list$x
  y <- coordinate_list$y
  
  coordinate_matrix <- cbind(x, y)
  
  return(coordinate_matrix)
  
}

refresh_image <- function(input_image, landmarks, curves, scale_bar) {
  
  rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
  
  points(convert_list_to_coordinate_matrix(scale_bar), col = "orange", pch = 1)
  lines(convert_list_to_coordinate_matrix(scale_bar), col = "orange", lty = "11", lwd = 2)
  
  if (NA %!in% landmarks) {
    
    if (!is.matrix(landmarks)) {
      points(landmarks[1], landmarks[2], pch = 19, col = "red")
    } else {
      points(landmarks, pch = 19, col = "red")
    }
    
  }
  
  if (length(curves) != 0) {
    for (plot_single_curve in 1:length(curves)) {
      points(curves[[plot_single_curve]], pch = 19, col="blue")
      lines(curves[[plot_single_curve]], lty = 1, col = "blue", lwd = 2)
    }
  }
  
}

refresh_image_v2 <- function(input_image, curves, scale_bar) {
  
  rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
  
  points(convert_list_to_coordinate_matrix(scale_bar), col = "orange", pch = 1)
  lines(convert_list_to_coordinate_matrix(scale_bar), col = "orange", lty = "11", lwd = 2)
  
  colour_list = c("blue", "red", "darkgreen", "orange")
  
  for (plot_single_curve in 1:length(curves)) {
    points(curves[[plot_single_curve]], pch = 19, col = colour_list[plot_single_curve])
    lines(curves[[plot_single_curve]], lty = 1, col = colour_list[plot_single_curve], lwd = 2)
  }
  
}


digitise_image <- function(image_path,
                           n_landmarks,
                           n_curves,
                           n_semilandmarks,
                           scale = NULL,
                           verbose = FALSE,
                           external_window = TRUE,
                           output = c("txt", "tps")) {
  
  file_name <- substr(image_path, 1, nchar(image_path) - 4)
  
  input_image <- jpeg::readJPEG(image_path)
  
  if (external_window == TRUE) {
    dev.new(
      width = dim(input_image)[2], height = dim(input_image)[1], noRStudioGD = TRUE
    )
  }
  
  plot(
    seq(0, dim(input_image)[2], length.out = 10),
    seq(0, dim(input_image)[1], length.out = 10),
    type = "n", xlab = "",
    ylab = "", asp = 1, tck = 0,
    xaxt = "n", yaxt = "n", main = file_name)
  
  editor_mode = TRUE
  
  while (editor_mode == TRUE) {
    rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
    
    scale_bar <- locator(2, type = "o", lwd = 2, col = "orange", lty = "11")
    
    if (is.null(scale)) {
      scale <- as.numeric(readline(prompt = "Enter a scale value: "))
    }
    
    if (!is.numeric(scale) | scale <= 0) {
      stop(
        "Scale must be a positive integer."
      )
    }
    
    scale_factor <- scale_value(scale_bar, scale)
    
    if (verbose == FALSE) {
      
      all_points <- c(NA, NA)
      all_curves <- list()
      
      if (n_landmarks != 0) {
        
        for (point in 1:n_landmarks) {
          single_point <- locator(1, type = "p", col = "red", pch = 19)
          all_points <- rbind(all_points, c(single_point$x, single_point$y))
        }
        all_points <- all_points[2:nrow(all_points),]
        
      }
      
      if (n_curves != 0) {
        
        for (curve in 1:n_curves) {
          single_curve <- locator(type = "o", pch = 19, lwd = 2, col = "blue", lty = 1)
          
          curve_x <- single_curve$x
          curve_y <- single_curve$y
          
          geometric_curve <- c(curve_x[1], curve_y[1])
          for(geometric_point in 2:length(curve_x)) {
            geometric_curve <- rbind(geometric_curve,
                                     c(curve_x[geometric_point],
                                       curve_y[geometric_point]))
          }
          
          # L.A. Courtenay Method for Semilandmark equidistant placement - works on x axis
          
          line_interpolation <- approxfun(geometric_curve, method = "linear")
          
          if (length(n_semilandmarks) != 1) {
            step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / (n_semilandmarks[curve] - 1)
          } else {
            step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / (n_semilandmarks - 1)
          }
          semi_x <- seq(range(geometric_curve[,1])[1], range(geometric_curve[,1])[2], step)
          semi_y <- line_interpolation(semi_x)
          
          semilandmark_curve <- cbind(semi_x, semi_y)
          
          all_curves[[length(all_curves) + 1]] <- semilandmark_curve
          
          refresh_image(input_image, all_points, all_curves, scale_bar)
          
        }
      }
      
    } else {
      
      all_points <- c(NA, NA)
      all_curves <- list()
      
      if (n_landmarks != 0) {
        
        for (point in 1:n_landmarks) {
          
          continue_digitising = FALSE
          
          while(continue_digitising == FALSE) {
            single_point <- locator(1, type = "p", col = "red", pch = 19)
            
            continue_bool <- TRUE
            
            while(continue_bool == TRUE) {
              
              continue_verbose <- readline(prompt = "Would you like to keep this landmark? (y/n): ")
              
              if (continue_verbose == "y") {
                continue_digitising = TRUE
                all_points <- rbind(all_points, c(single_point$x, single_point$y))
                
                if (NA %in% all_points) {
                  all_points <- all_points[2:nrow(all_points),]
                }
                continue_bool <- FALSE
              } else if (continue_verbose == "n") {
                refresh_image(input_image, all_points, all_curves, scale_bar)
                continue_bool <- FALSE
              } else {
                cat("\nIncorrect Input\nPlease state whether you want to try again (y) or not (n)\n\n")
              }  
            }
            
          }
          
        }
        
      }
      
      if (n_curves != 0) {
        
        for (curve in 1:n_curves) {
          
          continue_digitising = FALSE
          
          while(continue_digitising == FALSE) {
            single_curve <- locator(type = "o", pch = 19, lwd = 2, col = "blue", lty = 1)
            
            curve_x <- single_curve$x
            curve_y <- single_curve$y
            
            geometric_curve <- c(curve_x[1], curve_y[1])
            for(geometric_point in 2:length(curve_x)) {
              geometric_curve <- rbind(geometric_curve,
                                       c(curve_x[geometric_point],
                                         curve_y[geometric_point]))
            }
            
            # L.A. Courtenay Method for Semilandmark equidistant placement - works on x axis
            
            line_interpolation <- approxfun(geometric_curve, method = "linear")
            
            if (length(n_semilandmarks) != 1) {
              step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / n_semilandmarks[curve]
            } else {
              step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / n_semilandmarks
            }
            semi_x <- seq(range(geometric_curve[,1])[1], range(geometric_curve[,1])[2], step)
            semi_y <- line_interpolation(semi_x)
            
            semilandmark_curve <- cbind(semi_x, semi_y)
            
            all_curves[[length(all_curves) + 1]] <- semilandmark_curve
            
            refresh_image(input_image, all_points, all_curves, scale_bar)
            
            continue_bool_2 <- TRUE
            
            while(continue_bool_2 == TRUE) {
              
              continue_verbose <- readline(prompt = "Would you like to keep these semilandmarks? (y/n): ")
              
              if (continue_verbose == "y") {
                continue_digitising = TRUE
                continue_bool_2 <- FALSE
              } else if (continue_verbose == "n") {
                
                if (length(all_curves) == 1){
                  all_curves <- list()
                } else {
                  all_curves <- all_curves[-length(all_curves)]
                }
                
                refresh_image(input_image, all_points, all_curves, scale_bar)
                continue_bool_2 <- FALSE
                
              } else {
                cat("\nIncorrect Input\nPlease state whether you want to try again (y) or not (n)\n\n")
              }
            }
            
          }
          
        }
      }
      
    }
    
    if (output == "txt") {
      
      if (n_curves != 0) {
        semilandmarks <- all_curves[[1]]
        
        if (n_curves != 1) {
          for (single_curve in 2:length(all_curves)) {
            semilandmarks <- abind::abind(semilandmarks, all_curves[[single_curve]], along = 1)
          }
        }
        
        if (n_landmarks != 0) {
          all_landmarks <- abind::abind(all_points, semilandmarks, along = 1)
          all_landmarks <- all_landmarks * scale_factor
        } else {
          all_landmarks <- semilandmarks * scale_factor
        }
        
      } else {
        
        all_landmarks <- all_points * scale_factor
        
      }
      
      colnames(all_landmarks) <- c("x", "y")
      
      save_path <- paste0(file_name, ".txt")
      
    } else if (output == "tps") {
      
      tps_output <- paste0("LM=", n_landmarks, "\n")
      
      if(n_landmarks > 1) {
        for (fixed_lm in 1:n_landmarks) {
          tps_output <- paste0(tps_output, all_points[fixed_lm, 1], " ", all_points[fixed_lm, 2], "\n")
        }
      } else {
        tps_output <- paste0(tps_output, all_points[1], " ", all_points[2], "\n")
      }
      
      if (n_curves != 0) {
        
        tps_output <- paste0(tps_output, "CURVES=", n_curves, "\n")
        
        for (curve_output in 1:n_curves) {
          
          tps_output <- paste0(tps_output, "POINTS=", dim(all_curves[[curve_output]])[1], "\n")
          
          semilandmarks <- all_curves[[curve_output]]
          
          for (curve_output_smlms in 1:dim(semilandmarks)[1]) {
            tps_output <- paste0(
              tps_output, semilandmarks[curve_output_smlms, 1], " ", semilandmarks[curve_output_smlms, 2], "\n"
            )
          }
          
        }
        
        tps_output <- paste0(tps_output, "IMAGE=", image_path, "\n")
        tps_output <- paste0(tps_output, "SCALE=", scale_factor, "\n")
        
      }
      
      save_path <- paste0(file_name, ".tps")
      
    } else {
      stop("Error in the type of output")
    }
    
    save_bool_loop <- TRUE
    
    while (save_bool_loop == TRUE) {
      
      editing_bool <- readline(prompt = "Would you like to save progress? (y/n): ")
      
      if (editing_bool == "y") {
        if (output == "tps") {
          cat(tps_output, file = save_path)
        } else {
          write.table(all_landmarks, save_path, row.names = FALSE, col.names = FALSE, sep = ",")
        }
        editor_mode = FALSE
        save_bool_loop = FALSE
      } else if (editing_bool == "n") {
        
        edit_exit_bool_loop <- TRUE
        while(edit_exit_bool_loop == TRUE) {
          
          repeat_bool <- readline(prompt = "Would you like to try again? (y/n): ")
          
          if (repeat_bool == "y") {
            edit_exit_bool_loop <- FALSE
            save_bool_loop <- FALSE
            editor_mode <- TRUE
          } else if (repeat_bool == "n") {
            edit_exit_bool_loop <- FALSE
            save_bool_loop <- FALSE
            editor_mode <- FALSE
          } else {
            cat("\nIncorrect Input\nPlease state whether you want to try again (y) or not (n)\n\n")
          }
        }
        
      } else {
        cat("\nIncorrect Input\nPlease state whether you want to save (y) or not (n)\n\n")
      }
      
    }
    
  }
  
  dev.off()
  cat("\n\nEditing Finished.")
  
}

digitise_images_2d <- function(data_path,
                               n_landmarks,
                               n_curves,
                               n_semilandmarks,
                               scale = NULL,
                               verbose = FALSE,
                               external_window = TRUE,
                               output = c("txt", "tps")) {
  
  if(missing(data_path)) {
    stop("No image path has been provided")
  }
  
  if (n_landmarks < 0 | n_landmarks %% 1 != 0) {
    stop("Invalid n_landmarks parameter provided")
  }
  
  if (n_curves < 0 | n_curves %% 1 != 0) {
    stop("Invalid n_curves parameter provided")
  }
  
  if (n_curves != 0) {
    if (length(n_semilandmarks) > 1) {
      
      for (smlm_value in 1:length(n_semilandmarks)) {
        if (n_semilandmarks[smlm_value] < 0 | n_semilandmarks[smlm_value] %% 1 != 0) {
          stop("Invalid n_semilandmarks parameter provided")
        }
        
        if (n_semilandmarks[smlm_value] <= 3) {
          stop("The number of semilandmarks on a curve must be more than 3")
        }
        
      }
      
    } else {
      if (n_semilandmarks < 0 | n_semilandmarks %% 1 != 0) {
        stop("Invalid n_semilandmarks parameter provided")
      }
      
      if (n_semilandmarks <= 3) {
        stop("The number of semilandmarks on a curve must be more than 3")
      }
      
    }
  }
  
  if (!is.logical(external_window)) {
    stop("external_visor parameter must be boolean (TRUE or FALSE)")
  }
  
  if (n_curves > 0) {
    if (length(n_semilandmarks) != n_curves) {
      warning(
        paste0(
          "\nOnly one n_semilandmarks parameter has been provided.\n",
          "All plotted n_curves will have the same number of semilandmarks unless",
          " specified otherwise."
        )
      )
    }
  }
  
  output <- match.arg(output)
  
  #directory_contents <- list.files(data_path, pattern = ".jpg")
  
  if (length(data_path) == 1) {
    
    digitise_image(
      image_path = data_path,
      n_landmarks = n_landmarks,
      n_curves = n_curves,
      n_semilandmarks = n_semilandmarks,
      scale = scale,
      verbose = verbose,
      external_window = external_window,
      output = output
    )
    
  } else if (length(data_path) > 1) {
    
    continue_processing_folder <- TRUE
    
    while(continue_processing_folder == TRUE) {
      
      for (image_to_digitise in 1:length(data_path)) {
        
        digitise_image(data_path[image_to_digitise],
                       n_landmarks = n_landmarks,
                       n_curves = n_curves,
                       n_semilandmarks = n_semilandmarks,
                       scale = scale,
                       verbose = verbose,
                       external_window = external_window,
                       output = output
        )
        
        correct_prompt <- FALSE
        
        if (image_to_digitise != length(data_path)) {
          
          while(correct_prompt == FALSE) {
            
            continue_processing_folder_prompt <- readline(
              prompt = "Would you like to continue processing files in this folder? (y/n): "
            )
            
            if (continue_processing_folder_prompt == "y") {
              continue_processing_folder <- TRUE
              correct_prompt <- TRUE
            } else if (continue_processing_folder_prompt == "n") {
              continue_processing_folder <- FALSE
              correct_prompt <- TRUE
            } else {
              cat("\nIncorrect Input\nPlease state whether you want to continue (y) or not (n)\n\n")
              correct_prompt <- FALSE
            }
            
          }
        } else {
          continue_processing_folder <- FALSE
        }
        
        if (continue_processing_folder == FALSE) {
          break
        }
        
      }
      
      if (continue_processing_folder == FALSE) {
        break
      }
      
    }
    
  } else {
    stop("Invalid data input")
  }
  
  cat("\n\nDigitisation completed.")
  
}

#

digitise_bsm <- function(image_path,
                         n_semilandmarks,
                         scale = NULL,
                         verbose = FALSE,
                         external_window = TRUE) {
  
  file_name <- substr(image_path, 1, nchar(image_path) - 4)
  
  input_image <- jpeg::readJPEG(image_path)
  
  if (external_window == TRUE) {
    dev.new(
      width = dim(input_image)[2], height = dim(input_image)[1], noRStudioGD = TRUE
    )
  }
  
  plot(
    seq(0, dim(input_image)[2], length.out = 10),
    seq(0, dim(input_image)[1], length.out = 10),
    type = "n", xlab = "",
    ylab = "", asp = 1, tck = 0,
    xaxt = "n", yaxt = "n", main = file_name)
  
  editor_mode = TRUE
  
  while (editor_mode == TRUE) {
    rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
    
    scale_bar <- locator(2, type = "o", lwd = 2, col = "orange", lty = "11")
    
    if (is.null(scale)) {
      scale <- as.numeric(readline(prompt = "Enter a scale value: "))
    }
    
    if (!is.numeric(scale) | scale <= 0) {
      stop(
        "Scale must be a positive integer."
      )
    }
    
    scale_factor <- scale_value(scale_bar, scale)
    
    all_curves <- list()
    
    if (verbose == FALSE) {
      
      single_curve <- locator(type = "o", pch = 19, lwd = 2, col = "blue", lty = 1)
      
      curve_x <- single_curve$x
      curve_y <- single_curve$y
      
      geometric_curve <- c(curve_x[1], curve_y[1])
      for(geometric_point in 2:length(curve_x)) {
        geometric_curve <- rbind(geometric_curve,
                                 c(curve_x[geometric_point],
                                   curve_y[geometric_point]))
      }
      
      # L.A. Courtenay Method for Semilandmark equidistant placement - works on x axis
      
      line_interpolation <- approxfun(geometric_curve, method = "linear")
      
      for (iteration in 1:length(n_semilandmarks)) {
        step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / (n_semilandmarks[iteration] - 1)
        
        semi_x <- seq(range(geometric_curve[,1])[1], range(geometric_curve[,1])[2], step)
        semi_y <- line_interpolation(semi_x)
        
        semilandmark_curve <- cbind(semi_x, semi_y)
        
        all_curves[[length(all_curves) + 1]] <- semilandmark_curve
        
      }
      
      refresh_image_v2(input_image, all_curves, scale_bar)
      
    } else {
      
      all_curves <- list()
      
      for (curve in 1:n_curves) {
        
        continue_digitising = FALSE
        
        while(continue_digitising == FALSE) {
          
          single_curve <- locator(type = "o", pch = 19, lwd = 2, col = "blue", lty = 1)
          
          curve_x <- single_curve$x
          curve_y <- single_curve$y
          
          geometric_curve <- c(curve_x[1], curve_y[1])
          for(geometric_point in 2:length(curve_x)) {
            geometric_curve <- rbind(geometric_curve,
                                     c(curve_x[geometric_point],
                                       curve_y[geometric_point]))
          }
          
          # L.A. Courtenay Method for Semilandmark equidistant placement - works on x axis
          
          line_interpolation <- approxfun(geometric_curve, method = "linear")
          
          for (iteration in 1:length(n_semilandmarks)) {
            step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / (n_semilandmarks[iteration] - 1)
            
            semi_x <- seq(range(geometric_curve[,1])[1], range(geometric_curve[,1])[2], step)
            semi_y <- line_interpolation(semi_x)
            
            semilandmark_curve <- cbind(semi_x, semi_y)
            
            all_curves[[length(all_curves) + 1]] <- semilandmark_curve
            
          }
          
          refresh_image_v2(input_image, all_curves, scale_bar)
          
          continue_bool_2 <- TRUE
          
          while(continue_bool_2 == TRUE) {
            
            continue_verbose <- readline(prompt = "Would you like to keep these semilandmarks? (y/n): ")
            
            if (continue_verbose == "y") {
              continue_digitising = TRUE
              continue_bool_2 <- FALSE
            } else if (continue_verbose == "n") {
              
              all_curves <- list()
              
              refresh_image_v2(input_image, all_points, all_curves, scale_bar)
              continue_bool_2 <- FALSE
              
            } else {
              cat("\nIncorrect Input\nPlease state whether you want to try again (y) or not (n)\n\n")
            }
          }
          
        }
        
      }
      
    }
    
    for (single_curve in 1:length(all_curves)) {
      all_curves[[single_curve]] <- all_curves[[single_curve]] * scale_factor
    }
    
    save_bool_loop <- TRUE
    
    while (save_bool_loop == TRUE) {
      
      editing_bool <- readline(prompt = "Would you like to save progress? (y/n): ")
      
      if (editing_bool == "y") {
        
        for (single_curve in 1:length(all_curves)) {
          
          save_path <- paste0(file_name, "_", n_semilandmarks[single_curve], "smlm.txt")
          
          landmark_data <- all_curves[[single_curve]]
          
          write.table(landmark_data, save_path, row.names = FALSE, col.names = FALSE, sep = ",")
          
        }
        
        editor_mode = FALSE
        save_bool_loop = FALSE
        
      } else if (editing_bool == "n") {
        
        edit_exit_bool_loop <- TRUE
        while(edit_exit_bool_loop == TRUE) {
          
          repeat_bool <- readline(prompt = "Would you like to try again? (y/n): ")
          
          if (repeat_bool == "y") {
            edit_exit_bool_loop <- FALSE
            save_bool_loop <- FALSE
            editor_mode <- TRUE
          } else if (repeat_bool == "n") {
            edit_exit_bool_loop <- FALSE
            save_bool_loop <- FALSE
            editor_mode <- FALSE
          } else {
            cat("\nIncorrect Input\nPlease state whether you want to try again (y) or not (n)\n\n")
          }
        }
        
      } else {
        cat("\nIncorrect Input\nPlease state whether you want to save (y) or not (n)\n\n")
      }
      
    }
    
  }
  
  dev.off()
  cat("\n\nEditing Finished.")
  
}

extract_bsm_profile_smlm <- function(data_path,
                                     n_semilandmarks = c(7, 15, 30),
                                     scale = NULL,
                                     verbose = FALSE,
                                     external_window = TRUE) {
  
  if(missing(data_path)) {
    stop("No image path has been provided")
  }
  
  if (length(n_semilandmarks) > 1) {
    
    for (smlm_value in 1:length(n_semilandmarks)) {
      if (n_semilandmarks[smlm_value] < 0 | n_semilandmarks[smlm_value] %% 1 != 0) {
        stop("Invalid n_semilandmarks parameter provided")
      }
      
      if (n_semilandmarks[smlm_value] <= 3) {
        stop("The number of semilandmarks on a curve must be more than 3")
      }
      
    }
    
  } else {
    if (n_semilandmarks < 0 | n_semilandmarks %% 1 != 0) {
      stop("Invalid n_semilandmarks parameter provided")
    }
    
    if (n_semilandmarks <= 3) {
      stop("The number of semilandmarks on a curve must be more than 3")
    }
    
  }
  
  if (!is.logical(external_window)) {
    stop("external_visor parameter must be boolean (TRUE or FALSE)")
  }
  
  if (length(data_path) == 1) {
    
    digitise_bsm(
      image_path = data_path,
      n_semilandmarks = n_semilandmarks,
      scale = scale,
      verbose = verbose,
      external_window = external_window
    )
    
  } else if (length(data_path) > 1) {
    
    continue_processing_folder <- TRUE
    
    while(continue_processing_folder == TRUE) {
      
      for (image_to_digitise in 1:length(data_path)) {
        digitise_bsm(
          image_path = data_path[image_to_digitise],
          n_semilandmarks = n_semilandmarks,
          scale = scale,
          verbose = verbose,
          external_window = external_window
        )
        
        correct_prompt <- FALSE
        
        if (image_to_digitise != length(data_path)) {
          
          while(correct_prompt == FALSE) {
            
            continue_processing_folder_prompt <- readline(
              prompt = "Would you like to continue processing files in this folder? (y/n): "
            )
            
            if (continue_processing_folder_prompt == "y") {
              continue_processing_folder <- TRUE
              correct_prompt <- TRUE
            } else if (continue_processing_folder_prompt == "n") {
              continue_processing_folder <- FALSE
              correct_prompt <- TRUE
            } else {
              cat("\nIncorrect Input\nPlease state whether you want to continue (y) or not (n)\n\n")
              correct_prompt <- FALSE
            }
            
          }
        } else {
          continue_processing_folder <- FALSE
        }
        
        if (continue_processing_folder == FALSE) {
          break
        }
        
      }
      
      if (continue_processing_folder == FALSE) {
        break
      }
      
    }
    
  } else {
    stop("Invalid data input")
  }
  
  cat("\n\nDigitisation completed.")
  
}

new_lmsmlm_method <- function(data_path,
                              n_semilandmarks = c(7, 15, 30),
                              scale = NULL,
                              verbose = FALSE,
                              external_window = TRUE) {
  
  digitise_bsm <- function(image_path,
                           n_semilandmarks,
                           scale = NULL,
                           verbose = FALSE,
                           external_window = TRUE) {
    
    scale_value <- function(scale_bar, scale) {
      
      x <- scale_bar$x
      y <- scale_bar$y
      
      pixel_distance <- sqrt(
        (x[2] - x[1])^2 + (y[2] - y[1])^2
      )
      
      scale_factor <- scale / pixel_distance
      
      return(scale_factor)
      
    }
    
    convert_list_to_coordinate_matrix <- function(coordinate_list) {
      
      x <- coordinate_list$x
      y <- coordinate_list$y
      
      coordinate_matrix <- cbind(x, y)
      
      return(coordinate_matrix)
      
    }
    
    refresh_image_v2 <- function(input_image, curves, scale_bar) {
      
      rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
      
      points(convert_list_to_coordinate_matrix(scale_bar), col = "orange", pch = 1)
      lines(convert_list_to_coordinate_matrix(scale_bar), col = "orange", lty = "11", lwd = 2)
      
      colour_list = c("blue", "red", "darkgreen", "orange")
      
      for (plot_single_curve in 1:length(curves)) {
        points(curves[[plot_single_curve]], pch = 19, col = colour_list[plot_single_curve])
        lines(curves[[plot_single_curve]], lty = 1, col = colour_list[plot_single_curve], lwd = 2)
      }
      
    }
    
    refresh_image_v3 <- function(input_image, first_curves, second_curves, scale_bar) {
      
      rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
      
      points(convert_list_to_coordinate_matrix(scale_bar), col = "orange", pch = 1)
      lines(convert_list_to_coordinate_matrix(scale_bar), col = "orange", lty = "11", lwd = 2)
      
      colour_list = c("blue", "red", "darkgreen", "orange")
      
      for (plot_single_curve in 1:length(first_curves)) {
        points(first_curves[[plot_single_curve]], pch = 19, col = colour_list[plot_single_curve])
        lines(first_curves[[plot_single_curve]], lty = 1, col = colour_list[plot_single_curve], lwd = 2)
      }
      
      for (plot_single_curve in 1:length(second_curves)) {
        points(second_curves[[plot_single_curve]], pch = 19, col = colour_list[plot_single_curve])
        lines(second_curves[[plot_single_curve]], lty = 1, col = colour_list[plot_single_curve], lwd = 2)
      }
      
    }
    
    `%!in%` = Negate(`%in%`)
    
    file_name <- substr(image_path, 1, nchar(image_path) - 4)
    
    input_image <- jpeg::readJPEG(image_path)
    
    if (external_window == TRUE) {
      dev.new(
        width = dim(input_image)[2], height = dim(input_image)[1], noRStudioGD = TRUE
      )
    }
    
    plot(
      seq(0, dim(input_image)[2], length.out = 10),
      seq(0, dim(input_image)[1], length.out = 10),
      type = "n", xlab = "",
      ylab = "", asp = 1, tck = 0,
      xaxt = "n", yaxt = "n", main = file_name)
    
    editor_mode = TRUE
    
    while (editor_mode == TRUE) {
      rasterImage(input_image, 1, 1, dim(input_image)[2], dim(input_image)[1])
      
      scale_bar <- locator(2, type = "o", lwd = 2, col = "orange", lty = "11")
      
      if (is.null(scale)) {
        scale <- as.numeric(readline(prompt = "Enter a scale value: "))
      }
      
      if (!is.numeric(scale) | scale <= 0) {
        stop(
          "Scale must be a positive integer."
        )
      }
      
      scale_factor <- scale_value(scale_bar, scale)
      
      all_first_curves <- list()
      all_second_curves <- list()
      
      if (verbose == FALSE) {
        
        # first curve
        
        first_curve <- locator(type = "o", pch = 19, lwd = 2, col = "blue", lty = 1)
        
        curve_x <- first_curve$x
        curve_y <- first_curve$y
        
        geometric_curve <- c(curve_x[1], curve_y[1])
        for(geometric_point in 2:length(curve_x)) {
          geometric_curve <- rbind(geometric_curve,
                                   c(curve_x[geometric_point],
                                     curve_y[geometric_point]))
        }
        
        # L.A. Courtenay Method for Semilandmark equidistant placement - works on x axis
        
        line_interpolation <- approxfun(geometric_curve, method = "linear")
        
        for (iteration in 1:length(n_semilandmarks)) {
          step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / (n_semilandmarks[iteration] - 1)
          
          semi_x <- seq(range(geometric_curve[,1])[1], range(geometric_curve[,1])[2], step)
          semi_y <- line_interpolation(semi_x)
          
          semilandmark_curve <- cbind(semi_x, semi_y)
          
          all_first_curves[[length(all_first_curves) + 1]] <- semilandmark_curve
          
        }
        
        refresh_image_v2(input_image, all_first_curves, scale_bar)
        
        # second curve
        
        second_curve <- locator(type = "o", pch = 19, lwd = 2, col = "blue", lty = 1)
        
        curve_x <- second_curve$x
        curve_y <- second_curve$y
        
        geometric_curve <- c(curve_x[1], curve_y[1])
        for(geometric_point in 2:length(curve_x)) {
          geometric_curve <- rbind(geometric_curve,
                                   c(curve_x[geometric_point],
                                     curve_y[geometric_point]))
        }
        
        # L.A. Courtenay Method for Semilandmark equidistant placement - works on x axis
        
        line_interpolation <- approxfun(geometric_curve, method = "linear")
        
        for (iteration in 1:length(n_semilandmarks)) {
          step <- (range(geometric_curve[,1])[2] - range(geometric_curve[,1])[1]) / (n_semilandmarks[iteration] - 1)
          
          semi_x <- seq(range(geometric_curve[,1])[1], range(geometric_curve[,1])[2], step)
          semi_y <- line_interpolation(semi_x)
          
          semilandmark_curve <- cbind(semi_x, semi_y)
          
          all_second_curves[[length(all_second_curves) + 1]] <- semilandmark_curve
          
        }
        
        refresh_image_v3(input_image, all_first_curves, all_second_curves, scale_bar)
        
      } else {
        
        stop("No verbose option available any more for this function.")
        
      }
      
      for (single_curve in 1:length(all_first_curves)) {
        all_first_curves[[single_curve]] <- all_first_curves[[single_curve]] * scale_factor
      }
      
      for (single_curve in 1:length(all_second_curves)) {
        all_second_curves[[single_curve]] <- all_second_curves[[single_curve]] * scale_factor
      }
      
      save_bool_loop <- TRUE
      
      while (save_bool_loop == TRUE) {
        
        editing_bool <- readline(prompt = "Would you like to save progress? (y/n): ")
        
        if (editing_bool == "y") {
          
          for (single_curve in 1:length(all_first_curves)) {
            
            save_path <- paste0(file_name, "_", n_semilandmarks[single_curve], "smlm.txt")
            
            extract_first_curve <- all_first_curves[[single_curve]]
            extract_second_curve <- all_second_curves[[single_curve]]
            
            extract_first_curve <- extract_first_curve[-nrow(extract_first_curve),]
            
            landmark_data <- abind::abind(extract_first_curve, extract_second_curve, along = 1)
            
            write.table(landmark_data, save_path, row.names = FALSE, col.names = FALSE, sep = ",")
            
          }
          
          editor_mode = FALSE
          save_bool_loop = FALSE
          
        } else if (editing_bool == "n") {
          
          edit_exit_bool_loop <- TRUE
          while(edit_exit_bool_loop == TRUE) {
            
            repeat_bool <- readline(prompt = "Would you like to try again? (y/n): ")
            
            if (repeat_bool == "y") {
              edit_exit_bool_loop <- FALSE
              save_bool_loop <- FALSE
              editor_mode <- TRUE
            } else if (repeat_bool == "n") {
              edit_exit_bool_loop <- FALSE
              save_bool_loop <- FALSE
              editor_mode <- FALSE
            } else {
              cat("\nIncorrect Input\nPlease state whether you want to try again (y) or not (n)\n\n")
            }
          }
          
        } else {
          cat("\nIncorrect Input\nPlease state whether you want to save (y) or not (n)\n\n")
        }
        
      }
      
    }
    
    dev.off()
    cat("\n\nEditing Finished.")
    
  }
  
  if(missing(data_path)) {
    stop("No image path has been provided")
  }
  
  if (length(n_semilandmarks) > 1) {
    
    for (smlm_value in 1:length(n_semilandmarks)) {
      if (n_semilandmarks[smlm_value] < 0 | n_semilandmarks[smlm_value] %% 1 != 0) {
        stop("Invalid n_semilandmarks parameter provided")
      }
      
      if (n_semilandmarks[smlm_value] <= 3) {
        stop("The number of semilandmarks on a curve must be more than 3")
      }
      
    }
    
  } else {
    if (n_semilandmarks < 0 | n_semilandmarks %% 1 != 0) {
      stop("Invalid n_semilandmarks parameter provided")
    }
    
    if (n_semilandmarks <= 3) {
      stop("The number of semilandmarks on a curve must be more than 3")
    }
    
  }
  
  if (!is.logical(external_window)) {
    stop("external_visor parameter must be boolean (TRUE or FALSE)")
  }
  
  if (length(data_path) == 1) {
    
    digitise_bsm(
      image_path = data_path,
      n_semilandmarks = n_semilandmarks,
      scale = scale,
      verbose = verbose,
      external_window = external_window
    )
    
  } else if (length(data_path) > 1) {
    
    continue_processing_folder <- TRUE
    
    while(continue_processing_folder == TRUE) {
      
      for (image_to_digitise in 1:length(data_path)) {
        digitise_bsm(
          image_path = data_path[image_to_digitise],
          n_semilandmarks = n_semilandmarks,
          scale = scale,
          verbose = verbose,
          external_window = external_window
        )
        
        correct_prompt <- FALSE
        
        if (image_to_digitise != length(data_path)) {
          
          while(correct_prompt == FALSE) {
            
            continue_processing_folder_prompt <- readline(
              prompt = "Would you like to continue processing files in this folder? (y/n): "
            )
            
            if (continue_processing_folder_prompt == "y") {
              continue_processing_folder <- TRUE
              correct_prompt <- TRUE
            } else if (continue_processing_folder_prompt == "n") {
              continue_processing_folder <- FALSE
              correct_prompt <- TRUE
            } else {
              cat("\nIncorrect Input\nPlease state whether you want to continue (y) or not (n)\n\n")
              correct_prompt <- FALSE
            }
            
          }
        } else {
          continue_processing_folder <- FALSE
        }
        
        if (continue_processing_folder == FALSE) {
          break
        }
        
      }
      
      if (continue_processing_folder == FALSE) {
        break
      }
      
    }
    
  } else {
    stop("Invalid data input")
  }
  
  cat("\n\nDigitisation completed.")
  
}
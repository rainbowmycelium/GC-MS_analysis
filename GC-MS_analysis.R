## function matching two lists of compounds by rt##
## sample_df = dataframe with chromatogram data from sample ##
## control_df = dataframe with chromatogram data from control sample ##
## returns one dataframe with two chromatograms aligned ##

match_by_rt <- function(sample_df, control_df, rt_window) {
  matched_list <- list()
  for (i in 1:nrow(sample_df)) {
    sample_rt <- sample_df$rt[i]
    matched_control <- control_df %>% filter(abs(rt - sample_rt) <= rt_window)
    if (nrow(matched_control) > 0) {
      matched_list[[length(matched_list) + 1]] <- matched_control
    } else {
      matched_control[1,] <- c(sample_rt[1], "-", "-", "-")
      matched_list[[length(matched_list) + 1]] <- matched_control
    }
  }
  matched_df <- do.call(rbind, matched_list)
  missing <- setdiff(control_df$rt, matched_df$rt)
  additional_df <- control_df %>% filter(rt %in% missing)
  additional_df_1 <- data.frame(rt=additional_df$rt, area=rep("-", times=nrow(additional_df)), compound=rep("-", times=nrow(additional_df)), qual=rep("-", times=nrow(additional_df)))
  colnames(additional_df) <- paste(colnames(additional_df), c("kontrola"), sep=".")
  additional_df <- cbind(additional_df_1,additional_df)
  colnames(matched_df) <- paste(colnames(matched_df), c("kontrola"), sep=".")
  matched_df <- as.data.frame(cbind(sample_df, matched_df))
  matched_df <- rbind(matched_df, additional_df)
  return(matched_df)
}


## function alligning multiple chromatograms ##
## proba = dataframe with chromatogram data from a sample ##
## kontrole = list of dataframes with chromatogram data from control samples ##
## okno_rt = acceptable retention time bias between peeks from samples and controls ##
## returns one data frame with multiple chromatograms aligned ##

multiple_match <- function(proba, kontrole, okno_rt){
  jeden <- match_by_rt(sample_df = proba, control_df = kontrole[[1]], rt_window = okno_rt[1])
  jeden$rt <- as.numeric(jeden$rt)
  lista_wynikow <- list(jeden)
  if(length(kontrole)==1){
    jeden$rt <- as.numeric(jeden$rt)
    return(jeden)
  }else{
    for(g in 2:length(kontrole)){
      probsko <- lista_wynikow[[g-1]][,c(1,2,3,4)]
      wyniczyska <- match_by_rt(sample_df = probsko, control_df = kontrole[[g]], rt_window = okno_rt[1])
      roznica <- nrow(wyniczyska) - nrow(lista_wynikow[[g-1]])
      kontrolsko <- lista_wynikow[[g-1]][,-c(1,2,3,4)]
      probsko_1 <- wyniczyska[,c(1,2,3,4)]
      rzedy <- seq(from=(nrow(probsko_1)-roznica+1), to=nrow(probsko_1), by = 1)
      if(ncol(probsko_1)==ncol(kontrolsko)){
        colnames(probsko_1) <- colnames(kontrolsko)
        kontrolsko <- rbind(kontrolsko, probsko_1[rzedy,])
      }else{
        probsko_2 <- probsko_1
        for(h in 1:((ncol(kontrolsko)/ncol(probsko_1))-1)){
          probsko_1 <- cbind(probsko_1, probsko_2) 
        }
        colnames(probsko_1) <- colnames(kontrolsko)
        kontrolsko <- rbind(kontrolsko, probsko_1[rzedy,])
      }
      colnames(kontrolsko) <- paste(colnames(kontrolsko),".",sep="")
      wyniczyska <- cbind(wyniczyska, kontrolsko)
      lista_wynikow[[g]] <- wyniczyska
    }
    koniec <- lista_wynikow[[length(kontrole)]]
    koniec$rt <- as.numeric(koniec$rt)
    koniec <- arrange(koniec, rt)
    return(koniec)
  }
}


## a function matching chromatograms of samples with chromatograms of control samples##
## and excluding peeks that appear in the controls, provided they have relative surface not bigger by 0.1 when compared to control ##
## proby = list of dataframes with chromatogram data  from samples ##
## kontrole = list of dataframes with chromatogram data from cotnrol samples ##
## okno_rt = acceptable retention time bias between peeks from samples and controls ##
## returns a list of samples chromatogram data data frames without peeks excluded ##

exclusion <- function(proby, kontrole, okno_rt){
  koncowe <- list()
  for(w in 1:length(proby)){
    res <- multiple_match(proba = proby[[w]], kontrole = kontrole, okno_rt = okno_rt)
    area_cols <- grep("area", names(res))
    res[area_cols] <- lapply(res[area_cols], function(x) ifelse(x == "-", NA, as.numeric(x)))
    select_rows <- apply(res[area_cols], 1, function(x) {
      first_area <- x[1]
      other_areas <- x[-1]
      !is.na(first_area) && 
        all(is.na(other_areas) | first_area - other_areas >= 0.1)
    })
    res1 <- res[select_rows, ]
    res1 <- res1[,1:4]
    koncowe[[w]] <- res1
    names(koncowe)[w] <- names(proby)[w]
  }
  return(koncowe)
}

## function that discharges rows of aligned chromatogram dataframes with id quality < 80 % ##
## unless there is a peek at the same rt time on annother chromatogram aligned which has id qual >= 80 ##
## df = dataframe of aligned chromatograms ##
## returns the same dataframe but without excluded rows ##

discharge_rows <- function(df) {
  qual_cols <- grep("qual", names(df))
  df[qual_cols] <- lapply(df[qual_cols], function(x) ifelse(x == "-", NA, as.numeric(x)))
  df$x <- ifelse(rowSums(df[qual_cols] < 80, na.rm = TRUE) > 0, "x", "")
  df <- df[rowSums(df[qual_cols] >= 80, na.rm = TRUE) > 0, ]
  return(df)
}

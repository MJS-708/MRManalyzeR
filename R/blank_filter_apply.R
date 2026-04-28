# Function to apply blank filter

blank_filter_apply = function(df, blank_samples, blank_filter){

  df_blank = df[rownames(df) %in% blank_samples, ]

  blank_df = data.frame(feature = names(df_blank),
                        blank_val = colMeans(df_blank, na.rm=T))

  blank_df$blank_thresh = blank_df$blank_val * blank_filter

  for(col_idx in which(!is.na(blank_df$blank_thresh))){

    thresh = blank_df$blank_thresh[col_idx]

    temp_col = df[,col_idx]
    temp_col = ifelse(temp_col < thresh, NA, temp_col)

    df[,col_idx] = temp_col
  }

  return(df)
}

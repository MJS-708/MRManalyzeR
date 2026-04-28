# function to calculate conc in sample

# Calculate the concentration of lipid in each sample - TO BE COMPLETED
    # conc from cal curve (ng/mL) - based on known conc of cal mixes and adjustment factor for each std group
# To find conc in each sample
    # conc = conc * (std_vol / sample_vol) * (sample_IS / std_IS)

calculate_conc = function(df, metadata_batch){

  sample_IS_vol_uL = unique(metadata_batch$sample_IS_vol_uL)
  cal_IS_vol_uL = unique(metadata_batch$cal_IS_vol_uL)
  sample_volume_uL = unique(metadata_batch$sample_volume_uL)
  cal_vol_uL = unique(metadata_batch$cal_vol_uL)

  if( all(c(length(sample_IS_vol_uL), length(cal_IS_vol_uL), length(sample_volume_uL), length(cal_vol_uL)) == 1)){

    IS_fact = sample_IS_vol_uL / cal_IS_vol_uL
    vol_fact = cal_vol_uL / sample_volume_uL

    df = df * vol_fact * IS_fact

  } else{

    if(nrow(df) != nrow(metadata_batch)){
      print("FAIL - variable volumes, implement sample wise. But metadata and lcms data lengths differ")
    }

    for(i in 1:nrow(df)){

      IS_fact = metadata_batch$sample_IS_vol_uL[i] / metadata_batch$cal_IS_vol_uL[i]
      vol_fact = metadata_batch$cal_vol_uL[i] / metadata_batch$sample_volume_uL[i]

      df[i,] = df[i,] * vol_fact * IS_fact

    }
  }
  return(df)
}

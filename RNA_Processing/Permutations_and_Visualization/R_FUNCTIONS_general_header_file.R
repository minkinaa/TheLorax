
split_on_delimiter = function(one_col, delim, pos_to_return){
  val_to_split = as.character(one_col[1])
  return(strsplit(val_to_split, delim)[[1]][pos_to_return])
}

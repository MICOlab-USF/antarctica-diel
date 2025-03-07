kegg_API <- function(kegg,WhichPath = "Last"){
  require(stringi)
  
  pathway <- system(paste0("curl https://rest.kegg.jp/link/pathway/", kegg), intern=TRUE)
  n.paths <- grep("\tpath:ko",pathway)
  
  if(isEmpty(n.paths)){
    return("-")
  } else if(WhichPath %in% c("First","Last")){
    
    ifelse(WhichPath == "First",
           stri_sub(pathway[n.paths[1]],16),
           stri_sub(pathway[n.paths[length(n.paths)]],16))
    
  } else{
    paste(stri_sub(pathway[n.paths],16),collapse = ";")
  }
}

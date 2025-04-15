pathway_API <- function(input,ReturnFull = 0){
  pathway <- system(paste0("curl https://rest.kegg.jp/get/", input), intern=TRUE)
  
  Name_idx <- grep("NAME",pathway)
  Class_idx <- grep("CLASS",pathway)
  Module_idx <- grep("^MODULE|^\\s+M\\d{5}", pathway)
  Definition_idx <- grep("DEFINITION",pathway)
  
  RelPaths_idx <- grep("REL_PATHWAY",pathway)
  RelPaths_end <- grep("///",pathway)-1
  
  name <- ifelse(isEmpty(Name_idx),
                 NA,
                 sub("^\\S+\\s+", "", pathway[Name_idx]))
  
  class <- ifelse(isEmpty(Class_idx),
                  NA,
                  sub("^\\S+\\s+", "", pathway[Class_idx]))
  
  if(length(Module_idx) == 0){
    module <- NA
  }else{
    module <- unlist(regmatches(pathway[Module_idx],gregexpr("\\bM\\d{5}\\b", pathway[Module_idx])))
  }
  
  if(isEmpty(Definition_idx)){
    definition <- NA
  }else{
    definition <- unique(unlist(regmatches(pathway[Definition_idx],gregexpr("K\\d{5}\\b", pathway[Definition_idx]))))
  }
  
  
  if(!isEmpty(RelPaths_idx)){
    RelPaths <- sub(".*\\sko", "ko", pathway[RelPaths_idx:RelPaths_end])
  } else {
    RelPaths <- NA
  }
  
  output <- list("ko" = input,
                 "name" = name,
                 "class" = class,
                 "Rel_paths" = RelPaths,
                 "module" = module,
                 "definition" = definition)
  
  if(ReturnFull == 1){
    output$Full <- pathway
  }
  
  return(output)
}
add_model_id <- function(res,i){
    model_id <- NULL;
    mid <- gsub("file",paste0("model",i,"_"),basename(tempfile()))
    res[,model_id:=mid,]
}

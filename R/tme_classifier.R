




#' TME classifier for gastric cancer
#'
#' @param eset expression set with genes names as row and sample id as column
#' @param tme_deconvolution default is FALSE, if set to TRUE, this process will be time consuming.
#' @param tme_data user can provide TME data estimated by CIBERSORT and MCPcounter, please modify the names of cell types, a reference data could be found by `eset_example`
#' @param adjust_eset default is FALSE, it means that variables with missing value will to be preserved and will be replaced by mean value
#' @param replace_na default is TRUE, variables with missing value will to be preserved and will be replaced by mean value of observations
#' @param log2trans default is FALSE,
#' @param scale default is FALSE
#' @param array default is FALSE, if expression set was derived from micro array, this parameter should be set to TRUE, which will affect the TME deconvolution process
#' @param save_data default is FALSE, if TRUE, processing data will be save to local path
#' @param perm permutation of CIBERSORT
#' @param method default is `ensemble`, other options: `svm`, `rf`, `nnet`, `knn`, `dt`, `xgboost`
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#' @examples
#' tme_classifier(eset = eset_example1, method = "ensemble", scale = T)
tme_classifier<-function(eset,
                         method            = "ensemble",
                         tme_data          = NULL,
                         tme_deconvolution = FALSE,
                         array             = FALSE,
                         save_data         = FALSE,
                         adjust_eset       = FALSE,
                         replace_na        = TRUE,
                         log2trans         = FALSE,
                         scale             = FALSE,
                         perm              = 100){


  cat(crayon::green("Step-1: Expression data preprocessing...\n"))

  # message("Step-1: Expression data preprocessing...")

  if(sum(is.na(eset)>0)) message(">>> Parameter `adjust_eset` must be FALSE, if variables with missing value wanted to be preserved")

  if(adjust_eset){
    message(">>> Variables with NA will be deleted when `adjust_eset` is TRUE")
    feas<-IOBR::feature_manipulation(data=eset, is_matrix = T)
    eset<-eset[rownames(eset)%in%feas,]
  }


  if(replace_na){

    if(sum(is.na(eset)>0)){
      message(">>> Retain NA variables, replaced by mean value of all observations")
      teset<-as.data.frame(t(eset))
      for(i in c(1:ncol(teset))){
        teset[is.na(teset[,i]),i]<-mean(teset[!is.na(teset[,i]),i])
      }
      eset<-t(teset)

    }else{
      message(">>> There are no missing values")
    }

    if(log2trans){

      if(log2trans) message(">>> If paramater `log2trans` is TRUE, `replace_na` must be TRUE to replace NA before log2 transformation.")

      feas<-IOBR::feature_manipulation(data=eset,is_matrix = T)
      eset<-eset[rownames(eset)%in%feas,]
      eset<- log2(eset)
      if(ncol(eset) <5000) IOBR::check_eset(eset)

      message(">>> log2 transforamtion was finished.")
    }

  }else{
    message(">>>-- Missing value were preserved.")
  }

  message(" ")

  cat(crayon::green("Step-2: TME deconvolution...\n"))
  # message("Step-2: TME deconvolution...")

  if(is.null(tme_data)){
    if(tme_deconvolution){

      message(">>>-- This process is time-consuming. Please wait patiently.")


      cibersort<-IOBR:: deconvo_tme(eset = eset,method = "cibersort", arrays = array, perm = perm )
      mcp      <-IOBR:: deconvo_tme(eset = eset,method = "mcpcounter")



      print(head(cibersort))
      print(head(mcp))

      cibersort<-as.data.frame(cibersort)
      mcp<-as.data.frame(mcp)
      tme_data<-merge(cibersort, mcp, by= "ID", all = FALSE)

      if(save_data) save(tme_data, file = "TME_data.RData")

    }
  }

  if(!tme_deconvolution) message(">>> This step was skipped, user can set parameter `tme_deconvolution` to TRUE or provide TME data to realize prediction.")

  cat(crayon::green(">>>--  More TME deconvolution algorithms can be reached from:  \n"))
  cat(crayon::green(">>>--  https://github.com/IOBR/IOBR  \n"))

  eset2<-eset[rownames(eset)%in%c(feas1, feas2), ]
  ####################################################################
  if(scale){

    message(">>>-- Scaling data...")
    if(sum(is.na(eset2))>0){
      message(">>> Before scaling, NA will be repleased by mean value")

      teset<-as.data.frame(t(eset2))
      for(i in c(1:ncol(teset))){
        teset[is.na(teset[,i]),i]<-mean(teset[!is.na(teset[,i]),i])
      }
      eset2<-t(teset)
      message(paste0(">>> Counts of NA after replacement = "), sum(is.na(eset2)))

      feas<-IOBR::feature_manipulation(data=eset2,is_matrix = T)
      eset2<-eset2[rownames(eset2)%in%feas,]

    }
    eset2<-scale(t(eset2), center = T,scale = T)
  }else{
    eset2<-as.data.frame(t(eset))
  }


  if(!is.null(tme_data)){

    eset2<-tibble:: rownames_to_column(as.data.frame(eset2), var = "ID")
    eset2<- merge(eset2, tme_data, by = "ID")
    eset2<- tibble::column_to_rownames(eset2, var = "ID")

    data("eset_example2")
    eset2<- IOBR::assimilate_data(eset_example2, eset2)
    eset2[is.na(eset2)]<-0

  }else{

    data("eset_example1")
    eset2<- IOBR::assimilate_data(eset_example1, eset2)
    eset2[is.na(eset2)]<-0
  }
  ###########################

  if(save_data) save(eset2, file = "Processed_data.RData")
  ###########################

  message(" ")
  cat(crayon::green("Step-3: Predicting TME phenotypes...\n"))
  # message("Step-3: Predicting TME phenotypes...")


  if(is.null(tme_data)|length(colnames(eset2)[grep(colnames(eset2),pattern = "CIBERSORT")])==0){
    svmModel<-svmModel
    NNModel<-NNModel
    DecTreeModel<-DecTreeModel
    rfModel<-rfModel
    knnModel<-knnModel
    xgModel<-xgModel
  }else{
    svmModel<-svmModel2
    NNModel<-NNModel2
    DecTreeModel<-DecTreeModel2
    rfModel<-rfModel2
    knnModel<-knnModel2
    xgModel<-xgModel2
  }

  if(method=="ensemble"){

    message(">>>--- Ensemble Model was used to predict TME phenotypes... ")

    p1 <- predict(svmModel, eset2, preProcess = c("center", "scale"),type = "prob")
    p2 <- predict(NNModel, eset2, preProcess = c("center", "scale"),type = "prob")
    p3 <- predict(DecTreeModel, eset2, preProcess = c("center", "scale"),type = "prob")
    p4 <- predict(rfModel, eset2, preProcess = c("center", "scale"),type = "prob")
    p5 <- predict(knnModel, eset2, preProcess = c("center", "scale"),type = "prob")
    p6 <- predict(xgModel, eset2, preProcess = c("center", "scale"),type = "prob")
    # --------------------------------------------------------------------
    # create simple average ensemble // some are [0,1] some are probabilities
    # simple average, better would be glm to test best combination

    p =  ((p1+p2+p3+p4+p5+p6)/6)

  }else if(method=="nnet"){
    message(">>>--- The Neural Network Model was used to predict TME phenotypes... ")
    p <- predict(NNModel, eset2, preProcess = c("center", "scale"),type = "prob")
  }else if(method=="dt"){
    message(">>>--- The Decsion Tree: C5.0 Model was used to predict TME phenotypes... ")
    p <- predict(DecTreeModel, eset2, preProcess = c("center", "scale"),type = "prob")
  }else if(method=="rf"){
    message(">>>--- The Random Forrest Model was used to predict TME phenotypes... ")
    p <- predict(rfModel, eset2, preProcess = c("center", "scale"),type = "prob")
  }else if(method=="knn"){
    message(">>>--- The KNN Model was used to predict TME phenotypes... ")
    p <- predict(knnModel, eset2, preProcess = c("center", "scale"),type = "prob")
  }else if(method=="svm"){
    message(">>>--- The Support Vector Machines Model was used to predict TME phenotypes... ")
    p <- predict(svmModel, eset2, preProcess = c("center", "scale"),type = "prob")
  }else if(method=="xgboost"){
    message(">>>--- The xgboost Model was used to predict TME phenotypes... ")
    p <- predict(xgModel, eset2, preProcess = c("center", "scale"),type = "prob")
  }


  submission <- data.frame(rownames(eset2), p)
  colnames(submission)[1] = "ID"
  cluster<-colnames(submission)[2:4]
  submission$res_ensemble<-  cluster[apply(submission[,2:4], 1, which.max)]

  colnames(submission)[5]<-"TMEcluster"

  message(" ")
  ############################################################
  cat(crayon::red(">>>--- DONE!\n"))

  return(submission)
}



#' @title Epigenetic Modification Prediction
#' @description Predicting sequences with 6mA sites.
#' @param FastaData Sequence file (.fasta format)
#' @param Species Model organism
#' @return MethStatus: Sequences with their methylation state (methylated or non-methylated)
#' @import stats devtools tidyverse seqinr Biostrings splitstackshape entropy party stringr tibble doParallel parallel e1071 caret randomForest gbm foreach iterators ftrCOOL
#' @export
#'
#' @examples
#' \donttest{
#' library(EpiSemble)
#' data<-system.file("exdata/test.fasta", package = "EpiSemble")
#' pred<-epiPred(FastaData=data, Species="Rice")
#' }
#' @references Chen, W., Lv, H., Nie, F., & Lin, H. (2019). i6mA-Pred: identifying DNA N6-methyladenine 	sites in the rice genome. Bioinformatics, 35(16), 2796-2800.
requireNamespace("Biostrings","caret","devtools","doParallel","e1071","entropy","foreach","ftrCOOL","gbm","iterators","parallel","party","randomForest","seqinr","splitstackshape","stats","stringr","tibble","tidyverse")
epiPred<- function(FastaData,Species){
  `%dopar%` <- foreach::`%dopar%`
  ImpFeatures<-function(Fastafile){
    FastaToTabular <- function (filename){

      #read fasta file

      file1 <- readLines(filename)

      #find the genename location by grepping >

      location <- which((stringr::str_sub(file1,1,1))==">")

      #start an empty vector to collect name and sequence

      name=c()
      sequence =c()



      #number of genes= number of loops
      #extract name first
      for ( i in 1:length(location)){
        name_line = location[i]
        name1 = file1[name_line]
        name=c(name,name1)
        start= location[i]+1
        end = location[i+1]-1
        if ( i < length (location)){

          end=end

        } else {

          end=length(file1)
        }

        lines = start:end
        sequence1= as.character(paste(file1[lines],collapse = ""))
        sequence =c(sequence,sequence1)
      }


      data <- tibble(name,sequence)



      #finally export the file
      #before that remove preexisting file
      unlink(c("dna_table.csv"),force=TRUE)
      as.matrix(data,"dna_table.csv")

      #function ends
    }
    ########alphabetcheck##########
    alphabetCheck<-function (sequences, alphabet = "aa", label = c())
    {
      if (length(sequences) == 0) {
        stop("ERROR: sequence parameter is empty")
      }
      if (length(label) != 0 && length(label) != length(sequences)) {
        stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
      }
      if (alphabet == "rna") {
        alphabet <- c("A", "C", "G", "U")
      }
      else if (alphabet == "dna") {
        alphabet <- c("A", "C", "G", "T")
      }
      else if (alphabet == "aa") {
        alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                      "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                      "W", "Y")
      }
      else {
        stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
      }
      alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                                 split = "")[[1]] %in% alphabet))
      flag = 0
      if (length(label) == length(sequences)) {
        flag = 1
        label = label[alphabetCheck]
      }
      else if (length(label) > 0 && length(label) != length(sequences)) {
        stop("ERROR: The number of labels is not equal to the number of sequences!")
      }
      if (is.null(names(sequences))) {
        names(sequences) <- as.character(1:length(sequences))
      }
      nonstanSeq <- names(sequences)[!alphabetCheck]
      if (length(nonstanSeq) != 0) {
        nonstanSeq <- toString(nonstanSeq)
        warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
        message(warMessage)
      }
      sequences = sequences[alphabetCheck]
      if (length(sequences) == 0) {
        stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
      }
      if (flag == 1) {
        names(label) = names(sequences)
      }
      seq_lab <- list(sequences = sequences, Lab = label)
      return(seq_lab)
    }
    ##########NCP_DNA##########

    ncp_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                       label = c())
    {
      if (length(seqs) == 1 && file.exists(seqs)) {
        seqs <- fa.read(seqs, alphabet = "dna")
        seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
        seqs <- seqs_Lab[[1]]
        label <- seqs_Lab[[2]]
      }
      else if (is.vector(seqs)) {
        seqs <- sapply(seqs, toupper)
        seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
        seqs <- seqs_Lab[[1]]
        label <- seqs_Lab[[2]]
      }
      else {
        stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
      }
      lenSeqs <- sapply(seqs, nchar)
      nucs <- list(A = c(1, 1, 1), C = c(0, 0, 1), G = c(1, 0,
                                                         0), T = c(0, 1, 0), U = c(0, 1, 0))
      numSeqs <- length(seqs)
      if (outFormat == "mat") {
        if (length(unique(lenSeqs)) > 1) {
          stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
        }
        if (binaryType == "strBin") {
          nucs <- c(A = "111", C = "001", G = "100", T = "010",
                    U = "010")
          featureMatrix <- sapply(seqs, function(x) {
            charList <- unlist(strsplit(x, split = ""))
            cods <- nucs[charList]
            return(cods)
          })
          featureMatrix <- t(featureMatrix)
          colnames(featureMatrix) <- paste("pos", 1:lenSeqs[1],
                                           sep = "")
          row.names(featureMatrix) <- names(seqs)
        }
        else if (binaryType == "logicBin") {
          nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                      TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                     FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
          featureMatrix <- sapply(seqs, function(x) {
            charList <- unlist(strsplit(x, split = ""))
            cods <- nucs[charList]
            cods <- unlist(cods)
            return(cods)
          })
          featureMatrix <- t(featureMatrix)
          temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
          temp2 <- rep(1:lenSeqs[1], each = 3)
          colnames(featureMatrix) <- paste("pos", temp2, "-",
                                           temp1, sep = "")
          row.names(featureMatrix) <- names(seqs)
        }
        else if (binaryType == "numBin") {
          featureMatrix <- sapply(seqs, function(x) {
            charList <- unlist(strsplit(x, split = ""))
            cods <- nucs[charList]
            cods <- unlist(cods)
            return(cods)
          })
          featureMatrix <- t(featureMatrix)
          temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
          temp2 <- rep(1:lenSeqs[1], each = 3)
          colnames(featureMatrix) <- paste("pos", temp2, "-",
                                           temp1, sep = "")
          row.names(featureMatrix) <- names(seqs)
        }
        else {
          stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
        }
        return(featureMatrix)
      }
      else if (outFormat == "txt") {
        nucs <- c(A = "111", C = "001", G = "100", T = "010",
                  U = "010")
        counter <- 0
        namesSeqs <- names(seqs)
        codes <- lapply(seqs, function(x) {
          counter <- counter + 1
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          namecods <- namesSeqs[counter]
          cods <- unlist(cods)
          cods <- c(namecods, cods)
          temp <- paste(cods, collapse = "\t")
          write(temp, outputFileDist, append = TRUE)
        })
      }
      else {
        stop("ERROR: outFormat should be 'mat' or 'txt' ")
      }
    }



    ##########GC_Content#############
    GC.content <- function(fasta_file){
      x <- read.fasta(file=fasta_file)
      tt<-function(x){
        res<-GC(x)
        val=round(res,4)
        return(val)
      }

      f_res<-lapply(x,tt)
      s=data.frame(f_res)

      rownames(s) <- c("GC-content")

      w=t(s)
      return(w)
    }
    ###########ONF##########
    oligo.freq <- function(fasta_file,f){
      x<- readDNAStringSet(fasta_file)
      y <- oligonucleotideFrequency(x,width = f)
      z <- data.frame(y)
      rownames(z) <- names(x)

      return(z)
    }

    #########AMIP###########

    AMIP<-function(fasta_file,n1=1,n2=4){
      x=readDNAStringSet(fasta_file)
      #calculating frequency of occurence of nucleotides k bases apart
      AMI_fun<-function(x){
        y=oligonucleotideFrequency(x, width=1, step=1)
        z<-matrix(nrow=(n2-n1)+1,ncol=4)
        F<-list()
        length(F)<-(n2-n1)+1
        R<-numeric((n2-n1)+1)
        for (i in 1:((n2-n1)+1)){
          z[i,]=oligonucleotideFrequency(x, width=1, step=i+n1-1)
          F[[i]]=rbind(y,z[i,])
          R[i]=mi.plugin(F[[i]])
        }
        R=round(R,4)
        mean_AMI<-round(mean(R),4)
        AMI<-list( mean_AMI)
        return(AMI)
      }
      res<-lapply(x,AMI_fun)
      ress= data.frame(res)
      ress = t(ress)
      row.names(ress) <- names(x)
      colnames(ress) <- c("Mean Of AMIP")
      return(ress)
    }
    ############ train data ###########
    species=Species
    if (species=="Arabidopsis"){
      fasta_file<-system.file("exdata/arabidopsis.fasta", package = "EpiSemble")
    }else if  (species=="Rice") {
      fasta_file<-system.file("exdata/rice.fasta", package = "EpiSemble")
    } else {
      message("The species is invalid")
    }

    res<-FastaToTabular(fasta_file)
    data<-as.vector(res[,2])
    mat<-as.matrix(ncp_dna(seqs = data,binaryType="strBin",outFormat="mat"))
    sequence<-rownames(mat)
    seq_id<-res[,1]
    ncp<-cbind(seq_id,sequence,mat)
    rownames(ncp)<-seq_id
    ncp_temp<-data.frame(ncp[,-1], stringsAsFactors = FALSE)
    ncp_final<-as.data.frame(apply(ncp_temp[,-1], 2, as.numeric))
    log_gc_temp<-log((GC.content(fasta_file))*100, base = exp(1))
    log_gc<-as.data.frame(as.numeric((ifelse(log_gc_temp>0,log_gc_temp,'0'))))
    onf<-oligo.freq(fasta_file, 2)

    AMIP_final<- NULL
    for (i in 1:6) {
      for(j in 2:7){
        if(i<j){
          AMIP1<-AMIP(fasta_file,i,j)

          colnames(AMIP1) <- paste('Mean',i,j, sep="_")
          AMIP_final<-cbind(AMIP_final, AMIP1)
        }
      }
    }

    temp1<- cbind(onf, gcc =log_gc[,1], ncp_final, AMIP_final)
    binary1<- rownames(temp1)
    temp2<-cbind(temp1, binary1)
    temp3<-cSplit(temp2, 'binary1', sep="_", type.convert=FALSE)
    binary<-ifelse(temp3$binary1_1=="negative",1,2)
    my_data_temp<-cbind(binary, temp1)
    inputData <-as.data.frame(my_data_temp)

    ###########Feature_Selection############
    cf1 <- cforest(binary ~ . , data= inputData, controls=cforest_unbiased(mtry=2,ntree=50))
    rf_varimp<-varimp(cf1)# get variable importance, based on mean decrease in accuracy
    temp <- as.matrix(rf_varimp)
    final <- as.matrix(temp[order(temp[, 1], decreasing = TRUE),])
    rf_varimp_40<-as.matrix(final[1:40,])
    rf_signif<-rownames(rf_varimp_40)
    colname_rf<-append("binary", rf_signif)
    rf_dataset<-inputData[,colname_rf]

    ########SwR#########
    base.mod <- lm(binary ~ 1 , data= inputData)  # base intercept only model
    all.mod <- lm(binary ~ . , data= inputData) # full model with all predictors
    stepMod <- step(base.mod, scope = list(lower = base.mod, upper = all.mod), direction = "both", trace = 0, steps = 1000)  # perform step-wise algorithm
    shortlistedVars <- names(unlist(stepMod[[1]])) # get the shortlisted variable.
    Swr_signif <- shortlistedVars[!shortlistedVars %in% "(Intercept)"]  # remove intercept
    colname_swr<-append("binary", Swr_signif)
    swr_dataset<-inputData[,colname_swr]

    common_columns = intersect(colname_rf, colname_swr)
    data_temp <- inputData[,common_columns]
    train_data_input<-cbind(Sequence=ncp_temp$sequence, data_temp)
    ######### test data #########
    fasta_file_test<-Fastafile

    res_test<-FastaToTabular(fasta_file_test)
    data_test<-as.vector(res_test[,2])
    mat_test<-as.matrix(ncp_dna(seqs = data_test,binaryType="strBin",outFormat="mat"))
    sequence<-rownames(mat_test)
    seq_id<-res_test[,1]
    ncp<-cbind(seq_id,sequence,mat_test)
    rownames(ncp)<-seq_id
    ncp_temp_test<-data.frame(ncp[,-1], stringsAsFactors = FALSE)
    ncp_final_test<-as.data.frame(apply(ncp_temp_test[,-1], 2, as.numeric))
    ncp_final_test<-cbind(sequence=ncp_temp_test$sequence, ncp_final_test)
    log_gc_temp<-log((GC.content(fasta_file_test))*100, base = exp(1))
    log_gc_test<-as.data.frame(as.numeric((ifelse(log_gc_temp>0,log_gc_temp,'0'))))
    onf_test<-oligo.freq(fasta_file_test, 2)

    AMIP_final_test<- NULL
    for (i in 1:6) {
      for(j in 2:7){
        if(i<j){
          AMIP1<-AMIP(fasta_file_test,i,j)

          colnames(AMIP1) <- paste('Mean',i,j, sep="_")
          AMIP_final_test<-cbind(AMIP_final_test, AMIP1)
        }
      }
    }

    temp_test<- cbind(ncp_final_test,onf_test, gcc =log_gc_test[,1], AMIP_final_test)
    test_col<- colnames(train_data_input[,-c(1:2)])
    test_data_temp<- temp_test[,test_col]
    #colnames(test_data_temp)<-paste0('IF',1:ncol(test_data_temp))
    test_data_input<-cbind(Sequence=temp_test$sequence, test_data_temp)
    Iresults<-list(Train=train_data_input,Test=test_data_input)
    return(Iresults)
  }
  impfnc<-ImpFeatures(Fastafile=FastaData)
  train_data_input<-impfnc$Train
  test_data_input<-impfnc$Test

  ############### Modelling ###########
  cl <- makeCluster(16) # use 16 workers
  registerDoParallel(cl) # register the parallel backen
  ##############Supprot Vector Machine#########
  # Fit trees in parallel and compute predictions on the test set
  predict_svmbag <- foreach(
    icount(101),
    .packages = "e1071",
    .combine = cbind
  )%dopar%{
    # bootstrap copy of training data
    train_data_svm<-train_data_input[,-c(1:2)]
    index <- sample(ncol(train_data_svm), size=(ncol(train_data_svm)-2),replace = FALSE)
    train_data_boot_svm<-train_data_svm[, index]
    train_res_svm<-cbind(binary=train_data_input$binary, train_data_boot_svm)
    test_data_svm<-test_data_input[,-1]
    test_res_svm<-test_data_svm[,index]
    # fit tree to bootstrap copy

    bagged_svm <- svm(
      as.factor(binary) ~ .,
      #control = svm.control(minsplit = 2, cp = 0.1),
      data = train_res_svm, method="class", kernel = "radial"
    )

    predict(bagged_svm, newdata =test_res_svm , type = "class")

  }
  mode<-function(x){which.max(tabulate(x))}
  x_svm<-as.matrix(predict_svmbag)
  svm_bagged<-apply(x_svm,1,mode)

  ##############Random Forest#########
  # Fit trees in parallel and compute predictions on the test set
  predict_rf <- foreach(
    icount(101),
    .packages = c("randomForest"),
    .combine = cbind
  ) %dopar% {
    # bootstrap copy of training data
    train_data_rf<-train_data_input[,-c(1:2)]
    index <- sample(ncol(train_data_rf), size=(ncol(train_data_rf)-2),replace = FALSE)
    train_data_boot_rf<-train_data_rf[, index]
    train_res_rf<-cbind(binary=train_data_input$binary, train_data_boot_rf)
    test_data_rf<-test_data_input[,-1]
    test_res_rf<-test_data_rf[,index]
    # fit tree to bootstrap copy
    bagged_rf <-  randomForest(as.factor(binary)~.,data = train_res_rf, ntree=500)

    predict(bagged_rf,newdata = test_res_rf)
  }
  x_rf<-as.matrix(predict_rf)

  rf_bagged<-apply(x_rf,1,mode)

  ##############Gradient_Boost#############
  # Fit trees in parallel and compute predictions on the test set
  predict_gb <- foreach(
    icount(101),
    .packages = c("caret","gbm"),
    .combine = cbind
  ) %dopar% {
    # bootstrap copy of training data
    train_data_gb<-train_data_input[,-c(1:2)]
    index <- sample(ncol(train_data_gb), size=(ncol(train_data_gb)-2),replace = FALSE)
    train_data_boot_gb<-train_data_gb[, index]
    train_res_gb<-cbind(binary=train_data_input$binary, train_data_boot_gb)
    test_data_gb<-test_data_input[,-1]
    test_res_gb<-test_data_gb[,index]

    # fit tree to bootstrap copy
    bagged_gb <- train(as.factor(binary)~ ., train_res_gb, method = 'gbm')

    predict(bagged_gb,newdata = test_res_gb)
  }
  x_gb<-as.matrix(predict_gb)

  gb_bagged<-apply(x_gb,1,mode)

  ##########Prediction_Table##########
  bagged_all<- cbind(svm_bagged,
                     rf_bagged, gb_bagged)

  colnames(bagged_all)<- c("SVM",
                           "RF", "Gradient Boosting")

  m<-as.matrix(bagged_all)
  bagged_ensemble<-apply(m,1,mode)

  bagged_final<- cbind(svm_bagged,
                       rf_bagged, gb_bagged, bagged_ensemble)

  colnames(bagged_final)<- c("SVM",
                             "RF", "Gradient Boosting", "Bagged Ensemble")

  temp_final1<- as.matrix(bagged_final[,ncol(bagged_final)])
  temp_final2<-(ifelse(temp_final1==2,'Methylated','Non Methylated'))
  colnames(temp_final2)<- "Status"
  MethStatus<-cbind(Sequence=test_data_input$Sequence,temp_final2)
  return(MethStatus)

}

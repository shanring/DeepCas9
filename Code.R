#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DeepCas9
# 2018
# Jiesi Luo, Southwest Medical University, China
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' encodeOntargetSeq and DeepCas9
#'
#' encodeOntargetSeq: convert spCas9 30 bp target sequence (4 bp + 20 
#' bp protospacer + PAM + 3 bp) into a one-hot matrix with 4 rows and 
#' 30 columns. The four rows corresponds to the four nucleotide A, G, 
#' T and C. The output file ("one-hot.csv") of this function can be used
#' directly for DeepCas9 function.
#'
#'
#'
#'
#' @usage encodeOntargetSeq(FileName)
#'
#' @param FileName A character name of the 30 bp target sequence fasta 
#' format file, including full path of file if it is located outside the 
#' current working directory.
#'
#' @return one-hot.csv A matrix containing encoded features. Sequence
#' features are represented as binary numbers.
#'
#' @examples
#' encodeOntargetSeq("textsequence.txt")
#' @export one-hot.csv
#'
#'
#'
#'
#' DeepCas9: predict the CRISPR sgRNA activity using the deep convolutional
#' neural network.
#' 
#'
#'
#' @usage DeepCas9(FileName)
#'
#' @param FileName A character name of the "one-hot" matrix file
#'
#' @return result_ontarget.csv A file containing four activity scores, 
#' including DeepCas9_HL60_score, DeepCas9_mEL4_score, DeepCas9_
#' 293T_score and DeepCas9_score.
#' 
#' @examples
#' DeepCas9("one-hot.csv")
#' @export DeepCas9_result.csv
#'
#' @author Jiesi Luo

###########################################################################

encodeOntargetSeq<-function(FileName) {
      
     DNAsequence<-read.fasta(FileName,seqtype=c("DNA"),as.string=TRUE,
                          forceDNAtolower=FALSE,seqonly=TRUE)

     one_hot<-diag(1,4,4)
  
     nuc_type<-c("A","T","C","G")

     test<-vector()

     for (num in 1:length(DNAsequence)) {

          mononucleotide<-substring(DNAsequence[num],1:30,1:30)

            for(nuc in 1:4) {
             
              mononucleotide<-str_replace_all(mononucleotide,nuc_type[nuc],str_c(as.character(one_hot[nuc,]),collapse=""))
               }
     test[num]<-str_split(str_c(mononucleotide,collapse =""),"")
       }
    test_array <-t(do.call(rbind,test))

    output<-rbind(DNAsequence,test_array)
    
    write.csv(t(output),file="one-hot.csv",row.names=FALSE)

}

###########################################################################


DeepCas9_scores<-function(FileName) {
  
        test<-read.csv(FileName)

        test_onehot<-t(test[,-1])
  
        test_sequences<-as.character(test[,1])

        dim(test_onehot)<-c(4,30,1,ncol(test_onehot))

        DeepCas9_HL60_score<-predict(DeepCas9_HL60_model,test_onehot)
 
        DeepCas9_mEL4_score<-predict(DeepCas9_mEL4_model,test_onehot)

        DeepCas9_293T_score<-predict(DeepCas9_293T_model,test_onehot)

        DeepCas9_score<-DeepCas9_mEL4_score*0.2+DeepCas9_293T_score*0.3+DeepCas9_HL60_score*0.5

       result_DeepCas9<-rbind(test_sequences,DeepCas9_HL60_score, DeepCas9_mEL4_score, DeepCas9_293T_score,DeepCas9_score)

       rownames(result_DeepCas9)<-c("test_sequences","DeepCas9_HL60_score", "DeepCas9_mEL4_score","DeepCas9_293T_score","DeepCas9_score")

 write.csv(t(result_DeepCas9),file="DeepCas9_result.csv")
    
 }











          


     

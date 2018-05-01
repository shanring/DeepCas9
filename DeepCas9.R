library(mxnet)
library(seqinr)
library(stringr)
DeepCas9_293T_model<-mx.model.load("DeepCas9_293T",iteration = 250)
DeepCas9_HL60_model<-mx.model.load("DeepCas9_HL60",iteration = 250)
DeepCas9_mEL4_model<-mx.model.load("DeepCas9_mEL4",iteration = 250)
source("Code.R")
encodeOntargetSeq("input_example.txt")
DeepCas9_scores("one-hot.csv")



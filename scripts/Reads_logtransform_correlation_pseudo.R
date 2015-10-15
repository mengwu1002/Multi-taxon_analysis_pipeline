Log_corr<-function(ref,rawreads,meta_file){
 Ref<-read.table(ref,sep="\t",header=T,row.names=1)
 Rawreads<-read.table(rawreads,sep="\t",header=T,row.names=1)
 meta<-read.table(meta_file,sep="\t",header=T)
 Combine<-cbind(Ref,Rawreads)
 Filter<-Combine[Combine[1]>0,]
 Pseduo<-Filter+0.1
 samples<-ncol(Pseduo)
 genes<-nrow(Pseduo)
 gene_names<-row.names(Pseduo)
 transform<-matrix(data=0,nrow=genes,dimnames=list(gene_names))
 for (x in 2:samples){
    ratio<-Pseduo[,x]/Pseduo[,1]
    log<-log(ratio,10)
    newcol<-matrix(t(log))
    colnames(newcol)<-names(Pseduo[x])
    transform<-cbind(transform,newcol)
 }
 Final<-transform[,-1]
 Filter_file<-paste(rawreads,"_filtered.txt",sep="")
 write.table(Filter,file=Filter_file,sep="\t")
 outputname<-paste(rawreads,"_logratio.txt",sep="")
 write.table(Final,file=outputname,sep="\t")

 sample_type<-levels(as.factor(meta$Group))
 type_length<-length(sample_type)
 correlation_data<-list()
 for (x in 1:type_length){
  dpg_name<-paste("Group_",sample_type[x],".txt",sep="")
  sub<-Rawreads[,which(meta$Group==sample_type[x])]
  dpg_n<-ncol(sub)
  index=0
  sum=0
  dpg_corr<-matrix(data="NA",nrow=dpg_n,ncol=dpg_n,dimnames=list(colnames(sub),colnames(sub)))
  for (i in 1:dpg_n){
    m<-i+1
    if (m<=dpg_n){
      for (j in m:dpg_n){
        Corr<-round(cor(sub[,i],sub[,j])**2,digits=3)
        sum<-sum+Corr
        dpg_corr[i,j]<-Corr
        index<-index+1;  
      }
    }
  }
#  correlation_data[[dpg_name]]<-dpg_corr
  average_Corr<-round(sum/index,digits=3)
  dpg_corr[dpg_n,1]<-average_Corr
  write.table(dpg_corr,file=dpg_name,sep="\t")
 }
 
 for (x in 1:type_length){
  dpg_name<-paste("Ref_",sample_type[x],".txt",sep="")
  dpg_pdf<-paste("Ref_",sample_type[x],".pdf",sep="")
  sub<-Rawreads[,which(meta$Group==sample_type[x])]
  dpg_n<-ncol(sub)
  pdf(dpg_pdf,width=4*dpg_n,height=4)
  par(mfrow=c(1,dpg_n))
  index=0
  sum=0
  colnames_string<-c(colnames(sub),"average")
  dpg_corr<-matrix(data="NA",nrow=1,ncol=dpg_n+1,dimnames=list("R2",colnames_string))
  for (i in 1:dpg_n){
        Corr<-round(cor(Ref[,1],sub[,i])**2,digits=3)
        sum<-sum+Corr
        dpg_corr[1,i]<-Corr
        index<-index+1;
        Corr_name<-paste("R2=",Corr,sep="")
        plot(log(Ref[,1],10),log(sub[,i],10),cex=0.3,pch=19,xlim=range(-1,6),ylim=range(-1,6),main=colnames(sub)[i])
        legend("bottomright",legend=Corr_name,bty="n")
  }
#sample_sub1<ste("Ref_",sample_type[x],".pdf",sep="")  correlation_data[[dpg_name]]<-dpg_corr

  average_Corr<-round(sum/index,digits=3)
  dpg_corr[1,dpg_n+1]<-average_Corr
   write.table(dpg_corr,file=dpg_name,sep="\t")
  dev.off()
 }

}

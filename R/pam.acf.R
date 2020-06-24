pam.acf <-
function (object) 
    {
        indiv <- unique(object$index_diffdata)
        max_col<-max(table(object$index_diffdata))
        id_acf<-matrix(NA,length(indiv),max_col-1)
        id_pacf<-matrix(NA,length(indiv),max_col-1)
        empty_index<-c()
        for (j in 1:length(indiv)) {
            length_j<-sum(object$index_diffdata == indiv[j])
            if (length_j>1)
            {
                id_acf[j,1:(length_j-1)]<-(acf(object$resid[object$index_diffdata == indiv[j]],plot=F,lag.max = length_j))$acf[-1]
                id_pacf[j,1:(length_j-1)]<-(pacf(object$resid[object$index_diffdata == indiv[j]],plot=F,lag.max = length_j))$acf
            } else
                empty_index<-c(empty_index,j)
        }
        if (!is.null(empty_index))
        {
            id_acf<-id_acf[-empty_index,]
            id_pacf<-id_pacf[-empty_index,]
            row_names<-unique(object$index_diffdata)[-empty_index]
        } else  row_names<-unique(object$index_diffdata)
        id_acf<-round(id_acf,4)
        colnames(id_acf)<-as.character(1:(dim(id_acf)[2]))
        colnames(id_acf) <- paste("AR",colnames(id_acf),sep = "_")
        rownames(id_acf)<- row_names
        rownames(id_acf) <- paste(names(object$index_data)[1],rownames(id_acf),sep = "_")
        
        id_pacf<-round(id_pacf,4)
        colnames(id_pacf)<-as.character(1:(dim(id_pacf)[2]))
        colnames(id_pacf) <- paste("MA",colnames(id_pacf),sep = "_")
        rownames(id_pacf)<-row_names
        rownames(id_pacf) <- paste(names(object$index_data)[1],rownames(id_pacf),sep = "_")
        
        if (any(is.na(id_acf))) {
            summary_acf<-suppressWarnings(do.call(rbind,apply(id_acf,2,function (x) summary(x,na.rm=T))))
            summary_acf<-t(round(summary_acf[,1:6],4))
        } else {
            summary_acf<-apply(id_acf,2,function (x) summary(x,na.rm=T))
            summary_acf<-round(summary_acf[1:6,],4)
        }
        if (any(is.na(id_pacf))) {
            summary_pacf<-suppressWarnings(do.call(rbind,apply(id_pacf,2,function (x) summary(x,na.rm=T))))
            summary_pacf<-t(round(summary_pacf[,1:6],4))
        } else {
            summary_pacf<-apply(id_pacf,2,function (x) summary(x,na.rm=T))
            summary_pacf<-round(summary_pacf[1:6,],4)
        }
        ret <- list(id_acf = id_acf, id_pacf = id_pacf, summary_acf = summary_acf,
                    summary_pacf = summary_pacf)
        ret
    }

#Omnivory_Petchey

Fraction.omnivores <- function(web) {
    TLs <- TLs(web)[,1]
    non.int.TL <- web[,TLs %% 1 != 0]
    if(is.matrix(non.int.TL))
        frac.omniv <- sum(apply(non.int.TL, 2, sum) > 1)  / length(web[,1])
    if(is.vector(non.int.TL)) 
        frac.omniv <- (sum(non.int.TL) > 1)  / length(web[,1])
    frac.omniv
}


Level.omnivory <- function(web) {

    TLs <- TLs(web)[,1]

    if( sum(is.na(TLs)) == length(TLs) )
       rr <- NA

    if( sum(is.na(TLs)) != length(TLs) ) {

        
        web.TLs <- matrix(rep(TLs, length(web[,1])), length(web[,1]), length(web[,1]))
        lo.pc <- numeric(length=length(web[,1]))
        for(i in 1:length(web[,1])) {
            tt <- web.TLs[web[,i]==1,i]
            if(length(tt)==0 | sum(!is.na(tt))==0 )
                lo.pc[i] = NA
            if(length(tt)>0 & sum(!is.na(tt))!=0)
                lo.pc[i] <- sd(tt)
        }
        rr <- mean(lo.pc, na.rm=T)
    }
    rr
}


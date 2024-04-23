#' Title
#'
#' @param X 
#' @param vecd.or.p 
#' @param N 
#' @param B 
#' @param alpha 
#' @param display 
#' @param graphics 
#' @param nbclus 
#'
#' @return
#' @export
#'
#' @examples
dependogram <- function(X,vecd.or.p,N=10,B=2000,alpha=0.05,display=TRUE,graphics=TRUE,nbclus=1) {


  if (nbclus>1) {

#    suppressWarnings(parallel.pkg.present <- require(parallel))
      parallel.pkg.present <- "package:parallel" %in% search()
      Rmpi.pkg.present <- "package:Rmpi" %in% search()
    if (all(!c(parallel.pkg.present,Rmpi.pkg.present))) stop("Either package parallel or Rmpi should be installed!")
#    suppressWarnings(rsprng.pkg.present <- require(rsprng))
#    if (!rsprng.pkg.present) stop("Package rsprng is not installed!")
    cluster.type <- if (parallel.pkg.present) "PSOCK" else "MPI" # We prefer to use "PSOCK" (i.e. parallel) because it's easier.
    
  }


  X <- as.matrix(X)

# if length(vecd.or.p)>1 then serial case, otherwise non serial case

  if (length(vecd.or.p) > 1) {
# we perform the serial case
    seriel <- 0
    vecd <- vecd.or.p
    p <- length(vecd)
    taille <- 2^p-p-1
#


    
    RnAs <- rep(0,2^p-p-1)
    Rn <- 0

# We start the cluster
    if (nbclus > 1) {

      B <- round(B/nbclus)*nbclus
      

      cl <- parallel::makeCluster(nbclus, type = cluster.type) 
                                        
      
      myfunc <- function(B,p) {
        RnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
        Rnstar <- rep(0,B)
        require(IndependenceTests)
# We call the C function dependogram
        .C("dependogramC",
           as.integer(N),
           as.integer(vecd),
           as.integer(length(vecd)),
           as.integer(p),
           as.numeric(X),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           B=as.integer(B),
           as.numeric(alpha),
           RnAs=as.numeric(RnAs),
           RnAsstar=as.numeric(RnAsstar),
           Rn=as.numeric(Rn),
           Rnstar=as.numeric(Rnstar),
           as.integer(seriel),PACKAGE="IndependenceTests")
        
      }
      
      out2 <- parallel::clusterCall(cl, myfunc, round(B/nbclus),p)
      
                   
      
      # We stop the cluster
      parallel::stopCluster(cl)

      out <- list(RnAs=c(),RnAsstar=c(),Rn=c(),Rnstar=c())


      for (clus in 1:nbclus) {
        
        out$Rnstar <- c(out$Rnstar,out2[[clus]]$Rnstar)
        out$RnAsstar <- c(out$RnAsstar,out2[[clus]]$RnAsstar)
        
      }

      out$Rn <- out2[[clus]]$Rn
      out$RnAs <- out2[[clus]]$RnAs
            
      
    } else {

      RnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
      Rnstar <- rep(0,B)

# We call the C function dependogram
      out <- .C("dependogramC",
		as.integer(N),
		as.integer(vecd),
		as.integer(length(vecd)),
		as.integer(p),
		as.numeric(X),
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		as.integer(B),
		as.numeric(alpha),
		RnAs=as.numeric(RnAs),
		RnAsstar=as.numeric(RnAsstar),
		Rn=as.numeric(Rn),
		Rnstar=as.numeric(Rnstar),
		as.integer(seriel),PACKAGE="IndependenceTests")
      
    }






    if (display) {
# require(combinat)
      RES <- as.list(1:(p-1))
      for (cardA in 2:p) {RES[[cardA]] <- as.matrix(combn(p,cardA))}
      nb <- 0
      for (cardA in 2:p) {
        for (j in 1:(choose(p,cardA))) {
          nb <- nb+1
          
          cat(c(nb,": A=",RES[[cardA]][,j],": ||RnA||=",round(out$RnAs[nb],3),"\n"))
          
          
        }
      }
    }
    

# Sort the elements on each line
    RnAsstar <- matrix(out$RnAsstar,nrow=(2^p-p-1),ncol=B,byrow=FALSE)
    matseuils <- t(apply(RnAsstar,FUN=sort,MARGIN=1))


    beta <- (1-alpha)^(1/taille)
# Contains the beta-quantiles of the stats RnA for each one of the 2^p-p-1 set A
    AllThresholds <- matseuils[,round(beta*B)]

    if (graphics) {
# We draw a vertical bar for each A with height ||RnA||
      plot(out$RnAs,type="h",ylim=c(0,max(c(max(AllThresholds),max(out$RnAs),max(out$Rnstar)))+0.1),xlim=c(0,2^p-p),main="Dependogram",xlab="Subsets",ylab="||RnA||")
    }

    Rnstar <- sort(out$Rnstar)
    GlobalThreshold <- Rnstar[round((1-alpha)*B)]
# abline(h=GlobalThreshold,col="red")


    if (graphics) { 
# We put a star for each beta-quantile of ||R_A||
      points((1:(2^p-p-1)),AllThresholds,pch="*")
    }

    res <- list(norm.RnA=out$RnAs,Rn=out$Rn,rA=AllThresholds,r=GlobalThreshold,RnAsstar=RnAsstar)
    return(res)


  }



  if (length(vecd.or.p) == 1) {
# We perform the serial case
    seriel <- 1
    p <- vecd.or.p
    vecd <- rep(ncol(X),p)
    taille <- 2^(p-1)-1




    SnAs <- rep(0,taille)
    Sn <- 0


# We start the cluster
    if (nbclus > 1) {
      
      B <- round(B/nbclus)*nbclus

      SnAsstar <- matrix(0,nrow=taille,ncol=B)
      Snstar <- rep(0,B)
      
      cl <- parallel::makeCluster(nbclus, type = cluster.type) 
                                        
      
      myfunc <- function(B,p) {
        SnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
        Snstar <- rep(0,B)
        require(IndependenceTests)
# We call the C function dependogram
        .C("dependogramC",
           as.integer(N),
           as.integer(vecd),
           as.integer(length(vecd)),
           as.integer(p),
           as.numeric(X),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           B=as.integer(B),
           as.numeric(alpha),
           SnAs=as.numeric(SnAs),
           SnAsstar=as.numeric(SnAsstar),
           Sn=as.numeric(Sn),
           Snstar=as.numeric(Snstar),
           as.integer(seriel),PACKAGE="IndependenceTests")
      }
      
      out2 <- parallel::clusterCall(cl, myfunc, round(B/nbclus),p)
      
                                        
      
      # We stop the cluster
      parallel::stopCluster(cl)

      out <- list(SnAs=c(),SnAsstar=c(),Sn=c(),Snstar=c())
      

      for (clus in 1:nbclus) {
        
        out$Snstar <- c(out$Snstar,out2[[clus]]$Snstar)
        out$SnAsstar <- c(out$SnAsstar,out2[[clus]]$SnAsstar)
        
      }

      out$Sn <- out2[[clus]]$Sn
      out$SnAs <- out2[[clus]]$SnAs

      
    } else {

      SnAsstar <- matrix(0,nrow=taille,ncol=B)
      Snstar <- rep(0,B)
# We call the C function dependogram
      out <- .C("dependogramC",
		as.integer(N),
		as.integer(vecd),
		as.integer(length(vecd)),
		as.integer(p),
		as.numeric(X),
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		as.integer(B),
		as.numeric(alpha),
		SnAs=as.numeric(SnAs),
		SnAsstar=as.numeric(SnAsstar),
		Sn=as.numeric(Sn),
		Snstar=as.numeric(Snstar),
		as.integer(seriel),PACKAGE="IndependenceTests")
      
    }







    if (display) {

# require(combinat)
      RES <- as.list(1:(p-1))
      for (cardA in 2:p) {RES[[cardA]] <- as.matrix(rbind(rep(1,choose(p-1,cardA-1)),as.matrix(combn(p-1,cardA-1)+1)))}
      nb <- 0
      for (cardA in 2:p) {
        for (j in 1:(choose(p-1,cardA-1))) {
          nb <- nb+1

          cat(c(nb, ": A=",RES[[cardA]][,j],": ||SnA||=",round(out$SnAs[nb],3),"\n"))

        }
      }
    }


    SnAsstar <- matrix(out$SnAsstar,nrow=taille,ncol=B,byrow=FALSE)
    Sn <- max(out$SnAs)


    beta <- (1-alpha)^(1/taille)


# The beta-quantile of S_n,A is computed by grouping all the values 
# S_n,A^* with |A|=k as in the article
# There are choose(p-1,|A|-1) sets A of size |A| containing 1
# We therefore must take in the matrix SnAsstar groups of choose(p-1,|A|-1) rows, for |A|=2 to p.
# There will be p-1 groups.
# For each one of these groups, we create a vector vecA by taking all the elements in the group, 
# then we compute AllThresholds[|A|]<-vecA[round(beta*B*choose(p-1,|A|-1))] for |A|=2 to p
# This vector AllThresholds will thus contain the heights of the horizontal bars to put on the dependogram (one horizontal bar 
# for each |A|)


    AllThresholds <- rep(0,p-1)
    begin <- 1
    end <- 0
    for (cardA in 2:p) {
      end <- end+choose(p-1,cardA-1)
      vecA <- as.vector(SnAsstar[begin:end,])
      vecA <- sort(vecA)
      AllThresholds[cardA-1] <- vecA[round(beta*B*choose(p-1,cardA-1))]
      begin <- end+1
    }

    if (graphics) {
# We draw a vertical bar for each A of height ||SnA||
      plot(out$SnAs,type="h",ylim=c(0,max(c(max(AllThresholds),max(out$SnAs),max(out$Snstar)))+0.1),xlim=c(0,2^(p-1)),main="Dependogram",xlab="Subsets",ylab="||SnA||")
    }
    
    Snstar <- sort(out$Snstar)
    GlobalThreshold <- Snstar[round((1-alpha)*B)]
# abline(h=GlobalThreshold,col="red")


# It remains to put the horizontal lines of the threshold values for each |A|


    if (graphics) {
      begin <- 1
      end <- 0
      for (cardA in 2:p) {
        end <- end+choose(p-1,cardA-1)

        segments(begin-0.5,AllThresholds[cardA-1],end+0.5,AllThresholds[cardA-1],lty=4)

        begin <- end+1
      }
    }





    res <- list(norm.SnA=out$SnAs,Sn=out$Sn,sA=AllThresholds,s=GlobalThreshold,SnAsstar=SnAsstar)
    return(res)
  }
  


}




library(shiny)
library(TMB)
library(tmbstan)

loadModel <- function( name )
{
  if(name %in% names(getLoadedDLLs()))
    dyn.unload(dynlib(name))          # unlink the C++ code if already linked
  compile(paste0(name,".cpp"),flags = "")
  dyn.load(dynlib(name))          # Dynamically link the C++ code
}

function(input, output, session) {


  fitMod <- reactive({

    loadModel("swnsSeal")

    years <- seq(1980,input$finalYear,by=1)
    raw <- c(0.204,0.417,1.849,2.246)
    counts <- c(raw[1]*1.5,raw[2]*1.5,2.105,2.42)
    time <- c(2007,2010,2016,2021)
    CVs <- c(0.5,0.5,0.07,0.08)
    SDs <- sqrt( log(CVs^2+1) )
  
    nT <- length(years)
  
    obsI_t <- rep(NA,nT)
    obsI_t[ years %in% time ] <- counts
  
    sdI_t <- numeric(nT)
    sdI_t[ years %in% time ] <- SDs
  
    rawI_t <- rep(NA,nT)
    rawI_t[ years %in% time ] <- raw
    rawI_t[ rawI_t==0 ] <- NA

    data <- list( obsI_t=obsI_t,
                  sdI_t=sdI_t,
                  rprior=c(log(input$rMu),input$rSD),
                  Kprior=c(log(input$KMu),input$KSD) )
  
    pars <- list( I0=input$I0*1e-3, r=0.1, K=40 )
  
    map <- list( I0=factor(NA) )
  
    obj <- MakeADFun( data       = data,
                      parameters = pars,
                      map        = map,
                      DLL        = "swnsSeal" )
  
    lower <- rep(0,3)
    upper <- c(Inf,0.4,Inf)
  
    opt <- nlminb( start     = obj$par,
                   objective = obj$fn,
                   gradient  = obj$gr,
                   #lower     = lower,
                   #upper     = upper,
                   control   = list("trace"=1, "eval.max"=5e4, "iter.max"=5e4) )  
  
    rpt <- obj$report()
    rpt$opt <- opt
    rpt$years <- years
    rpt$data  <- data
  
    fit <- data.frame()

    while( nrow(as.data.frame(fit))==0 )
    {
      mcinit <- list()
        for( i in 1:input$nChain )
          mcinit[[i]] <- rnorm(n=length(opt$par),mean=opt$par,sd=0.1)
    
      fit <- try( tmbstan( obj = obj,
                      chains = input$nChain,
                      iter = input$nIter,
                      init = mcinit,
                      control = list( adapt_delta=input$adapt_delta,
                                      max_treedepth=input$max_treedepth ) ) )
    }
    
      

    df <- as.data.frame(fit)
    df <- df[ ,-ncol(df)]
  
  
    # Get posterior I
  
    I_it <- matrix( data=NA, nrow=nrow(df), ncol=length(years) )
    for( i in 1:nrow(df) )
      I_it[i, ] <- obj$report(df[i, ])$I_t
  
  
    list( fit=fit,
          I_it=I_it,
          years=years,
          obsI_t=obsI_t,
          rawI_t=rawI_t )


  })


  output$timeSeries <- renderPlot({


    #layout( rbind( c(1,1), c(2,2), c(3,4) ) ) )


    par( mfrow=c(2,1), mar=c(3,6,0,0), oma=c(2,0,1,1) )

    out <- fitMod()

    fit    <- out$fit
    I_it   <- out$I_it
    years  <- out$years
    obsI_t <- out$obsI_t
    rawI_t <- out$rawI_t

  Ilo_t <- apply( I_it, 2, quantile, 0.002 )
  Imo_t <- apply( I_it, 2, median )
  Ihi_t <- apply( I_it, 2, quantile, 0.998 )


  x <- years
  plot( x    = range(years),
        y    = range(obsI_t),
        ylim = c(0,max( obsI_t, Ihi_t, na.rm=1 )),
        axes = FALSE,
        xlab = "",
        ylab = "" )
  mtext( side = 2,
         text = "SWNS pup production ('000)",
         line = 4,
         cex  = 1.5 )
  grid()
  box()
  axis( side=1, cex.axis=1.5 )
  axis( side=2, las=1, cex.axis=1.5 )
  polygon( x=c(x,rev(x)), y=c(Ilo_t,rev(Ihi_t)),
          col=rgb(0.8,0.8,0.8,0.5), border=NA )
  lines( x, Imo_t, lwd=1.5 )
  points( years, rawI_t, pch=16, cex=1.3 )
  points( years, obsI_t, cex=1.3 )
  abline( v=2021.5, lty=2 )
  par(font=2)
  #legend( x="topleft", bty="n", legend="(a)", cex=1.2)
  par(font=1)
  legend( x="topleft", bty="n", pch=c(NA,1,16), lty=c(1,0,0), cex=1.5,
          legend=c("Logistic model","Adjusted count","Raw count") )


  N_it <- I_it / input$pupAdultRatio

  N_qt <- apply( N_it, 2, quantile, c(0.025,0.975) )
  N_t <- apply( N_it, 2, median )
  plot( x    = range(years),
        y    = range(N_qt),
        type = "n",
        las  = 1,
        ylab = "",
        axes = FALSE )
  mtext( side = 2,
         text = "Total SWNS abundance ('000)",
         line = 4,
         cex  = 1.5 )
  grid()
  box()
  axis( side=1, cex.axis=1.5 )
  axis( side=2, las=1, cex.axis=1.5 )
  axis( side=3, labels=NA )
  mtext( side=4, text="Pup N : Total N", cex=1.5, line=3.3 )
  polygon( x=c(years,rev(years)), y=c(N_qt[1, ],rev(N_qt[2, ])),
           col=rgb(0.8,0.8,0.8,0.5), border=NA )  
  lines( x=years, y=N_t, lwd=1.2 )
  abline( v=2021.5, lty=2 )


  }, width=650, height=800)


  output$priorPost <- renderPlot({


    par( mfrow=c(1,2), mar=c(5,3,0,0), oma=c(0,3,1,1) )

    out <- fitMod()

    fit    <- out$fit
    I_it   <- out$I_it
    years  <- out$years
    obsI_t <- out$obsI_t
    rawI_t <- out$rawI_t

  # Posterior plots

  df <- as.data.frame(fit)

  # r
  rprior <- rnorm(1e5,log(input$rMu),input$rSD)
  plot( x=density( df$r ),
        lwd=2,
        las=1,
        xlab="",
        ylab="",
        main="",
        xlim=c(0,1),
        xaxs="i" )
  par(font=2)
  mtext( side = 1,
         text = "r",
         line = 4,
         cex  = 1.5 )
  par(font=1)
  mtext( side = 2,
         text = "Density",
         line = 4,
         cex  = 1.5 )   
  grid()
  box()
  lines( density( exp(rprior) ), lwd=2, lty=3, col="red" )
  lines( density( df$r ), lwd=2 )
  legend( x      = "topright",
          bty    = "n",
          legend = c("Prior","Posterior"),
          lwd    = 2,
          col    = c("red","black"),
          lty    = c(3,1),
          cex    = 1.5 )

  # K
  Kprior <- rnorm(1e5,log(input$KMu),input$KSD)
  plot( x=density( df$K ),
        lwd=2,
        las=1,
        xlab="",
        ylab="",
        main="",
        xlim=c(0,20),
        xaxs="i" )
  par(font=2)
  mtext( side = 1,
         text = "K",
         line = 4,
         cex  = 1.5 )
  par(font=1)     
  grid()
  box()
  lines( density( exp(Kprior) ), lwd=2, lty=3, col="red" )
  lines( density( df$K ), lwd=2 )
  legend( x      = "topright",
          bty    = "n",
          legend = c("Prior","Posterior"),
          lwd    = 2,
          col    = c("red","black"),
          lty    = c(3,1),
          cex    = 1.5 )

  }, width=700, height=350)

  output$diagnos <- renderTable({
    
    out <- fitMod()

    fit <- out$fit

    d <- as.data.frame(fit)

    diags <- c("nSamps",
               "r_Rhat",
               "r_BulkESS",
               "r_TailESS",
               "K_Rhat",
               "K_BulkESS",
               "K_TailESS",
               "nDiverg",
               "nSMT")

    for( i in 1:input$nChain )
      diags <- c( diags, paste0("EBFMI_chain",i) )

    df <- data.frame( Diagnostic = diags,
                      Value      = rep(NA,length(diags)),
                      Target     = rep(NA,length(diags)) )

    df[diags=="nSamps","Value"]    <- nrow(d)
    df[diags=="r_Rhat","Value"]    <- format(round(summary(fit)[[1]]["r","Rhat"],3),nsmall=3)
    df[diags=="r_BulkESS","Value"] <- round(ess_bulk( d$r ))
    df[diags=="r_TailESS","Value"] <- round(ess_tail( d$r ))
    df[diags=="K_Rhat","Value"]    <- format(round(summary(fit)[[1]]["K","Rhat"],3),nsmall=3)
    df[diags=="K_BulkESS","Value"] <- round(ess_bulk( d$K ))
    df[diags=="K_TailESS","Value"] <- round(ess_tail( d$K ))
    df[diags=="nDiverg","Value"]   <- sum(get_divergent_iterations(fit))
    df[diags=="nSMT","Value"]      <- sum(get_max_treedepth_iterations(fit))
    
    ebfmi <- get_bfmi(fit)
    
    for( i in 1:input$nChain )
      df[diags==paste0("EBFMI_chain",i),"Value"] <- format(round(ebfmi[i],3),nsmall=3)


    df[diags=="nSamps","Target"]    <- "-"
    df[diags=="r_Rhat","Target"]    <- "< 1.010"
    df[diags=="r_BulkESS","Target"] <- paste(">",input$nChain*100)
    df[diags=="r_TailESS","Target"] <- paste(">",input$nChain*100)
    df[diags=="K_Rhat","Target"]    <- "< 1.010"
    df[diags=="K_BulkESS","Target"] <- paste(">",input$nChain*100)
    df[diags=="K_TailESS","Target"] <- paste(">",input$nChain*100)
    df[diags=="nDiverg","Target"]   <- "0"
    df[diags=="nSMT","Target"]      <- "0"

    for( i in 1:input$nChain )
      df[diags==paste0("EBFMI_chain",i),"Target"] <- "> 0.3"

    df

  })



}



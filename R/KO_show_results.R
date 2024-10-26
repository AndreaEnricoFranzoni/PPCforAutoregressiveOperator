#' @import ggplot2
#' @import gridExtra
#' @importFrom grid textGrob
#' @importFrom grid gpar

KO_show_results <-function(results_ko){
  
  if (length(results_ko) == 0) {
    stop("First parameter (Result) must not be NULL")
  } else{
    cat("\n")
    cat("Results of KO algorithm")
    cat("\n")
  }
  
  
  #printing parameters found
  cat("\n")
  cat("Alpha:", results_ko$Alpha, "\n")
  cat("\n")
  cat("Number of PPCs retained:", results_ko$`Number of PPCs retained`, "\n")
  cat("\n")
  
  exp_pw = data.frame(
    m = results_ko$`Explanatory power PPCs`
  )
  colnames(exp_pw) = "Explanatory power"
  print(exp_pw)
  cat("\n")
  
  
  
  #One step ahead prediction plot
  quartz()
  plot(x = results_ko$`Function discrete evaluations points`,
       y = results_ko$`One-step ahead predictio`,
       type="l", 
       col='black',
       lwd = 2,
       xlab = "Units domain",
       ylab = "Value",
       main = "1-step ahead prediction",
       ylim = c(min(s$`One-step ahead predictio`,s$f_n),max(s$`One-step ahead predictio`,s$f_n)))
  points(x = results_ko$`Function discrete evaluations points`,
         y = results_ko$f_n,
         type="l",
         col = "blue",
         lwd = 2
  )
  legend("topleft",                    # Posizione: "topleft" (in alto a sinistra)
         legend = c(expression(f[n]), expression(f[n+1])),  # Etichette
         col = c("blue", "black"),
         lty = 1,
         lwd = 2)            # Titolo opzionale
  
  
  #CV errors to be plotted
  if(length(results_ko) == 15){
    
    #CV alpha
    if(results_ko$CV=="CV_alpha"){
      quartz()
      # Plotting them equispaced in order to avoid problems
      x_ax <- seq(1, length(results_ko$`Validation errors`))
      
      plot(x = x_ax, 
           results_ko$`Validation errors`, 
           type = "o", 
           col = "blue", 
           pch = 16,
           xlab = "Alpha", 
           ylab = "Validation error", 
           main = "Validation error for cv on regularization parameter", 
           xaxt = "n")
      
      # Names of the real values of alpha
      axis(1, at = x_ax, labels = as.character(results_ko$Alphas),las = 2)
      grid()
    }
    
    if(results_ko$CV=="CV_k"){
      quartz()
      plot(x = 1:(length(results_ko$`Validation errors`)), 
           y = results_ko$`Validation errors`, 
           type = "o", 
           col = "blue", 
           pch = 16,
           xlab = "k", 
           ylab = "Validation error", 
           main = "Validation error for cv on the number of PPCs") 
    }
    
    if(results_ko$CV=="CV"){
      
    }
    
    
  }
  
  #plot of directions and weights
  for (i in 1:results_ko$`Number of PPCs retained`) {
    
    quartz()
    title_ppc = paste("Principal Predictive Component ",i)
    # Crea i grafici con ggplot2
    grafico1 <- ggplot(data = data.frame(results_ko$`Function discrete evaluations points`, results_ko$`Directions of PPCs`[,i]), aes(x = results_ko$`Function discrete evaluations points`, y = results_ko$`Directions of PPCs`[,i])) +
      geom_line() +
      #ggtitle("Direction") +
      labs(subtitle = "Direction", x = "Units domain", y = "y")
    
    grafico2 <- ggplot(data = data.frame(results_ko$`Function discrete evaluations points`, results_ko$`Weights of PPCs`[,i]), aes(x = results_ko$`Function discrete evaluations points`, y = results_ko$`Weights of PPCs`[,i])) +
      geom_line() +
      #ggtitle("Weight") +
      labs(subtitle = "Weight", x = "Units domain", y = "y")
    
    # Unisci i grafici con un titolo generale
    grid.arrange(grafico1, grafico2, ncol = 2, top = textGrob(title_ppc, gp = gpar(fontsize = 20)))
  }
}
#' @import ggplot2
#' @import gridExtra
#' @importFrom grid textGrob
#' @importFrom grid gpar

KO_show_results <- function(results_ko,hp_ko=NULL,true_alphas=FALSE){
  
  lenght_res_with_er = 15
  
  #checking correct input
  if (length(results_ko) == 0) {
    stop("First parameter (results_ko) must not be NULL")
  } else{
    if(length(results_ko)!=lenght_res_with_er && length(results_ko)!=(lenght_res_with_er-1)){
      stop("results_ko could not come from KO")
    }
    cat("\n")
    cat("Results of KO algorithm")
    cat("\n")
  }
  
  
  
  #Print parameters found: alpha, number of PPCs and Explanatory power
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
  
  
  
  #Plot pvalues of ADF
  if(!(is.null(hp_ko))){
    
    if(length(hp_ko)!=1){ stop("Second parameter (hp_ko) could be not from hp-check KO")}
    
    data_hp <- data.frame(x = results_ko$`Function discrete evaluations points`,y = hp_ko$`Pvalues ADF`)
    
    plot_pvalues <- ggplot(data_hp, aes(x = x, y = y, fill = y)) +
      geom_bar(stat = "identity", color = "black") +  # 'stat = "identity"' per utilizzare i valori di y
      scale_fill_gradient2(low = "blue", mid="white", high = "red", midpoint = 0.1) +
      labs(title = "Pointwise ADF test", x = "x", y = "P-values") +
      geom_hline(yintercept = c(0.01,0.05,0.1), color = "red") + 
      theme_minimal() + 
      theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5),
            panel.grid.major.x = element_line( color = "gray", linewidth = 0.1))
    
    quartz()
    print(plot_pvalues)}
  
  
  
  #Plot ine step ahead prediction plot
  plot_pred <- ggplot() +
    geom_line(data = data.frame(x=results_ko$`Function discrete evaluations points`,y=results_ko$`One-step ahead prediction`), aes(x = x, y = y, color = "Current")) +
    geom_line(data = data.frame(x=results_ko$`Function discrete evaluations points`,y=results_ko$f_n),                         aes(x = x, y = y, color = "Next")) +
    labs(title = "1-step ahead prediction", x = "x", y = "y") +
    scale_color_manual( name = "Time instants", values = c("Current" = "black", "Next" = "blue")) +
    theme_minimal() +
    theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5),
          panel.grid.major = element_line(color = "gray", linewidth = 0.5),
          panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25))
  
  quartz()
  print(plot_pred)
  
  
  
  #CV errors to be plotted
  if(length(results_ko) == lenght_res_with_er){
    
    #CV alpha
    if(results_ko$CV=="CV_alpha"){
      
      
      if(true_alphas==FALSE){  #put the alphas equispaced
        opt_par       <- data.frame(x=which(results_ko$Alphas==results_ko$Alpha),y=min(results_ko$`Validation errors`))
        plot_cv_alpha <- ggplot( data = data.frame(x = 1:length(results_ko$`Validation errors`), y = results_ko$`Validation errors`), aes(x = x, y = y)) +
          geom_point( color = "blue", size = 1) +  
          geom_line( color = "blue") +  
          geom_point(data = opt_par, aes(x = x, y = y), color = "red", size = 4, shape = 8) +
          labs(title = "Validation error for CV on regularization parameter", x = "alpha", y = "Validation error") +
          theme_minimal() +
          scale_x_continuous( breaks = seq(1,length(results_ko$`Validation errors`), by = 1), labels = as.character(results_ko$Alphas)) + 
          theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5),
                panel.grid.major.x = element_line( color = "gray", linewidth = 0.5))}
      
      if(true_alphas==TRUE){  #use the true alphas 
        opt_par       <- data.frame(x=results_ko$Alpha,y=min(results_ko$`Validation errors`))
        plot_cv_alpha <- ggplot( data = data.frame(x = results_ko$Alphas, y = results_ko$`Validation errors`), aes(x = x, y = y)) +
          geom_point( color = "blue", size = 1) +  
          geom_line( color = "blue") +
          geom_point(data = opt_par, aes(x = x, y = y), color = "red", size = 4, shape = 8) +
          labs(title = "Validation error for CV on regularization parameter", x = "alpha", y = "Validation error") +
          theme_minimal() +
          scale_x_continuous( breaks = results_ko$Alphas) + 
          theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5),
                panel.grid.major.x = element_line( color = "gray", linewidth = 0.5))}
      
      quartz()
      print(plot_cv_alpha)}
    
    
    #CV k
    if(results_ko$CV=="CV_k"){
      opt_par   <- data.frame(x=results_ko$`Number of PPCs retained`,y=min(results_ko$`Validation errors`))
      plot_cv_k <- ggplot(data=data.frame(x = 1:length(results_ko$`Validation errors`), y = results_ko$`Validation errors`), aes(x = x, y = y)) +
        geom_point(color = "blue", size = 1) +  
        geom_line(color = "blue") +  
        geom_point(data = opt_par, aes(x = x, y = y), color = "red", size = 4, shape = 8) +
        labs(title = "Validation error for CV on number of PPCs", x = "k", y = "Validation error") +
        theme_minimal() +
        scale_x_continuous(breaks = seq(1,length(results_ko$`Validation errors`), by = 1)) + 
        theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5),
              panel.grid.major.x = element_line( color = "gray", linewidth = 0.5)) 
      
      quartz()
      print(plot_cv_k)}
    
    
    #CV alpha-k
    if(results_ko$CV=="CV"){
      
      #preparing the data for the heatmap
      data <- expand.grid(x = results_ko$Alphas, y = results_ko$K_s)
      tot_alpha = length(results_ko$Alphas)
      tot_k = length(results_ko$K_s)
      for(i in 1:tot_alpha){
        err_alpha_fix = unlist(results_ko$`Validation errors`[i]) 
        if(length(err_alpha_fix)<tot_k){  err_alpha_fix = c(err_alpha_fix,rep(max(err_alpha_fix),tot_k-length(err_alpha_fix)))}
        data[which(data$x==results_ko$Alphas[i]),3]=err_alpha_fix
        if(true_alphas==FALSE){  data[which(data$x==results_ko$Alphas[i]),1]=i}
        if(true_alphas==TRUE) {  data[which(data$x==results_ko$Alphas[i]),1]=results_ko$Alphas[i]}}
      colnames(data) = c("x","y","z")
      
      if(true_alphas==FALSE){
        #best_pair = data.frame(x=which(results_ko$Alphas==results_ko$Alpha),y=results_ko$`Number of PPCs retained`)
        plot_cv_alpha_k <- ggplot(data, aes(x = x, y = y, fill = z)) +
          geom_tile() +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(data$z)) +
          #geom_point(data = data.frame(x = results_ko$Alpha, y = results_ko$`Number of PPCs retained`), aes(x = x, y = y), color = "red", size = 4, shape = 8) + 
          labs(title = "Validation error for CV on reg param and number PPCs", x = "alpha", y = "k", fill = "Validation error") +
          theme_minimal() + 
          scale_x_continuous( breaks = seq(1,length(results_ko$Alphas), by = 1), labels = as.character(results_ko$Alphas)) + 
          theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5))}
      
      if(true_alphas==TRUE){
        plot_cv_alpha_k <- ggplot(data, aes(x = x, y = y, fill = z)) +
          geom_tile() +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(data$z)) +
          #geom_point(data = data.frame(x = results_ko$Alpha, y = results_ko$`Number of PPCs retained`), aes(x = x, y = y), color = "red", size = 4, shape = 8) + 
          labs(title = "Validation error for CV on reg param and number PPCs", x = "alpha", y = "k", fill = "Validation error") +
          theme_minimal() + 
          theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 16, hjust = 0.5))}
      
      quartz()
      print(plot_cv_alpha_k)}
  }
  
  
  
  #plot of directions and weights
  for (i in 1:results_ko$`Number of PPCs retained`) {
    
    quartz()
    # Crea i grafici con ggplot2
    p1 <- ggplot(data = data.frame(results_ko$`Function discrete evaluations points`, results_ko$`Directions of PPCs`[,i]), aes(x = results_ko$`Function discrete evaluations points`, y = results_ko$`Directions of PPCs`[,i])) +
      geom_line() +
      ggtitle("Direction") +
      labs(x = "x", y = "y") + 
      theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 10, hjust = 0.5))
    
    p2 <- ggplot(data = data.frame(results_ko$`Function discrete evaluations points`, results_ko$`Weights of PPCs`[,i]), aes(x = results_ko$`Function discrete evaluations points`, y = results_ko$`Weights of PPCs`[,i])) +
      geom_line() +
      ggtitle("Weight") +
      labs(x = "x", y = "y") + 
      theme(plot.title = element_text( family = "Arial", face = "bold", color = "black", size = 10, hjust = 0.5))
    
    grid.arrange(p1, p2, ncol = 2, top = textGrob(paste("Principal Predictive Component ",i), gp = gpar(fontsize = 20)))
  }
}
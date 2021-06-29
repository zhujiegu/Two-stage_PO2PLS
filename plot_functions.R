pca_plot <- function(data, xPC, yPC, group, name, sdev){
  var <- sdev^2
  exp_var <- 100*var/sum(var)
  ggplot(data, aes_string(x=paste('PC',xPC,sep = ''), y=paste('PC',yPC,sep = ''))) + 
    geom_point(size = 2, aes(shape = group, color = group)) + 
    stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group)) +
    ggtitle(name) +
    xlab(paste('PC', xPC, '(', round(exp_var[xPC], digits = 2), '%)')) + 
    ylab(paste('PC', yPC, '(', round(exp_var[yPC], digits = 2), '%)')) + 
    theme_bw()+
    theme(plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size=15, face = "bold"),
          legend.title = element_text(size=15, face = "bold"), axis.title = element_text(size=15),
          legend.position = 'bottom',
          legend.direction = "horizontal")
}

pls_plot <- function(data, group, name, var){
  ggplot(data, aes(x=comp.1, y=comp.2, color = group)) + geom_point(size = 2, aes(shape = group)) + 
    stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group)) +
    ggtitle(name) +
    xlab(paste('comp 1 (', round(100*var[1], digits = 2), '%)')) + 
    ylab(paste('comp 2 (', round(100*var[2], digits = 2), '%)')) + 
    theme_bw()+
    theme(plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size=15, face = "bold"),
          legend.title = element_text(size=15, face = "bold"), axis.title = element_text(size=15),
          legend.position = 'bottom',
          legend.direction = "horizontal")
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
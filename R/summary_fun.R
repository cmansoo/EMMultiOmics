#' summary graphs
#' 
#' @param box_tables result from `cv_helper`
#' @export
summary_plot <- function(box_tables){
  g1 <- ggplot2::ggplot(data=box_tables, ggplot2::aes(x=as.factor(group),y=AGE))+
    ggplot2::geom_boxplot()+
    ggplot2::xlab(NULL)+
    ggplot2::ylab('beta age')+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15))
  
  g2 <- ggplot2::ggplot(data=box_tables,ggplot2::aes(x=as.factor(group),y=PRIOR_GLIOMA))+
    ggplot2::geom_boxplot()+
    ggplot2::xlab(NULL)+
    ggplot2::ylab('beta prior glioma')+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15))
  
  g3 <- ggplot2::ggplot(data=box_tables,ggplot2::aes(x=as.factor(group),y=SEX))+
    ggplot2::geom_boxplot()+
    ggplot2::xlab(NULL)+
    ggplot2::ylab('beta sex')+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15))
  
  g4 <- ggplot2::ggplot(data=box_tables,ggplot2::aes(x=as.factor(group),y=PRETREATMENT_HISTORY))+
    ggplot2::geom_boxplot()+
    ggplot2::xlab(NULL)+
    ggplot2::ylab('beta pretreatment history')+
    ggplot2::theme(axis.text.x = element_text(angle = 15))
  
  gridExtra::grid.arrange(g1,g2,g3,g4, ncol = 2)
}








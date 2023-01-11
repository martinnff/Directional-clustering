
# Sizer para conxuntos de nivel

# Parámetros
## x = mostra circular
## h = vector cos parámetros de concentración a representar
## grid.length = tamaño do vector de bw a graficar

# Saidas
## Gráfico

Ccluster <- function(x,h=NULL,grid_length=50,cut=1E-4){
  
  dens_plot <- ggplot2::ggplot() + ggplot2::scale_x_continuous(name = '', 
                                             breaks = c(0,3*pi/2,pi,pi/2),
                                             labels = c("0",expression(pi/2),expression(pi),expression(3*pi/2))) 
  
  
  eje_bw=data.frame(rad=c(0+pi/4,pi/2+pi/4,pi+pi/4,3*pi/2+pi/4),
                    h=c(1/NPCirc::bw.rt(x),
                        NPCirc::bw.CV(x,upper=100), 
                        NPCirc::bw.pi(x), 
                        1/NPCirc::bw.boot(x)),
                    name=c("bw.rt", "bw.CV", "bw.pi", "bw.boot"))
  
  if(is.null(h)){
    bw.l=seq(0.001,max(eje_bw$h)+((max(eje_bw$h)-min(eje_bw$h))/2),length=grid_length)
  }
  if(length(h)>1){
    bw.l=seq(min(h),max(h),length=grid_length)
  }
  eje_by= round(max(bw.l))/5
  
  eje=data.frame(rad=rep(0,length(seq(0,max(bw.l),by=eje_by)))[-1],
                 h=seq(0,max(bw.l),by=eje_by)[-1],
                 name=as.character(round(seq(0,max(bw.l),by=eje_by)))[-1])
  if(length(h)<=1 & !is.null(h)){
    stop("h must ne numeric whith length > 1")
  }
  if(class(x)[1]!="circular"){
    stop("x must be circular data")
  }
  
  count1=0
  count2=0
  points = data.frame(x=x,h=rep(max(bw.l)+1,length(x)))
  points1=data.frame(matrix(ncol=3))
  points2=data.frame(matrix(ncol=3))
  colnames(points1)<-c('rad','h','f')
  colnames(points2)<-c('rad','h','f')
  for(i in seq_len(length(bw.l))){
    
    
    y=NPCirc::kern.den.circ(x, t=NULL, bw=bw.l[i], from=circular(0), to=circular(2*pi), len=250)
    i=sort(bw.l)[i]
    dx=y$x;dy=y$y
    
    
    for(j in seq_len(length(dx))){
      
      if(dy[j]>cut){
        count1=count1+1
        points1[count1,1]=dx[j]
        points1[count1,2]=i
        points1[count1,3]=dy[j]
        
      }else{
        count2=count2+1
        points2[count2,1]=dx[j]
        points2[count2,2]=i
        points2[count2,3]=0
      }
    }
    if(!any(is.na(points1))){
      dens_plot <- dens_plot + aes(x=rad,y=h,col=f) +
        geom_point(data = points1,shape=15, 
                   size = log(points1$h+1)+1) 
    }
    if(!any(is.na(points2))){
      dens_plot <- dens_plot + aes(x=rad,y=h) +
        geom_point(data = points2, shape=15,
                   colour = "#999999", 
                   size = log(points2$h+1)+1)
    }
    
  }
  dens_plot <- dens_plot +
    ggplot2::geom_point(aes(x=points$x,
                            y=rep(max(bw.l)+(max(bw.l)/20),
                            length(points$x))),
                        shape=16, 
                        colour = 1, 
                        size = 2,
                        na.rm = T)
  
  dens_plot <- dens_plot + ggplot2::geom_text(data=eje_bw,
                                     label=c("bw.rt", "bw.CV", "bw.pi", "bw.boot"),
                                     check_overlap = F,
                                     color='black'
  ) + 
    ggplot2::geom_text(data=eje,
              label=as.character(round(seq(0,max(bw.l),by=eje_by)))[-1],
              check_overlap = F,
              color='white'
    ) +
    ggplot2::scale_y_continuous(name="",
                       breaks=round(seq(0,max(bw.l),length=5)),
                       labels=rep('',length=5),
                       limits=c(min(bw.l)-max(bw.l)/4.2, max(bw.l)+(max(bw.l)/15))) +
    
    ggplot2::coord_polar(theta='x',start=pi/2)
  
  
  
  dens_plot
}
pdf('sizer.pdf')
Ccluster(x)
dev.off()





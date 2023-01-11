#library(HDiR)


set.seed(1)
x=NPCirc::rcircmix(100,model=12) # 250


x=circular::circular(x)
pdf('density_circ.pdf')
d<-NPCirc::kern.den.circ(x, 
                         t=NULL, 
                         bw=NPCirc::bw.CV(x), 
                         from=circular::circular(0), 
                         to=circular::circular(2*pi), 
                         len=500)
plot(d,main='',xlab='',ylab='',xlim=c(-1.2,1.2),ylim=c(-1.2,1.2))
circular::points.circular(x)
dev.off()


# argumentos:
## x = directional sample
## len.grid = grid steps number
## bw = concentration parameter
## plot.mode (bool)
## plot.tree (bool)

# Saidas:
## Grafics
## theshold tau to identify the core-clusters
## kernel density estimation
## number of groups




modeAnalysis <- function(x,bw=NULL,plot.mode=T,plot.tree=T,len.grid=1000){
  if(is.null(bw)){
    bw=NPCirc::bw.CV(x)
  }
  y=NPCirc::kern.den.circ(x,
                          t=NULL,
                          bw=bw, 
                          from=circular::circular(0),
                          to=circular::circular(2*pi),
                          len=len.grid)
  f.est = y
# density to probability conversion
  {
  t=0
  for(i in seq_len(length(y$y))){
    l=y$y[which(y$y>y$y[i])]
    t[i]=length(l)/length(y$y)
  }
  y$y=t
  
}

# Empirical mode function
  mod_e <-ggplot2::ggplot() + 
    ggplot2::scale_x_continuous(name=expression(tau)) + 
    ggplot2::scale_y_continuous(name='Nº modes')
  {
  modes=0
  cuts<-sort(y$y)
  for(i in seq_len(length(y$y))){
    cuts_y=0
    last=y$y[1]>=cuts[i]
    for(j in seq_len(length(y$y))){
      actual=0
      if(y$y[j]>=cuts[i]){
        actual=T
      }else{
        actual=F
      }
      if(last!=actual){
        
        cuts_y=cuts_y+1
        last=actual
      }
    }
    modes[i]<-cuts_y
  }
  modes=modes/2
  datos<-data.frame(x=sort(y$y),y=modes)
  mod_e <- mod_e + ggplot2::aes(x=x,y=y) +
    ggplot2::geom_line(data = datos) 
  t=max(cuts[which(modes==max(modes))])

}
  
# Tree diagram
  dens_plot <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(name = '', 
                       breaks = c(0,1,2,3,4,5,6),
                       labels = rep('',7)) +
    ggplot2::scale_y_continuous(name = expression(tau))
  {
    
    xc=0
    yc=0
    m=0
    last=y$y[1]
    for(i in c(seq_len(length(y$x)-1),1:20)){
      if((y$y[i] < last ) & (y$y[i] < y$y[i+1])){
        xc=c(xc,y$x[i])
        yc=c(yc,y$y[i])
        m=c(m,1)
        last=y$y[i]
      }
      if((y$y[i] > last ) & (y$y[i] > y$y[i+1])){
        xc=c(xc,y$x[i])
        yc=c(yc,y$y[i])
        m=c(m,0)
        last=y$y[i]
      }
    }
    
    max(y$y)
    
    dend=data.frame(x=xc,y=yc,m=m)
    
    dend=data.frame(x=xc,y=yc,m=m)
    joins=ifelse(dend[,3]==0,0,1)
    dend=cbind(dend,joins)
    dend=dend[order(dend[,2],decreasing=F),]
    dend=unique(dend)
    
    xmin=dend[which.max(dend[,2]),1]
    for(i in seq_len(nrow(dend))){
      
      if(dend[i,1]<xmin){
        dend[i,1]=2*pi-(xmin-dend[i,1])
      }else{
        dend[i,1]=dend[i,1]-xmin
      }
    }
    leafs=dend[which(dend[,3]==1),]
    nodes=dend[which(dend[,3]==0),]
    grid=seq(0,1,length=nrow(dend))
    
    r=0
    while(r<(length(nodes)*2)){
      r=r+1
      ai=0
      for(i in seq_len(nrow(leafs))){
        
        i=i-ai
        
        if((nrow(nodes) == 0) || (nrow(leafs) == 0)){
          break
        }
        
        for(j in seq_len(nrow(nodes))){
          
          if((leafs[i,4]<=2 & nodes[j,4]<2)){
            if(leafs[i,1]>2*pi-xmin){
              x1=xmin-(2*pi-leafs[i,1])
            }else{
              x1=leafs[i,1]+xmin
            }
            if(nodes[j,1]>2*pi-xmin){
              x2=xmin-(2*pi-nodes[j,1])
            }else{
              x2=nodes[j,1]+xmin
            }
            node<-unlist(nodes[j,])
            ind2<-(x2%/%(2*pi/length(y$y)))+1
            ind1<-as.numeric(x1%/%(2*pi/length(y$y))+1)
            inds<-sort(c(ind1,ind2))
            m<-sort(c(y$x[ind1],y$x[ind2]))
            rev_i<-c(inds[1]:1,length(y$y):inds[2])
            check1<-!any(c(y$y[inds[1]:inds[2]])>node[2])
            check2<-!any(c(y$y[rev_i])>(node[2]))
            
            if(check1 || check2){
              
              check=T
              if(leafs[i,3]!=1 & leafs[i,4]==1){
                if(leafs[i,2]<nodes[j,2]){
                  check=F
                }
              }
              if(leafs[i,2]!=nodes[j,2] & check){
                ys=c(leafs[i,2],nodes[j,2])
                x1=leafs[i,1]
                x2=nodes[j,1]
                y1=max(ys)
                y2=min(ys)
                coords_h=data.frame(x=c(x1,x2),y=c(y1,y1))
                coords_v=data.frame(x=c(x1,x1),y=c(y1,y2))
                dens_plot <- dens_plot + ggplot2::aes(x=x,y=y) +
                  ggplot2::geom_line(data = coords_h) +
                  ggplot2::geom_line(data = coords_v)
                if(leafs[i,3]==1){
                  dens_plot <- dens_plot +
                    ggplot2::aes(x=x,y=y) +
                    ggplot2::geom_text(data=leafs[i,]-c(0,0.02),
                              label=as.character(i),
                              check_overlap = T
                    )
                }
                
                leafs[i,4]=leafs[i,4]+1
                nodes[j,4]=nodes[j,4]+1
                if(any(nodes[j,1]==leafs[,1])){
                  ind=which(leafs[,1] == nodes[j,1])
                  if(length(leafs[ind,4])>1){
                    leafs=leafs[-ind,]
                    i=i-length(ind)
                    ai=ai+length(ind)
                  }
                }
                break
              }
            }
          }
          
          if((nrow(nodes) == 0) || (nrow(leafs) == 0)){
            break
          }
        }
      }
      
      ind1=-which((leafs[,3]==1 & leafs[,4]==2))
      if(length(ind1)>0){
        leafs=leafs[ind1,]}
      leafs=leafs[-which(leafs[,4]==3),]
      leafs=rbind(leafs[which(leafs[,4]>=1),],nodes[which(nodes[,4]>=1),])
      nodes=nodes[which(nodes[,4]<=1),]
      ind=0
      count=0
      for(i in seq_len(nrow(leafs))){
        
        if(any(leafs[i,1]==leafs[-i,1])){
          count=count+1
          ind[count]=-which.min(leafs[which(leafs[,1]==leafs[i,1]),4])
        }
      }
      if(length(ind)>=1){
        ind<-unique(ind)
        leafs=leafs[ind,]
      }  
      if(nrow(nodes)==1 & nodes[1,4]==0){
        nodes[1,1]=leafs[1,1]
      }
      

    }


    }
  if(plot.mode){
  plot(mod_e + ggplot2::theme_bw() + ggplot2::theme(panel.border = ggplot2::element_blank(), 
                                                    panel.grid.major = ggplot2::element_blank(),
                                  panel.grid.minor = ggplot2::element_blank(),
                                  axis.line = ggplot2::element_line(colour = "black"))
  )
  }
  if(plot.tree){
  plot(dens_plot + ggplot2::theme_bw() + ggplot2::theme(panel.border = ggplot2::element_blank(), 
                                               panel.grid.major = ggplot2::element_blank(),
                                      panel.grid.minor = ggplot2::element_blank(), 
                                      axis.line = ggplot2::element_line(colour = "black"))
  )
  }
  out=list()
  out$tau=t
  out$bw=bw
  out$nclust=max(modes)
  out$f.est=f.est
  cat(paste(" tau: ",round(out$tau,4)),"\n",
      paste("bw: ",round(out$bw,4)),"\n",
      paste("nclust: ",out$nclust),"\n")
  return(out)


}

pdf('tree_circ.pdf')
a<-modeAnalysis(x,len.grid=200,plot.mode= F)
dev.off()
pdf('modes_circ.pdf')
a<-modeAnalysis(x,len.grid=200,plot.tree= F)
dev.off()


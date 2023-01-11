# función 2

# 1ª Etapa: clasifica os puntos pertencentes a cada core-cluster
# 2ª Etapa: clasifica os puntos restantes mediante o método
#           proposto por Azzalini

# parametros
## x = mostra
## f.est = densidade estimada
## tau = punto de corte dos core-clusters
## K número de etapas para clasificar os puntos 
## non pertencentes aos core clusters

# Saidas
## groups = dataframe etiquetando cada punto a cada grupo e indicando en que etapa foi clasificado
## K número de etapas empregadas para a clasificación
## unclassified = puntos sin clasificar si hay

kernelClust <- function(x,f.est,tau, bw = NPCirc::bw.CV(x),K=NULL){
 
  y=list()
  y$y=f.est$y
  y$x=f.est$x

  
  # reordenar x para que o mínimo de densidad coincida no 0
  {
    ind<-which.min(y$y)
    xmin<-y$x[ind]
    for(i in seq_len(length(y$x))){
      if(y$x[i]<xmin){
        y$x[i]=2*pi-(xmin-y$x[i])
      }else{
        y$x[i]=y$x[i]-xmin
      }
    }
    for(i in seq_len(length(x))){
      if(x[i]<xmin){
        x[i]=2*pi-(xmin-x[i])
      }else{
        x[i]=x[i]-xmin
      }
    }
  }
  
  # paso de densidade a probabilidade
  {
    t=0
    for(i in seq_len(length(y$y))){
      l=y$y[which(y$y<=y$y[i])]
      t[i]=length(l)/length(y$y)
    }
    y$y=t
  }
  x<-sort(x)
  
  # cluster cores
  {
    
    x<-sort(x)
    t=round(1-tau,3)
    groups=list()
    clasified_ix=0
    group=data.frame(matrix(ncol=3))
    count1=0
    count2=0
    indx<-order(y$x,decreasing=F)
    y$x<-y$x[indx]
    y$y<-y$y[indx]
    last=y$y[1]>t
    
    for(i in seq_len(length(x))){
      
      index=(x[i]%/%((2*pi)/length(y$y)))+1
      if(i>1){
        index2 = (x[i-1]%/%((2*pi)/length(y$y)))+1
        last=!any(y$y[index2:index]<=t)
      }
      actual=y$y[index]>t
      
      if(actual==T & last == F){
        count1=count1+1
      }
      if(actual==T){
        count2=count2+1
        group[count2,1]=x[i]
        group[count2,2]=count1
        group[count2,3]=1
        clasified_ix=c(clasified_ix,i)
      }
      if(actual==F & last==T){
        groups[[count1]]=group
        group=data.frame(matrix(ncol=4))
        count2=0
        
      }
      last=actual
      if(i<length(x) & actual==T){
        index3= (x[i+1]%/%((2*pi)/length(y$y)))+1
        nxt=any(y$y[index:index3]<=t)
        if(nxt){
          groups[[count1]]=group
          group=data.frame(matrix(ncol=2))
          count2=0
          last=F
        }
      }
      if(i==length(x) & actual==T){
        index3 = (x[1]%/%((2*pi)/length(y$y)))+1
        nxt=any(y$y[c(index3:0,index:length(y$y))]<=t)
        if(nxt){
          groups[[count1]]=group
          group=data.frame(matrix(ncol=4))
          count2=0
          last=F
        }
      }
      

      if(i == (length(x)) & count1==1){
        groups[[count1]]=group
        group=data.frame(matrix(ncol=4))
        
      }
    }
    
    m=0
    for(i in seq_len(length(groups))){
      m[i]<-max(groups[[i]][[2]])
    }
    ind<-order(m,decreasing=T)
    count=0
    
    {
      for(j in seq_len(length(groups))){
        for(k in seq_len(length(groups[[j]][[1]]))){
          if(groups[[j]][[1]][[k]]>(2*pi-xmin)){
            groups[[j]][[1]][[k]]=groups[[j]][[1]][[k]]-(2*pi-xmin)
          }else{
            groups[[j]][[1]][[k]]=groups[[j]][[1]][[k]]+xmin
          }
        }
      }
      
      for(k in seq_len(length(x))){
        if(x[k]>(2*pi-xmin)){
          x[k]=x[k]-(2*pi-xmin)
        }else{
          x[k]=x[k]+xmin
        }
      }
    }
    
    
    
}

  # k etapas
  unclasified<-x[-clasified_ix]
  if(is.null(K)){
    K=round(length(unclasified)*1.5)
  }

  plot_c <- ggplot() + scale_x_continuous(name = '', 
                                          breaks = c(0,pi/2,pi,3*pi/2),
                                          labels = c("0",expression(pi/2),
                                                     expression(pi),
                                                     expression(3*pi/2)),
                                          limits = c(0,2*pi)) 
  
  {
  if(length(unclasified)>0){
    if(K>0){
    for(r in 1:K){
 
      dens_m<-list()
      for(i in seq_len(length(groups))){
        dens_m[[i]]<-NPCirc::kern.den.circ(circular(groups[[i]][[1]]), 
                                   t=NULL, bw=bw, 
                                   from=circular(0), 
                                   to=circular(2*pi),
                                   len=1000)
      }
      log_ratios_k=data.frame(matrix(ncol=length(unclasified),
                                     nrow=length(groups)))
      for(k in seq_len(length(unclasified))){
        
        count=0
        log_ratios=0
        for(i in 1:length(groups)){
          
          
          count=count+1
          ind<-unclasified[k]%/%((2*pi)/length(dens_m[[i]]$x))+1

          l_m<-dens_m[[i]]$y[ind]
          l_d=0
          for(s in 1:length(groups)){
            l_d<-c(l_d, dens_m[[s]]$y[ind])
          }
          l_d<-l_d[-i]
          l_max=max(l_d)
          log_ratios[count]=log(l_m/l_max)
          
        }
        log_ratios_k[,k]=log_ratios
      }
      group=apply(log_ratios_k,2,which.max)
      values=apply(log_ratios_k,2,max)
      log_ratios=data.frame(group=group,
                            values=values)
      cut<-quantile(log_ratios[,2],1-(i/100))
      ind<-which(log_ratios[,2]>=cut)
      
      for(l in ind){
        group=log_ratios[l,1]
        n=nrow(groups[[group]])+1
        groups[[group]][n,1]=unclasified[l]
        groups[[group]][n,2]=group
        groups[[group]][n,3]=r+1
      }
      unclasified=unclasified[-ind]

      if(length(unclasified)<1){
        K=r
        break
      }
    }
    }
    }
    m=0
    for(i in seq_len(length(groups))){
      m[i]<-max(groups[[i]][[2]])
    }
    ind<-order(m,decreasing=T)

    
    count=0
    cols<-palette(RColorBrewer::brewer.pal(n=length(groups),
                                           name='Set2'))

    for(i in ind){
    
      count=count+1
      points<-data.frame(rad=circular(groups[[i]][[1]]),
                         h=rep(2,length(groups[[i]][[1]])))
      plot_c <- plot_c + aes(x=rad,y=h) + 
        geom_point(data = points,
                   colour=cols[count],
                   alpha = 6/10)
    }
    points<-data.frame(rad=circular(unclasified),
                       h=rep(2,length(unclasified)))
    plot_c <- plot_c + aes(x=rad,y=h) +
      
      geom_point(data = points,
                 colour='black',
                 alpha = 6/10) +
      
      scale_y_continuous(name="",
                         breaks=c(0,1,2),
                         labels=c('','',''),
                         limits=c(0, 2.2)) + 
      
      coord_polar(theta='x',
                  start =3/2*pi,
                  direction=-1) 
    
    df<-data.frame(matrix(ncol=3))
    colnames(df)=rep('',3)

    for(g in seq_len(length(groups))){
      colnames(groups[[g]])=rep('',3)
      df<- rbind(df,groups[[g]])
    }
    df=df[-1,]
    colnames(df)<-c('x','Group','Step')
    out=list(groups=df,
             unclassified=unclasified,
             K=K)
  }
  
  plot(plot_c)
  cat(paste(" Unclassified: ",length(out$unclassified)),"\n",
      paste("Stages: ",out$K),"\n",
      paste("n_clust: ",length(groups)),"\n")
  return(out)
}
s<-kernelKlust(x,f.est=a$f.est,tau=a$tau)
names(s)


s$groups


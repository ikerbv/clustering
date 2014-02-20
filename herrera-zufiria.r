#HERRERA-ZUFIRIA MODEL

library(igraph);
library(ggplot2);
library(reshape);

nvertices <- 20;
s_exponent <- 1;
cc <- 40;


iker.harmonic <- function(n=1, m=1){
  vs <- seq(1, n);
  vvs <- sapply(vs, FUN=function(x){ return (1/x^m)});
  return (sum(vvs));
}

iker.zipf.CDF <- function(N=1, k=1, s=1){
  return (iker.harmonic(k,s)/iker.harmonic(N,s));
}

iker.zipf.pdf <- function(N=1, k=1, s=1){
  return (k^(-s)/iker.harmonic(N,s));
}

iker.zipf <- function(n=10, N=100, s=1){
  #GENERACION DE LAS UNIFORMES X,Y
  x <- ceiling(runif(1e4*n, min=0, max=N));
  y <- runif(1e4*n);
  z <- x[iker.zipf.pdf(N=N, k=x, s=s)>=y];
  return(sample(x=z, size=n, replace=F))
}

iker.v0 <- function(graph=G){
  v0 <- floor(runif(1, min=1, max=vcount(graph=graph)));
  return(v0);
}

iker.v0.with_neighbours <- function(graph=G, mode='out'){
  v0 <- iker.v0(graph);
  while(neighborhood.size(graph, order=1, nodes=v0, mode=mode) == 1){
    v0 <- iker.v0(graph=graph);
  }
  return(v0);
}

iker.rw <- function(graph=G, v0=NA, L=7, mode='out', echo=F){
  
  if(is.na(v0))
    v0 <- iker.v0.with_neighbours(graph=graph, mode=mode);
  
  nei <- neighborhood(graph, order=1, nodes=v0, mode=mode)[[1]][-1];
  V(graph)$color <- 'white';
  V(graph)[v0]$color <- 'light blue';
  E(graph)$color <- 'grey';
  v1 <- v0
  
  for(i in seq(1, L)){
    if(echo)
      print(paste('Salto: ', i));
    vv <- floor(runif(1, min=1, max=length(nei)));
    v1 <- nei[vv];
    if(echo)
      print(paste('V1 es: ', v1));
    
    e01 <- get.edge.ids(graph=graph, vp=c(v0, v1))
    V(graph)[v1]$color <- 'light green';
    E(graph)[e01]$color <- 'red';    
    nei <- neighborhood(graph, order=1, nodes=v1, mode=mode)[[1]][-1]
    if(length(nei)==0){
      v1 <- iker.v0.with_neighbours(graph, mode);
      nei <- neighborhood(graph, order=1, nodes=v1, mode=mode)[[1]][-1]
    }
    v0  <- v1;
  }
  
  if(echo)
    plot(graph);
  
  return(v1)
  
}

herrera.zufiria <- function(n=1000, s=1, cc=40, isdirected=T, directedRW=T){
  #ncero <- max(11,m);
  ncero <- 11;
  mode_rw <- ifelse(isdirected, 'out', 'all');
  mode_rw <- ifelse(directedRW, mode_rw, 'all');
  
  G <- graph.ring(n=ncero, directed=isdirected);
  pcc <- rbinom(n=n, size=1, prob=as.double(cc/100));
  degout <- iker.zipf(n=n, N=n, s=s);
  G <- G + vertices(names=seq(ncero+1,n))
  V(G)$color <- 'white';
  V(G)$pcc <- pcc;
  V(G)$deg <- degout;
  
  
  for(v in V(G)[ncero+1:n]){
    
    #m1 <- rpois(n=1, lambda=m)
    m1 <- V(G)[v]$deg;
    vv <- iker.rw(graph=G, mode=mode_rw);
    G <- add.edges(G, edges=c(v,vv), directed=isdirected);
    
    #AJOUTER LES ARÈTES
    for(l in seq(2,m1)){
      L_ <- ifelse(V(G)[vv]$pcc == 1, 1, 2)
      vv1 <- iker.rw(graph=G, v0=vv, mode=mode_rw, L=L_);
      if(get.edge.ids(graph=G, vp=c(v,vv1)) == 0)
        G <- add.edges(graph=G, edges=c(v, vv1), directed=T);
      vv <- vv1;
    }
  }
  
  
  plot(G)
  
  return(G)
}

G <- herrera.zufiria(n=nvertices, s=s_exponent, cc=cc, directedRW=F);
ds <- data.frame(din=degree(G, mode='in'), dout=degree(G, mode='out'), tr=transitivity(G, type='local'))
iker <- melt(ds);
iker <- iker[iker$variable != 'tr',];
p1 <- ggplot(data=ds, aes(x=din, y=dout, colour=tr)) + geom_point() + ggtitle('Directed Herrera-Zufiria');
p2 <- ggplot(data=iker, aes(x=value, fill=variable)) + geom_bar(position='dodge') + ggtitle('In-degree and Out-degree \n Directed Herrera-Zufiria');
p3 <- ggplot(data=ds, aes(x=tr)) + geom_bar(fill='blue', ) + geom_density();


save(c(G,ds), file='herrera_zufiria.rda');

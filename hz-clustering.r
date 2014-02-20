#HERRERA-ZUFIRIA MODEL

source('lib_hz.r');

nvertices <- 2000;
s_exponent <- 1;
new_edges <- 2;
cc <- 100;

# G <- herrera.zufiria(n=nvertices, s=s_exponent, cc=cc, directedRW=F);
# ds <- data.frame(din=degree(G, mode='in'), dout=degree(G, mode='out'), tr=transitivity(G, type='local'))
# iker <- melt(ds);
# iker <- iker[iker$variable != 'tr',];
# p1 <- ggplot(data=ds, aes(x=din, y=dout, colour=tr)) + geom_point() + ggtitle('Directed Herrera-Zufiria');
# p2 <- ggplot(data=iker, aes(x=value, fill=variable)) + geom_bar(position='dodge') + ggtitle('In-degree and Out-degree \n Directed Herrera-Zufiria');
# p3 <- ggplot(data=ds, aes(x=tr)) + geom_bar(fill='blue', ) + geom_density();
# 
# save(list=c('G','ds'), file='herrera_zufiria22.rda');

HZ <- hz.clustering(n=nvertices, s=new_edges, cc=cc, isdirected=F, directedRW=F, log=F);
tt <- transitivity(HZ$G, type='local');
hist(tt)
log <- HZ$log;

name <- paste('HZ',nvertices, new_edges, format(Sys.time(), "%Y%m%d-%H:%M"), sep='_');

save(HZ, file=paste(name,'rda',sep='.'));

p1.HZ <- ggplot(data=log[log$v %in% c(5, 50, 100, 500, 1000),], aes(x=iter, y=tr, colour=as.factor(v))) + geom_line() + ggtitle('Evolution of vertex clustering');
ggsave(p1.HZ, filename='HZ1.jpg');
transi <- aggregate(data=log, tr ~ iter, FUN=mean);

p2.HZ <- ggplot(data=transi[transi$iter %% 5 == 0, ], aes(x=iter, y=tr)) + geom_line() + ggtitle('Global Clustering Coefficient');
ggsave(p2.HZ, filename='HZ2.jpg');

p3.HZ <- ggplot(data=log[log$iter %% 5 == 0,], aes(x=factor(iter), y=tr)) + geom_boxplot(fill='light blue') + ggtitle('Clustering coefficient dispersion');
ggsave(p3.HZ, filename='HZ3.jpg');

p4.HZ <- ggplot(log[log$iter %in% c(50, 80, 100, 200, 500, 780, 900),], aes(x=tr, color=factor(iter))) + geom_density(na.rm=T)
ggsave(p4.HZ, filename='HZ4.jpg');

gini.degree <- aggregate(data=log[!is.na(log$tr),], din ~ iter, FUN=ineq);
gini.tr <- aggregate(data=log[!is.na(log$tr),], tr ~ iter, FUN=ineq);

p5.HZ <- ggplot(data=gini.degree[gini.degree$iter %% 5 == 0, ], aes(x=iter, y=din)) + geom_line() + ggtitle('Gini Index \n Vertex degree')
ggsave(p5.HZ, filename='HZ5.jpg');

p6.HZ <- ggplot(data=gini.tr[gini.tr$iter %% 5 == 0, ], aes(x=iter, y=tr)) + geom_line() + ggtitle('Gini Index \n Vertex Clustering')
ggsave(p6.HZ, filename='HZ6.jpg');

p7.HZ <- ggplot(data=log[ log$iter == max(log$iter), ], aes(x=v, y=din, colour=factor(round(tr,digits=1)))) + geom_point()
ggsave(p7.HZ, filename='HZ7.jpg');

ggplot(data=log[ log$iter == max(log$iter), ], aes(x=factor(round(tr,digits=1)), y=din)) + geom_point()



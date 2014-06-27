library(igraph)

set.seed(20)
x = sample(1:20, 50, replace=TRUE)
g = graph(x, directed=FALSE)
plot(g)

# remove loops and add small components
x[x==2] = 1:sum(x==2) 
x[x==20] = 1:sum(x==20) + 5
x = c(x, 
      21,23, 25,21, 23,25, 24,25,
      28,30, 26,29, 20,27, 1,20, 
      29,27, 27, 26, 25,3, 31,32,
      32,33, 33,31, 31,5, 34,35,
      14,1)

ug = graph(x, directed=FALSE)
dg = graph(x, directed=TRUE)

pdf("/Users/kokrah/Dropbox/large_scale_networks/figures/ug.pdf")
plot(ug, vertex.color=clusters(ug)$membership+1)
dev.off()

clusters(dg, mode="weak")
pdf("/Users/kokrah/Dropbox/large_scale_networks/figures/weak_dg.pdf")
plot(dg, edge.arrow.size=0.4, 
     vertex.color=clusters(dg, mode="weak")$membership+1)
dev.off()

clusters(dg, mode="strong")
pdf("/Users/kokrah/Dropbox/large_scale_networks/figures/strong_dg.pdf")
plot(dg, edge.arrow.size=0.4, 
     vertex.color=ifelse(clusters(dg, mode="strong")$membership==8, "blue",
                         ifelse(clusters(dg, mode="strong")$membership==13, "red",
                                ifelse(clusters(dg, mode="strong")$membership==15, "green",
                                       ifelse(clusters(dg, mode="strong")$membership==3, "yellow",
                                              "white")))))
dev.off()

pdf("/Users/kokrah/Dropbox/large_scale_networks/figures/erdos_renyi.pdf",
    height=8, width=12)
par(mfrow=c(1,2))
erdos.renyi = erdos.renyi.game(1000, .003)
clusters(erdos.renyi)
plot(erdos.renyi, vertex.size=2, vertex.label=NA, layout=layout.fruchterman.reingold,
     main="Erdos & Renyi (n=1000, p=0.002)")
hist(degree(erdos.renyi), col="light blue")
dev.off()


pdf("/Users/kokrah/Dropbox/large_scale_networks/figures/barabasi.pdf",
    height=8, width=12)
par(mfrow=c(1,2))
barabasi = barabasi.game(1000, directed=FALSE)
clusters(barabasi)
plot(barabasi, vertex.size=2, vertex.label=NA, layout=layout.fruchterman.reingold,
     main="Barabasi & Albert (n=1000)")
hist(degree(barabasi), breaks=30, col="light blue")
dev.off()

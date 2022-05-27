library(igraph)
setwd("~/epistasis3/3cell/dataLUAD")
links0=read.csv("GPairsPvalues",sep=" ",header=FALSE)
links=links0[1:10,]
links[,13] = -log10(links[,13]+ 10^{-16})
colnames(links)[13] = "weight"
net <- graph_from_data_frame(d=links[,c(1:2,13)],  directed=FALSE) 
e= eigen_centrality(net)
sort(e$vector,decreasing=TRUE)[1:12]



n=gorder(net)

library(Matrix)
cols=hcl.colors(n)
perm=order(e$vector,decreasing = FALSE)
invperm=invPerm(perm)
net1=permute(net,invperm)
l=layout_in_circle(net1) #,order=V(net)[perm])
plot(net1,layout=l,vertex.color=cols)


plot(net)
plot(net,layout=layout_with_fr)

l=layout.grid(net)
l=layout.mds(net)
l=layout.norm(l,xmin=-1,xmax=1,ymin=-2,ymax=2)
l=layout_with_kk(net)
l=layout_with_lgl(net,repulserad = 0.02)

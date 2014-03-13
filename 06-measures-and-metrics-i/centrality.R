edge.list = read.csv('/Users/daviddarmon/Documents/R/journal-club/edge_list.dat', header=FALSE)

# Get the number of edges in the network.
V = dim(edge.list)[1]

# Get the number of vertices in the network.
n = length(unique(c(edge.list[,1], edge.list[,2])))

# Form the adjacency matrix.
adj.mat = matrix(0, ncol = n, nrow = n)

for (edge.ind in 1:V){
  # The file is of the format
  #   from_node (j), to_node (i)
  # which by Newman's notation gives a_{ij}
  
  j = edge.list[edge.ind, 1]; i = edge.list[edge.ind, 2];
  
  # Use this if the network is undirected.
  # adj.mat[i, j] = 1; adj.mat[j, i] = 1;
  
  # Use this if the network is directed.
  adj.mat[i, j] = 1
}

degree.centrality = adj.mat%*%matrix(1, nrow=n, ncol = 1)

eig.decomp = eigen(adj.mat)

# Get out the largest eigenvalue of A.

lambda.1 = max(Mod(eig.decomp$values))

max.alpha = 1/lambda.1

eig.centrality = abs(eig.decomp$vectors[, 1])

alpha = 0.6

if (alpha > max.alpha){
  cat(sprintf('Warning: You\'ve chosen too large of an alpha. The Katz centrality returns nonsense!\n' ))
}

katz.centrality = solve(diag(n) - alpha*adj.mat)%*%rep(1, n)

show(cbind(degree.centrality, eig.centrality, katz.centrality))
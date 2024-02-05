
# Null model
# Curveball algorithm
# The function is described in a paper by Strona G. et al. 2014.
# Strona, G. et al. 2014. A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals. -Nat. Comm. 5: 4114. doi: 10.1038/ncomms5114 
# The function takes a matrix (e.g. polygon by species matrix), makes a swap (process described in their paper) and returns a new matrix

# The algorithm ignores spatial environmental gradients as a structuring factor

# Null model constraints
# 1. keeps the number of taxa per polygon fixed => fixed row sum
# 2. keeps the number of time a species occurrs in polygons fixed => fixed column sum
# 3. only constructs connected graphs

curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}


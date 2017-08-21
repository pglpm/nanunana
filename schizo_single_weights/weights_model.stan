data {
  int<lower=0> numquantities;
  int<lower=0> numdata;
  int<lower=2> numcategories;
  for (ic in 1:numcategories)
    vector[numquantities] x[ic,numdata];
    real p[ic];
}
parameters {
  real mu[numcategories];
  
}


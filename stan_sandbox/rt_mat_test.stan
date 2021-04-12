functions{
    matrix rt_mat_1(
    int r, // number of rows
    int c, // number of columns
    int[,] ind, // indices
    vector cov_u // unique covariate values
  ){
    int a_flat[r*c] = to_array_1d(ind);
    matrix[r,c] out = to_matrix(cov_u[a_flat], r, c, 0);
    return(out);
  }
  

    matrix rt_mat_2(
    int r, // number of rows
    int c, // number of columns
    int[,] ind, //
    vector cov_u
  ){
    matrix[r,c] out;
    for(i in 1:c){
      out[,i] = cov_u[ind[,i]];
    }
    return(out);
  }
}
data {
  int r;
  int c;
  int ind[r,c];
  vector[2] cov_u;
}

transformed data {
  matrix[r,c] a = rt_mat_1(r,c,ind,cov_u);
  matrix[r,c] b = rt_mat_2(r,c,ind,cov_u);
  print(a);
  print(b);
  print(to_array_1d(ind));
}

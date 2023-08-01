#' @export
CompTestER = function(a,b,cdfX = NULL){
  ab = abs(a*b)
  var_vec = c(var(a),var(b))
  df_min = min(ab)/sqrt(max(var_vec))
  if(is.null(cdfX)){cdfX = CDF_p2(df_min, max(ab))}
  xout = ab/sqrt(sum(var_vec)-1)
  p_value = approx(x = cdfX$x, y = cdfX$y, xout = xout, yleft = 1, yright = 0, method = 'linear')$y
  return(list(pp=p_value,zz=safe_z(p_value)*sign(a*b)))
}

CDF_p2 = function(df_min, df_max){
  len_x = 1e5
  max_df_prod = df_max * 2
  grid_x = poly_space(0, max_df_prod, len_x, order = 10)
  grid_y = 2 * besselK(grid_x, 0) / pi
  integrand = (grid_y[1:(len_x - 1)] + grid_y[2:len_x])/2 * diff(grid_x)
  int_grid_y = rev(cumsum(rev(integrand))); int_grid_y[1] = 1
  cdfX = list(x = grid_x[1:(len_x - 1)], y = int_grid_y)
  return(cdfX)
}

poly_space = function(a, b, n, order = 1){
  k = (b - a)^(1 - order)
  linsp = pracma::linspace(a, b, n)
  return(k * (linsp - a) ^ order + a)
}

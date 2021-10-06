kfold_ksmm = function(x, y, x_dim, cost_range, sigma_range, kernel, nFold = 5, nCores = 1, ...) 
{
  
  folds_ind = createFolds(y, k = nFold, list = FALSE)
  valid_err_mat = matrix(0, nrow = length(cost_range), ncol = length(sigma_range))
  # rownames(valid_err_mat) = sigma_range
  # colnames(valid_err_mat) = cost_range
  
  for (i in 1:nFold) {
    cat(i, "th fold fit \n")
    train_x = x[folds_ind != i, ]
    train_y = y[folds_ind != i]
    valid_x = x[folds_ind == i, ]
    valid_y = y[folds_ind == i]
    
    cost_valid_err_mat = matrix(NA, nrow = length(cost_range), ncol = length(sigma_range))
    
    for (j in 1:length(sigma_range)) {
      sigma = sigma_range[j]
      K = make_ksmm_kernel(train_x, train_x, kernel = kernel, sigma = sigma, dim = x_dim, nCores = nCores)
      res = mclapply(cost_range, function(cost) {
        ksmm_fit = ksmm(x = train_x, y = train_y, K = K, x_dim = x_dim, cost = cost, kernel = kernel, sigma = sigma, nCores = nCores, ...) 
        pred_y = predict.ksmm(ksmm_fit, new_x = valid_x)
        tt = table(valid_y, pred_y)
        err = 1 - sum(diag(tt)) / sum(tt)
        return(err)
      }, mc.cores = nCores)
      cost_valid_err = unlist(res)
      cost_valid_err_mat[, j] = cost_valid_err
    }
    valid_err_mat = valid_err_mat + cost_valid_err_mat / nFold
  }
  
  valid_err_vec = as.vector(valid_err_mat)
  opt_ind = max(which(valid_err_vec == min(valid_err_vec)))
  
  params = expand.grid(cost = cost_range, sigma = sigma_range)
  opt_params = params[opt_ind, ]
  opt_params = c(cost = opt_params$cost, sigma = opt_params$sigma)
  out = list()
  out$opt_params = opt_params
  out$opt_valid_err = min(valid_err_vec)
  out$valid_err = valid_err_mat
  return(out)
}

ksmm = function(x = NULL, y, K = NULL, x_dim, cost, kernel = c("linear", "rbf"), sigma = 1,
                maxit = 1e+4, epsilon = 5e-2, init_alpha = 1e-6, nCores = 1){
  
  
  call = match.call()
  kernel = match.arg(kernel)
  # Make kernel
  cat("Proceeding: make kernel matrix", "\n")
  if (is.null(K)) {
    K = make_ksmm_kernel(x, x, kernel = kernel, sigma = sigma, dim = x_dim, nCores = nCores)
  }
  cat('done', "\n")
  
  N = nrow(x)
  
  C = cost
  # a = rep(min(cost / 1000, 0.1), N)
  # a = runif(N, 0, min(cost / 1000, 0.1))
  # a = runif(N, 0, cost / 1000)
  a = rep(1e-5, N)
  
  n = dim(K[[1, 1]])[1]
  p = dim(K[[1, 1]])[2]
  temp = matrix(list(), nrow = N, ncol = 1)
  
  for(i in 1:N){
    temp[[i, 1]] = diag(a[i]*y[i], nrow = n, ncol = p)
  }
  
  cell2mat = function(m) 
  {
    return(do.call(rbind, lapply(1:nrow(m), function(x) do.call(cbind, m[x, ]))))
  }
  
  ttemp = do.call(rbind, lapply(1:nrow(temp), function(x) temp[[x, 1]]))
  KK = cell2mat(K)
  
  W = crossprod(ttemp, KK) %*% ttemp
  
  sd_tmp = crossprod(ttemp, KK)
  sd = lapply(1:(ncol(sd_tmp) / p), function(i) sd_tmp[, (p * (i - 1) + 1):(p * i)])
  
  a_old = a
  # a_old = rep(0, N)
  # aa = c()
  # k = 1
  # for (i in 1:80) {
  #   for (j in 1:80) {
  #     aa[k] = sum(is.infinite(K[[i, j]]))
  #     k = k + 1
  #   }
  # }
  
  cat("Proceeding: fit ksmm", "\n")
  
  sol = smo_smm(K, sd, W, a, as.vector(y), as.double(C), nrow(W), ncol(W), as.double(epsilon), as.integer(maxit))
  a = sol[[1]]
  b = sol[[2]]
  W = sol[[3]]
  
  fit_tmp = rep(0, nrow(K))
  for (k in 1:nrow(K)) {
    temp = matrix(0, nrow(K[[1, 1]]), ncol(K[[1, 1]]))
    for (i in 1:length(a)) {
      temp = temp + a[i] * y[i] * K[[i, k]]
    }
    fit_tmp[k] = sum(temp * W / norm(W, type = "f"))
  }
  fitted = fit_tmp + b
  
  cat('done', "\n")
  
  out = list()
  out$alpha = a
  out$b = b
  out$W = W
  out$fitted = fitted
  out$fitted_y = sign(fitted)
  out$x = x
  out$x_dim = x_dim
  out$y = y
  out$kernel = kernel
  out$sigma = sigma
  out$K = K
  out$call = call
  class(out) = "ksmm"
  return(out)
}

predict.ksmm = function(object, new_x, nCores = 1)
{
  x = object$x
  x_dim = object$x_dim
  y = object$y
  alpha = object$alpha
  b = object$b
  W = object$W
  kernel = object$kernel
  sigma = object$sigma
  K = object$K
  
  ker_test = make_ksmm_kernel(x = x, y = new_x, kernel = kernel, sigma = sigma, dim = x_dim, nCores = nCores)
  
  pred_tmp = rep(0, nrow(new_x))
  for (k in 1:nrow(new_x)) {
    temp = matrix(0, nrow(ker_test[[1, 1]]), ncol(ker_test[[1, 1]]))
    for (i in 1:length(alpha)) {
      temp = temp + alpha[i] * y[i] * ker_test[[i, k]]
    }
    pred_tmp[k] = sum(temp * W / norm(W, type = "f"))
  }
  
  pred_y = ifelse(pred_tmp + b > 0, 1, -1)
  return(pred_y)
}


make_ksmm_kernel = function(x, y, kernel, sigma, dim, nCores = 1)
{
  n1 = nrow(x)
  n2 = nrow(y)
  
  ker_train = matrix(list(), n1, n2)
  
  
  kernel_fun = function(x, y, type = "linear", sigma = NULL) {
    dot_fun = switch(type,
                     linear = vanilladot(),
                     rbf = rbfdot(sigma = sigma))
    
    K = kernelMatrix(dot_fun, x, y)
    if (type == "linear") {
      K = K 
      # + diag(sigma, ncol(x))
    }
    return(K)
  }
  
  ij = expand.grid(n1 = 1:n1, n2 = 1:n2)
  
  kern_list = mclapply(1:nrow(ij), FUN = function(tt) {
    i = ij[tt, 1]
    j = ij[tt, 2]
    res = kernel_fun(x = matrix(x[i, ], dim[1], dim[2]),
                     y = matrix(y[j, ], dim[1], dim[2]),
                     type = kernel,
                     sigma = sigma)
    return(res)
  }, mc.cores = nCores)
  
  kern_mat = matrix(kern_list, nrow = n1, ncol = n2)
  return(kern_mat)
}

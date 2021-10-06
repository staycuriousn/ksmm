
sim_gen = function(n, sd, seed, type = c("type1", "type2")) 
{
  sim_type = match.arg(type)
  out = list()
  if (sim_type == "type1") {
    set.seed(seed)
    v_vec = u_vec = list()
    for (i in 1:10) {
      v_vec[[i]] = runif(n = n, min = -1, max = 1)
    }
    v_mat = matrix(unlist(v_vec), nrow = n)
    qr = qr(v_mat)
    u_mat = qr.Q(qr)
    simu_data = list()
    for (i in 1:nrow(u_mat)) {
      mat = u_mat[i, 1] + matrix(0, nrow = 80, ncol = 10)
      for (j in 2:ncol(u_mat)) {
        mat = cbind(mat, (u_mat[i, j] + matrix(0, nrow = 80, ncol = 10)))
      }
      simu_data[[i]] = mat + rnorm(n = 8000, mean = 0, sd = 1e-3)
    }
    
    rank = 10
    nn = 10
    x = matrix(0, nrow = 80, ncol = 10 * nn)
    for (i in 1:10) {
      x[i, i * 10 - 1] = 1
    }
    
    W = x
    
    y = c()
    for (j in 1:n) {
      y[j] = ifelse(sum(diag(t(W) %*% simu_data[[j]])) >= 0, 1, -1)
    }
    
    simulation_dat = as.vector(simu_data[[1]])
    for (i in 2:n) {
      simulation_dat = rbind(simulation_dat, as.vector(simu_data[[i]]))
    }
    simulation_dat = simulation_dat + matrix(rnorm((nrow(simulation_dat) * ncol(simulation_dat)), mean = 0, sd = sd), 
                                             nrow = nrow(simulation_dat), ncol = ncol(simulation_dat))
    out$X = simulation_dat
    out$y = y
  }
  
  if (sim_type == "type2") {
    set.seed(seed)
    v_vec = u_vec = list()
    for (i in 1:10) {
      v_vec[[i]] = runif(n = n, min = -1, max = 1)
    }
    v_mat = matrix(unlist(v_vec), nrow = n)
    qr = qr(v_mat)
    u_mat = qr.Q(qr)
    simu_data = list()
    for (i in 1:nrow(u_mat)) {
      mat = u_mat[i, 1] + matrix(0, nrow = 80, ncol = 10)
      for (j in 2:ncol(u_mat)) {
        mat = cbind(mat, (u_mat[i, j] + matrix(0, nrow = 80, ncol = 10)))
      }
      simu_data[[i]] = mat + rnorm(n = 8000, mean = 0, sd = 1e-3)
    }
    
    rank = 10
    nn = 10
    set.seed(1)
    A_temp = matrix(rnorm(80 * 80, 10), 80, 80)
    A = t(A_temp) + A_temp
    ed = eigen(A)
    WW = NULL
    for (i in 1:10) {
      WW = cbind(WW, ed$vectors[, 1:10])
    }
    # W = x
    W = WW
    
    y_temp = c()
    xx = c()
    # for (j in 1:n) {
    #   wx = sum(diag(t(W) %*% simu_data[[j]]))
    #   y_temp[j] = -3 * wx^2 + 2 * wx + 8
    #   xx[j] = wx
    # }
    
    for (j in 1:n) {
      # wx = sum(diag(t(W) %*% simu_data[[j]]))
      wx = 2*sum(diag(t(W) %*% simu_data[[j]]^2)) + 1.8
      # y_temp[j] = sin(pi * wx)
      xx[j] = wx
    }
    
    # plot(xx, y_temp)
    # abline(h = 0)
    
    # for (j in 1:n) {
    #   y_temp[j] = sin(pi * sum(diag(t(W) %*% simu_data[[j]])))
    # }
    
    # hist(y_temp)
    y = ifelse(xx > 0, 1, -1)
    table(y)
    # y = c()
    # for (j in 1:n) {
    #   y[j] = sum(diag(t(W) %*% simu_data[[j]]))
    # }
    
    simulation_dat = as.vector(simu_data[[1]])
    for (i in 2:n) {
      simulation_dat = rbind(simulation_dat, as.vector(simu_data[[i]]))
    }
    simulation_dat = simulation_dat + matrix(rnorm((nrow(simulation_dat) * ncol(simulation_dat)), mean = 0, sd = sd), 
                                             nrow = nrow(simulation_dat), ncol = ncol(simulation_dat))
    out$X = simulation_dat
    out$y = y
  }
  return(out)
}


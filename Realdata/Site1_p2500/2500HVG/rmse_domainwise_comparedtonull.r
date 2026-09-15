##### Site 1: domain-wise relative MSE vs null model
##### Relative MSE < 1 means better prediction than null.
##### No log transformation.

fname_vec = c(
  "MTL_transfer_cluster",
  "MTL_transfer_debiased",
  "MTL_notransfer",
  "MTL_notransfer_merged",
  "STL_transfer_cluster",
  "STL_transfer_debiased"
)

legendvec = c(
  "TMTL(Fused)",
  "TMTL(Debiased)",
  "MTL(Target)",
  "MTL(Full)",
  "TSTL(Fused)",
  "TSTL(Debiased)"
)

donorvec = 1:12

#### Read summary files without overwriting existing variables

read_summary = function(folder, domain) {
  
  path = file.path(
    folder,
    paste0("summary_logtransform_target", domain, ".Rdata")
  )
  
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  
  e = new.env(parent = emptyenv())
  load(path, envir = e)
  
  if (!exists("prediction_mse_vec", envir = e, inherits = FALSE)) {
    stop("prediction_mse_vec not found in: ", path)
  }
  
  as.list(e)
}

#### Number of response tasks

data_env = new.env(parent = emptyenv())
load(
  "clr_data_donor_nofiltered_2500hvg1.Rdata",
  envir = data_env
)

K = ncol(data_env$Y)

stopifnot(
  length(K) == 1L,
  is.finite(K),
  K > 0
)

#### Domain labels

if (file.exists("domain_index.csv")) {
  domain_names = as.character(
    read.csv("domain_index.csv")[[2]][donorvec]
  )
} else {
  domain_names = paste("domain", donorvec)
}

stopifnot(
  length(domain_names) == length(donorvec),
  !anyNA(domain_names)
)

#### Storage

n_donor_vec = numeric(length(donorvec))
nrep_vec = integer(length(donorvec))

mse_mat_list = vector("list", length(donorvec))
null_mse_list = vector("list", length(donorvec))
relative_mse_mat_list = vector("list", length(donorvec))
summary_list = vector("list", length(donorvec))

#### Domain-wise calculation

for (l in seq_along(donorvec)) {
  
  domain = donorvec[l]
  
  target = read_summary("MTL_notransfer", domain)
  n_test = target$test_sample_size
  nrep = length(target$prediction_mse_vec)
  
  if (!is.numeric(n_test) ||
      length(n_test) != 1L ||
      !is.finite(n_test) ||
      n_test <= 0) {
    stop("Invalid test_sample_size for domain ", domain)
  }
  
  if (nrep != 20L) {
    stop(
      "Expected 20 replications for domain ", domain,
      ", but found ", nrep
    )
  }
  
  n_donor_vec[l] = n_test
  nrep_vec[l] = nrep
  
  # Check prediction values and, when available, test sample size.
  check_summary = function(x, folder) {
    
    v = x$prediction_mse_vec
    
    if (!is.numeric(v) ||
        length(v) != nrep ||
        any(!is.finite(v)) ||
        any(v < 0)) {
      stop(
        "Invalid prediction_mse_vec: ",
        folder, ", domain ", domain
      )
    }
    
    if (!is.null(x$test_sample_size) &&
        !isTRUE(all.equal(
          as.numeric(x$test_sample_size),
          as.numeric(n_test)
        ))) {
      stop(
        "Test sample size mismatch: ",
        folder, ", domain ", domain
      )
    }
    
    v
  }
  
  #### Null model MSE
  # Following rmse.r, nullmodel stores SSE over cells and responses.
  
  null = read_summary("nullmodel", domain)
  null_sse = check_summary(null, "nullmodel")
  
  null_mse = null_sse / (K * n_test)
  
  if (any(null_mse <= 0)) {
    stop("Null MSE must be positive for domain ", domain)
  }
  
  #### Method MSE
  
  mse_mat_l = matrix(
    NA_real_,
    nrow = length(fname_vec),
    ncol = nrep,
    dimnames = list(
      legendvec,
      paste0("rep", seq_len(nrep))
    )
  )
  
  for (a in seq_along(fname_vec)) {
    
    x = read_summary(fname_vec[a], domain)
    v = check_summary(x, fname_vec[a])
    
    if (fname_vec[a] == "STL_transfer_cluster") {
      # TSTL(Fused) stores RMSE in the original code.
      mse_mat_l[a, ] = v^2
    } else {
      # Other methods store SSE over cells and responses.
      mse_mat_l[a, ] = v / (K * n_test)
    }
  }
  
  #### Relative MSE, paired by replication
  # Column r is divided by the null MSE from replication r.
  # Replication order must represent the same splits across methods.
  
  relative_mse_mat_l = sweep(
    mse_mat_l,
    MARGIN = 2,
    STATS = null_mse,
    FUN = "/"
  )
  
  mse_mat_list[[l]] = mse_mat_l
  null_mse_list[[l]] = null_mse
  relative_mse_mat_list[[l]] = relative_mse_mat_l
  
  #### Mean and SD of replication-level relative MSE
  
  tab_l = data.frame(
    method = legendvec,
    mean = rowMeans(relative_mse_mat_l),
    sd = apply(relative_mse_mat_l, 1, sd),
    row.names = NULL
  )
  
  # Keep full precision in summary_list.
  summary_list[[l]] = tab_l
  
  # Round only for display and CSV output.
  tab_print = tab_l
  tab_print$mean = round(tab_print$mean, 3)
  tab_print$sd = round(tab_print$sd, 3)
  
  cat("\n==========================================================\n")
  cat(
    "Site 1 | domain ", domain, " : ", domain_names[l],
    "   n_test = ", n_test,
    "   nrep = ", nrep, "\n",
    sep = ""
  )
  cat("Relative MSE = method MSE / null MSE\n")
  cat("==========================================================\n")
  print(tab_print)
  
  write.csv(
    tab_print,
    paste0(
      "relative_mse_null_site1_domain", domain, "_",
      gsub("[^A-Za-z0-9]", "", domain_names[l]),
      ".csv"
    ),
    row.names = FALSE
  )
}

#### Combined table: one row per domain and method

summary_all = do.call(
  rbind,
  lapply(seq_along(donorvec), function(l) {
    data.frame(
      site = 1L,
      domain = donorvec[l],
      domain_name = domain_names[l],
      n_test = n_donor_vec[l],
      nrep = nrep_vec[l],
      summary_list[[l]],
      row.names = NULL
    )
  })
)

summary_all_print = summary_all
summary_all_print$mean = round(summary_all_print$mean, 3)
summary_all_print$sd = round(summary_all_print$sd, 3)

# write.csv(
#   summary_all_print,
#   "relative_mse_null_site1_all_domains.csv",
#   row.names = FALSE
# )

#### Wide tables: methods in rows, domains in columns

mean_wide = sapply(
  relative_mse_mat_list,
  function(x) rowMeans(x)
)

sd_wide = sapply(
  relative_mse_mat_list,
  function(x) apply(x, 1, sd)
)

colnames(mean_wide) = colnames(sd_wide) = domain_names
rownames(mean_wide) = rownames(sd_wide) = legendvec

cat("\n--- Site 1: relative MSE mean by domain ---\n")
print(round(mean_wide, 3))

cat("\n--- Site 1: relative MSE SD by domain ---\n")
print(round(sd_wide, 3))

# write.csv(
#   data.frame(
#     method = legendvec,
#     round(mean_wide, 3),
#     check.names = FALSE,
#     row.names = NULL
#   ),
#   "relative_mse_null_site1_mean_wide.csv",
#   row.names = FALSE
# )

# write.csv(
#   data.frame(
#     method = legendvec,
#     round(sd_wide, 3),
#     check.names = FALSE,
#     row.names = NULL
#   ),
#   "relative_mse_null_site1_sd_wide.csv",
#   row.names = FALSE
# )
# 
# #### Save full-precision results for subsequent plotting
# 
# saveRDS(
#   list(
#     site = 1L,
#     domains = donorvec,
#     domain_names = domain_names,
#     methods = legendvec,
#     n_test = n_donor_vec,
#     nrep = nrep_vec,
#     method_mse = mse_mat_list,
#     null_mse = null_mse_list,
#     relative_mse = relative_mse_mat_list,
#     summary = summary_all
#   ),
#   "relative_mse_null_site1_results.rds"
# )










































##### Site 1: domain-wise relative MSE vs null model
##### Relative MSE < 1 means better prediction than null.
##### No log transformation.

fname_vec = c(
  "MTL_transfer_cluster",
  "MTL_transfer_debiased",
  "MTL_notransfer",
  "MTL_notransfer_merged",
  "STL_transfer_cluster",
  "STL_transfer_debiased"
)

legendvec = c(
  "TMTL(Fused)",
  "TMTL(Debiased)",
  "MTL(Target)",
  "MTL(Full)",
  "TSTL(Fused)",
  "TSTL(Debiased)"
)

donorvec = 1:12

#### Read summary files without overwriting existing variables

read_summary = function(folder, domain) {
  
  path = file.path(
    folder,
    paste0("summary_logtransform_target", domain, ".Rdata")
  )
  
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  
  e = new.env(parent = emptyenv())
  load(path, envir = e)
  
  if (!exists("prediction_mse_vec", envir = e, inherits = FALSE)) {
    stop("prediction_mse_vec not found in: ", path)
  }
  
  as.list(e)
}

#### Number of response tasks

data_env = new.env(parent = emptyenv())
load(
  "clr_data_donor_nofiltered_2500hvg1.Rdata",
  envir = data_env
)

K = ncol(data_env$Y)

stopifnot(
  length(K) == 1L,
  is.finite(K),
  K > 0
)

#### Domain labels

if (file.exists("domain_index.csv")) {
  domain_names = as.character(
    read.csv("domain_index.csv")[[2]][donorvec]
  )
} else {
  domain_names = paste("domain", donorvec)
}

stopifnot(
  length(domain_names) == length(donorvec),
  !anyNA(domain_names)
)

#### Storage

n_donor_vec = numeric(length(donorvec))
nrep_vec = integer(length(donorvec))

mse_mat_list = vector("list", length(donorvec))
null_mse_list = vector("list", length(donorvec))
relative_mse_mat_list = vector("list", length(donorvec))
summary_list = vector("list", length(donorvec))

#### Domain-wise calculation

for (l in seq_along(donorvec)) {
  
  domain = donorvec[l]
  
  target = read_summary("MTL_notransfer", domain)
  n_test = target$test_sample_size
  nrep = length(target$prediction_mse_vec)
  
  if (!is.numeric(n_test) ||
      length(n_test) != 1L ||
      !is.finite(n_test) ||
      n_test <= 0) {
    stop("Invalid test_sample_size for domain ", domain)
  }
  
  if (nrep != 20L) {
    stop(
      "Expected 20 replications for domain ", domain,
      ", but found ", nrep
    )
  }
  
  n_donor_vec[l] = n_test
  nrep_vec[l] = nrep
  
  # Check prediction values and, when available, test sample size.
  check_summary = function(x, folder) {
    
    v = x$prediction_mse_vec
    
    if (!is.numeric(v) ||
        length(v) != nrep ||
        any(!is.finite(v)) ||
        any(v < 0)) {
      stop(
        "Invalid prediction_mse_vec: ",
        folder, ", domain ", domain
      )
    }
    
    if (!is.null(x$test_sample_size) &&
        !isTRUE(all.equal(
          as.numeric(x$test_sample_size),
          as.numeric(n_test)
        ))) {
      stop(
        "Test sample size mismatch: ",
        folder, ", domain ", domain
      )
    }
    
    v
  }
  
  #### Null model MSE
  # Following rmse.r, nullmodel stores SSE over cells and responses.
  
  null = read_summary("nullmodel", domain)
  null_sse = check_summary(null, "nullmodel")
  
  null_mse = null_sse / (K * n_test)
  
  if (any(null_mse <= 0)) {
    stop("Null MSE must be positive for domain ", domain)
  }
  
  #### Method MSE
  
  mse_mat_l = matrix(
    NA_real_,
    nrow = length(fname_vec),
    ncol = nrep,
    dimnames = list(
      legendvec,
      paste0("rep", seq_len(nrep))
    )
  )
  
  for (a in seq_along(fname_vec)) {
    
    x = read_summary(fname_vec[a], domain)
    v = check_summary(x, fname_vec[a])
    
    if (fname_vec[a] == "STL_transfer_cluster") {
      # TSTL(Fused) stores RMSE in the original code.
      mse_mat_l[a, ] = v^2
    } else {
      # Other methods store SSE over cells and responses.
      mse_mat_l[a, ] = v / (K * n_test)
    }
  }
  
  #### Relative MSE, paired by replication
  # Column r is divided by the null MSE from replication r.
  # Replication order must represent the same splits across methods.
  
  relative_mse_mat_l = sweep(
    mse_mat_l,
    MARGIN = 2,
    STATS = null_mse,
    FUN = "/"
  )
  
  mse_mat_list[[l]] = mse_mat_l
  null_mse_list[[l]] = null_mse
  relative_mse_mat_list[[l]] = relative_mse_mat_l
  
  #### Mean and SD of replication-level relative MSE
  
  tab_l = data.frame(
    method = legendvec,
    mean = rowMeans(relative_mse_mat_l),
    sd = apply(relative_mse_mat_l, 1, sd),
    row.names = NULL
  )
  
  # Keep full precision in summary_list.
  summary_list[[l]] = tab_l
  
  # Round only for display and CSV output.
  tab_print = tab_l
  tab_print$mean = round(tab_print$mean, 3)
  tab_print$sd = round(tab_print$sd, 3)
  
  cat("\n==========================================================\n")
  cat(
    "Site 1 | domain ", domain, " : ", domain_names[l],
    "   n_test = ", n_test,
    "   nrep = ", nrep, "\n",
    sep = ""
  )
  cat("Relative MSE = method MSE / null MSE\n")
  cat("==========================================================\n")
  print(tab_print)
  
  write.csv(
    tab_print,
    paste0(
      "relative_mse_null_site1_domain", domain, "_",
      gsub("[^A-Za-z0-9]", "", domain_names[l]),
      ".csv"
    ),
    row.names = FALSE
  )
}

#### Combined table: one row per domain and method

summary_all = do.call(
  rbind,
  lapply(seq_along(donorvec), function(l) {
    data.frame(
      site = 1L,
      domain = donorvec[l],
      domain_name = domain_names[l],
      n_test = n_donor_vec[l],
      nrep = nrep_vec[l],
      summary_list[[l]],
      row.names = NULL
    )
  })
)

summary_all_print = summary_all
summary_all_print$mean = round(summary_all_print$mean, 3)
summary_all_print$sd = round(summary_all_print$sd, 3)

# write.csv(
#   summary_all_print,
#   "relative_mse_null_site1_all_domains.csv",
#   row.names = FALSE
# )

#### Wide tables: methods in rows, domains in columns

mean_wide = sapply(
  relative_mse_mat_list,
  function(x) rowMeans(x)
)

sd_wide = sapply(
  relative_mse_mat_list,
  function(x) apply(x, 1, sd)
)

colnames(mean_wide) = colnames(sd_wide) = domain_names
rownames(mean_wide) = rownames(sd_wide) = legendvec

cat("\n--- Site 1: relative MSE mean by domain ---\n")
print(round(mean_wide, 3))

cat("\n--- Site 1: relative MSE SD by domain ---\n")
print(round(sd_wide, 3))

# write.csv(
#   data.frame(
#     method = legendvec,
#     round(mean_wide, 3),
#     check.names = FALSE,
#     row.names = NULL
#   ),
#   "relative_mse_null_site1_mean_wide.csv",
#   row.names = FALSE
# )

# write.csv(
#   data.frame(
#     method = legendvec,
#     round(sd_wide, 3),
#     check.names = FALSE,
#     row.names = NULL
#   ),
#   "relative_mse_null_site1_sd_wide.csv",
#   row.names = FALSE
# )
# 
# #### Save full-precision results for subsequent plotting
# 
# saveRDS(
#   list(
#     site = 1L,
#     domains = donorvec,
#     domain_names = domain_names,
#     methods = legendvec,
#     n_test = n_donor_vec,
#     nrep = nrep_vec,
#     method_mse = mse_mat_list,
#     null_mse = null_mse_list,
#     relative_mse = relative_mse_mat_list,
#     summary = summary_all
#   ),
#   "relative_mse_null_site1_results.rds"
# )

##### Site 1: grouped barplot of mean relative MSE (base R)
##### Rows of mean_wide = 6 methods; columns = 12 domains.
##### Bars show means over replications; no error bars, as in the reference.

colvec = c("red1", "pink2", "green1", "green4", "blue", "skyblue3")

stopifnot(
  nrow(mean_wide) == 6L,
  ncol(mean_wide) == length(donorvec),
  length(domain_names) == ncol(mean_wide),
  length(legendvec) == nrow(mean_wide),
  all(is.finite(mean_wide)),
  all(mean_wide >= 0)
)

plot_relative_mse_site1 = function() {
  
  # Reserve a separate bottom panel for the legend.
  layout(matrix(c(1, 2), ncol = 1), heights = c(6.2, 0.8))
  on.exit(layout(1), add = TRUE)
  
  par(
    mar = c(11.5, 5, 2, 1),
    mgp = c(3, 0.7, 0),
    las = 1,
    family = "sans",
    xaxs = "i",
    yaxs = "i"
  )
  
  # Use zero as the baseline and always include the null reference at 1.
  ymax = max(1.2, ceiling(max(mean_wide) * 1.08 / 0.2) * 0.2)
  yticks = seq(0, ymax, by = 0.2)
  
  # Obtain bar positions without drawing.
  bar_x = barplot(
    mean_wide,
    beside = TRUE,
    space = c(0, 1),
    plot = FALSE
  )
  domain_x = colMeans(bar_x)
  xlim = c(min(bar_x) - 1, max(bar_x) + 1)
  
  plot.new()
  plot.window(xlim = xlim, ylim = c(0, ymax))
  
  # Light grid behind the bars.
  abline(h = yticks, col = "gray90", lwd = 0.7)
  abline(v = domain_x, col = "gray90", lwd = 0.7)
  
  barplot(
    mean_wide,
    beside = TRUE,
    space = c(0, 1),
    col = colvec,
    border = NA,
    axes = FALSE,
    axisnames = FALSE,
    add = TRUE
  )
  
  # Null model has relative MSE exactly 1.
  abline(h = 1, col = "black", lwd = 1.1)
  
  axis(
    2, at = yticks,
    labels = formatC(yticks, format = "f", digits = 1),
    col = "gray60", col.axis = "gray30", cex.axis = 1
  )
  box(col = "gray65", lwd = 0.7)
  mtext("Relative MSE (vs null model)", side = 2, line = 3.2,
        font = 2, cex = 1.05, las = 0)
  title(main = "Site 1", cex.main = 1.1)
  
  # Full domain names, rotated to avoid overlap.
  text(
    x = domain_x,
    y = -0.025 * ymax,
    labels = domain_names,
    srt = 45,
    adj = c(1, 1),
    xpd = NA,
    cex = 0.85,
    font = 2,
    col = "gray30"
  )
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = legendvec,
    fill = colvec,
    border = NA,
    ncol = 6,
    bty = "o",
    box.col = "black",
    cex = 0.95,
    x.intersp = 0.6,
    y.intersp = 1.1
  )
}

# Save a vector PDF. The full labels need a sufficiently wide device.
pdf("relative_mse_null_site1_barplot.pdf", width = 15, height = 7.5)
tryCatch(plot_relative_mse_site1(), finally = dev.off())

# Optional: draw in the current graphics window.
# plot_relative_mse_site1()

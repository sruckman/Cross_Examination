library(shiny)
library(ggplot2)
library(gridExtra)

# Chromosome positions: 0 to 100 cM in 0.5 cM steps = 202 positions
POS_SEQ <- seq(0, 100, length.out = 202L)
N_POS   <- length(POS_SEQ)

# Background colors for dark mode
BG_PAGE       <- "#000000"
BG_PANEL      <- "#141414"
BG_PLOT       <- "#000000"
BG_PLOT_LIGHT <- "#111111"
COL_GRID      <- "#2a2a2a"

# F2 colors (global so plot helpers can access them)
F2_COLS <- c("1" = "#CC0000", "2" = "#1E90FF")
F2_LBLS <- c("1" = "P1",     "2" = "P2")

# Teaching-mode fixed QTL positions (same every run when teaching mode is on)
TEACH_POS <- list(
  "1qtl"      = 50L,
  "2qtl"      = c(30L, 70L),
  "3qtl"      = c(20L, 50L, 80L),
  "polygenic" = as.integer(round(seq(10, 90, length.out = 10L)))
)

# Teaching-mode fixed GWAS causal SNP indices (out of 1000 SNPs)
TEACH_SNP <- list(
  "1qtl"      = 500L,
  "2qtl"      = c(250L, 750L),
  "3qtl"      = c(200L, 500L, 800L),
  "polygenic" = as.integer(round(seq(0.1, 0.9, length.out = 10L) * 1000L))
)

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

get_founder_colors <- function(n) {
  base <- c(
    "#CC0000",  # dark red
    "#66BB6A",  # medium green
    "#80D8FF",  # light sky blue
    "#FFA726",  # amber
    "#AB47BC",  # purple
    "#26C6DA",  # cyan
    "#E040FB",  # vivid magenta
    "#9CCC65",  # lime green
    "#5C6BC0",  # indigo
    "#FFEE58",  # bright yellow
    "#EC407A",  # hot pink
    "#A1887F"   # rosy brown
  )
  if (n <= 12L) return(base[seq_len(n)])
  colorRampPalette(base)(n)
}

theme_dark_mpp <- function(base_size = 16) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.background   = element_rect(fill = BG_PLOT,  color = NA),
      panel.background  = element_rect(fill = BG_PLOT,  color = NA),
      panel.grid.major  = element_line(color = COL_GRID, linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = "#ffffff"),
      axis.title        = element_text(color = "#ffffff"),
      strip.background  = element_rect(fill = "#111111", color = NA),
      strip.text        = element_text(color = "#ffffff", size = 11),
      legend.background = element_rect(fill = BG_PLOT, color = NA),
      legend.text       = element_text(color = "#ffffff"),
      legend.title      = element_text(color = "#ffffff"),
      legend.key        = element_rect(fill = BG_PLOT, color = NA)
    )
}

smooth_lod <- function(lod, window = 5L) {
  n    <- length(lod)
  half <- window %/% 2L
  lod[!is.finite(lod)] <- 0
  vapply(seq_len(n), function(i) {
    mean(lod[max(1L, i - half):min(n, i + half)])
  }, numeric(1L))
}

chisq_lod <- function(statistic, df) {
  log_p <- pchisq(as.numeric(statistic), df = as.numeric(df),
                  lower.tail = FALSE, log.p = TRUE)
  -log_p / log(10)
}

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================

simulate_rils <- function(n_founders, n_gen, n_total, design) {
  haps <- matrix(1L, nrow = n_total, ncol = N_POS)
  for (r in seq_len(n_total)) {
    if (design == "hub") {
      hub   <- 1L
      other <- if (n_founders >= 2L) sample.int(n_founders - 1L, 1L) + 1L else 1L
      pool  <- c(hub, other)
    } else {
      pool <- seq_len(n_founders)
    }
    state       <- rep(sample(pool, 1L), N_POS)
    cross_rate  <- if (design == "hub") 1.0 else 0.1
    total_cross <- rpois(1L, n_gen * cross_rate)
    if (total_cross > 0L) {
      cpts <- sort(runif(total_cross, 0, 100))
      for (cp in cpts) {
        idx <- which(POS_SEQ >= cp)[1L]
        if (is.na(idx) || idx > N_POS) next
        state[idx:N_POS] <- sample(pool, 1L)
      }
    }
    haps[r, ] <- state
  }
  haps
}

# teaching = TRUE uses fixed QTL positions so peaks stay in the same place
# across repeated Simulate clicks — useful for classroom demonstrations.
get_qtl_config <- function(model, n_founders, hub_design = FALSE, teaching = FALSE) {
  avail    <- if (hub_design && n_founders > 1L) seq(2L, n_founders) else seq_len(n_founders)
  n_effect <- max(1L, min(10L, round(0.1 * length(avail))))

  hub_scale <- if (hub_design && n_founders > 2L) {
    f_hub   <- 0.5 / (n_founders - 1L)
    f_magic <- 1.0 / n_founders
    sqrt(f_magic * (1 - f_magic) / (f_hub * (1 - f_hub)))
  } else 1.0

  make_eff <- function(effect_size) {
    eff      <- numeric(n_founders)
    who      <- sample(avail, min(n_effect, length(avail)))
    scaled   <- effect_size * hub_scale
    eff[who] <- runif(length(who), scaled * 0.75, scaled * 1.25)
    eff
  }

  switch(model,
    "null"      = list(pos = integer(0), eff = list()),
    "1qtl"      = { p <- if (teaching) TEACH_POS[["1qtl"]] else sample(10L:90L, 1L)
                    list(pos = p, eff = list(make_eff(1.5))) },
    "2qtl"      = { p <- if (teaching) TEACH_POS[["2qtl"]] else sort(sample(10L:90L, 2L))
                    list(pos = p, eff = list(make_eff(1.2), make_eff(1.2))) },
    "3qtl"      = { p <- if (teaching) TEACH_POS[["3qtl"]] else sort(sample(10L:90L, 3L))
                    list(pos = p, eff = list(make_eff(1.0), make_eff(1.0), make_eff(1.0))) },
    "polygenic" = { p <- if (teaching) TEACH_POS[["polygenic"]] else sort(sample(5L:95L, 10L))
                    list(pos = p, eff = lapply(seq_len(10L), function(i) make_eff(0.7))) }
  )
}

simulate_phenotype <- function(true_haps, qtl_cfg) {
  pheno <- rnorm(nrow(true_haps), 0, 1)
  for (i in seq_along(qtl_cfg$pos)) {
    col_idx <- which.min(abs(POS_SEQ - qtl_cfg$pos[i]))
    pheno   <- pheno + qtl_cfg$eff[[i]][true_haps[, col_idx]]
  }
  pheno
}

qtl_scan <- function(haps_list, case_idx_list, control_idx_list) {
  n_reps <- length(haps_list)
  n_pos  <- ncol(haps_list[[1L]])
  lod    <- numeric(n_pos)

  for (p in seq_len(n_pos)) {
    lvls <- sort(unique(unlist(lapply(seq_len(n_reps), function(ri) {
      c(haps_list[[ri]][case_idx_list[[ri]], p],
        haps_list[[ri]][control_idx_list[[ri]], p])
    }))))
    if (length(lvls) < 2L) next

    seq_depth <- 200L

    if (n_reps == 1L) {
      cf     <- haps_list[[1L]][case_idx_list[[1L]], p]
      kf     <- haps_list[[1L]][control_idx_list[[1L]], p]
      cf_raw <- as.integer(table(factor(cf, levels = lvls)))
      kf_raw <- as.integer(table(factor(kf, levels = lvls)))
      n_cf   <- sum(cf_raw); n_kf <- sum(kf_raw)
      cf_obs <- if (n_cf > seq_depth) as.integer(rmultinom(1L, seq_depth, cf_raw / n_cf)) else cf_raw
      kf_obs <- if (n_kf > seq_depth) as.integer(rmultinom(1L, seq_depth, kf_raw / n_kf)) else kf_raw
      tab    <- rbind(cf_obs, kf_obs)
      tryCatch({
        total <- sum(tab)
        row_s <- rowSums(tab)
        col_s <- colSums(tab)
        exp_v <- as.numeric(outer(row_s, col_s) / total)
        obs_v <- as.numeric(tab)
        G     <- 2 * sum(ifelse(obs_v > 0L, obs_v * log(obs_v / exp_v), 0))
        df    <- (nrow(tab) - 1L) * (ncol(tab) - 1L)
        if (is.finite(G) && G > 0 && df > 0) lod[p] <- chisq_lod(G, df)
      }, error = function(e) NULL)

    } else {
      k       <- length(lvls)
      tab_arr <- array(0L, dim = c(2L, k, n_reps))
      for (ri in seq_len(n_reps)) {
        cf     <- haps_list[[ri]][case_idx_list[[ri]], p]
        kf     <- haps_list[[ri]][control_idx_list[[ri]], p]
        cf_raw <- as.integer(table(factor(cf, levels = lvls)))
        kf_raw <- as.integer(table(factor(kf, levels = lvls)))
        n_cf   <- sum(cf_raw); n_kf <- sum(kf_raw)
        tab_arr[1L, , ri] <- if (n_cf > seq_depth) as.integer(rmultinom(1L, seq_depth, cf_raw / n_cf)) else cf_raw
        tab_arr[2L, , ri] <- if (n_kf > seq_depth) as.integer(rmultinom(1L, seq_depth, kf_raw / n_kf)) else kf_raw
      }
      tryCatch({
        test <- suppressWarnings(mantelhaen.test(tab_arr, correct = FALSE))
        stat <- as.numeric(test$statistic)
        df   <- as.numeric(test$parameter)
        if (is.finite(stat) && is.finite(df) && stat > 0 && df > 0) {
          lod[p] <- chisq_lod(stat, df)
        }
      }, error = function(e) NULL)
    }
  }
  lod
}

# ============================================================
# BIPARENTAL F2 SIMULATION FUNCTIONS
# ============================================================

sim_f2_gamete <- function() {
  state <- sample(1L:2L, 1L)
  g     <- rep(state, N_POS)
  n_x   <- rpois(1L, 1.0)
  if (n_x > 0L) {
    cpts <- sort(runif(n_x, 0, 100))
    for (cp in cpts) {
      idx <- which(POS_SEQ >= cp)[1L]
      if (!is.na(idx) && idx <= N_POS) {
        state        <- 3L - state
        g[idx:N_POS] <- state
      }
    }
  }
  g
}

simulate_f2 <- function(n_individuals) {
  geno <- matrix(0L, nrow = n_individuals, ncol = N_POS)
  for (i in seq_len(n_individuals))
    geno[i, ] <- sim_f2_gamete()
  geno
}

simulate_phenotype_f2 <- function(geno, qtl_pos, qtl_eff) {
  pheno <- rnorm(nrow(geno), 0, 1)
  for (i in seq_along(qtl_pos)) {
    col_idx <- which.min(abs(POS_SEQ - qtl_pos[i]))
    pheno   <- pheno + (geno[, col_idx] - 1.5) * 2 * qtl_eff[i]
  }
  pheno
}

qtl_scan_biparental <- function(geno, pheno) {
  n   <- length(pheno)
  lod <- numeric(N_POS)
  ss0 <- sum((pheno - mean(pheno))^2)
  for (p in seq_len(N_POS)) {
    g <- geno[, p]
    if (length(unique(g)) < 2L) next
    tryCatch({
      fit   <- lm.fit(cbind(1, g), pheno)
      ss1   <- sum(fit$residuals^2)
      f_val <- ((ss0 - ss1) / 1) / (ss1 / (n - 2L))
      if (is.finite(f_val) && f_val > 0) {
        pv <- pf(f_val, 1L, n - 2L, lower.tail = FALSE)
        if (is.finite(pv) && pv > 0) lod[p] <- -log10(pv)
      }
    }, error = function(e) NULL)
  }
  lod
}

# ============================================================
# GWAS SIMULATION FUNCTIONS
# ============================================================

# Vectorized MAGIC-style mosaic simulation.
# Uses a Markov chain over SNP positions: at each step, each individual
# switches founder with probability p_switch = n_cross / (n_snps - 1).
# This is ~50x faster than the per-individual loop for large n_ind,
# making sample sizes up to 50 000 feasible.
simulate_gwas_human <- function(n_ind, n_snps, n_founders = 12L,
                                n_cross = 8L, maf_min = 0.05) {
  positions    <- seq(0, 100, length.out = n_snps)
  mafs         <- runif(n_snps, maf_min, 0.45)
  founder_haps <- matrix(
    rbinom(n_founders * n_snps, 1L, rep(mafs, each = n_founders)),
    nrow = n_founders, ncol = n_snps
  )

  p_switch <- n_cross / max(n_snps - 1L, 1L)

  geno      <- matrix(0L, nrow = n_ind, ncol = n_snps)
  block_vis <- matrix(0L, nrow = n_ind, ncol = n_snps)

  h1 <- sample.int(n_founders, n_ind, replace = TRUE)
  h2 <- sample.int(n_founders, n_ind, replace = TRUE)

  for (pos in seq_len(n_snps)) {
    if (pos > 1L) {
      sw1 <- which(runif(n_ind) < p_switch)
      sw2 <- which(runif(n_ind) < p_switch)
      if (length(sw1) > 0L) h1[sw1] <- sample.int(n_founders, length(sw1), replace = TRUE)
      if (length(sw2) > 0L) h2[sw2] <- sample.int(n_founders, length(sw2), replace = TRUE)
    }
    geno[, pos]      <- as.integer(founder_haps[h1, pos]) + as.integer(founder_haps[h2, pos])
    block_vis[, pos] <- h1
  }

  block_idx <- rep(seq_len(n_founders), each = ceiling(n_snps / n_founders))[seq_len(n_snps)]

  list(geno = geno, positions = positions, n_snps = n_snps,
       block_idx = block_idx, n_blocks = n_founders,
       n_founders = n_founders, block_vis = block_vis)
}

qtl_scan_gwas <- function(geno, pheno) {
  n      <- length(pheno)
  g_c    <- scale(geno)
  p_c    <- as.vector(scale(pheno))
  r      <- as.vector(crossprod(g_c, p_c)) / (n - 1L)
  r      <- pmax(pmin(r, 1 - 1e-10), -(1 - 1e-10))
  t_stat <- r * sqrt((n - 2L) / (1 - r^2))
  pv     <- 2 * pt(abs(t_stat), df = n - 2L, lower.tail = FALSE)
  ifelse(is.finite(pv) & pv > 0, -log10(pv), 0)
}

# ============================================================
# PLOT HELPER FUNCTIONS
# (Used by both renderPlot and downloadHandler so saved plots
#  always match what is displayed on screen.)
# ============================================================

make_mpp_cases_plot <- function(r) {
  nf           <- r$n_founders
  cols         <- get_founder_colors(nf)
  case_haps    <- r$obs_haps[r$disp_case, , drop = FALSE]
  control_haps <- r$obs_haps[r$disp_ctrl, , drop = FALSE]
  if (length(r$qtl_pos) > 0L && length(r$qtl_eff) > 0L) {
    qtl_col         <- which.min(abs(POS_SEQ - r$qtl_pos[1L]))
    effect_founders <- which(r$qtl_eff[[1L]] > 0)
    at_qtl          <- case_haps[, qtl_col]
    is_effect       <- at_qtl %in% effect_founders
    case_haps       <- case_haps[order(!is_effect, at_qtl), , drop = FALSE]
  }
  nc       <- nrow(case_haps)
  nk       <- nrow(control_haps)
  df_cases <- data.frame(
    position = rep(POS_SEQ, each = nc), ril = rep(seq_len(nc), times = N_POS),
    founder  = factor(as.vector(case_haps), levels = seq_len(nf)), group = "Cases")
  df_ctrl  <- data.frame(
    position = rep(POS_SEQ, each = nk), ril = rep(seq_len(nk), times = N_POS),
    founder  = factor(as.vector(control_haps), levels = seq_len(nf)), group = "Controls")
  df       <- rbind(df_cases, df_ctrl)
  df$group <- factor(df$group, levels = c("Cases", "Controls"))
  show_legend <- nf <= 16
  p <- ggplot(df, aes(x = position, y = ril, fill = founder)) +
    geom_raster() +
    facet_wrap(~group, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = cols, drop = FALSE,
      guide = if (show_legend) guide_legend(title = "Founder", nrow = ceiling(nf / 8)) else "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Position (cM)", y = "RIL") +
    theme_dark_mpp() +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_text(color = "#ffffff", face = "bold", size = 22),
          legend.position = if (show_legend) "bottom" else "none")
  if (length(r$qtl_pos) > 0)
    p <- p + geom_vline(xintercept = r$qtl_pos, color = "#ffffff", linewidth = 1.6)
  p
}

make_mpp_lod_plot <- function(r) {
  lod_df <- data.frame(position = POS_SEQ, lod = r$lod)
  y_max  <- max(r$lod_thresh * 1.5, max(r$lod) * 1.15, na.rm = TRUE)
  p <- ggplot(lod_df, aes(x = position, y = lod)) +
    geom_line(color = "#D0E8FF", linewidth = 1) +
    annotate("segment", x = 0, xend = 100, y = r$lod_thresh, yend = r$lod_thresh,
             linetype = "dashed", color = "#64B5F6", linewidth = 0.7) +
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(x = "Position (cM)", y = "LOD") +
    theme_dark_mpp() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#222222", linewidth = 0.4),
          plot.background  = element_rect(fill = BG_PLOT_LIGHT, color = NA),
          panel.background = element_rect(fill = BG_PLOT_LIGHT, color = NA))
  if (length(r$qtl_pos) > 0)
    p <- p + geom_vline(xintercept = r$qtl_pos, color = "#ffffff", linewidth = 1.6)
  p
}

make_mpp_freq_plot <- function(r) {
  nf           <- r$n_founders
  cols         <- get_founder_colors(nf)
  nc           <- length(r$case_idx)
  nk           <- length(r$control_idx)
  case_haps    <- r$obs_haps[r$case_idx,    , drop = FALSE]
  control_haps <- r$obs_haps[r$control_idx, , drop = FALSE]
  rows <- vector("list", nf * N_POS)
  k    <- 1L
  for (f in seq_len(nf)) {
    for (pos in seq_len(N_POS)) {
      rows[[k]] <- data.frame(
        position = POS_SEQ[pos],
        founder  = factor(f, levels = seq_len(nf)),
        diff     = (if (nc > 0) sum(case_haps[, pos] == f) / nc else 0) -
                   (if (nk > 0) sum(control_haps[, pos] == f) / nk else 0)
      )
      k <- k + 1L
    }
  }
  df <- do.call(rbind, rows)
  p <- ggplot(df, aes(x = position, y = diff, color = founder)) +
    geom_line(linewidth = 1.1, alpha = 1) +
    annotate("segment", x = 0, xend = 100, y = 0, yend = 0,
             color = "#666666", linewidth = 0.7) +
    scale_color_manual(values = cols, drop = FALSE,
      guide = if (nf <= 16) guide_legend(title = "Founder", nrow = ceiling(nf / 8)) else "none") +
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.01))) +
    labs(x = "Position (cM)", y = "Freq diff") +
    theme_dark_mpp() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#222222", linewidth = 0.4),
          legend.position  = if (nf <= 16) "bottom" else "none",
          plot.background  = element_rect(fill = BG_PLOT_LIGHT, color = NA),
          panel.background = element_rect(fill = BG_PLOT_LIGHT, color = NA))
  if (length(r$qtl_pos) > 0)
    p <- p + geom_vline(xintercept = r$qtl_pos, color = "#ffffff", linewidth = 1.6)
  p
}

make_f2_hap_plot <- function(r) {
  idx <- r$hap_idx
  nd  <- length(idx)
  df  <- data.frame(
    position = rep(POS_SEQ, each = nd),
    ind      = rep(seq_len(nd), times = N_POS),
    geno     = factor(as.vector(r$geno[idx, ]), levels = 1L:2L)
  )
  p <- ggplot(df, aes(x = position, y = ind, fill = geno)) +
    geom_raster() +
    scale_fill_manual(values = F2_COLS, labels = F2_LBLS, name = "Allele",
                      guide = guide_legend(direction = "horizontal")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Position (cM)", y = "F2 individual") +
    theme_dark_mpp() +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), legend.position = "bottom")
  if (length(r$qtl_pos) > 0)
    p <- p + geom_vline(xintercept = r$qtl_pos, color = "#ffffff", linewidth = 1.6)
  p
}

make_f2_lod_plot <- function(r) {
  y_max <- max(r$lod_thresh * 1.5, max(r$lod) * 1.15, na.rm = TRUE)
  df    <- data.frame(position = POS_SEQ, lod = r$lod)
  p <- ggplot(df, aes(x = position, y = lod)) +
    geom_line(color = "#D0E8FF", linewidth = 1) +
    annotate("segment", x = 0, xend = 100, y = r$lod_thresh, yend = r$lod_thresh,
             linetype = "dashed", color = "#64B5F6", linewidth = 0.7) +
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(x = "Position (cM)", y = "LOD") +
    theme_dark_mpp() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#222222", linewidth = 0.4),
          plot.background  = element_rect(fill = BG_PLOT_LIGHT, color = NA),
          panel.background = element_rect(fill = BG_PLOT_LIGHT, color = NA))
  if (length(r$qtl_pos) > 0)
    p <- p + geom_vline(xintercept = r$qtl_pos, color = "#ffffff", linewidth = 1.6)
  p
}

make_gwas_hap_plot <- function(r) {
  cols      <- get_founder_colors(r$n_founders)
  n_disp    <- min(50L, length(r$case_idx))
  disp_case <- r$case_idx[seq_len(n_disp)]
  disp_ctrl <- r$ctrl_idx[seq_len(n_disp)]
  if (!anyNA(r$qtl_snp) && !is.null(r$qtl_snp)) {
    causal_col <- r$qtl_snp[1L]
    disp_case  <- disp_case[order(r$block_vis[disp_case, causal_col])]
    disp_ctrl  <- disp_ctrl[order(r$block_vis[disp_ctrl, causal_col])]
  }
  nc       <- length(disp_case); nk <- length(disp_ctrl)
  df_cases <- data.frame(
    position = rep(r$positions, each = nc),
    ind      = rep(seq_len(nc), times = r$n_snps),
    founder  = factor(as.vector(r$block_vis[disp_case, ]), levels = seq_len(r$n_founders)),
    group    = "Cases")
  df_ctrl  <- data.frame(
    position = rep(r$positions, each = nk),
    ind      = rep(seq_len(nk), times = r$n_snps),
    founder  = factor(as.vector(r$block_vis[disp_ctrl, ]), levels = seq_len(r$n_founders)),
    group    = "Controls")
  df       <- rbind(df_cases, df_ctrl)
  df$group <- factor(df$group, levels = c("Cases", "Controls"))
  p <- ggplot(df, aes(x = position, y = ind, fill = founder)) +
    geom_raster() +
    facet_wrap(~group, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = cols, drop = FALSE, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Position (cM)", y = "Individual") +
    theme_dark_mpp() +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_text(color = "#ffffff", size = 22))
  if (!anyNA(r$qtl_snp) && !is.null(r$qtl_snp))
    p <- p + geom_vline(xintercept = r$positions[r$qtl_snp],
                        color = "#ffffff", linewidth = 0.6, alpha = 0.9)
  p
}

make_gwas_manhattan <- function(r) {
  block_col <- ifelse(r$block_idx %% 2L == 0L, "#80D8FF", "#4FC3F7")
  df        <- data.frame(position = r$positions, lod = r$lod, col = block_col)
  y_max     <- max(r$lod_thresh * 1.5, max(r$lod) * 1.15, na.rm = TRUE)
  p <- ggplot(df, aes(x = position, y = lod, color = col)) +
    geom_point(size = 2.0, alpha = 0.9) +
    annotate("segment", x = 0, xend = 100, y = r$lod_thresh, yend = r$lod_thresh,
             linetype = "dashed", color = "#64B5F6", linewidth = 0.7) +
    scale_color_identity(guide = "none") +
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(x = "Position (cM)", y = "LOD") +
    theme_dark_mpp() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#222222", linewidth = 0.4),
          plot.background  = element_rect(fill = BG_PLOT_LIGHT, color = NA),
          panel.background = element_rect(fill = BG_PLOT_LIGHT, color = NA))
  if (!anyNA(r$qtl_snp) && !is.null(r$qtl_snp))
    p <- p + geom_vline(xintercept = r$positions[r$qtl_snp],
                        color = "#ffffff", linewidth = 0.6, alpha = 0.9)
  p
}

# ============================================================
# UI
# ============================================================

ui <- navbarPage(
  title = "Cross Examination",
  tags$head(tags$style(HTML(paste0("

    /* Navbar */
    .navbar, .navbar-default {
      background-color: ", BG_PANEL, " !important;
      border-color: #253050 !important;
    }
    .navbar-default .navbar-brand { color: #90CAF9 !important; font-size: 17px !important; }
    .navbar-default .navbar-nav > li > a { color: #ffffff !important; font-size: 15px !important; }
    .navbar-default .navbar-nav > .active > a,
    .navbar-default .navbar-nav > .active > a:focus,
    .navbar-default .navbar-nav > .active > a:hover {
      background-color: #2E75B6 !important; color: #ffffff !important;
    }
    .navbar-default .navbar-nav > li > a:hover { background-color: #253050 !important; }
    .tab-content { background-color: ", BG_PAGE, " !important; }

    /* ===================================================
       FULL DARK MODE
       =================================================== */

    html, body, .container-fluid {
      background-color: ", BG_PAGE, " !important;
      color: #ffffff !important;
      font-family: Arial, sans-serif;
      font-size: 15px;
    }

    .container-fluid > h2, h2, h3 { color: #ffffff !important; }
    h4 { color: #ffffff; margin-bottom: 4px; margin-top: 12px; font-weight: normal !important; }
    hr { margin: 8px 0; border-color: #253050; }

    .well {
      background: ", BG_PANEL, " !important;
      border: 1px solid #253050 !important;
      box-shadow: none;
      color: #ffffff !important;
    }

    label, .control-label { color: #ffffff !important; font-size: 15px !important; }
    p small { color: #90CAF9; font-size: 13px; }

    .well { font-size: 15px !important; }
    .well h4 { font-size: 17px !important; }

    .form-control, input, textarea, select {
      background-color: #0d0d0d !important;
      color: #ffffff !important;
      border: 1px solid #253050 !important;
      box-shadow: none !important;
    }
    .form-control:focus, input:focus {
      border-color: #4FC3F7 !important;
      box-shadow: 0 0 0 2px rgba(79,195,247,0.2) !important;
    }

    .selectize-input {
      background: #141414 !important;
      color: #ffffff !important;
      border: 2px solid #333333 !important;
      box-shadow: none !important;
      font-size: 15px !important;
      padding: 7px 10px !important;
    }
    .selectize-input.focus {
      border-color: #4FC3F7 !important;
      box-shadow: 0 0 0 2px rgba(79,195,247,0.2) !important;
    }
    .selectize-input > input { color: #ffffff !important; font-size: 15px !important; }
    .selectize-dropdown,
    .selectize-dropdown-content {
      background: #141414 !important;
      border: 2px solid #333333 !important;
      color: #ffffff !important;
      font-size: 15px !important;
    }
    .selectize-dropdown .option { color: #ffffff !important; padding: 7px 10px; font-size: 15px !important; }
    .selectize-dropdown .option:hover,
    .selectize-dropdown .option.active,
    .selectize-dropdown .option.selected {
      background: #2E75B6 !important;
      color: #ffffff !important;
    }
    .selectize-control.single .selectize-input:after {
      border-top-color: #666666 !important;
    }

    .irs--shiny .irs-line  { background: #333333; border-color: #333333; height: 8px; border-radius: 4px; }
    .irs--shiny .irs-bar   { background: #2E75B6; border-color: #2E75B6; height: 8px; }
    .irs--shiny .irs-from,
    .irs--shiny .irs-to,
    .irs--shiny .irs-single { background: #2E75B6; color: #ffffff !important; font-size: 13px !important; }
    .irs--shiny .irs-handle { background: #4FC3F7 !important; border: 2px solid #ffffff !important; width: 20px !important; height: 20px !important; top: 21px !important; }
    .irs--shiny .irs-min,
    .irs--shiny .irs-max   { color: #aaaaaa !important; font-size: 12px !important; background: transparent !important; }
    .irs--shiny .irs-grid-text { display: none !important; }
    .irs--shiny .irs-grid-pol { display: none !important; }

    #simulate {
      background: #2E75B6 !important;
      color: #fff !important;
      border: none !important;
    }
    #simulate:hover { background: #4194D4 !important; }
    #simulate:active { background: #1a5a96 !important; }

    .col-sm-9, .col-sm-3 { background: ", BG_PAGE, " !important; }

    .shiny-plot-output {
      background: ", BG_PLOT, " !important;
      border-radius: 4px;
    }

    .shiny-output-error, .shiny-output-error-message { color: #EF9A9A !important; }
    .shiny-notification {
      background: #16213e !important;
      color: #ffffff !important;
      border: 1px solid #253050 !important;
    }

    ::-webkit-scrollbar { width: 6px; height: 6px; }
    ::-webkit-scrollbar-track { background: ", BG_PAGE, "; }
    ::-webkit-scrollbar-thumb { background: #2E75B6; border-radius: 3px; }
    ::-webkit-scrollbar-thumb:hover { background: #4FC3F7; }

    /* ── Teaching mode box ── */
    .teach-box {
      background: rgba(46, 117, 182, 0.12);
      border: 1px solid #2E75B6;
      border-radius: 4px;
      padding: 7px 10px 5px 10px;
      margin-bottom: 6px;
    }
    .teach-box label { color: #80CFFF !important; font-size: 14px !important; }

    /* ── Download button ── */
    .shiny-download-link {
      display: block !important;
      width: 100% !important;
      text-align: center !important;
      background: #1a5a96 !important;
      color: #fff !important;
      border: none !important;
      padding: 8px !important;
      border-radius: 4px !important;
      font-size: 14px !important;
      margin-top: 6px !important;
      text-decoration: none !important;
      box-sizing: border-box !important;
    }
    .shiny-download-link:hover {
      background: #2E75B6 !important;
      color: #fff !important;
      text-decoration: none !important;
    }
  ")))),

  # ── Tab 1: MPP / BSA ─────────────────────────────────────
  tabPanel("MPP / BSA",
  sidebarLayout(
    sidebarPanel(
      width = 3,

      div(class = "teach-box",
        checkboxInput("teaching_mode",
          "Teaching Mode (fixed QTL positions)", value = FALSE)
      ),

      h4("Population Design"),
      selectInput("design", "Crossing design",
        choices = c("Fully intercrossed" = "magic",
                    "Hub-and-spoke"       = "hub"),
        width = "100%"),
      sliderInput("n_founders", "Founders",
                  min=2, max=30, value=8, step=2, width="100%"),
      sliderInput("n_gen", "Generations of recombination / selfing",
                  min=1, max=50, value=10, step=1, width="100%"),
      conditionalPanel(
        condition = "input.design == 'hub'",
        p(tags$small(style="color:#FFCC80;",
          "Hub-and-spoke crosses one inbred hub line to each spoke founder separately. ",
          "The F1 is selfed for several generations to produce inbred RILs. ",
          "Each RIL carries a mosaic of hub and one spoke's alleles. ",
          "Phenotyping is done on the RILs directly. 5-8 generations is a typical range. ",
          tags$br(),
          "Power note: each spoke allele appears in only ~1/(founders-1) of all RILs. ",
          "This dilutes BSA signal compared to fully intercrossed. ",
          "More generations = smaller haplotype blocks = sharper QTL peaks. ",
          "Family structure (RILs from the same spoke share more ancestry) is not simulated here."
        ))
      ),

      hr(),
      h4("Pool"),
      p(tags$small("Cases = top N by phenotype. Controls = random N from unselected.")),
      sliderInput("pool_size", "Individuals per pool (N)",
                  min=100, max=2000, value=200, step=100, width="100%"),
      sliderInput("n_reps", "Replicates (fresh phenotype draw per rep)",
                  min=1, max=10, value=1, step=1, width="100%"),
      p(tags$small(
        "Each replicate re-draws phenotype noise and re-selects pools from the same RILs. ",
        "CMH combines replicates: the QTL signal is consistent across reps (same genetics) so it stacks up, ",
        "while non-QTL noise is random across reps and averages out. ",
        "LOD at true QTL peaks rises and scan baseline variance falls as replicates increase."
      )),

      hr(),
      h4("Genetic Model"),
      selectInput("genetic_model", "QTL architecture",
        choices = c("Null - no QTL" = "null",
                    "1 QTL"         = "1qtl",
                    "2 QTLs"        = "2qtl",
                    "3 QTLs"        = "3qtl",
                    "Polygenic"     = "polygenic"),
        width = "100%"),

      hr(),
      h4("Significance Threshold"),
      p(tags$small("Bonferroni (202 tests, α=0.05): LOD 3.6")),

      br(),
      actionButton("simulate", "Simulate",
                   width = "100%",
                   style = "background:#2E75B6; color:#fff; border:none;
                            padding:8px; border-radius:4px; font-size:14px;"),
      downloadButton("download_mpp", "Save Plot")
    ),

    mainPanel(
      width = 9,

      h4("Founder Haplotypes"),
      p(tags$small("Each row is one founder: pure inbred, no recombination.")),
      plotOutput("founder_plot", height = "160px"),

      h4("RIL Haplotype Mosaic"),
      p(tags$small("A random sample of RILs. Colors = founder ancestry across 100 cM.")),
      plotOutput("hap_plot", height = "300px"),

      h4("Cases vs. Controls"),
      p(tags$small(
        "Top N = cases (top); random N from unselected = controls (bottom). ",
        "Under a QTL peak, one founder color dominates in cases."
      )),
      plotOutput("cases_plot", height = "360px"),

      h4("QTL Scan"),
      p(tags$small(
        "Fisher's exact test (1 rep) or CMH (multiple reps) on founder counts, smoothed with 5-position rolling average. ",
        tags$span(style="color:#ffffff; font-weight:bold;", "White line"),
        " = true QTL.  ",
        tags$span(style="color:#64B5F6;", "Blue dashed"),
        " = LOD threshold."
      )),
      plotOutput("lod_plot", height = "200px"),

      h4("Founder Frequency: Cases − Controls"),
      p(tags$small("Difference in founder frequency between pools at each position.")),
      plotOutput("freq_plot", height = "280px")
    )
  )
  ), # end tabPanel MPP / BSA

  # ── Tab 2: Biparental QTL (F2) ───────────────────────────
  tabPanel("Biparental QTL",
    sidebarLayout(
      sidebarPanel(
        width = 3,

        div(class = "teach-box",
          checkboxInput("f2_teaching",
            "Teaching Mode (fixed QTL positions)", value = FALSE)
        ),

        h4("Cross Design"),
        p(tags$small(
          "Two inbred or distinct parents (P1 × P2) are crossed to make F1 offspring. ",
          "F1 individuals are crossed to produce F2 offspring. ",
          "Each F2 chromosome is a mosaic of P1 and P2 segments created by recombination."
        )),
        sliderInput("f2_n_ind", "F2 individuals",
                    min=100, max=2000, value=200, step=100, width="100%"),

        hr(),
        h4("Genetic Model"),
        selectInput("f2_model", "QTL architecture",
          choices = c("Null - no QTL" = "null",
                      "1 QTL"         = "1qtl",
                      "2 QTLs"        = "2qtl",
                      "3 QTLs"        = "3qtl",
                      "Polygenic"     = "polygenic"),
          width = "100%"),

        hr(),
        h4("Significance Threshold"),
        p(tags$small("Bonferroni (202 tests, α=0.05): LOD 3.6")),
        br(),
        actionButton("f2_simulate", "Simulate",
                     width="100%",
                     style="background:#2E75B6; color:#fff; border:none;
                            padding:8px; border-radius:4px; font-size:14px;"),
        downloadButton("download_f2", "Save Plot")
      ),
      mainPanel(
        width = 9,
        h4("Parents"),
        p(tags$small("Parent 1 (red) and Parent 2 (blue). Both are fully inbred — every locus is fixed.")),
        plotOutput("f2_parents_plot", height="70px"),

        h4("F2 Haplotype Mosaic"),
        p(tags$small(
          "Each row is one F2 individual. ",
          "Colors show which parent's allele is present at each position. ",
          "Recombination during meiosis creates the alternating red/blue blocks."
        )),
        plotOutput("f2_hap_plot", height="280px"),

        h4("QTL Scan"),
        p(tags$small(
          "Additive regression (F-test) at each position, smoothed with 5-position rolling average. ",
          tags$span(style="color:#ffffff; font-weight:bold;", "White line"),
          " = true QTL.  ",
          tags$span(style="color:#64B5F6;", "Blue dashed"),
          " = LOD threshold."
        )),
        plotOutput("f2_lod_plot", height="200px"),
      )
    )
  ), # end tabPanel Biparental QTL

  # ── Tab 3: GWAS ──────────────────────────────────────────
  tabPanel("GWAS",
    sidebarLayout(
      sidebarPanel(
        width = 3,

        div(class = "teach-box",
          checkboxInput("gwas_teaching",
            "Teaching Mode (fixed QTL positions)", value = FALSE)
        ),

        h4("Population"),
        p(tags$small(
          "Outbred population. 1000 SNPs across 100 cM (10 SNPs/cM). ",
          "Each chromosome is divided into large LD blocks. ",
          "All SNPs within a block are co-inherited — a causal SNP pulls up every SNP in its block."
        )),
        sliderInput("gwas_n_ind", "Individuals",
                    min=500, max=100000, value=10000, step=500, width="100%"),
        sliderInput("gwas_block_cM", "Avg. LD block size (cM)",
                    min=5, max=50, value=12, step=1, width="100%"),
        p(tags$small(
          "Small blocks = diverse outbred population with long recombination history. ",
          "Large blocks = bottlenecked, domesticated, or recently hybridized species."
        )),

        hr(),
        h4("Genetic Model"),
        selectInput("gwas_model", "Genetic architecture",
          choices = c("Null - no QTL" = "null",
                      "1 causal SNP"  = "1qtl",
                      "2 causal SNPs" = "2qtl",
                      "3 causal SNPs" = "3qtl",
                      "Polygenic"     = "polygenic"),
          width="100%"),

        hr(),
        h4("Significance Threshold"),
        p(tags$small("Bonferroni (α = 0.05, 1000 tests): LOD 4.3")),
        br(),
        actionButton("gwas_simulate", "Simulate",
                     width="100%",
                     style="background:#2E75B6; color:#fff; border:none;
                            padding:8px; border-radius:4px; font-size:14px;"),
        downloadButton("download_gwas", "Save Plot")
      ),
      mainPanel(
        width = 9,
        h4("Cases vs. Controls"),
        p(tags$small(
          "Top N by phenotype = cases (top); random N from unselected = controls (bottom). ",
          "Colors = LD block haplotype. Under a causal block, one color dominates in cases."
        )),
        plotOutput("gwas_hap_plot", height="360px"),

        h4("Manhattan Plot"),
        p(tags$small(
          "Chi-square test of allele frequency difference (cases vs. controls) at each of 1000 SNPs. LOD = −log₁₀(p). ",
          tags$span(style="color:#ffffff; font-weight:bold;", "White line"),
          " = true causal SNP.  ",
          tags$span(style="color:#64B5F6;", "Blue dashed"),
          " = LOD threshold. ",
          "All SNPs in a causal block rise together due to LD."
        )),
        plotOutput("gwas_manhattan", height="220px")
      )
    )
  ) # end tabPanel GWAS

) # end navbarPage

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  results <- eventReactive(input$simulate, {
    pool_size <- input$pool_size
    n_total   <- 3L * pool_size
    n_display <- min(50L, pool_size)
    hub       <- input$design == "hub"
    n_reps    <- input$n_reps
    teach     <- isTRUE(input$teaching_mode)
    if (teach) set.seed(42L)
    qtl_cfg   <- get_qtl_config(input$genetic_model, input$n_founders,
                                hub_design = hub, teaching = teach)

    founder_haps <- matrix(rep(seq_len(input$n_founders), each = N_POS),
                           nrow = input$n_founders, ncol = N_POS)

    obs_haps <- simulate_rils(input$n_founders, input$n_gen, n_total, input$design)

    haps_list     <- vector("list", n_reps)
    case_idx_list <- vector("list", n_reps)
    ctrl_idx_list <- vector("list", n_reps)

    for (rep_i in seq_len(n_reps)) {
      pheno      <- simulate_phenotype(obs_haps, qtl_cfg)
      ranked     <- order(pheno, decreasing = TRUE)
      ci         <- ranked[seq_len(pool_size)]
      unselected <- ranked[(pool_size + 1L):n_total]
      ki         <- sample(unselected, pool_size)
      haps_list[[rep_i]]     <- obs_haps
      case_idx_list[[rep_i]] <- ci
      ctrl_idx_list[[rep_i]] <- ki
    }

    raw_lod <- qtl_scan(haps_list, case_idx_list, ctrl_idx_list)
    lod     <- smooth_lod(raw_lod, window = 5L)

    ci1 <- case_idx_list[[1L]]
    ki1 <- ctrl_idx_list[[1L]]
    disp_all  <- sample.int(n_total, min(n_display * 2L, n_total))
    disp_case <- ci1[seq_len(n_display)]
    disp_ctrl <- ki1[seq_len(n_display)]

    list(
      obs_haps     = obs_haps,
      founder_haps = founder_haps,
      lod          = lod,
      qtl_pos      = qtl_cfg$pos,
      qtl_eff      = qtl_cfg$eff,
      n_founders   = input$n_founders,
      pool_size    = pool_size,
      n_reps       = n_reps,
      n_display    = n_display,
      lod_thresh   = 3.6,
      case_idx     = ci1,
      control_idx  = ki1,
      disp_all     = disp_all,
      disp_case    = disp_case,
      disp_ctrl    = disp_ctrl,
      teaching     = teach,
      design       = input$design,
      n_founders_in = input$n_founders,
      n_gen        = input$n_gen,
      genetic_model = input$genetic_model
    )
  })

  output$founder_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)
    df <- data.frame(
      position = rep(POS_SEQ, each = nf),
      founder  = rep(seq_len(nf), times = N_POS),
      fill_f   = factor(rep(seq_len(nf), times = N_POS), levels = seq_len(nf))
    )
    ggplot(df, aes(x = position, y = founder, fill = fill_f)) +
      geom_raster() +
      scale_fill_manual(values = cols, drop = FALSE, guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
      scale_y_continuous(expand = c(0, 0),
        breaks = if (nf <= 20) seq_len(nf) else c(1, seq(10, nf, by=10)),
        labels = if (nf <= 20) paste0("F", seq_len(nf)) else
                   c("F1", paste0("F", seq(10, nf, by=10)))) +
      labs(x = "Position (cM)", y = "Founder") +
      theme_dark_mpp() +
      theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9))
  }, bg = BG_PLOT)

  output$hap_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)
    disp_haps <- r$obs_haps[r$disp_all, , drop = FALSE]
    nd <- nrow(disp_haps)
    df <- data.frame(
      position = rep(POS_SEQ, each = nd),
      ril      = rep(seq_len(nd), times = N_POS),
      founder  = factor(as.vector(disp_haps), levels = seq_len(nf))
    )
    show_legend <- nf <= 16
    ggplot(df, aes(x = position, y = ril, fill = founder)) +
      geom_raster() +
      scale_fill_manual(values = cols, drop = FALSE,
        guide = if (show_legend) guide_legend(title="Founder", nrow=ceiling(nf/8)) else "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Position (cM)", y = "RIL") +
      theme_dark_mpp() +
      theme(panel.grid = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = if (show_legend) "bottom" else "none")
  }, bg = BG_PLOT)

  output$cases_plot <- renderPlot({
    req(results()); make_mpp_cases_plot(results())
  }, bg = BG_PLOT)

  output$lod_plot <- renderPlot({
    req(results()); make_mpp_lod_plot(results())
  }, bg = BG_PLOT_LIGHT)

  output$freq_plot <- renderPlot({
    req(results()); make_mpp_freq_plot(results())
  }, bg = BG_PLOT_LIGHT)

  # MPP download handler
  output$download_mpp <- downloadHandler(
    filename = function() paste0("mpp_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      req(results())
      r  <- results()
      p1 <- make_mpp_cases_plot(r)
      p2 <- make_mpp_lod_plot(r)
      p3 <- make_mpp_freq_plot(r)
      tm <- if (isTRUE(r$teaching)) "On" else "Off"
      design_label <- switch(r$design,
        "magic" = "Fully intercrossed", "hub" = "Hub-and-spoke", r$design)
      caption <- paste0(
        "Design: ", design_label,
        "  |  Founders: ", r$n_founders_in,
        "  |  Generations: ", r$n_gen,
        "  |  Pool: ", r$pool_size,
        "  |  Replicates: ", r$n_reps,
        "  |  Model: ", r$genetic_model,
        "  |  LOD threshold: 3.6 (Bonferroni)",
        "  |  Teaching mode: ", tm
      )
      png(file, width = 2000, height = 2200, res = 150, bg = BG_PLOT)
      gridExtra::grid.arrange(
        p1, p2, p3, ncol = 1,
        bottom = grid::textGrob(caption,
          gp = grid::gpar(fontsize = 10, col = "#aaaaaa"), hjust = 0.5)
      )
      dev.off()
    }
  )

  # ── F2 / Biparental QTL server ─────────────────────────

  f2_results <- eventReactive(input$f2_simulate, {
    teach   <- isTRUE(input$f2_teaching)
    if (teach) set.seed(42L)
    qtl_pos <- switch(input$f2_model,
      "null"      = integer(0),
      "1qtl"      = if (teach) TEACH_POS[["1qtl"]] else sample(10L:90L, 1L),
      "2qtl"      = if (teach) TEACH_POS[["2qtl"]] else sort(sample(10L:90L, 2L)),
      "3qtl"      = if (teach) TEACH_POS[["3qtl"]] else sort(sample(10L:90L, 3L)),
      "polygenic" = if (teach) TEACH_POS[["polygenic"]] else sort(sample(5L:95L, 10L))
    )
    qtl_eff <- switch(input$f2_model,
      "null"      = numeric(0),
      "1qtl"      = 0.15,
      "2qtl"      = c(0.12, 0.12),
      "3qtl"      = c(0.10, 0.10, 0.10),
      "polygenic" = rep(0.06, 10)
    )
    geno    <- simulate_f2(input$f2_n_ind)
    nd      <- min(100L, nrow(geno))
    hap_idx <- sample.int(nrow(geno), nd)
    pheno   <- simulate_phenotype_f2(geno, qtl_pos, qtl_eff)
    raw     <- qtl_scan_biparental(geno, pheno)
    lod     <- smooth_lod(raw, window = 5L)
    list(geno = geno, pheno = pheno, lod = lod, qtl_pos = qtl_pos,
         lod_thresh = 3.6, hap_idx = hap_idx,
         teaching = teach, n_ind = input$f2_n_ind, model = input$f2_model)
  })

  output$f2_parents_plot <- renderPlot({
    req(f2_results())
    df <- data.frame(
      position = rep(POS_SEQ, 2L),
      ind      = rep(c(2L, 1L), each = N_POS),
      geno     = factor(rep(c(1L, 2L), each = N_POS), levels = 1L:2L)
    )
    ggplot(df, aes(x = position, y = ind, fill = geno)) +
      geom_raster() +
      scale_fill_manual(values = F2_COLS, labels = F2_LBLS, name = NULL,
                        guide = guide_legend(direction = "horizontal")) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
      scale_y_continuous(breaks = c(1L, 2L), labels = c("Parent 2", "Parent 1"),
                         expand = c(0.2, 0)) +
      labs(x = NULL, y = NULL) +
      theme_dark_mpp(base_size = 13) +
      theme(panel.grid = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), legend.position = "right",
            legend.text = element_text(size = 11))
  }, bg = BG_PLOT)

  output$f2_hap_plot <- renderPlot({
    req(f2_results()); make_f2_hap_plot(f2_results())
  }, bg = BG_PLOT)

  output$f2_lod_plot <- renderPlot({
    req(f2_results()); make_f2_lod_plot(f2_results())
  }, bg = BG_PLOT_LIGHT)

  # F2 download handler
  output$download_f2 <- downloadHandler(
    filename = function() paste0("biparental_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      req(f2_results())
      r  <- f2_results()
      p1 <- make_f2_hap_plot(r)
      p2 <- make_f2_lod_plot(r)
      tm <- if (isTRUE(r$teaching)) "On" else "Off"
      caption <- paste0(
        "F2 individuals: ", r$n_ind,
        "  |  Model: ", r$model,
        "  |  LOD threshold: 3.6 (Bonferroni)",
        "  |  Teaching mode: ", tm
      )
      png(file, width = 2000, height = 1400, res = 150, bg = BG_PLOT)
      gridExtra::grid.arrange(
        p1, p2, ncol = 1,
        bottom = grid::textGrob(caption,
          gp = grid::gpar(fontsize = 10, col = "#aaaaaa"), hjust = 0.5)
      )
      dev.off()
    }
  )

  # ── GWAS server ────────────────────────────────────────

  gwas_results <- eventReactive(input$gwas_simulate, {
    withProgress(message = "Simulating GWAS...", value = 0, {
      n_snps   <- 1000L
      gwas_ef  <- 0.5
      n_ind    <- input$gwas_n_ind
      teach    <- isTRUE(input$gwas_teaching)
      if (teach) set.seed(42L)
      block_cM <- input$gwas_block_cM
      n_cross  <- max(1L, round(100L / block_cM))

      setProgress(0.2, detail = "Generating population")
      pop  <- simulate_gwas_human(n_ind, n_snps, n_cross = n_cross)
      geno <- pop$geno
      nb   <- pop$n_blocks

      pick_snp_in_block <- function(b) sample(which(pop$block_idx == b), 1L)

      qtl_snp <- NA_integer_
      pheno   <- rnorm(n_ind)

      setProgress(0.5, detail = "Applying genetic model")

      if (input$gwas_model == "1qtl") {
        qtl_snp <- if (teach) TEACH_SNP[["1qtl"]] else pick_snp_in_block(sample.int(nb, 1L))
        g       <- geno[, qtl_snp]
        pheno   <- rnorm(n_ind) + (g - mean(g)) * gwas_ef

      } else if (input$gwas_model == "2qtl") {
        qtl_snp <- if (teach) TEACH_SNP[["2qtl"]] else sapply(sample.int(nb, 2L), pick_snp_in_block)
        for (qs in qtl_snp) { g <- geno[, qs]; pheno <- pheno + (g - mean(g)) * gwas_ef * 0.8 }

      } else if (input$gwas_model == "3qtl") {
        qtl_snp <- if (teach) TEACH_SNP[["3qtl"]] else sapply(sample.int(nb, min(3L, nb)), pick_snp_in_block)
        for (qs in qtl_snp) { g <- geno[, qs]; pheno <- pheno + (g - mean(g)) * gwas_ef * 0.7 }

      } else if (input$gwas_model == "polygenic") {
        qtl_snp <- if (teach) TEACH_SNP[["polygenic"]] else sort(sample(round(0.05*n_snps):round(0.95*n_snps), 10L))
        for (cs in qtl_snp) { g <- geno[, cs]; pheno <- pheno + (g - mean(g)) * gwas_ef * 0.3 }
      }

      pool_size <- min(200L, floor(n_ind / 3L))
      ranked    <- order(pheno, decreasing = TRUE)
      case_idx  <- ranked[seq_len(pool_size)]
      ctrl_idx  <- sample(ranked[(pool_size + 1L):n_ind], pool_size)

      pheno_cc            <- rep(0L, n_ind)
      pheno_cc[case_idx]  <- 1L

      setProgress(0.75, detail = "Running association scan")
      lod <- qtl_scan_gwas(geno, pheno_cc)
      setProgress(1.0)

      list(lod = lod, positions = pop$positions, qtl_snp = qtl_snp,
           n_snps = n_snps, block_idx = pop$block_idx, n_blocks = nb,
           n_founders = pop$n_founders, block_vis = pop$block_vis,
           case_idx = case_idx, ctrl_idx = ctrl_idx,
           lod_thresh = 4.3, teaching = teach,
           n_ind = n_ind, model = input$gwas_model,
           block_cM = block_cM)
    })
  })

  output$gwas_hap_plot <- renderPlot({
    req(gwas_results()); make_gwas_hap_plot(gwas_results())
  }, bg = BG_PLOT)

  output$gwas_manhattan <- renderPlot({
    req(gwas_results()); make_gwas_manhattan(gwas_results())
  }, bg = BG_PLOT_LIGHT)

  # GWAS download handler
  output$download_gwas <- downloadHandler(
    filename = function() paste0("gwas_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      req(gwas_results())
      r  <- gwas_results()
      p1 <- make_gwas_hap_plot(r)
      p2 <- make_gwas_manhattan(r)
      tm <- if (isTRUE(r$teaching)) "On" else "Off"
      caption <- paste0(
        "Individuals: ", r$n_ind,
        "  |  LD block: ", r$block_cM, " cM",
        "  |  Model: ", r$model,
        "  |  LOD threshold: ", r$lod_thresh,
        "  |  Teaching mode: ", tm
      )
      png(file, width = 2000, height = 1600, res = 150, bg = BG_PLOT)
      gridExtra::grid.arrange(
        p1, p2, ncol = 1,
        bottom = grid::textGrob(caption,
          gp = grid::gpar(fontsize = 10, col = "#aaaaaa"), hjust = 0.5)
      )
      dev.off()
    }
  )

}

shinyApp(ui, server)

#' Run all statistics defined by the YAML `stats_report:` block
#'
#' Iterates over `comparisons:`, `correlations:`, and `linear_models:` and
#' returns a list of three tidy data frames suitable for writing to xlsx
#' tabs.
#'
#' Method auto-selection for comparisons:
#' \itemize{
#'   \item 2 levels  → t-test + Wilcoxon
#'   \item 3+ levels → ANOVA + Kruskal-Wallis (+ Tukey HSD + pairwise Wilcoxon post-hoc)
#' }
#'
#' @param de A `struct::DatasetExperiment` (from `run_MRManalyzeR()` or `readRDS()`).
#' @param st_params The `stats_report:` block from the YAML (a list).
#' @return A list with elements `stats`, `correlations`, `linear_models`,
#'   each a data frame (possibly with 0 rows if the corresponding YAML
#'   block was empty).
#' @export
run_stats = function(de, st_params){

  comparisons   = st_params$comparisons   %||% list()
  correlations  = st_params$correlations  %||% list()
  linear_models = st_params$linear_models %||% list()

  # p-value adjustment is hard-coded to BH (= FDR). Bonferroni / others
  # were dropped from the YAML for simplicity.
  p_adjust      = "BH"
  sig_threshold = st_params$sig_threshold %||% 0.05

  list(
    stats         = .run_comparisons(de, comparisons,   p_adjust, sig_threshold),
    correlations  = .run_correlations(de, correlations, p_adjust, sig_threshold),
    linear_models = .run_linear_models(de, linear_models, p_adjust, sig_threshold)
  )
}


# ---- Comparisons -----------------------------------------------------------

#' @keywords internal
#' @noRd
.run_comparisons = function(de, comparisons, p_adjust, sig_threshold){

  comparisons = .filter_enabled(comparisons)
  if(length(comparisons) == 0)
    return(.empty_stats_df())

  rows = lapply(comparisons, function(comp){
    .run_one_comparison(de, comp, p_adjust, sig_threshold)
  })

  out = dplyr::bind_rows(rows)
  if(nrow(out) > 0){
    out$p_adj       = stats::p.adjust(out$p_value, method = p_adjust)
    out$significant = out$p_adj < sig_threshold
  }
  out
}

#' @keywords internal
#' @noRd
.run_one_comparison = function(de, comp, p_adjust, sig_threshold){

  comp_name = comp$name %||% "unnamed_comparison"
  factor_name = comp$compare$factor
  levels_keep = comp$compare$levels
  if(is.null(factor_name) || is.null(levels_keep))
    stop(sprintf("Comparison '%s': compare.factor and compare.levels required.", comp_name))

  # Build a single combined subset: user-supplied conditions + the level
  # restriction. de_subset() supports vector values via %in%, so we can pass
  # `levels_keep` directly. Doing it in one call avoids mutating the DE in
  # place (which trips SummarizedExperiment's assay-replacement check).
  cond = comp$subset
  cond = if(is.null(cond)) list() else as.list(cond)

  if(!factor_name %in% colnames(de$sample_meta))
    stop(sprintf("Comparison '%s': factor '%s' not in sample_meta.",
                 comp_name, factor_name))
  cond[[factor_name]] = levels_keep
  de_sub = de_subset(de, conditions = cond)

  group = factor(de_sub$sample_meta[[factor_name]], levels = levels_keep)
  if(nlevels(droplevels(group)) < 2){
    warning(sprintf("Comparison '%s': fewer than 2 non-empty levels — skipped.",
                    comp_name))
    return(.empty_stats_df())
  }

  data = de_sub$data
  feats = colnames(data)
  ngrp  = nlevels(droplevels(group))

  if(ngrp == 2){
    .compare_two(comp_name, data, group, feats, levels_keep)
  } else {
    .compare_many(comp_name, data, group, feats, levels_keep)
  }
}

#' @keywords internal
#' @noRd
.compare_two = function(comp_name, data, group, feats, levels_keep){
  g1 = levels_keep[1]; g2 = levels_keep[2]
  i1 = which(group == g1); i2 = which(group == g2)

  rows = lapply(feats, function(f){
    x = data[[f]]
    x1 = x[i1]; x2 = x[i2]
    m1 = mean(x1, na.rm = TRUE); m2 = mean(x2, na.rm = TRUE)
    fc     = suppressWarnings(m2 / m1)
    log2fc = suppressWarnings(log2(fc))

    tt = tryCatch(stats::t.test(x1, x2),     error = function(e) NULL)
    wt = tryCatch(stats::wilcox.test(x1, x2, exact = FALSE),
                  error = function(e) NULL, warning = function(w) NULL)

    list(
      .row(comp_name, "t_test", f, 2, g1, g2, m1, m2, fc, log2fc,
           statistic = if(is.null(tt)) NA_real_ else unname(tt$statistic),
           df        = if(is.null(tt)) NA_character_ else as.character(round(tt$parameter, 2)),
           p_value   = if(is.null(tt)) NA_real_ else tt$p.value),
      .row(comp_name, "wilcoxon", f, 2, g1, g2, m1, m2, fc, log2fc,
           statistic = if(is.null(wt)) NA_real_ else unname(wt$statistic),
           df        = NA_character_,
           p_value   = if(is.null(wt)) NA_real_ else wt$p.value)
    )
  })

  dplyr::bind_rows(unlist(rows, recursive = FALSE))
}

#' @keywords internal
#' @noRd
.compare_many = function(comp_name, data, group, feats, levels_keep){
  rows = lapply(feats, function(f){
    x = data[[f]]
    df_in = data.frame(value = x, grp = group)

    # Global tests
    av = tryCatch(stats::aov(value ~ grp, data = df_in),  error = function(e) NULL)
    av_p = if(is.null(av)) NA_real_ else summary(av)[[1]][["Pr(>F)"]][1]
    av_f = if(is.null(av)) NA_real_ else summary(av)[[1]][["F value"]][1]
    av_df = if(is.null(av)) NA_character_ else
      sprintf("%d,%d", summary(av)[[1]][["Df"]][1], summary(av)[[1]][["Df"]][2])

    kw = tryCatch(stats::kruskal.test(value ~ grp, data = df_in),
                  error = function(e) NULL)

    global = list(
      .row(comp_name, "anova",    f, length(levels_keep), NA, NA, NA, NA, NA, NA,
           statistic = av_f, df = av_df, p_value = av_p),
      .row(comp_name, "kruskal",  f, length(levels_keep), NA, NA, NA, NA, NA, NA,
           statistic = if(is.null(kw)) NA_real_ else unname(kw$statistic),
           df        = if(is.null(kw)) NA_character_ else as.character(kw$parameter),
           p_value   = if(is.null(kw)) NA_real_ else kw$p.value)
    )

    # Tukey HSD post-hoc on aov fit
    tk_rows = list()
    if(!is.null(av)){
      tk = tryCatch(stats::TukeyHSD(av), error = function(e) NULL)
      if(!is.null(tk)){
        tk_df = as.data.frame(tk$grp)
        tk_df$pair = rownames(tk_df)
        tk_rows = lapply(seq_len(nrow(tk_df)), function(k){
          pr = strsplit(tk_df$pair[k], "-", fixed = TRUE)[[1]]
          g_a = pr[2]; g_b = pr[1]                  # TukeyHSD label = "B-A"
          m_a = mean(x[group == g_a], na.rm = TRUE)
          m_b = mean(x[group == g_b], na.rm = TRUE)
          fc_pair = suppressWarnings(m_b / m_a)
          .row(comp_name, "tukey", f, 2, g_a, g_b, m_a, m_b,
               fc_pair, suppressWarnings(log2(fc_pair)),
               statistic = tk_df$diff[k],
               df = NA_character_,
               p_value = tk_df$`p adj`[k])
        })
      }
    }

    c(global, tk_rows)
  })

  dplyr::bind_rows(unlist(rows, recursive = FALSE))
}


# ---- Correlations (pairwise feature x feature) -----------------------------

#' @keywords internal
#' @noRd
.run_correlations = function(de, correlations, p_adjust, sig_threshold){
  correlations = .filter_enabled(correlations)
  if(length(correlations) == 0)
    return(.empty_corr_df())

  rows = lapply(correlations, function(corr){
    name    = corr$name    %||% "unnamed_correlation"
    methods = corr$methods %||% c("pearson", "spearman")

    # Accept either `subsets:` (plural list) or `subset:` (singular, legacy).
    subsets = corr$subsets
    if(is.null(subsets)){
      if(is.null(corr$subset)) subsets = list(list())   # one analysis on all samples
      else                     subsets = list(as.list(corr$subset))
    }

    include_self = isTRUE(corr$include_self %||% FALSE)

    per_subset = lapply(subsets, function(sub){
      sub_lab = if(length(sub) == 0) "all"
                else paste(sprintf("%s=%s", names(sub), unlist(sub)),
                           collapse = ", ")

      de_sub = de_subset(de, conditions = if(length(sub)) as.list(sub) else NULL)
      X = as.matrix(de_sub$data)

      if(ncol(X) < 2 || nrow(X) < 3){
        warning(sprintf("Correlation '%s' / %s: insufficient data (n=%d, features=%d).",
                        name, sub_lab, nrow(X), ncol(X)))
        return(.empty_corr_df())
      }

      per_method = lapply(methods, function(m){
        .pairwise_cor_long(X, method = m, name = name, sub_lab = sub_lab,
                           include_self = include_self)
      })
      dplyr::bind_rows(per_method)
    })

    dplyr::bind_rows(per_subset)
  })

  out = dplyr::bind_rows(rows)
  if(nrow(out) > 0){
    out$p_adj       = stats::p.adjust(out$p_value, method = p_adjust)
    out$significant = out$p_adj < sig_threshold
  }
  out
}

#' Compute pairwise correlations + vectorized p-values; return long format
#' (lower triangle; diagonal/self-pairs included only when `include_self = TRUE`).
#' @keywords internal
#' @noRd
.pairwise_cor_long = function(X, method, name, sub_lab,
                              include_self = FALSE){

  # Pairwise complete-obs correlation matrix
  r_mat = suppressWarnings(stats::cor(X, use = "pairwise.complete.obs",
                                      method = method))

  # Per-pair sample count (handles NA pattern)
  ok    = !is.na(X)
  n_mat = crossprod(ok)                               # features x features

  # Vectorized t-statistic + p-value (valid for Pearson; an asymptotic
  # approximation for Spearman that matches what cor.test() reports for
  # large n with `exact = FALSE`).
  denom = 1 - r_mat^2
  denom[denom <= 0] = NA_real_
  t_mat = r_mat * sqrt((n_mat - 2) / denom)
  p_mat = 2 * stats::pt(-abs(t_mat), df = n_mat - 2)
  diag(p_mat) = NA_real_                              # self-correlation

  # Lower triangle → long (omit diagonal / self-pairs unless requested)
  feats = colnames(X)
  idx = which(lower.tri(r_mat, diag = include_self), arr.ind = TRUE)
  data.frame(
    correlation = name,
    subset      = sub_lab,
    method      = method,
    feature_a   = feats[idx[, 1]],
    feature_b   = feats[idx[, 2]],
    n           = n_mat[idx],
    estimate    = r_mat[idx],
    statistic   = t_mat[idx],
    p_value     = p_mat[idx],
    stringsAsFactors = FALSE
  )
}


# ---- Linear models ---------------------------------------------------------

#' @keywords internal
#' @noRd
.run_linear_models = function(de, linear_models, p_adjust, sig_threshold){
  linear_models = .filter_enabled(linear_models)
  if(length(linear_models) == 0)
    return(.empty_lm_df())

  rows = lapply(linear_models, function(mod){
    name = mod$name %||% "unnamed_model"
    rhs  = mod$formula %||% stop(sprintf("Linear model '%s': `formula` required.", name))
    cond = mod$subset; if(!is.null(cond)) cond = as.list(cond)
    de_sub = de_subset(de, conditions = cond)

    feats = colnames(de_sub$data)
    smeta = de_sub$sample_meta

    inner = lapply(feats, function(f){
      df_in = cbind(.feature = de_sub$data[[f]], smeta)
      fm = stats::as.formula(paste(".feature", rhs))
      fit = tryCatch(stats::lm(fm, data = df_in), error = function(e) NULL)
      if(is.null(fit)) return(NULL)
      co = tryCatch(stats::coef(summary(fit)), error = function(e) NULL)
      if(is.null(co)) return(NULL)

      # ---- friendlier term labels (S9) ----
      # Default lm() names terms like "(Intercept)" and "TreatmentHDM_Curdlan".
      # Rename to make the reference level explicit:
      #   "(Intercept)"            -> "Intercept (reference: Treatment=Ctrl, Sex=Female)"
      #   "TreatmentHDM_Curdlan"   -> "Treatment: HDM_Curdlan (vs Ctrl)"
      raw_terms = rownames(co)
      term_labels = .pretty_lm_terms(raw_terms, fit, smeta)

      data.frame(
        model      = name,
        formula    = rhs,
        feature    = f,
        term       = term_labels,
        term_raw   = raw_terms,
        estimate   = co[, "Estimate"],
        std_error  = co[, "Std. Error"],
        statistic  = co[, "t value"],
        p_value    = co[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(inner)
  })

  out = dplyr::bind_rows(rows)
  if(nrow(out) > 0){
    out$p_adj       = stats::p.adjust(out$p_value, method = p_adjust)
    out$significant = out$p_adj < sig_threshold
  }
  out
}


# ---- Empty-frame templates (so xlsx writer always sees the schema) --------

#' @keywords internal
#' @noRd
.empty_stats_df = function() data.frame(
  comparison = character(), method = character(), feature = character(),
  n_groups = integer(), group1 = character(), group2 = character(),
  mean_g1 = numeric(), mean_g2 = numeric(),
  FC = numeric(), log2FC = numeric(),
  statistic = numeric(), df = character(),
  p_value = numeric(), p_adj = numeric(), significant = logical(),
  stringsAsFactors = FALSE
)

#' @keywords internal
#' @noRd
.empty_corr_df = function() data.frame(
  correlation = character(), subset = character(), method = character(),
  feature_a = character(), feature_b = character(), n = integer(),
  estimate = numeric(), statistic = numeric(),
  p_value = numeric(), p_adj = numeric(), significant = logical(),
  stringsAsFactors = FALSE
)

#' @keywords internal
#' @noRd
.empty_lm_df = function() data.frame(
  model = character(), formula = character(), feature = character(),
  term = character(), term_raw = character(),
  estimate = numeric(), std_error = numeric(), statistic = numeric(),
  p_value = numeric(), p_adj = numeric(), significant = logical(),
  stringsAsFactors = FALSE
)

#' Translate raw lm() coefficient names into something readable.
#'
#' Examples (with Treatment ∈ {Ctrl, HDM_Curdlan}, Sex ∈ {Female, Male}):
#'  "(Intercept)"           -> "Intercept (ref: Treatment=Ctrl; Sex=Female)"
#'  "TreatmentHDM_Curdlan"  -> "Treatment: HDM_Curdlan (vs Ctrl)"
#'  "SexMale:TreatmentHDM_Curdlan" -> kept as-is (interaction); doc'd elsewhere
#' @keywords internal
#' @noRd
.pretty_lm_terms = function(raw_terms, fit, smeta){
  # Pull factor reference levels from the model frame xlevels
  xlev = tryCatch(fit$xlevels, error = function(e) list())
  ref_levels = vapply(names(xlev), function(v) xlev[[v]][1], character(1))
  names(ref_levels) = names(xlev)

  vapply(raw_terms, function(t){
    if(t == "(Intercept)"){
      if(length(ref_levels) == 0) return("Intercept")
      ref_str = paste(sprintf("%s=%s", names(ref_levels), unname(ref_levels)),
                      collapse = "; ")
      return(sprintf("Intercept (ref: %s)", ref_str))
    }
    # Try to peel "<factor><level>" off the start of the term
    for(v in names(ref_levels)){
      if(startsWith(t, v)){
        lvl = sub(paste0("^", v), "", t)
        # Skip pure numeric coefficients (continuous predictors)
        if(nzchar(lvl) && lvl %in% xlev[[v]])
          return(sprintf("%s: %s (vs %s)", v, lvl, ref_levels[[v]]))
      }
    }
    t
  }, character(1), USE.NAMES = FALSE)
}

#' Drop entries with `enabled: False`. Default = TRUE if absent.
#' @keywords internal
#' @noRd
.filter_enabled = function(entries){
  if(length(entries) == 0) return(entries)
  keep = vapply(entries, function(e){
    en = e$enabled
    if(is.null(en)) TRUE else isTRUE(en)
  }, logical(1))
  entries[keep]
}

#' @keywords internal
#' @noRd
.row = function(comp_name, method, feature, n_groups,
                group1, group2, mean_g1, mean_g2, FC, log2FC,
                statistic, df, p_value){
  data.frame(
    comparison = comp_name, method = method, feature = feature,
    n_groups = as.integer(n_groups),
    group1 = as.character(group1), group2 = as.character(group2),
    mean_g1 = as.numeric(mean_g1), mean_g2 = as.numeric(mean_g2),
    FC = as.numeric(FC),
    log2FC = as.numeric(log2FC),
    statistic = as.numeric(statistic),
    df = as.character(df),
    p_value = as.numeric(p_value),
    p_adj = NA_real_, significant = NA,
    stringsAsFactors = FALSE
  )
}

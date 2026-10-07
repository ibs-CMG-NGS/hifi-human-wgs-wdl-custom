# aggregate_6mA.py 출력(샘플당 1개)을 모아 특징(feature) 단위로 그룹 효과 검정.
#
# 모델: logit(6mA 비율) ~ group + run
#   group: A=무처치, B=급성 H2O2, C=급성 H2O2+약물, D=만성 H2O2+약물 (2026-10-07 사용자 확인)
#   대비(contrast): H2O2 = B-A, drug = C-B (급성 H2O2에서 약물 효과), chronic = D-C (약물 하 만성 vs 급성)
#   F-test = 세 대비 동시 검정 (4군 간 차이 여부)
#   - 샘플마다 SMRT cell이 달라 cell 단위 기술적 변동이 그룹과 겹치므로 두 척도로 본다:
#     abs : 특징의 6mA 비율 자체
#     rel : 특징 비율 / 그 샘플의 상염색체 전체 비율 (농축도; cell 전체 수준 차이 제거)
#   - limma eBayes(trend=TRUE)로 특징 간 분산 공유
#   - run 안 그룹 라벨 permutation으로 "p<α 특징 수"의 null 분포를 만들어 실제와 비교
# 사용: Rscript stats_6mA.R samples.tsv AGG_DIR OUT_DIR [N_PERM]
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
samples <- read.delim(args[1], stringsAsFactors = FALSE)
agg_dir <- args[2]
out_dir <- args[3]
n_perm <- if (length(args) >= 4) as.integer(args[4]) else 200
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# 층별 최소 N_valid (모든 샘플에서). 전역 비율 ~0.5%를 기준으로 기대 호출 수 ≥ ~50 정도
MIN_VALID <- c(element = 1e5, cgi = 1e5, rep_class = 1e5, rep_family = 1e5, rep_name = 2e4,
               gene_body = 2e4, gene_promoter = 1e4, bin100k = 2e4, ctx5 = 1e5, chrom = 1e5)
ALPHAS <- c(1e-3, 1e-4)
AUTOSOMES <- paste0("chr", 1:19)

agg <- lapply(seq_len(nrow(samples)), function(i) {
  f <- file.path(agg_dir, paste0(samples$sample_id[i], ".6mA.agg.tsv.gz"))
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$sample_id <- samples$sample_id[i]
  d
})
agg <- do.call(rbind, agg)

grp <- samples$group
run <- factor(samples$run)
COEFS <- c("H2O2", "drug", "chronic")
# 계수 순서: (Intercept), gB, gC, gD, run2, run3
CONTRASTS <- cbind(H2O2 = c(0, 1, 0, 0, 0, 0), drug = c(0, -1, 1, 0, 0, 0), chronic = c(0, 0, -1, 1, 0, 0))

# ---- 1. 샘플별 전역 비율 ----
glob <- subset(agg, layer == "chrom" & feature %in% AUTOSOMES)
glob <- aggregate(cbind(N_valid, N_mod, n_sites, n_hit_sites) ~ sample_id, glob, sum)
glob <- glob[match(samples$sample_id, glob$sample_id), ]
glob$rate <- glob$N_mod / glob$N_valid
chrm <- subset(agg, layer == "chrom" & feature == "chrM")
glob$chrM_rate <- (chrm$N_mod / chrm$N_valid)[match(glob$sample_id, chrm$sample_id)]
glob$group <- grp
glob$run <- samples$run
write.table(glob, file.path(out_dir, "global_rate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
gfit <- lm(qlogis(rate) ~ g + run, data = data.frame(rate = glob$rate, g = factor(grp, levels = c("A", "B", "C", "D")), run))
sink(file.path(out_dir, "global_rate_model.txt"))
print(glob[, c("sample_id", "group", "run", "rate", "chrM_rate")], row.names = FALSE)
cat("\n## logit(global rate) ~ group + run (reference = A)\n")
print(summary(gfit))
sink()

# ---- 2. 층별 특징 검정 ----
perm_labels <- function() {
  g <- grp
  for (r in unique(samples$run)) {
    idx <- which(samples$run == r)
    g[idx] <- sample(grp[idx])
  }
  g
}
fit_p <- function(Y, g) {
  X <- model.matrix(~ factor(g, levels = c("A", "B", "C", "D")) + run)
  fit <- eBayes(contrasts.fit(lmFit(Y, X), CONTRASTS), trend = TRUE)
  list(fit = fit, F_p = topTable(fit, coef = COEFS, number = Inf, sort.by = "none")$P.Value)
}

perm_rows <- list()
for (L in names(MIN_VALID)) {
  d <- subset(agg, layer == L)
  if (L == "chrom") d <- subset(d, feature %in% c(AUTOSOMES, "chrX", "chrY", "chrM"))
  nv <- xtabs(N_valid ~ feature + sample_id, d)[, samples$sample_id, drop = FALSE]
  nm <- xtabs(N_mod ~ feature + sample_id, d)[, samples$sample_id, drop = FALSE]
  keep <- apply(nv, 1, min) >= MIN_VALID[[L]]
  nv <- unclass(nv[keep, , drop = FALSE])
  nm <- unclass(nm[keep, , drop = FALSE])
  if (nrow(nv) < 3) next
  abs_y <- qlogis((nm + 0.5) / (nv + 1))
  rel_y <- abs_y - matrix(qlogis(glob$rate), nrow(abs_y), ncol(abs_y), byrow = TRUE)
  for (scale in c("abs", "rel")) {
    Y <- if (scale == "abs") abs_y else rel_y
    r <- fit_p(Y, grp)
    res <- data.frame(feature = rownames(Y), mean_rate = rowSums(nm) / rowSums(nv),
                      min_N_valid = apply(nv, 1, min),
                      F_p = r$F_p, F_fdr = p.adjust(r$F_p, "BH"))
    for (cf in COEFS) {
      tt <- topTable(r$fit, coef = cf, number = Inf, sort.by = "none")
      res[[paste0(cf, "_logFC")]] <- tt$logFC
      res[[paste0(cf, "_p")]] <- tt$P.Value
      res[[paste0(cf, "_fdr")]] <- tt$adj.P.Val
    }
    # 그룹별 평균 비율 (해석용)
    for (g in c("A", "B", "C", "D")) res[[paste0("rate_", g)]] <- rowSums(nm[, grp == g]) / rowSums(nv[, grp == g])
    res <- res[order(res$F_p), ]
    write.table(res, gzfile(file.path(out_dir, sprintf("%s.%s.tsv.gz", L, scale))),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # permutation (특징 수가 충분할 때만 의미 있음)
    if (nrow(Y) >= 200) {
      obs <- sapply(ALPHAS, function(a) sum(r$F_p < a))
      perm <- t(replicate(n_perm, {
        p <- fit_p(Y, perm_labels())$F_p
        sapply(ALPHAS, function(a) sum(p < a))
      }))
      for (j in seq_along(ALPHAS)) {
        perm_rows[[length(perm_rows) + 1]] <- data.frame(
          layer = L, scale = scale, n_features = nrow(Y), alpha = ALPHAS[j], observed = obs[j],
          perm_median = median(perm[, j]), perm_95 = quantile(perm[, j], 0.95), perm_max = max(perm[, j]),
          empirical_p = (1 + sum(perm[, j] >= obs[j])) / (1 + n_perm))
      }
    }
    cat(sprintf("%-14s %-3s features=%6d  F_fdr<0.05: %d\n", L, scale, nrow(Y), sum(res$F_fdr < 0.05)))
  }
}
perm_tab <- do.call(rbind, perm_rows)
write.table(perm_tab, file.path(out_dir, "permutation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(perm_tab, row.names = FALSE)

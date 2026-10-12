# Fig 4e: 영역 단위 5mC 차등 검정 (condition + run), run 안 permutation으로 보정, 기술적 교란 점검.
#
# 입력: dmr_region_counts.py 출력 (영역 × 샘플 메틸화 read 수 / coverage)
# 반응: M = log2((meth + 1) / (cov - meth + 1)); 영역 coverage ≥ 50이 보장되므로 가중치 없이 limma-trend로 평균-분산 관계를 처리
#   (영역별 가중치를 주면 lmFit이 영역마다 따로 적합해 permutation이 수백 배 느려짐)
# 모형 1 (주): M ~ condition + run, limma eBayes(trend=TRUE)
#   대비: H2O2 = B-A, drug = C-B, chronic = D-C, DvsA = D-A, F = 4군 차이
# 모형 2 (교란 점검): M ~ condition + run + tech,  tech = logit(STR 수축 지수, 2 bp) — 라이브러리 artifact 지표
# 보정: run 안 조건 라벨 permutation(13,824 배치 중 N_PERM개 무작위)으로 "p < α 영역 수"와 최소 p의 null 분포
# 패턴: F p < 1e-3 영역마다 run 효과를 뺀 샘플값과 tech의 상관(|r|) — 교란을 따라가는 영역 비율을 배경과 비교
# 사용: Rscript dmr_limma.R samples.tsv tech.tsv COUNT_DIR OUT_DIR SET [N_PERM]
suppressPackageStartupMessages({library(limma); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
samples <- fread(args[1]); tech <- fread(args[2])
cdir <- args[3]; odir <- args[4]; set <- args[5]
n_perm <- if (length(args) >= 6) as.integer(args[6]) else 200
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

meth <- as.matrix(fread(file.path(cdir, paste0(set, ".meth.tsv.gz"))), rownames = 1)
cov <- as.matrix(fread(file.path(cdir, paste0(set, ".cov.tsv.gz"))), rownames = 1)
reg <- fread(file.path(cdir, paste0(set, ".regions.tsv.gz")))
sids <- samples$sample_id
meth <- meth[, sids]; cov <- cov[, sids]
M <- log2((meth + 1) / (cov - meth + 1))
grp <- factor(samples$group, levels = c("A", "B", "C", "D"))
run <- factor(samples$run)
tc <- tech[match(sids, tech$sample), contraction_idx_2bp]
tc <- log(tc / (1 - tc))
CON <- cbind(H2O2 = c(0, 1, 0, 0, 0, 0), drug = c(0, -1, 1, 0, 0, 0), chronic = c(0, 0, -1, 1, 0, 0),
             DvsA = c(0, 0, 0, 1, 0, 0))
ALPHAS <- c(1e-3, 1e-4)

fit_all <- function(g, extra = NULL) {
  X <- model.matrix(~ g + run)
  if (!is.null(extra)) X <- cbind(X, tech = extra)
  con <- if (is.null(extra)) CON else rbind(CON, 0)
  f <- lmFit(M, X)
  f2 <- eBayes(contrasts.fit(f, con), trend = TRUE)
  list(fit = f2, p = cbind(f2$p.value, F = f2$F.p.value))
}

# 주 모형과 교란 점검 모형
main <- fit_all(grp)
adj <- fit_all(grp, tc)
cnt <- function(p) sapply(ALPHAS, function(a) colSums(p < a))

# permutation: run 안에서 조건 라벨을 섞음
runs <- split(seq_along(sids), samples$run)
perm_lab <- function() {
  g <- as.character(grp)
  for (r in runs) g[r] <- g[r][sample(length(r))]
  factor(g, levels = levels(grp))
}
null_cnt <- array(NA, c(n_perm, ncol(main$p), length(ALPHAS)), dimnames = list(NULL, colnames(main$p), ALPHAS))
null_minp <- matrix(NA, n_perm, ncol(main$p), dimnames = list(NULL, colnames(main$p)))
seen <- character(0)
i <- 0
while (i < n_perm) {
  g <- perm_lab()
  key <- paste(g, collapse = "")
  if (key %in% seen || key == paste(grp, collapse = "")) next
  seen <- c(seen, key); i <- i + 1
  p <- fit_all(g)$p
  null_cnt[i, , ] <- cnt(p)
  null_minp[i, ] <- apply(p, 2, min)
}

obs <- cnt(main$p); obs_adj <- cnt(adj$p)
summ <- rbindlist(lapply(colnames(main$p), function(k) rbindlist(lapply(seq_along(ALPHAS), function(j) data.table(
  set = set, n_regions = nrow(M), contrast = k, alpha = ALPHAS[j], n_hits = obs[k, j],
  null_median = median(null_cnt[, k, j]), null_q95 = quantile(null_cnt[, k, j], 0.95),
  perm_p_count = (1 + sum(null_cnt[, k, j] >= obs[k, j])) / (n_perm + 1),
  n_hits_tech_adjusted = obs_adj[k, j],
  min_p = min(main$p[, k]), perm_p_minp = (1 + sum(null_minp[, k] <= min(main$p[, k]))) / (n_perm + 1))))))
fwrite(summ, file.path(odir, paste0(set, ".summary.tsv")), sep = "\t")

# 영역별 결과 + 교란 상관
Xr <- model.matrix(~ run)
resid <- M - t(apply(M, 1, function(y) Xr %*% qr.solve(Xr, y)))
tc_r <- tc - Xr %*% qr.solve(Xr, tc)
r_tech <- as.vector(cor(t(resid), tc_r))
gm <- sapply(levels(grp), function(l) rowMeans(M[, grp == l, drop = FALSE]))
res <- data.table(region_id = rownames(M), reg[match(rownames(M), reg$region_id), .(chrom, start, end, name, ncpg)],
                  meth_frac_mean = rowSums(meth) / rowSums(cov),
                  mA = gm[, "A"], mB = gm[, "B"], mC = gm[, "C"], mD = gm[, "D"],
                  logFC_H2O2 = main$fit$coefficients[, "H2O2"], logFC_drug = main$fit$coefficients[, "drug"],
                  logFC_chronic = main$fit$coefficients[, "chronic"], logFC_DvsA = main$fit$coefficients[, "DvsA"],
                  p_H2O2 = main$p[, "H2O2"], p_drug = main$p[, "drug"], p_chronic = main$p[, "chronic"],
                  p_DvsA = main$p[, "DvsA"], p_F = main$p[, "F"], p_F_tech_adj = adj$p[, "F"], r_tech = r_tech)
res[, fdr_F := p.adjust(p_F, "BH")]
setorder(res, p_F)
fwrite(res, file.path(odir, paste0(set, ".regions.tsv.gz")), sep = "\t")

hit <- res$p_F < 1e-3
cat(sprintf("## %s: %d regions\n", set, nrow(M)))
print(summ[, .(contrast, alpha, n_hits, null_median, null_q95, perm_p_count, n_hits_tech_adjusted, min_p, perm_p_minp)])
cat(sprintf("F p<1e-3 hits: %d; |r_tech|>0.7 in hits %.3f vs background %.3f; median |r_tech| hits %.3f vs background %.3f\n",
            sum(hit), mean(abs(res$r_tech[hit]) > 0.7), mean(abs(res$r_tech) > 0.7),
            median(abs(res$r_tech[hit])), median(abs(res$r_tech))))
cat("top 15 by F:\n")
print(head(res[, .(chrom, start, end, name, meth_frac_mean, mA, mB, mC, mD, p_F, fdr_F, p_F_tech_adj, r_tech)], 15))

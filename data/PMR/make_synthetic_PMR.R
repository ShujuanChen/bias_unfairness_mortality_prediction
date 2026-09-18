#!/usr/bin/env Rscript
#
# Synthetic Census linked mortality records.

set.seed(6089)

N_IN <- 3750
N_OUT <- 8250
n <- N_IN + N_OUT

out_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
})

N_LSOA_E <- 1800
N_LSOA_W <- 200
lsoa11_e <- sprintf("E01%06d", seq_len(N_LSOA_E))
lsoa11_w <- sprintf("W01%06d", seq_len(N_LSOA_W))
lsoa11 <- c(lsoa11_e, lsoa11_w)
lsoa01_e <- sprintf("E01%06d", seq_len(N_LSOA_E))
lsoa01_w <- sprintf("W01%06d", seq_len(N_LSOA_W))
lsoa01 <- c(lsoa01_e, lsoa01_w)

draw_cat <- function(n, p, missing = 0, numeric = FALSE, integer = FALSE) {
  x <- sample(names(p), n, replace = TRUE, prob = as.numeric(p))
  if (missing > 0) x[runif(n) < missing] <- NA_character_
  if (integer) as.integer(x) else if (numeric) as.numeric(x) else x
}
draw_num <- function(n, mean, sd, lo, hi, missing = 0, digits = NULL) {
  x <- pmin(pmax(rnorm(n, mean, sd), lo), hi)
  if (!is.null(digits)) x <- round(x, digits)
  if (missing > 0) x[runif(n) < missing] <- NA_real_
  x
}
draw_date <- function(n, from, to, missing = 0) {
  a <- as.Date(from); b <- as.Date(to)
  x <- format(a + floor(runif(n) * as.numeric(b - a)), "%Y-%m-%d")
  if (missing > 0) x[runif(n) < missing] <- NA_character_
  x
}
write_csv_file <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "NA")
  message("wrote ", path, "  (", nrow(df), " rows, ", ncol(df), " columns)")
}

P_ICD <- c("C349" = 0.12799, "I219" = 0.062529, "I259" = 0.049509,
           "I251" = 0.049321, "C509" = 0.042132, "U071" = 0.041821,
           "C259" = 0.035507, "C159" = 0.030578, "J449" = 0.030546,
           "J440" = 0.030486, "C61" = 0.025204, "I64" = 0.018988,
           "C719" = 0.018268, "F03" = 0.017696, "C189" = 0.01761,
           "C56" = 0.01668, "J189" = 0.016616, "C800" = 0.014708,
           "C64" = 0.014704, "G309" = 0.013763, "J441" = 0.013521,
           "C20" = 0.01325, "C679" = 0.012805, "J841" = 0.011661,
           "I619" = 0.011448, "J180" = 0.010572, "C221" = 0.010309,
           "C169" = 0.010299, "G20" = 0.010287, "C19" = 0.010158,
           "C809" = 0.00986, "G122" = 0.009594, "K709" = 0.009302,
           "C900" = 0.009273, "K703" = 0.009002, "C920" = 0.008876,
           "C439" = 0.008835, "I802" = 0.00881, "C220" = 0.008696,
           "K704" = 0.008462, "K746" = 0.008337, "G35" = 0.007906,
           "C260" = 0.007809, "I269" = 0.007508, "F019" = 0.007085,
           "C459" = 0.006694, "C541" = 0.006679, "I609" = 0.006591,
           "C859" = 0.006175, "I639" = 0.005837, "N390" = 0.005669,
           "I119" = 0.005491, "I713" = 0.005141, "I710" = 0.005099,
           "C80" = 0.005052, "I489" = 0.004831, "I110" = 0.004665,
           "R99" = 0.004665, "X700" = 0.004576, "K559" = 0.004512)

pmr <- data.frame(.row = seq_len(n), stringsAsFactors = FALSE)
p_sampling_weight <- c("1" = 0.655458, "20" = 0.344542)
pmr$sampling_weight <- draw_cat(n, p_sampling_weight, missing = 0, numeric = TRUE, integer = TRUE)
p_age_census <- c("64" = 0.054512, "68" = 0.054024, "67" = 0.052573,
                  "66" = 0.05226, "69" = 0.051137, "63" = 0.050503,
                  "65" = 0.045949, "62" = 0.043679, "61" = 0.040235,
                  "60" = 0.036545, "59" = 0.034208, "58" = 0.032723,
                  "57" = 0.031626, "56" = 0.029573, "55" = 0.028604,
                  "54" = 0.028112, "53" = 0.027773, "52" = 0.026856,
                  "51" = 0.026003, "50" = 0.025829, "49" = 0.025484,
                  "48" = 0.025033, "47" = 0.024558, "46" = 0.024278,
                  "45" = 0.023255, "44" = 0.022289, "43" = 0.021492,
                  "42" = 0.021003, "41" = 0.019968, "40" = 0.019915)
pmr$age_census <- draw_cat(n, p_age_census, missing = 0, numeric = TRUE, integer = TRUE)
p_sex_census <- c("1" = 0.541338, "2" = 0.458662)
pmr$sex_census <- draw_cat(n, p_sex_census, missing = 0, numeric = TRUE, integer = TRUE)
p_ethpuk11_census <- c("01" = 0.890996, "04" = 0.021705, "09" = 0.017774,
                       "02" = 0.014922, "15" = 0.009947, "10" = 0.009116,
                       "14" = 0.007725, "13" = 0.007371, "12" = 0.003101,
                       "11" = 0.003072, "18" = 0.002872, "05" = 0.002546,
                       "16" = 0.00204, "08" = 0.001822, "07" = 0.00175,
                       "17" = 0.001387, "03" = 0.000888, "06" = 0.000884,
                       "XX" = 0.000081)
pmr$ethpuk11_census <- draw_cat(n, p_ethpuk11_census, missing = 0)
p_health_census <- c("2" = 0.361805, "3" = 0.248809, "1" = 0.207728,
                     "4" = 0.130881, "5" = 0.050696, "X" = 0.000081)
pmr$health_census <- draw_cat(n, p_health_census, missing = 0)
p_hlqpuk11_census <- c("10" = 0.341589, "15" = 0.217679, "11" = 0.128242,
                       "12" = 0.11834, "14" = 0.080997, "16" = 0.058703,
                       "13" = 0.054369, "XX" = 0.000081)
pmr$hlqpuk11_census <- draw_cat(n, p_hlqpuk11_census, missing = 0)
p_ecocatpuk11_census <- c("1" = 0.337503, "7" = 0.315699, "8" = 0.165418,
                          "3" = 0.147542, "5" = 0.031502, "2" = 0.001172,
                          "4" = 0.000835, "6" = 0.000248, "X" = 0.000081)
pmr$ecocatpuk11_census <- draw_cat(n, p_ecocatpuk11_census, missing = 0)
p_tenhuk11_census <- c("0" = 0.38025, "1" = 0.304654, "3" = 0.110912,
                       "4" = 0.088895, "5" = 0.077511, "NA" = 0.015075,
                       "9" = 0.010329, "2" = 0.005002, "7" = 0.004185,
                       "8" = 0.001756, "6" = 0.001305, "X" = 0.000126)
pmr$tenhuk11_census <- draw_cat(n, p_tenhuk11_census, missing = 0)
p_hhchuk11_census <- c("04" = 0.236133, "02" = 0.145503, "07" = 0.112575,
                       "03" = 0.093809, "01" = 0.062395, "06" = 0.061513,
                       "05" = 0.05174, "26" = 0.049782, "12" = 0.040587,
                       "21" = 0.033426, "22" = 0.01644, "NA" = 0.015075,
                       "23" = 0.014134, "19" = 0.012636, "18" = 0.010178,
                       "13" = 0.009749, "14" = 0.008823, "15" = 0.008635,
                       "20" = 0.007057, "16" = 0.003779, "25" = 0.003323,
                       "08" = 0.001364, "17" = 0.001052, "XX" = 0.000126,
                       "09" = 0.00005, "24" = 0.000045, "10" = 0.000037,
                       "11" = 0.000034)
pmr$hhchuk11_census <- draw_cat(n, p_hhchuk11_census, missing = 0)
p_residence_type_census <- c("H" = 0.984925, "C" = 0.015075)
pmr$residence_type_census <- draw_cat(n, p_residence_type_census, missing = 0)
p_ruralurban_code_census <- c("C1" = 0.440263, "A1" = 0.321328,
                              "D1" = 0.094023, "E1" = 0.059325,
                              "A2" = 0.037417, "F1" = 0.034382,
                              "D2" = 0.004148, "E2" = 0.003747,
                              "F2" = 0.003169, "C2" = 0.002199)
pmr$ruralurban_code_census <- draw_cat(n, p_ruralurban_code_census, missing = 0)
p_nssec_census <- c("4.1" = 0.086361, "9.1" = 0.08574, "7.1" = 0.082197,
                    "13.4" = 0.071473, "12.2" = 0.064804, "13.3" = 0.056162,
                    "5" = 0.049231, "10" = 0.0456, "14.1" = 0.043283,
                    "12.1" = 0.042703, "3.1" = 0.032697, "11.1" = 0.03134,
                    "8.1" = 0.029571, "6" = 0.026656, "7.2" = 0.025171,
                    "12.4" = 0.024575, "2" = 0.023947, "3.2" = 0.01925,
                    "13.2" = 0.018517, "14.2" = 0.017649, "13.1" = 0.015548,
                    "12.3" = 0.014752, "4.3" = 0.014416, "12.6" = 0.014019,
                    "3.3" = 0.011402, "4.2" = 0.010731, "7.3" = 0.009962,
                    "11.2" = 0.005849, "7.4" = 0.004851, "15" = 0.00412,
                    "9.2" = 0.004034, "12.7" = 0.00334, "4.4" = 0.002902,
                    "3.4" = 0.002085, "12.5" = 0.001955, "8.2" = 0.001554,
                    "13.5" = 0.001191, "1" = 0.00028, "XXXX" = 0.000081)
pmr$nssec_census <- draw_cat(n, p_nssec_census, missing = 0)

pmr$LSOA11CD <- sample(lsoa11_e, n, replace = TRUE)

p_excluded_age <- c("80" = 0.03976, "79" = 0.03862, "81" = 0.03839,
                    "78" = 0.03774, "82" = 0.03695, "77" = 0.0361,
                    "76" = 0.03609, "75" = 0.03505, "83" = 0.0348,
                    "74" = 0.03408, "84" = 0.03363, "73" = 0.03305,
                    "72" = 0.03163, "85" = 0.03138, "71" = 0.03002,
                    "86" = 0.02855, "70" = 0.02685, "87" = 0.02609,
                    "88" = 0.02338, "89" = 0.02155, "90" = 0.0192,
                    "91" = 0.01448, "39" = 0.0094, "38" = 0.00879,
                    "92" = 0.00877, "20" = 0.00829, "37" = 0.00819,
                    "19" = 0.00803, "21" = 0.00787, "36" = 0.00778,
                    "30" = 0.00763, "31" = 0.00761, "22" = 0.00751,
                    "35" = 0.00742, "32" = 0.0074, "29" = 0.00736,
                    "23" = 0.00733, "28" = 0.00724, "26" = 0.00719,
                    "18" = 0.00716, "34" = 0.00715, "25" = 0.0071,
                    "27" = 0.00709, "33" = 0.00704, "24" = 0.00702,
                    "93" = 0.00665, "17" = 0.00665, "15" = 0.00658,
                    "16" = 0.00643, "14" = 0.00637, "13" = 0.00621,
                    "0" = 0.00612, "3" = 0.00607, "1" = 0.00606,
                    "12" = 0.00604, "94" = 0.00602, "2" = 0.006,
                    "4" = 0.00599, "11" = 0.00591, "5" = 0.00581,
                    "10" = 0.00562, "6" = 0.0056, "7" = 0.00557,
                    "9" = 0.0054, "8" = 0.00539, "95" = 0.00466,
                    "96" = 0.0037, "97" = 0.00264, "98" = 0.0018,
                    "99" = 0.00115, "100" = 0.00075, "101" = 0.00044,
                    "102" = 0.00026, "103" = 0.00016, "104" = 0.00009,
                    "105" = 0.00004, "106" = 0.00003, "107" = 0.00002,
                    "108" = 0.00001, "109" = 0.00001, "110" = 0, "111" = 0,
                    "113" = 0, "114" = 0)
excluded <- sample(n, N_OUT)
by_age <- excluded[runif(N_OUT) < 0.9701]
pmr$age_census[by_age] <- draw_cat(length(by_age), p_excluded_age,
                                   numeric = TRUE, integer = TRUE)
pmr$LSOA11CD[setdiff(excluded, by_age)] <-
  sample(lsoa11_w, length(setdiff(excluded, by_age)), replace = TRUE)

pmr$dod_deaths <- NA_character_
is_death <- pmr$sampling_weight == 1
pmr$dod_deaths[is_death] <- draw_date(sum(is_death), "2011-03-27", "2023-02-15")
pmr$fic10und_deaths <- NA_character_
pmr$fic10und_deaths[is_death] <- draw_cat(sum(is_death), P_ICD)

pmr$.row <- NULL
pmr <- pmr[sample(n), , drop = FALSE]
write_csv_file(pmr, file.path(out_dir, "PMR.csv"))

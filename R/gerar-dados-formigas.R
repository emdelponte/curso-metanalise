# Gera dados SIMULADOS no desenho de Farji-Brener & Werenkraut (2015),
# Ecological Entomology 40:150-158 (formigas cortadeiras, fertilidade do solo
# e desempenho de plantas). Os valores NÃO são os dados do artigo: o número de
# estudos e de efeitos por variável-foco segue o artigo, e os efeitos médios
# verdadeiros foram ajustados para ficar próximos da Figura 1.
#
# Saída: dados/formigas_cortadeiras.csv (uma linha = um efeito, isto é, uma
# medida de um estudo, com média, DP e n de ninho e controle).

library(tidyverse)
set.seed(2015)

# Estudos: S01-S11 só solo; S12-S21 solo e planta; S22-S34 só planta
solo   <- sprintf("S%02d", 1:21)
planta <- sprintf("S%02d", 12:34)

# Desenho por variável-foco: estudos, medidas, n de efeitos, efeito médio e
# desvio entre estudos (tau), valores aproximados da Figura 1 do artigo.
desenho <- tribble(
  ~foco,         ~medidas,                  ~n_estudos, ~n_efeitos, ~mu,   ~tau, ~fonte, ~cv,
  "nutrientes",  c("C", "N", "P", "K"),     19,         66,         1.5,   1.8,  "solo",   0.35,
  "cations",     c("Al", "Ca", "Mg", "Na"), 13,         29,         1.6,   1.4,  "solo",   0.35,
  "pH",          c("pH_agua", "pH_CaCl2"),                    6,          7,         -1.1,   2.0,  "solo",   0.06,
  "crescimento", c("biomassa", "altura", "diametro"), 13, 26,      1.6,   2.0,  "planta", 0.20,
  "densidade",   c("cobertura", "individuos"), 11,      12,        -0.15,  2.0,  "planta", 0.20,
  "riqueza",     c("especies_total", "especies_nativas"),        8,          10,        -0.5,   1.1,  "planta", 0.20
)

# Distribui n_efeitos entre os estudos de cada foco (cada estudo tem >= 1)
distribui <- function(n_estudos, n_efeitos, n_medidas) {
  extra <- n_efeitos - n_estudos
  m <- rep(1, n_estudos)
  while (extra > 0) {
    i <- sample(which(m < n_medidas), 1)
    m[i] <- m[i] + 1
    extra <- extra - 1
  }
  m
}

simula_foco <- function(foco, medidas, n_estudos, n_efeitos, mu, tau, fonte, cv) {
  ids <- sample(if (fonte == "solo") solo else planta, n_estudos) |> sort()
  m <- distribui(n_estudos, n_efeitos, length(medidas))
  map2_dfr(ids, m, function(id, k) {
    n_ninho <- sample(4:30, 1)         # mesmo n para todas as medidas do estudo
    n_ctrl <- n_ninho
    d_estudo <- rnorm(1, mu, tau)      # efeito verdadeiro do estudo
    med <- sample(medidas, k)
    map_dfr(med, function(medida) {
      d_ij <- d_estudo + rnorm(1, 0, 0.3)       # variação entre medidas
      escala <- if (foco == "pH") 5.5 else runif(1, 5, 80)
      mc <- escala
      dp <- cv * mc
      mn <- max(mc + d_ij * dp, 0.05 * mc)
      # estatísticas amostrais observadas
      obs_media <- function(m0, n) rnorm(1, m0, dp / sqrt(n))
      obs_dp <- function(n) dp * sqrt(rchisq(1, n - 1) / (n - 1))
      tibble(
        estudo = id, foco = foco, medida = medida,
        media_ninho = obs_media(mn, n_ninho), dp_ninho = obs_dp(n_ninho), n_ninho = n_ninho,
        media_controle = obs_media(mc, n_ctrl), dp_controle = obs_dp(n_ctrl), n_controle = n_ctrl
      )
    })
  })
}

formigas <- pmap_dfr(desenho, simula_foco) |>
  mutate(across(starts_with(c("media", "dp")), \(x) round(x, 3)))

dir.create("dados", showWarnings = FALSE)
write_csv(formigas, "dados/formigas_cortadeiras.csv")

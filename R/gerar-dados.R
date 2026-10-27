# Gera os conjuntos de dados SIMULADOS usados no site (pasta dados/).
# Substitua cada CSV pelos dados reais quando estiverem disponíveis,
# mantendo os mesmos nomes de colunas.

library(tidyverse)
set.seed(2026)

# 1. Fungicida × produtividade da soja (mofo branco) -------------------------
# Um ensaio por estudo; usado nos capítulos 3, 4, 7, 8, 9 e 10.
k <- 30
fungicida_soja <- tibble(
  estudo = sprintf("E%02d", 1:k),
  ano = sample(2010:2023, k, replace = TRUE),
  regiao = sample(c("Norte", "Sul"), k, replace = TRUE, prob = c(0.4, 0.6)),
  latitude = if_else(regiao == "Norte", runif(k, 12, 20), runif(k, 20, 30)),
  pressao = sample(c("Baixa", "Alta"), k, replace = TRUE),
  n = sample(3:5, k, replace = TRUE),
  cv = round(runif(k, 4, 15), 1)
) |>
  mutate(
    prod_ctrl = round(rnorm(k, 3000, 400) - if_else(pressao == "Alta", 400, 0)),
    efeito = 150 + 15 * (latitude - 20) + if_else(pressao == "Alta", 250, 0) +
      rnorm(k, 0, 120),
    prod_trat = round(prod_ctrl + efeito + rnorm(k, 0, cv / 100 * prod_ctrl / sqrt(n))),
    media_geral = round((prod_trat + prod_ctrl) / 2, 1),
    latitude = round(latitude, 1)
  ) |>
  select(-efeito)

# 2. Tebuconazol × severidade da giberela (eficácia de controle) -------------
k <- 20
giberela <- tibble(
  estudo = sprintf("G%02d", 1:k),
  n = sample(4:6, k, replace = TRUE),
  sev_ctrl = round(runif(k, 10, 60), 1)
) |>
  mutate(
    sev_trat = round(sev_ctrl * exp(rnorm(k, log(0.45), 0.25)), 1),
    sd_ctrl = round(sev_ctrl * runif(k, 0.15, 0.35), 1),
    sd_trat = round(sev_trat * runif(k, 0.20, 0.45), 1)
  ) |>
  relocate(sd_ctrl, .after = sev_ctrl)

# 3. Tratamento de sementes × incidência de tombamento (binário) --------------
k <- 15
tombamento <- tibble(
  estudo = sprintf("T%02d", 1:k),
  n_trat = sample(c(100, 200, 400), k, replace = TRUE),
  n_ctrl = n_trat,
  p_ctrl = runif(k, 0.08, 0.40)
) |>
  mutate(
    p_trat = plogis(qlogis(p_ctrl) + rnorm(k, log(0.45), 0.3)),
    doentes_trat = rbinom(k, n_trat, p_trat),
    doentes_ctrl = rbinom(k, n_ctrl, p_ctrl),
    sadias_trat = n_trat - doentes_trat,
    sadias_ctrl = n_ctrl - doentes_ctrl
  ) |>
  select(estudo, doentes_trat, sadias_trat, doentes_ctrl, sadias_ctrl)

# 4. Ferrugem asiática × produtividade (parcelas por estudo) -----------------
ferrugem_parcelas <- map_dfr(1:15, function(i) {
  n <- sample(8:20, 1)
  altitude <- sample(c("Baixa", "Alta"), 1)
  inclinacao <- rnorm(1, if_else(altitude == "Baixa", 28, 18), 5)
  sev <- runif(n, 0, 80)
  tibble(
    estudo = sprintf("F%02d", i),
    altitude = altitude,
    parcela = 1:n,
    severidade = round(sev, 1),
    produtividade = round(3600 - inclinacao * sev + rnorm(n, 0, 300))
  )
})

# 5. Vários ensaios por estudo (modelo multinível) ---------------------------
n_estudos <- 12
efeito_estudo <- rnorm(n_estudos, 350, 150)
ensaios_multinivel <- map_dfr(1:n_estudos, function(i) {
  m <- sample(2:6, 1)
  n <- 4
  prod_ctrl <- round(rnorm(m, 3000, 350))
  cv <- round(runif(m, 5, 14), 1)
  tibble(
    estudo = sprintf("M%02d", i),
    ensaio = sprintf("M%02d-%d", i, 1:m),
    n = n,
    cv = cv,
    prod_ctrl = prod_ctrl,
    prod_trat = round(prod_ctrl + efeito_estudo[i] + rnorm(m, 0, 250) +
                        rnorm(m, 0, cv / 100 * prod_ctrl / sqrt(n)))
  )
}) |>
  mutate(media_geral = round((prod_trat + prod_ctrl) / 2, 1))

# 6. Viés de publicação (estudos pequenos "não significativos" faltando) ------
# Semente própria para que a assimetria do funil fique visível no exemplo.
set.seed(3)
k_total <- 150
vies <- tibble(
  n = round(exp(runif(k_total, log(3), log(200)))),
  vi = 2 * (0.10 * 3000)^2 / n
) |>
  mutate(
    yi = rnorm(k_total, 60, 20) + rnorm(k_total, 0, sqrt(vi)),
    p = 2 * pnorm(-abs(yi / sqrt(vi))),
    publicado = p < 0.05 | runif(k_total) < 0.10
  ) |>
  filter(publicado) |>
  mutate(estudo = sprintf("P%02d", row_number()),
         yi = round(yi, 1), vi = round(vi, 1)) |>
  select(estudo, n, yi, vi)

# Salva -----------------------------------------------------------------------
dir.create("dados", showWarnings = FALSE)
write_csv(fungicida_soja, "dados/fungicida_soja.csv")
write_csv(giberela, "dados/giberela_tebuconazol.csv")
write_csv(tombamento, "dados/tombamento_sementes.csv")
write_csv(ferrugem_parcelas, "dados/ferrugem_parcelas.csv")
write_csv(ensaios_multinivel, "dados/ensaios_multinivel.csv")
write_csv(vies, "dados/vies_publicacao.csv")

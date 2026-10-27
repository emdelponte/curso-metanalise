# Gera dados SIMULADOS no desenho de Nkebiwe et al. (2016), Field Crops
# Research 196:389-401 (posicionamento de fertilizantes e produtividade).
# Os valores NÃO são os dados do artigo (o apêndice online não estava
# acessível): o número de estudos (39), de conjuntos de dados (772) e a
# distribuição por cultura, técnica de posicionamento (e, portanto,
# profundidade) e tipo de fertilizante seguem o artigo, e os efeitos médios
# foram ajustados para ficar próximos das Figuras 2 a 5.
#
# Saída: dados/fertilizantes_posicionamento.csv (uma linha = um "conjunto de
# dados": médias, DP e n do fertilizante posicionado e da aplicação a lanço).

library(tidyverse)
set.seed(2016)

# Culturas: n de estudos, n de conjuntos de dados, razão de resposta
# (posicionado / lanço) aproximada da Figura 3 e produtividade-base (kg/ha)
culturas <- tribble(
  ~cultura,           ~n_estudos, ~n_dados, ~rr,   ~base,
  "milho",            12,         408,      1.040, 9000,
  "trigo de inverno",  7,         112,      1.090, 7000,
  "soja",              4,         134,      1.000, 3200,
  "nabo forrageiro",   2,          40,      1.000, 3000,
  "beterraba",         2,          20,      1.035, 70000,
  "centeio",           2,          16,      1.000, 6000,
  "colza",             2,          12,      1.250, 3500,
  "batata",            1,           6,      1.220, 40000,
  "arroz",             1,           6,      1.130, 7000,
  "sorgo",             1,           4,      1.120, 5000,
  "repolho chines",    1,           4,      1.090, 50000,
  "gramineas",         1,           3,      1.000, 8000,
  "alface",            1,           2,      1.000, 30000,
  "trigo de primavera",1,           3,      1.220, 6000,
  "couve-flor",        1,           2,      1.000, 25000
)

# Técnicas (a profundidade decorre da técnica, como no artigo)
tecnicas <- tribble(
  ~tecnica,                    ~profundidade, ~n_dados, ~rr,
  "banda superficial",         "0 cm",        101,      1.040,
  "semente",                   "5-10 cm",      42,      1.000,
  "abaixo da semente",         "5-10 cm",       3,      1.070,
  "banda rasa",                "5-10 cm",     215,      1.030,
  "sulco rasa (faca)",         "5-10 cm",      10,      1.060,
  "ponto raso",                "5-10 cm",       2,      1.300,
  "injecao rasa",              "5-10 cm",      37,      1.010,
  "ponto profundo",            "5-10 cm",       8,      1.190,
  "banda profunda",            ">10 cm",      202,      1.020,
  "sulco profunda (faca)",     ">10 cm",       38,      1.190,
  "injecao profunda",          ">10 cm",      114,      1.060
)

# Tipos de fertilizante (Figura 5 e texto)
fertilizantes <- tribble(
  ~fertilizante,         ~n_dados, ~rr,
  "ureia + P soluvel",       12,   1.273,
  "amonio + P soluvel",     163,   1.147,
  "esterco liquido",         24,   1.116,
  "ureia",                   64,   1.110,
  "amonio",                 185,   1.038,
  "potassio",               134,   1.016,
  "P soluvel",              136,   1.000,
  "N sem descricao",         48,   1.014,
  "esterco solido",           6,   1.079
)

stopifnot(sum(culturas$n_dados) == 772, sum(tecnicas$n_dados) == 772,
          sum(fertilizantes$n_dados) == 772, sum(culturas$n_estudos) == 39)

mu_geral <- log(1.037)
# Efeito aditivo (em ln R) de cada nível, centrado na média ponderada
centra <- function(tab) log(tab$rr) - weighted.mean(log(tab$rr), tab$n_dados)
culturas$delta      <- centra(culturas)
tecnicas$delta      <- centra(tecnicas)
fertilizantes$delta <- centra(fertilizantes)

# Estudos: cultura, ano e número de conjuntos de dados
estudos <- culturas |>
  uncount(n_estudos, .id = "k") |>
  group_by(cultura) |>
  mutate(n_estudo = {
    # divide os conjuntos da cultura entre seus estudos (mínimo 1 cada)
    m <- n(); tot <- first(n_dados)
    if (m == 1) tot else 1 + as.vector(rmultinom(1, tot - m, rgamma(m, 0.7)))
  }) |>
  ungroup() |>
  mutate(estudo = sprintf("N%02d", sample(n())),
         ano = c(sample(1982:1999, 6, replace = TRUE),
                 sample(2000:2015, 33, replace = TRUE))[sample(n())],
         b_estudo = rnorm(n(), 0, 0.03),
         dp_reportado = FALSE) |>
  select(estudo, cultura, delta_cult = delta, base, ano, n_estudo, b_estudo)

# Estudos sem DP reportado (o DP será imputado na análise): ~42% dos dados
sem_dp <- character(0); total <- 0
for (i in sample(nrow(estudos))) {
  if (total + estudos$n_estudo[i] <= 328) {
    sem_dp <- c(sem_dp, estudos$estudo[i]); total <- total + estudos$n_estudo[i]
  }
}

dados <- estudos |>
  uncount(n_estudo) |>
  mutate(tecnica = sample(rep(tecnicas$tecnica, tecnicas$n_dados)),
         fertilizante = sample(rep(fertilizantes$fertilizante, fertilizantes$n_dados))) |>
  left_join(select(tecnicas, tecnica, profundidade, delta_tec = delta), by = "tecnica") |>
  left_join(select(fertilizantes, fertilizante, delta_fert = delta), by = "fertilizante") |>
  mutate(
    n_trat = sample(3:8, n(), replace = TRUE, prob = c(.15, .35, .25, .15, .06, .04)),
    n_ctrl = n_trat,
    cv = pmin(pmax(rlnorm(n(), log(0.11), 0.35), 0.03), 0.40),
    lnR = mu_geral + delta_cult + delta_tec + delta_fert + b_estudo + rnorm(n(), 0, 0.02),
    # recentra para que o efeito médio ponderado (pesos de efeitos aleatórios
    # aproximados) seja o do artigo, apesar do ruído dos estudos grandes
    peso = 1 / (2 * cv^2 / n_trat + 0.01),
    lnR = lnR - weighted.mean(lnR, peso) + log(1.029),
    media_ctrl_verd = base * exp(rnorm(n(), 0, 0.15)),
    media_trat_verd = media_ctrl_verd * exp(lnR),
    media_lanco = rnorm(n(), media_ctrl_verd, cv * media_ctrl_verd / sqrt(n_ctrl)),
    media_posic = rnorm(n(), media_trat_verd, cv * media_trat_verd / sqrt(n_trat)),
    dp_lanco = cv * media_ctrl_verd * sqrt(rchisq(n(), n_ctrl - 1) / (n_ctrl - 1)),
    dp_posic = cv * media_trat_verd * sqrt(rchisq(n(), n_trat - 1) / (n_trat - 1)),
    sem = estudo %in% sem_dp,
    dp_lanco = if_else(sem, NA_real_, dp_lanco),
    dp_posic = if_else(sem, NA_real_, dp_posic)
  ) |>
  arrange(estudo) |>
  group_by(estudo) |>
  mutate(conjunto = sprintf("%s-%03d", estudo, row_number())) |>
  ungroup() |>
  transmute(estudo, conjunto, ano, cultura, fertilizante, tecnica, profundidade,
            n_posic = n_trat, media_posic = round(media_posic, 1),
            dp_posic = round(dp_posic, 1),
            n_lanco = n_ctrl, media_lanco = round(media_lanco, 1),
            dp_lanco = round(dp_lanco, 1))

dir.create("dados", showWarnings = FALSE)
write_csv(dados, "dados/fertilizantes_posicionamento.csv", na = "")

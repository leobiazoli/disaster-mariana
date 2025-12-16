# Generalized synthetic control for birth rate

# Packages -----------------------------------------------------------------
library(stringr)
library(dplyr)
library(MatchIt)
library(tidyr)
library(gsynth)
library(MatchIt)
library(ggplot2)
library(tidyverse)
library(fect)
library(boot)

# Reading data -----------------------------------------------------------------

df_cs <- read.csv("dados_biazoli.csv", header = T, sep = ',', encoding = 'UTF-8')
df_psm <- read.csv("dados_biazoli_psm.csv", header = T, sep = ',', encoding = 'UTF-8')

# merging data

sim <- read_delim(
  "sim_mg_2000a2020.csv",
  delim = ";",
  skip = 3,
  locale = locale(encoding = "LATIN1")
) |>
  slice(1:(n() - 10)) |>
  select(-Total) |>
  slice(-1) 

sim <- sim %>%
  mutate(across(2:22, ~ ifelse(.x == "-", 0, .x))) %>%
  mutate(across(2:22, as.numeric))

sinasc <- read_delim(
  "sinasc_mg_2000a2020.csv",
  delim = ";",
  skip = 3,
  locale = locale(encoding = "LATIN1")
) |>
  slice(1:(n() - 12)) |>
  select(-Total) |>
  slice(-1) 

sim_long <- sim %>%
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "ano",
    values_to = "mortes_infantis"
  ) 

sinasc_long <- sinasc %>%
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "ano",
    values_to = "nascimentos"
  ) 

dados_juntos <- left_join(sim_long, sinasc_long, by=c("Município", "ano")) |>
  separate(Município, into = c("cod", "mun"), sep = " ", extra = "merge")

dados_juntos <- dados_juntos |>
  mutate(cod = as.character(cod), ano = as.numeric(ano)) |>
  mutate(mun = str_to_lower(mun))

glimpse(dados_juntos)

data_final <- left_join(df_psm, dados_juntos, by=c("ano", "mun"))
  
# GSC analysis ----------------------------------------

mg3 = data_final |>
  drop_na(part_agro, part_industria, part_servicos,
          part_administracao, taxa_estab, mortes_infantis, nascimentos, dens_pop) |>
  mutate(tx_nasc = (nascimentos/pop)*1000,
         tx_mort = (mortes_infantis/pop)*1000) |> # dividido por nascimentos
  mutate(treat = ifelse(ano >= 2015 & afetadas == 1, 1, 0))

# PSM year 2015 (tx_nasc)
mg3_15 <- mg3 |>
  filter(ano == 2015)

model_psm = matchit(afetadas ~ part_industria + part_servicos + pop + PIBreal_pc + tx_nasc,
                 method = "nearest", data = mg3_15, ratio = 3)

summary(model_psm)

pards = match.data(model_psm)

mun_psm <- pards |>
  pull(mun)

data_psm <- mg3 |>
  filter(mun %in% mun_psm)

# para tx_mort
model_psm2 = matchit(afetadas ~ pop + PIBreal_pc + tx_mort + trab_analfabeto + trab_sup + mortes_infantis,
                    method = "nearest", data = mg3_15, ratio = 3)

pards2 = match.data(model_psm2)

mun_psm2 <- pards2 |>
  pull(mun)

data_psm2 <- mg3 |>
  filter(mun %in% mun_psm2)

# Model 1
out_psm2 <- gsynth(tx_nasc ~ treat + part_industria + PIBreal_pc + part_servicos, 
                   data = data_psm2,  index = c("mun","ano"),
                 estimator = "mc", nlambda = 5,
                 se = T, inference = "nonparametric", # se = TRUE inf = parametric
                 r = c(0, 5), CV = T, force = "unit", # CV = T
                 parallel = TRUE, min.T0 = 5, 
                 nboots = 100, seed = 02139) # nboots = 1000

out_psm3 <- gsynth(tx_mort ~ treat + part_industria + PIBreal_pc + part_servicos, 
                   data = data_psm,  index = c("mun","ano"),
                   estimator = "mc", nlambda = 5,
                   se = T, inference = "nonparametric", # se = TRUE inf = parametric
                   r = c(0, 5), CV = T, force = "unit", # CV = T
                   parallel = TRUE, min.T0 = 5, 
                   nboots = 100, seed = 02139) # nboots = 1000

plot(out_psm2, type = "gap")

plot(out_psm2, type = "counterfactual", raw = "none", main="")

plot(out_psm2, type = "counterfactual", raw = "all")

plot(out_psm2, type = "counterfactual", id = "mariana")

out_ub <- fect(tx_nasc ~ treat + part_industria + PIBreal_pc + part_servicos, 
               data = data_psm_t,  index = c("mun","ano"),
               se = T, method = "mc", nlambda = 5, 
               CV = T, force = "unit", vartype = "bootstrap",
               parallel = T, min.T0 = 5, 
               nboots = 100, seed = 02139) 

plot(out_ub, show.points = FALSE)

plot(out_ub, type = "equiv", ylim = c(-1, 1), show.stats =  FALSE, cex.legend = 0.6,
     cex.text = 0.5, main = "", ylab = "Effect")

#plot(out_psm2, type = "counterfactual", raw = "band")

plot(out_ub, type = "box")

out_pc <- fect(tx_nasc ~ treat + part_industria + PIBreal_pc + part_servicos, 
               data = data_psm_t,  index = c("mun","ano"),
               se = T, method = "mc", nlambda = 5, 
               CV = 0, force = "unit", placeboTest = T, placebo.period = c(-2, 0),
               parallel = T, min.T0 = 5, vartype = "bootstrap",
               nboots = 10, seed = 02139) 

plot(out_pc, show.points = FALSE, main = "", ylab = "Effect", cex.legend = 0.6, cex.text = 0.65)

round(out_psm2$est.att, 2)

# cummulative effect
cumu2 <- cumuEff(out_psm2, cumu = TRUE, id = NULL)
cumu2$est.catt

# Model 2 - Best model-------------
model_psm_teste= matchit(afetadas ~ part_industria + part_servicos + part_agro + part_administracao + 
                      pop + PIBreal_pc + tx_nasc + trab_analfabeto + trab_sup,
                    method = "nearest", data = mg3_15, ratio = 3)

summary(model_psm_teste)

pards_t = match.data(model_psm_teste)

mun_psm_t <- pards_t |>
  pull(mun)

data_psm_t <- mg3 |>
  filter(mun %in% mun_psm_t)


out_psm2_t <- gsynth(tx_nasc ~ treat + part_industria + PIBreal_pc + part_servicos, 
                   data = data_psm_t,  index = c("mun","ano"),
                   estimator = "mc", nlambda = 5,
                   se = T, inference = "nonparametric", # se = TRUE inf = parametric
                   r = c(0, 5), CV = T, force = "unit", # CV = T
                   parallel = TRUE, min.T0 = 5, 
                   nboots = 100, seed = 02139) # nboots = 1000

plot(out_psm2_t, type = "gap")#, id = "mariana", main="")

plot(out_psm2_t, type = "counterfactual", raw = "none", main="")

plot(out_psm2_t, type = "counterfactual", raw = "all")

plot(out_psm2_t, type = "counterfactual", id = "mariana", main="")

out_psm2_t$est.att

mean(out_psm2_t$Y.bar[,1])

cumu2 <- cumuEff(out_psm2_t, cumu = TRUE, id = NULL)
cumu2$est.catt

# Descriptive analysis --------------------------------------------------
data_psm_t %>% 
  filter(afetadas == 1) %>% 
  ggplot(aes(x = ano, y = tx_nasc, group = mun,
             color = ifelse(mun == "mariana", "Mariana", "Outros"))) +
  geom_line() +
  geom_vline(xintercept = 2015, color = "black", 
             linetype = "dashed") +
  scale_color_manual(values = c("Mariana" = "red", "Outros" = "gray90")) +
  labs(
    title = "Rompimento da barragem do Fundão em 2015",
    y = "Taxa de nascimento",
    x = NULL,
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = "none")

data_psm_t %>% 
  ggplot(aes(x = ano, y = tx_nasc, group = mun,
             color = ifelse(afetadas == 1, "Afetadas", "Não afetadas"))) +
  geom_line() +
  geom_vline(xintercept = 2015, color = "orange", 
             linetype = "dashed") +
  scale_color_manual(values = c("Afetadas" = "red", "Não afetadas" = "gray90")) +
  labs(
    title = "Evolução da taxa de nascimento — municípios afetados e não afetados",
    subtitle = "Rompimento da barragem do Fundão em 2015",
    y = "Taxa de nascimento",
    x = NULL,
    color = NULL
  ) +
  theme_classic()

data_psm_t %>% 
  group_by(ano, afetadas) %>% 
  summarise(media_tx_nasc = mean(tx_nasc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = ano, y = media_tx_nasc, 
             color = factor(afetadas, labels = c("Não afetadas", "Afetadas")))) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 2015, color = "orange", linetype = "dashed") +
  scale_color_manual(values = c("Não afetadas" = "gray60", "Afetadas" = "red")) +
  labs(
    #title = "Trajetória média da taxa de nascimento",
    #subtitle = "Municípios afetados e não afetados — rompimento da barragem de Fundão (2015)",
    y = "Taxa média de nascimento",
    x = NULL,
    color = NULL
  ) +
  theme_classic() 


# Calcular a média dos municípios afetados
media_afetadas <- data_psm_t %>%
  filter(afetadas == 1) %>%
  group_by(ano) %>%
  summarise(media_tx_nasc = mean(tx_nasc, na.rm = TRUE)) %>%
  ungroup()

# Gráfico
data_psm_t %>% 
  filter(afetadas == 1) %>% 
  ggplot() +
  # Linhas individuais
  geom_line(aes(x = ano, y = tx_nasc,
                group = mun,
                color = ifelse(mun == "mariana", "Mariana", "Outros"))) +
  # Linha média
  geom_line(data = media_afetadas, aes(x = ano, y = media_tx_nasc),
            color = "blue", size = 1.2) +
  geom_vline(xintercept = 2015, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("Mariana" = "red", "Outros" = "gray90")) +
  labs(
    #title = "Rompimento da barragem do Fundão em 2015",
    y = "crude birth rate",
    x = "year",
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = "none")

# grafico do artigo:

# Médias por grupo
media_afetadas <- data_psm_t %>%
  filter(afetadas == 1) %>%
  group_by(ano) %>%
  summarise(media_tx_nasc = mean(tx_nasc, na.rm = TRUE)) %>%
  mutate(grupo = "Affected")

media_nao_afetadas <- data_psm_t %>%
  filter(afetadas == 0) %>%
  group_by(ano) %>%
  summarise(media_tx_nasc = mean(tx_nasc, na.rm = TRUE)) %>%
  mutate(grupo = "Non-affected")

# Gráfico
data_psm_t %>%
  filter(afetadas == 1) %>%
  ggplot() +
  # Linhas individuais dos afetados
  geom_line(aes(x = ano, y = tx_nasc, group = mun,
                color = ifelse(mun == "mariana", "Mariana", "Other affected")),
            size = 0.7, alpha = 0.7) +
  # Linha média afetadas
  geom_line(data = media_afetadas,
            aes(x = ano, y = media_tx_nasc, color = grupo),
            size = 1.2) +
  # Linha média não afetadas
  geom_line(data = media_nao_afetadas,
            aes(x = ano, y = media_tx_nasc, color = grupo),
            size = 1.2) +
  # Linha vertical (evento)
  geom_vline(xintercept = 2015, color = "black", linetype = "dashed") +
  # Escala de cores
  scale_color_manual(values = c(
    "Mariana" = "red",
    "Other affected" = "gray80",
    "Affected" = "blue",
    "Non-affected" = "black"
  )) +
  # Rótulos e tema
  labs(
    y = "Crude birth rate",
    x = "Year",
    color = "Trajectory type"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )


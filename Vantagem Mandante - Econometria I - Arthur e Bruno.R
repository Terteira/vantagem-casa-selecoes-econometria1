# =============================================================================
# Econometria I - FGV EESP - 2º semestre de 2026
# Professores: Gabriela do Couto Fonseca e Leon Labre
#
# Trabalho: Artigo Empírico
# Integrantes: Arthur Souza de Almeida Taveira
#              Bruno Lopes Peixoto
#
# Pergunta de pesquisa: Jogar em casa aumenta a probabilidade de vitória e o
# saldo de gols em partidas internacionais de futebol entre seleções nacionais?
#
# Fontes dos dados:
#  - International football results from 1872 to 2026, de Mart Jürisoo.
#    https://github.com/martj42/international_results
#  - GeoNames (cities500 e countryInfo), licença CC BY 4.0, para a altitude
#    das cidades. https://www.geonames.org
#
# Uso de IA: conforme as Orientações da disciplina; link da conversa no relatório.
# =============================================================================

# Limpeza do Ambiente
rm(list=ls())
gc()
graphics.off()
cat("\014")

# Define o Diretório a Ser Usado
setwd("D:/FGV/4o Semestre/Econometria I/Trabalho")   # ajustar em outra máquina

# Pacotes (instalar uma vez):
# install.packages(c("sandwich", "lmtest", "marginaleffects"))
library(sandwich)          # erros-padrão robustos / cluster
library(lmtest)            # coeftest()
library(marginaleffects)   # efeito marginal médio do probit

# -----------------------------------------------------------------------------
# 1. Download dos dados (só baixa o que falta; TRUE baixa tudo de novo)
# -----------------------------------------------------------------------------
atualizar_dados <- FALSE
options(timeout = max(600, getOption("timeout")))   # até 10 min por arquivo

gh  <- "https://raw.githubusercontent.com/martj42/international_results/master/"
geo <- "https://download.geonames.org/export/dump/"
urls <- c("results.csv"      = paste0(gh, "results.csv"),
          "goalscorers.csv"  = paste0(gh, "goalscorers.csv"),
          "shootouts.csv"    = paste0(gh, "shootouts.csv"),
          "former_names.csv" = paste0(gh, "former_names.csv"),
          "cities500.zip"    = paste0(geo, "cities500.zip"),     # ~10 MB
          "countryInfo.txt"  = paste0(geo, "countryInfo.txt"))
pasta <- "dados/brutos"
dir.create(pasta, recursive = TRUE, showWarnings = FALSE)

# Registro do download: uma linha por arquivo, com data e tamanho
log_arq <- file.path(pasta, "log_download.csv")
log_dl  <- if (file.exists(log_arq)) read.csv(log_arq) else
  data.frame(arquivo = character(), data = character(), bytes = numeric())

for (arq in names(urls)) {
  destino <- file.path(pasta, arq)
  if (atualizar_dados || !file.exists(destino)) {
    tmp <- paste0(destino, ".part")                  # baixa num arquivo temporário
    ok  <- tryCatch(download.file(urls[arq], tmp, mode = "wb", quiet = TRUE) == 0,
                    error = function(e) FALSE)
    if (!ok || !file.exists(tmp) || file.size(tmp) == 0) {
      unlink(tmp)
      stop("Falha no download de ", arq, ". Verifique a conexão e rode de novo.")
    }
    file.copy(tmp, destino, overwrite = TRUE)        # só vira o arquivo final
    unlink(tmp)                                      # se o download terminou
    log_dl <- rbind(log_dl[log_dl$arquivo != arq, ],
                    data.frame(arquivo = arq, data = as.character(Sys.Date()),
                               bytes = file.size(destino)))
    write.csv(log_dl, log_arq, row.names = FALSE)
  }
}

# Descompacta se o .txt não existe ou se o .zip é mais novo que ele
zip_cid <- file.path(pasta, "cities500.zip")
txt_cid <- file.path(pasta, "cities500.txt")
if (!file.exists(txt_cid) || file.mtime(zip_cid) > file.mtime(txt_cid))
  unzip(zip_cid, exdir = pasta)


# -----------------------------------------------------------------------------
# 2. Importação
# -----------------------------------------------------------------------------
ler <- function(arq) read.csv(file.path(pasta, arq), fileEncoding = "UTF-8")
jogos    <- ler("results.csv")
gols     <- ler("goalscorers.csv")
disputas <- ler("shootouts.csv")
nomes    <- ler("former_names.csv")

# Se o encoding falhar, read.csv para no 1º acento e lê só parte da base
stopifnot(nrow(jogos) > 49000)

# GeoNames: arquivos separados por tabulação, sem cabeçalho. Lemos só as
# colunas usadas ("NULL" pula as demais). na.strings = "" é essencial: o
# código ISO da Namíbia é "NA" e seria lido como ausente.
# Coluna 17 = "dem" (modelo digital de elevação, SRTM/GTOPO30), preferida à
# coluna 16 ("elevation"), que tem muitos vazios. -9999 = sem dado (oceano).
col <- rep("NULL", 19)
col[c(2, 3, 4, 9, 15, 17)] <- c("character", "character", "character",
                                "character", "numeric", "numeric")
cidades <- read.delim(txt_cid, header = FALSE, quote = "", comment.char = "",
                      na.strings = "", colClasses = col, encoding = "UTF-8")
names(cidades) <- c("nome", "ascii", "alternativos", "iso", "pop", "altitude")

# Altitude -9999 = a célula do modelo de relevo caiu no mar (cidade costeira).
# Não é falta de dado: a cidade fica ao nível do mar, por isso vira 0.
cidades$altitude[cidades$altitude < -1000] <- 0

# countryInfo: coluna 1 = ISO, 5 = nome do país, 6 = capital
paises <- read.delim(file.path(pasta, "countryInfo.txt"), header = FALSE, quote = "",
                     comment.char = "#", na.strings = "", colClasses = "character",
                     encoding = "UTF-8")[, c(1, 5, 6)]
names(paises) <- c("iso", "pais", "capital")


# -----------------------------------------------------------------------------
# 3. Tratamento das bases
# -----------------------------------------------------------------------------
# Corrige tipos, padroniza nomes, liga as tabelas e cria marcadores de
# qualidade. Não cria variáveis do modelo nem exclui jogos da amostra.
diag3 <- list(jogos_na_base_original = nrow(jogos))

# 3.1 Funções auxiliares (usadas nas seções 3 e 4) -----------------------------
# Altitude de uma cidade dentro do(s) país(es) indicado(s) (códigos ISO
# separados por espaço). Procura pelo nome, pelo nome sem acento e pelos nomes
# alternativos, sem diferenciar maiúsculas; entre homônimos, fica a mais populosa.
achar_altitude <- function(cidade, isos, por_pais) {
  if (is.na(isos)) return(NA_real_)
  c0 <- do.call(rbind, por_pais[strsplit(isos, " ")[[1]]])
  if (is.null(c0)) return(NA_real_)
  c0 <- c0[order(-c0$pop), ]
  h <- which(tolower(c0$nome) == tolower(cidade) | tolower(c0$ascii) == tolower(cidade))
  if (length(h) == 0)
    h <- which(grepl(paste0(",", tolower(cidade), ","),
                     paste0(",", tolower(c0$alternativos), ","), fixed = TRUE))
  if (length(h) > 0) c0$altitude[h[1]] else NA_real_
}

# Quantas vezes cada chave de 'alvo' aparece em 'chaves' (0 se nenhuma)
contar <- function(chaves, alvo) {
  x <- as.integer(table(chaves)[alvo])
  ifelse(is.na(x), 0L, x)
}

# Diferença de Elo antes de cada jogo (mandante - visitante, sem o bônus).
# d precisa estar em ordem cronológica e ter home_team, away_team, saldo, casa
# e torneio; K é o fator por tipo de torneio.
calc_elo <- function(d, bonus, K) {
  times <- unique(c(d$home_team, d$away_team))
  elo <- rep(1500, length(times)); nj <- rep(0, length(times))
  im <- match(d$home_team, times); iv <- match(d$away_team, times)
  K <- K[as.character(d$torneio)]
  n <- abs(d$saldo); G <- ifelse(n <= 1, 1, ifelse(n == 2, 1.5, (11 + n) / 8))
  res <- (sign(d$saldo) + 1) / 2
  dif <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    m <- im[i]; v <- iv[i]
    dif[i] <- elo[m] - elo[v]                                  # antes do jogo
    E <- 1 / (10^(-(dif[i] + bonus * d$casa[i]) / 400) + 1)
    elo[m] <- elo[m] + K[i] * G[i] * (res[i] - E) * ifelse(nj[m] < 20, 2, 1)
    elo[v] <- elo[v] - K[i] * G[i] * (res[i] - E) * ifelse(nj[v] < 20, 2, 1)
    nj[c(m, v)] <- nj[c(m, v)] + 1
  }
  dif
}

# Média dos últimos j valores ANTERIORES (NA na primeira posição)
media_ant <- function(x, j) {
  s <- c(0, cumsum(x)); k <- seq_along(x) - 1; ini <- pmax(k - j, 0)
  ifelse(k == 0, NA, (s[k + 1] - s[ini + 1]) / (k - ini))
}

# Mediana dos valores não ausentes ANTERIORES (NA enquanto não há nenhum)
mediana_ant <- function(x) {
  out <- rep(NA_real_, length(x)); vistos <- numeric(0)
  for (i in seq_along(x)) {
    if (length(vistos) > 0) out[i] <- median(vistos)
    if (!is.na(x[i])) vistos <- c(vistos, x[i])
  }
  out
}

# Aplica f a x separadamente em cada grupo, preservando a ordem
por_time <- function(x, grupo, f, ...) ave(x, grupo, FUN = function(v) f(v, ...))

# Diferença mandante - visitante de uma variável com sufixos _m e _v
dif_mv <- function(base, var) base[[paste0(var, "_m")]] - base[[paste0(var, "_v")]]

# 3.2 Tipos e ordem cronológica -------------------------------------------------
jogos$data    <- as.Date(jogos$date)
gols$data     <- as.Date(gols$date)
disputas$data <- as.Date(disputas$date)
stopifnot(!anyNA(jogos$data), is.logical(jogos$neutral),
          is.numeric(jogos$home_score), is.numeric(jogos$away_score))

sem_placar <- is.na(jogos$home_score) | is.na(jogos$away_score)   # jogos ainda não disputados
diag3$jogos_sem_placar_removidos <- sum(sem_placar)
jogos <- jogos[!sem_placar, ]
jogos <- jogos[order(jogos$data, seq_len(nrow(jogos))), ]           # empate na data: ordem do arquivo
jogos$id_jogo <- seq_len(nrow(jogos))
rownames(jogos) <- NULL
disputas$first_shooter[disputas$first_shooter == ""] <- NA

# 3.3 Chave de ligação (data + mandante + visitante) ---------------------------
# Chaves duplicadas são mantidas e marcadas; a ligação só é segura se elas não
# aparecem nos outros arquivos.
jogos$chave    <- paste(jogos$date, jogos$home_team, jogos$away_team)
gols$chave     <- paste(gols$date, gols$home_team, gols$away_team)
disputas$chave <- paste(disputas$date, disputas$home_team, disputas$away_team)
jogos$chave_duplicada <- duplicated(jogos$chave) | duplicated(jogos$chave, fromLast = TRUE)
stopifnot(!any(gols$chave %in% jogos$chave[jogos$chave_duplicada]),
          !any(disputas$chave %in% jogos$chave[jogos$chave_duplicada]))
diag3$jogos_com_chave_duplicada <- sum(jogos$chave_duplicada)

# 3.4 País-sede com o nome atual ------------------------------------------------
# former_names traduz nomes antigos de seleções; usamos só o nome, sem as datas.
jogos$pais_sede <- ifelse(jogos$country %in% nomes$former,
                          nomes$current[match(jogos$country, nomes$former)], jogos$country)

# Estados que se dividiram: em jogo não neutro, o país-sede passa a ser o mandante.
sucessores <- list(
  "Soviet Union"   = c("Russia", "Ukraine", "Belarus", "Georgia", "Armenia", "Azerbaijan",
                       "Kazakhstan", "Uzbekistan", "Turkmenistan", "Kyrgyzstan",
                       "Tajikistan", "Moldova", "Lithuania", "Latvia", "Estonia"),
  "Yugoslavia"     = c("Serbia", "Croatia", "Slovenia", "Bosnia and Herzegovina",
                       "North Macedonia", "Montenegro", "Kosovo"),
  "Czech Republic" = "Czechoslovakia",
  "Slovakia"       = "Czechoslovakia")
for (pais in names(sucessores)) {
  ajusta <- jogos$country == pais & !jogos$neutral & jogos$home_team %in% sucessores[[pais]]
  jogos$pais_sede[ajusta] <- jogos$home_team[ajusta]
}
diag3$nao_neutros_mandante_fora_do_pais <- sum(!jogos$neutral & jogos$home_team != jogos$pais_sede)
diag3$neutros_mandante_no_proprio_pais  <- sum(jogos$neutral & jogos$home_team == jogos$pais_sede)

# 3.5 País-sede -> código ISO ---------------------------------------------------
# Tabela manual para o que não casa com o countryInfo; Estados que se dividiram
# recebem todos os códigos sucessores.
iso_manual <- c(
  # Reino Unido e Irlanda
  "England" = "GB", "Scotland" = "GB", "Wales" = "GB", "Northern Ireland" = "GB",
  "Ireland" = "GB IE", "Alderney" = "GG", "Republic of Ireland" = "IE", "Éire" = "IE",
  "Irish Free State" = "IE",
  # Estados que se dividiram ou mudaram de nome
  "Soviet Union" = "RU UA BY GE AM AZ KZ UZ TM KG TJ MD LT LV EE",
  "Yugoslavia" = "RS HR SI BA MK ME XK", "FR Yugoslavia" = "RS ME XK",
  "Serbia and Montenegro" = "RS ME XK", "Czechoslovakia" = "CZ SK",
  "Bohemia" = "CZ", "Bohemia and Moravia" = "CZ", "German DR" = "DE", "Saarland" = "DE",
  "United Arab Republic" = "EG SY", "Ethiopia" = "ET ER", "Pakistan" = "PK BD",
  "Malaya" = "MY SG", "Manchuria" = "CN", "Burma" = "MM", "Ceylon" = "LK",
  "North Vietnam" = "VN", "Vietnam Republic" = "VN", "Yemen AR" = "YE", "Yemen DPR" = "YE",
  "Mandatory Palestine" = "IL PS", "Northern Cyprus" = "CY", "Belgian Congo" = "CD",
  "Congo-Léopoldville" = "CD", "Congo-Kinshasa" = "CD", "Zaïre" = "CD", "Dahomey" = "BJ",
  "Upper Volta" = "BF", "Gold Coast" = "GH", "Portuguese Guinea" = "GW",
  "French Somaliland" = "DJ", "Tanganyika" = "TZ", "Zanzibar" = "TZ",
  "Northern Rhodesia" = "ZM", "Southern Rhodesia" = "ZW", "Nyasaland" = "MW",
  "Swaziland" = "SZ", "British Guiana" = "GY", "Dutch Guyana" = "SR",
  "Dutch East Indies" = "ID", "Netherlands Antilles" = "CW BQ SX AW", "New Hebrides" = "VU",
  "Western Samoa" = "WS", "Tahiti" = "PF", "Macedonia" = "MK",
  # Grafias que diferem do countryInfo
  "Czech Republic" = "CZ", "North Macedonia" = "MK", "Ivory Coast" = "CI",
  "Cape Verde" = "CV", "DR Congo" = "CD", "Congo" = "CG", "Eswatini" = "SZ",
  "Gambia" = "GM", "Macau" = "MO", "Palestine" = "PS", "Timor-Leste" = "TL",
  "São Tomé and Príncipe" = "ST", "Saint Kitts and Nevis" = "KN",
  "Saint Vincent and the Grenadines" = "VC", "United States Virgin Islands" = "VI",
  "Curaçao" = "CW", "Bonaire" = "BQ", "Saint Barthélemy" = "BL", "Réunion" = "RE",
  "Netherlands" = "NL", "Turkey" = "TR", "United States" = "US", "South Korea" = "KR",
  "North Korea" = "KP", "Russia" = "RU", "Iran" = "IR", "Syria" = "SY", "Laos" = "LA",
  "Brunei" = "BN", "Vietnam" = "VN", "Taiwan" = "TW", "Moldova" = "MD",
  "Kosovo" = "XK", "Tanzania" = "TZ", "Bolivia" = "BO", "Venezuela" = "VE",
  "Micronesia" = "FM")

locais <- unique(jogos[, c("city", "country")])
locais$iso <- ifelse(locais$country %in% names(iso_manual), iso_manual[locais$country],
                     paises$iso[match(locais$country, paises$pais)])
sem_iso <- sort(unique(locais$country[is.na(locais$iso)]))
diag3$paises_sede_sem_iso <- length(sem_iso)
if (length(sem_iso) > 0) print(sem_iso)   # acrescentar em iso_manual

# 3.6 Cidade -> altitude (GeoNames) ---------------------------------------------
cidade_manual <- c(
  # grafias e estádios
  "So Kon Po" = "Hong Kong", "Washington, D.C." = "Washington", "Ta' Qali" = "Attard",
  "Macoya" = "Tunapuna", "Paynesville" = "Monrovia", "North Sound" = "Saint John's",
  "Faro-Loulé" = "Loulé", "Dededo" = "Dededo Village", "Breslau" = "Wrocław",
  # bairros e estádios de uma região metropolitana
  "Kasarani" = "Nairobi", "Atteridgeville" = "Pretoria", "Phokeng" = "Rustenburg",
  "Saulspoort" = "Rustenburg", "KwaMashu" = "Durban", "Umlazi" = "Durban",
  "Clermont" = "Durban", "Puerto Ordaz" = "Ciudad Guayana", "Providence" = "Georgetown",
  "Luque" = "Asunción", "Nausori" = "Suva", "Siu Sai Wan" = "Hong Kong",
  "Leonora" = "Georgetown", "Arad" = "Manama", "East Riffa" = "Ar Rifā‘",
  "Paya Rumput" = "Malacca", "Sulaibikhat" = "Kuwait City", "Jurong Stadium" = "Singapore",
  "Palo Seco" = "San Fernando", "Bacolet" = "Scarborough", "Arnos Vale" = "Kingstown",
  "Addu City" = "Hithadhoo", "‘Atele" = "Nuku'alofa", "Yona" = "Hagåtña",
  "Boulari" = "Nouméa",
  "El Alto" = "La Paz",                                  # ausente da GeoNames: cidade vizinha
  "Thailand" = "Bangkok", "New Caledonia" = "Nouméa")    # base traz o país no lugar da cidade
locais$busca <- ifelse(locais$city %in% names(cidade_manual),
                       cidade_manual[locais$city], locais$city)
locais$iso[locais$city == "Breslau"] <- "PL"             # hoje na Polônia

por_pais <- split(cidades, cidades$iso)
locais$altitude <- unname(mapply(achar_altitude, locais$busca, locais$iso,
                                 MoreArgs = list(por_pais = por_pais)))
k_local <- match(paste(jogos$city, jogos$country), paste(locais$city, locais$country))
jogos$alt_sede       <- locais$altitude[k_local]
jogos$cidade_sem_alt <- is.na(jogos$alt_sede)                                  # marcador
diag3$pct_cidades_com_altitude <- round(100 * mean(!is.na(locais$altitude)), 1)
diag3$pct_jogos_com_altitude   <- round(100 * mean(!jogos$cidade_sem_alt), 1)

# Lista das cidades não encontradas (para as limitações)
dir.create("resultados", showWarnings = FALSE)
nao_casadas <- locais[is.na(locais$altitude), c("city", "country", "iso")]
nao_casadas$jogos <- contar(paste(jogos$city, jogos$country),
                            paste(nao_casadas$city, nao_casadas$country))
nao_casadas <- nao_casadas[order(-nao_casadas$jogos), ]
write.csv(nao_casadas, "resultados/cidades_sem_altitude.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

# Cidade não encontrada: mediana das sedes encontradas no mesmo país
med_pais <- tapply(jogos$alt_sede, jogos$country, median, na.rm = TRUE)
jogos$alt_sede[jogos$cidade_sem_alt] <- med_pais[jogos$country[jogos$cidade_sem_alt]]
diag3$jogos_ainda_sem_altitude <- sum(is.na(jogos$alt_sede))

# 3.7 Ligação com goalscorers ---------------------------------------------------
# Gols de pênalti de cada lado; só confiáveis nos jogos com todos os gols registrados.
stopifnot(all(gols$chave %in% jogos$chave))
gols_registrados    <- contar(gols$chave, jogos$chave)
jogos$pen_mandante  <- contar(gols$chave[gols$penalty & gols$team == gols$home_team], jogos$chave)
jogos$pen_visitante <- contar(gols$chave[gols$penalty & gols$team == gols$away_team], jogos$chave)
total_gols <- jogos$home_score + jogos$away_score
stopifnot(all(gols_registrados <= total_gols))

# 0 a 0 conta como completo só em torneio-década com cobertura >= 95%
celula_cobertura <- paste(jogos$tournament, 10 * (as.integer(format(jogos$data, "%Y")) %/% 10))
cobertura <- tapply((gols_registrados == total_gols)[total_gols > 0],
                    celula_cobertura[total_gols > 0], mean)
celulas_cobertas <- names(cobertura)[cobertura >= 0.95]
jogos$gols_completos <- (total_gols > 0 & gols_registrados == total_gols) |
  (total_gols == 0 & celula_cobertura %in% celulas_cobertas)   # marcador
diag3$torneio_decada_cobertos   <- length(celulas_cobertas)
diag3$jogos_com_gols_completos  <- sum(jogos$gols_completos)
diag3$dos_quais_0_a_0           <- sum(jogos$gols_completos & total_gols == 0)

# 3.8 Ligação com shootouts -----------------------------------------------------
n_antes <- nrow(jogos)
p <- match(jogos$chave, disputas$chave)
jogos$disputa          <- !is.na(p)
jogos$disputa_vencedor <- disputas$winner[p]
jogos$disputa_primeiro <- disputas$first_shooter[p]
jogos$ida_e_volta <- jogos$disputa & jogos$home_score != jogos$away_score   # marcador: decidido no agregado
stopifnot(nrow(jogos) == n_antes, sum(jogos$disputa) == sum(disputas$chave %in% jogos$chave))
diag3$disputas_sem_jogo       <- sum(!disputas$chave %in% jogos$chave)
diag3$jogos_com_disputa       <- sum(jogos$disputa)
diag3$disputas_ida_e_volta    <- sum(jogos$ida_e_volta)

# 3.9 Marcadores de torneio e de seleções da FIFA -------------------------------
# Só marcam; os filtros ficam para a seção 5. Membro da FIFA = já disputou Copa
# ou eliminatórias (regra fixa no tempo).
fora_fifa <- "CONIFA|ConIFA|Island Games|Muratti Vase|Viva World Cup|FIFI Wild Cup|ELF Cup|Inter Games|Tynwald Hill|Niamh Challenge|World Unity Cup|Tournament of Peoples|Coupe de l'Outre-Mer|Olympic Games"
jogos$torneio_excluido <- grepl(fora_fifa, jogos$tournament)
copa <- jogos$tournament %in% c("FIFA World Cup", "FIFA World Cup qualification")
membros <- unique(c(jogos$home_team[copa], jogos$away_team[copa]))
jogos$fifa_mandante  <- jogos$home_team %in% membros
jogos$fifa_visitante <- jogos$away_team %in% membros
diag3$selecoes_na_base     <- length(unique(c(jogos$home_team, jogos$away_team)))
diag3$membros_fifa         <- length(membros)
diag3$jogos_torneio_excluido <- sum(jogos$torneio_excluido)
diag3$jogos_selecao_fora_fifa <- sum(!jogos$fifa_mandante | !jogos$fifa_visitante)

# 3.10 Conferências e resumo ----------------------------------------------------
stopifnot(nrow(jogos) == diag3$jogos_na_base_original - sum(sem_placar),
          !anyDuplicated(jogos$id_jogo), !is.unsorted(jogos$data))
diag3$jogos_ao_fim <- nrow(jogos)
print(data.frame(diagnostico = names(diag3), valor = sapply(diag3, format)), row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Criação e Transformação de Variáveis
# -----------------------------------------------------------------------------
# Uma linha = um jogo, visto pelo mandante; diferenças = mandante - visitante.
# Nenhum jogo é excluído aqui.
n_jogos <- nrow(jogos)

# 4.1 Base para o histórico -----------------------------------------------------
# Só jogos entre seleções da FIFA e fora dos torneios excluídos alimentam o
# histórico; os demais ficam com NA nessas variáveis.
hist <- !jogos$torneio_excluido & jogos$fifa_mandante & jogos$fifa_visitante
stopifnot(!is.unsorted(jogos$id_jogo))
diag4 <- list(jogos_no_historico = sum(hist))

# 4.2 Variáveis dependentes, de interesse e tipo de torneio --------------------
# Vêm antes do histórico porque o Elo usa o saldo e o tipo de torneio.
jogos$saldo     <- jogos$home_score - jogos$away_score    # com prorrogação, sem a disputa
jogos$vitoria   <- as.integer(jogos$saldo > 0)            # decidido nos pênaltis = empate
jogos$casa      <- as.integer(!jogos$neutral)
jogos$resultado <- factor(sign(jogos$saldo), levels = c(-1, 0, 1),
                          labels = c("Derrota", "Empate", "Vitória"))   # descritiva

# Tipo de torneio: cada regra sobrescreve as anteriores (a ordem importa)
nome_torneio <- jogos$tournament
finais <- c("UEFA Euro", "Copa América", "African Cup of Nations", "AFC Asian Cup",
            "Gold Cup", "CONCACAF Championship", "Oceania Nations Cup", "Confederations Cup")
jogos$torneio <- "Outros torneios"
jogos$torneio[grepl("qualification", nome_torneio)]           <- "Eliminatórias continentais"
jogos$torneio[grepl("Nations League", nome_torneio)]          <- "Liga das Nações"
jogos$torneio[nome_torneio %in% finais]                       <- "Continental (final)"
jogos$torneio[nome_torneio == "FIFA World Cup qualification"] <- "Eliminatórias Copa"
jogos$torneio[nome_torneio == "FIFA World Cup"]               <- "Copa do Mundo"
jogos$torneio[nome_torneio == "Friendly"]                     <- "Amistoso"
jogos$torneio <- relevel(factor(jogos$torneio), ref = "Amistoso")
jogos$oficial <- as.integer(jogos$torneio != "Amistoso")

# 4.3 Variáveis de histórico (só jogos anteriores) ------------------------------
# Elo: método do eloratings.net, com rating inicial e K provisório nossos
# (detalhes no relatório). O bônus entra na atualização, não no controle.
K_tipo <- c("Copa do Mundo" = 60, "Continental (final)" = 50, "Eliminatórias Copa" = 40,
            "Eliminatórias continentais" = 40, "Liga das Nações" = 40,
            "Outros torneios" = 30, "Amistoso" = 20)
hb <- jogos[hist, ]
for (b in c(100, 0, 150)) {                               # 100 = principal; 0 e 150 = robustez
  nome <- if (b == 100) "dif_elo" else paste0("dif_elo_b", b)
  jogos[[nome]] <- NA
  jogos[[nome]][hist] <- calc_elo(hb, b, K_tipo) / 100
}

# Demais variáveis: cada jogo aparece uma vez por seleção, em ordem cronológica
longo <- rbind(
  data.frame(id = hb$id_jogo, data = hb$data, time = hb$home_team, lado = "m",
             sg = hb$saldo,  alt_casa = ifelse(hb$casa == 1, hb$alt_sede, NA)),
  data.frame(id = hb$id_jogo, data = hb$data, time = hb$away_team, lado = "v",
             sg = -hb$saldo, alt_casa = NA))
longo <- longo[order(longo$time, longo$id), ]
longo$pts      <- ifelse(longo$sg > 0, 3, ifelse(longo$sg == 0, 1, 0))
longo$n_ant    <- por_time(longo$id, longo$time, seq_along) - 1
longo$descanso <- pmin(por_time(as.numeric(longo$data), longo$time, function(v) c(NA, diff(v))), 60)
for (j in c(5, 10, 20)) {
  longo[[paste0("forma_pts_", j)]] <- por_time(longo$pts, longo$time, media_ant, j)
  longo[[paste0("forma_sg_", j)]]  <- por_time(longo$sg,  longo$time, media_ant, j)
}
longo$alt_hab <- por_time(longo$alt_casa, longo$time, mediana_ant)

# Por lado: traz o histórico, completa a altitude habitual e calcula subida e
# descida (km). Sem jogo em casa anterior, a altitude habitual é a mediana das
# sedes do país (exceção à regra de informação anterior; ver relatório).
vars_hist <- c("n_ant", "descanso", "alt_hab", paste0("forma_pts_", c(5, 10, 20)),
               paste0("forma_sg_", c(5, 10, 20)))
alt_pais <- tapply(hb$alt_sede[hb$casa == 1], hb$pais_sede[hb$casa == 1], median)
for (lado in c("m", "v")) {
  x <- longo[longo$lado == lado, ][match(jogos$id_jogo, longo$id[longo$lado == lado]), ]
  for (v in vars_hist) jogos[[paste0(v, "_", lado)]] <- x[[v]]
  time <- if (lado == "m") jogos$home_team else jogos$away_team
  a <- paste0("alt_hab_", lado)
  sem <- hist & is.na(jogos[[a]])
  jogos[[paste0("alt_hab_aprox_", lado)]] <- sem                       # marcador
  jogos[[a]][sem] <- alt_pais[time[sem]]
  d <- jogos$alt_sede - jogos[[a]]
  jogos[[paste0("subida_", lado)]]  <- pmax(d, 0) / 1000
  jogos[[paste0("descida_", lado)]] <- pmax(-d, 0) / 1000
}

# Diferenças mandante - visitante (janela de forma de 10 jogos = principal)
jogos$dif_experiencia <- log1p(jogos$n_ant_m) - log1p(jogos$n_ant_v)
jogos$dif_descanso    <- dif_mv(jogos, "descanso")
for (j in c(5, 10, 20)) {
  s <- if (j == 10) "" else paste0("_", j)
  jogos[[paste0("dif_forma_pts", s)]] <- dif_mv(jogos, paste0("forma_pts_", j))
  jogos[[paste0("dif_forma_sg", s)]]  <- dif_mv(jogos, paste0("forma_sg_", j))
}
# Altitude: visitante - mandante (positivo favorece o mandante)
jogos$dif_subida  <- jogos$subida_v  - jogos$subida_m
jogos$dif_descida <- jogos$descida_v - jogos$descida_m

# 4.4 Demais variáveis ----------------------------------------------------------
ano <- as.integer(format(jogos$data, "%Y"))
jogos$decada_jogo <- factor(10 * (ano %/% 10))
jogos$periodo <- cut(ano, c(-Inf, 1959, 1979, 1999, Inf),
                     labels = c("até 1959", "1960-79", "1980-99", "2000-26"))

# Anfitrião: jogo em casa numa edição (torneio + ano) com sede concentrada, ou
# seja, em que no máximo 1/3 das seleções jogou em casa. Exclui torneios em ida
# e volta e mantém co-sedes. Só vale 1 com casa = 1 (já é a interação).
edicao  <- paste(jogos$tournament, ano)
em_casa <- jogos$casa == 1 & jogos$home_team == jogos$pais_sede
n_times <- tapply(seq_len(n_jogos), edicao, function(i)
  length(unique(c(jogos$home_team[i], jogos$away_team[i]))))
n_casa  <- tapply(seq_len(n_jogos), edicao, function(i)
  length(unique(jogos$home_team[i][em_casa[i]])))
anf <- em_casa & (n_casa / n_times)[edicao] <= 1/3 &
  jogos$torneio %in% c("Copa do Mundo", "Continental (final)", "Outros torneios")
jogos$anf_copa        <- as.integer(anf & jogos$torneio == "Copa do Mundo")
jogos$anf_continental <- as.integer(anf & jogos$torneio == "Continental (final)")
jogos$anf_regional    <- as.integer(anf & jogos$torneio == "Outros torneios")
# Conferência: edições de Copa e de finais continentais sem sede concentrada
edicoes_sem_sede <- sort(unique(edicao[jogos$torneio %in% c("Copa do Mundo", "Continental (final)") &
                                         (n_casa / n_times)[edicao] > 1/3]))

# Confronto: par de seleções em ordem alfabética (efeitos fixos e cluster por par)
jogos$confronto <- ifelse(jogos$home_team < jogos$away_team,
                          paste(jogos$home_team, "-", jogos$away_team),
                          paste(jogos$away_team, "-", jogos$home_team))

# Extensões de pênaltis (ordem de cobrança sorteada até julho de 2003)
jogos$dif_penalti     <- ifelse(jogos$gols_completos, jogos$pen_mandante - jogos$pen_visitante, NA)
jogos$venceu_disputa  <- ifelse(jogos$disputa, as.integer(jogos$disputa_vencedor == jogos$home_team), NA)
jogos$cobrou_primeiro <- as.integer(jogos$disputa_primeiro == jogos$home_team)
jogos$ordem_sorteada  <- jogos$data < as.Date("2003-07-01")

# 4.5 Conferências e resumo -----------------------------------------------------
inv <- !jogos$neutral & jogos$away_team == jogos$pais_sede   # país-sede listado como visitante
stopifnot(
  nrow(jogos) == n_jogos, !anyDuplicated(jogos$id_jogo),
  all(jogos$descanso_m >= 0 & jogos$descanso_m <= 60, na.rm = TRUE),
  all(jogos$descanso_v >= 0 & jogos$descanso_v <= 60, na.rm = TRUE),
  all(c(jogos$anf_copa, jogos$anf_continental, jogos$anf_regional) %in% 0:1),
  all(jogos$casa %in% 0:1), all(jogos$vitoria %in% 0:1),
  !any((jogos$anf_copa | jogos$anf_continental | jogos$anf_regional) & jogos$casa == 0),
  !any(inv & hist),
  all(is.na(longo$forma_pts_10[longo$n_ant == 0])),     # primeiro jogo sem histórico
  all(is.na(longo$descanso[longo$n_ant == 0])))

vars_modelo <- c("dif_elo", "dif_forma_pts", "dif_forma_sg", "dif_experiencia",
                 "dif_descanso", "alt_hab_m", "alt_hab_v", "dif_subida", "dif_descida")
diag4 <- c(diag4,
           as.list(setNames(colSums(is.na(jogos[hist, vars_modelo])), paste0("ausentes_", vars_modelo))),
           list(pais_sede_como_visitante = sum(inv), dos_quais_no_historico = sum(inv & hist),
                dif_elo_min = round(min(jogos$dif_elo, na.rm = TRUE), 1),
                dif_elo_max = round(max(jogos$dif_elo, na.rm = TRUE), 1),
                anfitriao_copa = sum(jogos$anf_copa), anfitriao_continental = sum(jogos$anf_continental),
                anfitriao_regional = sum(jogos$anf_regional),
                altitude_habitual_aproximada = sum(jogos$alt_hab_aprox_m | jogos$alt_hab_aprox_v)))
print(data.frame(diagnostico = names(diag4), valor = sapply(diag4, format)), row.names = FALSE)
cat("Edições de Copa e finais continentais sem sede concentrada:\n"); print(edicoes_sem_sede)


# -----------------------------------------------------------------------------
# 5. Seleção da Amostra Final
# -----------------------------------------------------------------------------
# Filtros em sequência, com a tabela de fluxo. Subamostras de robustez ficam
# para a seção 9.

# 5.1 Amostra principal ---------------------------------------------------------
controles <- c("saldo", "vitoria", "casa", "dif_elo", "dif_forma_pts", "dif_forma_sg",
               "dif_experiencia", "dif_descanso", "torneio", "decada_jogo")
filtros <- list(
  "Base completa"                              = rep(TRUE, nrow(jogos)),
  "Sem torneios excluídos (CONIFA, Olimpíadas)" = !jogos$torneio_excluido,
  "As duas seleções filiadas à FIFA"           = jogos$fifa_mandante & jogos$fifa_visitante,
  "As duas com pelo menos 10 jogos anteriores" = jogos$n_ant_m >= 10 & jogos$n_ant_v >= 10 &
    !is.na(jogos$n_ant_m) & !is.na(jogos$n_ant_v),
  # incoerentes (neutro em casa ou mando fora do país): excluídos, não reclassificados
  "Casa/neutro coerente com o país-sede"       = !(jogos$neutral & jogos$home_team == jogos$pais_sede) &
    !(!jogos$neutral & jogos$home_team != jogos$pais_sede),
  "Variáveis do modelo sem ausentes"           = complete.cases(jogos[, controles]))
acumulado <- Reduce(`&`, filtros, accumulate = TRUE)

fluxo <- data.frame(
  etapa       = names(filtros),
  jogos       = sapply(acumulado, sum),
  removidos   = c(NA, -diff(sapply(acumulado, sum))),
  pct_neutros = sapply(acumulado, function(s) round(100 * mean(jogos$neutral[s]), 1)))
print(fluxo, row.names = FALSE)
write.csv(fluxo, "resultados/fluxo_amostra.csv", row.names = FALSE, fileEncoding = "UTF-8")

am <- droplevels(jogos[acumulado[[length(acumulado)]], ])
rownames(am) <- NULL

# 5.2 Subamostras das extensões -------------------------------------------------
# Marcadores dentro da amostra principal; disputas seguem a regra da amostra.
am$am_altitude <- !is.na(am$dif_subida) & !is.na(am$dif_descida)
am$am_penalti  <- !is.na(am$dif_penalti)
am$am_disputa  <- am$disputa

# 5.3 Conferências e resumo (subamostras, ausentes e outliers) -----------------
# Outliers ficam na amostra; a influência é testada na seção 9 (saldo +-5).
stopifnot(!anyDuplicated(am$id_jogo), !anyNA(am[, controles]),
          all(am$n_ant_m >= 10 & am$n_ant_v >= 10),
          !any(am$torneio_excluido), all(am$fifa_mandante & am$fifa_visitante))
diag5 <- list(
  amostra_principal       = nrow(am),
  pct_da_base             = round(100 * nrow(am) / nrow(jogos), 1),
  primeiro_jogo           = format(min(am$data)),
  ultimo_jogo             = format(max(am$data)),
  ext_altitude            = sum(am$am_altitude),
  altitude_aproximada     = sum(am$am_altitude & (am$cidade_sem_alt | am$alt_hab_aprox_m | am$alt_hab_aprox_v)),
  ext_penaltis            = sum(am$am_penalti),
  ext_disputas            = sum(am$am_disputa),
  disputas_com_ordem      = sum(!is.na(am$cobrou_primeiro)),
  disputas_ordem_sorteada = sum(!is.na(am$cobrou_primeiro) & am$ordem_sorteada),
  saldo_acima_de_5        = sum(abs(am$saldo) > 5),
  chave_duplicada         = sum(am$chave_duplicada),
  ida_e_volta             = sum(am$ida_e_volta),
  cidade_aproximada       = sum(am$cidade_sem_alt),
  sem_altitude            = sum(!am$am_altitude))
print(data.frame(diagnostico = names(diag5), valor = sapply(diag5, format)), row.names = FALSE)

# 5.4 Salvar as bases -----------------------------------------------------------
# Em CSV, fatores viram texto: ao reler, refazer torneio (referência "Amistoso"),
# decada_jogo e periodo.
dir.create("dados/finais", recursive = TRUE, showWarnings = FALSE)
write.csv(jogos, "dados/finais/results_tratado.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(am,    "dados/finais/results_final.csv",   row.names = FALSE, fileEncoding = "UTF-8")


# -----------------------------------------------------------------------------
# 6. Análise descritiva
# -----------------------------------------------------------------------------
# Duas tabelas e duas figuras da amostra principal, salvas em
# resultados/descritivas para o paper. Figuras sem título interno (título e
# fonte vão no texto, ABNT) e com vírgula decimal nos eixos.
pasta_desc <- "resultados/descritivas"
dir.create(pasta_desc, recursive = TRUE, showWarnings = FALSE)
em_casa_am <- am$casa == 1
eixo_virgula <- function(lado) {                     # eixo com vírgula decimal
  marcas <- axTicks(lado)
  axis(lado, at = marcas, labels = format(marcas, decimal.mark = ","), las = 1)
}

# 6.1 Tabela 1: médias em casa e em campo neutro --------------------------------
# Mostra a diferença bruta e o desequilíbrio dos controles entre os dois grupos.
vars_desc <- data.frame(
  "Saldo de gols"               = am$saldo,
  "Vitória do mandante"         = am$vitoria,
  "Empate"                      = as.integer(am$resultado == "Empate"),
  "Derrota do mandante"         = as.integer(am$resultado == "Derrota"),
  "Dif. de Elo (centenas)"      = am$dif_elo,
  "Dif. de forma (pontos)"      = am$dif_forma_pts,
  "Dif. de forma (saldo)"       = am$dif_forma_sg,
  "Dif. de experiência (log)"   = am$dif_experiencia,
  "Dif. de descanso (dias)"     = am$dif_descanso,
  check.names = FALSE)
tab1 <- data.frame(
  variavel     = c(names(vars_desc), "Jogos"),
  em_casa      = c(colMeans(vars_desc[em_casa_am, ]), sum(em_casa_am)),
  neutro       = c(colMeans(vars_desc[!em_casa_am, ]), sum(!em_casa_am)),
  desvio_total = c(sapply(vars_desc, sd), NA))
tab1$diferenca <- c(tab1$em_casa[-nrow(tab1)] - tab1$neutro[-nrow(tab1)], NA)
# Diferença padronizada: diferença / desvio combinado sqrt((s2_casa + s2_neutro) / 2)
desvio_comb <- sqrt((sapply(vars_desc[em_casa_am, ], var) + sapply(vars_desc[!em_casa_am, ], var)) / 2)
tab1$dif_padronizada <- tab1$diferenca / c(desvio_comb, NA)
tab1[, -1] <- round(tab1[, -1], 3)
write.csv(tab1, file.path(pasta_desc, "tab1_casa_neutro.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 6.2 Tabela 2: por tipo de torneio ---------------------------------------------
# Quanto cada tipo pesa na amostra, quanto é jogado em campo neutro e a
# diferença bruta de saldo entre casa e neutro. Grupo vazio fica NA ("-" no paper).
media_ou_na <- function(x) if (length(x) > 0) round(mean(x), 3) else NA
tab2 <- do.call(rbind, lapply(split(am, am$torneio), function(d) data.frame(
  jogos         = nrow(d),
  jogos_casa    = sum(d$casa == 1),
  jogos_neutros = sum(d$casa == 0),
  pct_neutros   = round(100 * mean(d$casa == 0), 1),
  saldo_casa    = media_ou_na(d$saldo[d$casa == 1]),
  saldo_neutro  = media_ou_na(d$saldo[d$casa == 0]))))
tab2 <- cbind(torneio = rownames(tab2), tab2)
tab2 <- tab2[order(-tab2$jogos), ]
write.csv(tab2, file.path(pasta_desc, "tab2_tipo_torneio.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 6.3 Figura 1: saldo médio por década, em casa e neutro -------------------------
# Só décadas com jogos neutros suficientes para uma média estável.
min_neutros <- 50
decadas <- names(which(table(am$decada_jogo[!em_casa_am]) >= min_neutros))
media_dec <- function(sel) tapply(am$saldo[sel], am$decada_jogo[sel], mean)[decadas]
saldo_casa_dec   <- media_dec(em_casa_am)
saldo_neutro_dec <- media_dec(!em_casa_am)

png(file.path(pasta_desc, "fig1_saldo_por_decada.png"), width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 1, 1))                         # sem espaço para título interno
x <- as.integer(decadas)
plot(x, saldo_casa_dec, type = "b", pch = 19, yaxt = "n",
     ylim = c(min(0, saldo_casa_dec, saldo_neutro_dec),          # não corta médias negativas
              1.3 * max(saldo_casa_dec, saldo_neutro_dec)),      # folga para a legenda
     xlab = "Década", ylab = "Saldo médio do mandante (gols)")
eixo_virgula(2)
lines(x, saldo_neutro_dec, type = "b", pch = 1, lty = 2)
abline(h = 0, col = "grey60")
legend("top", c("Em casa", "Campo neutro"), pch = c(19, 1), lty = c(1, 2), bty = "n", horiz = TRUE)
dev.off()

# 6.4 Figura 2: diferença de Elo em casa e neutro --------------------------------
# Mostra se o mandante em casa já é, em média, mais forte (por que controlar o Elo).
png(file.path(pasta_desc, "fig2_dif_elo_casa_neutro.png"), width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 1, 1))
dens_casa   <- density(am$dif_elo[em_casa_am])
dens_neutro <- density(am$dif_elo[!em_casa_am])
plot(dens_casa, main = "", xlab = "Diferença de Elo (mandante - visitante, centenas)",
     ylab = "Densidade", ylim = range(0, dens_casa$y, dens_neutro$y), yaxt = "n")
eixo_virgula(2)
lines(dens_neutro, lty = 2)
abline(v = 0, col = "grey60")
legend("topright", c("Em casa", "Campo neutro"), lty = c(1, 2), bty = "n")
dev.off()

# 6.5 Resumo --------------------------------------------------------------------
print(tab1, row.names = FALSE)
print(tab2, row.names = FALSE)
cat("Arquivos salvos em", pasta_desc, ":", paste(list.files(pasta_desc), collapse = ", "), "\n")


# -----------------------------------------------------------------------------
# 7. Estimação (MQO e Probit)
# -----------------------------------------------------------------------------
# Saldo de gols por MQO e vitória por probit (principal), mais as especificações
# alternativas. Erros-padrão com agrupamento duplo (mandante e visitante), com o
# ajuste de pequenas amostras padrão do sandwich. Nos modelos não lineares, os
# efeitos são efeitos marginais médios (EMM), com erro-padrão pelo método delta.
pasta_est <- "resultados/estimacao"
dir.create(pasta_est, recursive = TRUE, showWarnings = FALSE)

# 7.1 Modelos principais e alternativos -----------------------------------------
# Brutos (só casa) e completos; alternativas para a vitória: MPL e logit.
controles_mod <- c("anf_copa", "anf_continental", "anf_regional", "dif_elo", "dif_forma_pts",
                   "dif_forma_sg", "dif_experiencia", "dif_descanso", "torneio", "decada_jogo")
completo <- function(y) reformulate(c("casa", controles_mod), y)
modelos <- list(
  mqo_bruto    = lm(saldo ~ casa, data = am),
  mqo          = lm(completo("saldo"), data = am),
  probit_bruto = glm(vitoria ~ casa, family = binomial("probit"), data = am),
  probit       = glm(completo("vitoria"), family = binomial("probit"), data = am),
  mpl          = lm(completo("vitoria"), data = am),
  logit        = glm(completo("vitoria"), family = binomial("logit"), data = am))
vcovs <- lapply(modelos, vcovCL, cluster = ~ home_team + away_team)

# Efeitos fixos de par (só nos modelos lineares): par dirigido, mandante ->
# visitante, porque o saldo é sempre do ponto de vista do mandante. Estimado
# pela transformação "dentro do par" (tira a média de cada par). Cada par
# dirigido fica dentro de um cluster de mandante, então os efeitos fixos não
# entram no ajuste de graus de liberdade (mesma convenção de fixest e Stata).
par_dirigido <- paste(am$home_team, "->", am$away_team)
dentro_par <- function(y) {
  X  <- model.matrix(completo(y), am)[, -1]
  Xd <- apply(X, 2, function(v) v - ave(v, par_dirigido))
  Xd <- Xd[, colSums(Xd^2) > 1e-8]                        # tira o que não varia dentro do par
  yd <- am[[y]] - ave(am[[y]], par_dirigido)
  m  <- lm(yd ~ Xd - 1)
  names(m$coefficients) <- sub("^Xd", "", names(m$coefficients))
  list(m = m, V = vcovCL(m, cluster = am[c("home_team", "away_team")]))
}
fe_par <- list(mqo_par = dentro_par("saldo"), mpl_par = dentro_par("vitoria"))

# 7.2 Efeito de casa e vantagem total do anfitrião -----------------------------
# Casa: casa de 0 para 1 em todos os jogos, com as dummies de anfitrião em zero
# (efeito fora de torneio com sede, como beta1 no MQO). Anfitrião: casa e a
# dummy da sede mudam juntas, só nos jogos daquela sede (beta1 + beta_sede).
efeito <- function(k, dados, vars) {
  r <- avg_comparisons(modelos[[k]], variables = setNames(rep(list(0:1), length(vars)), vars),
                       cross = length(vars) > 1, newdata = dados, vcov = vcovs[[k]])
  c(estimativa = r$estimate, ep = r$std.error, p = r$p.value)
}
soma_coef <- function(m, V, nomes) {                     # soma de coeficientes (modelos lineares)
  w <- as.numeric(names(coef(m)) %in% nomes)
  est <- sum(w * coef(m)); ep <- sqrt(drop(w %*% V %*% w))
  c(estimativa = est, ep = ep, p = 2 * pnorm(-abs(est / ep)))
}
sedes   <- c(copa = "anf_copa", continental = "anf_continental", regional = "anf_regional")
sem_anf <- transform(am, anf_copa = 0, anf_continental = 0, anf_regional = 0)
efeitos <- list(mqo_bruto_casa    = efeito("mqo_bruto", am, "casa"),
                probit_bruto_casa = efeito("probit_bruto", am, "casa"))
for (k in c("mqo", "probit", "mpl", "logit")) {
  efeitos[[paste0(k, "_casa")]] <- efeito(k, sem_anf, "casa")
  for (s in names(sedes))
    efeitos[[paste(k, s, sep = "_")]] <- efeito(k, am[am[[sedes[s]]] == 1, ], c("casa", sedes[s]))
}
for (k in names(fe_par)) {
  efeitos[[paste0(k, "_casa")]] <- soma_coef(fe_par[[k]]$m, fe_par[[k]]$V, "casa")
  for (s in names(sedes))
    efeitos[[paste(k, s, sep = "_")]] <- soma_coef(fe_par[[k]]$m, fe_par[[k]]$V, c("casa", sedes[s]))
}

# 7.3 Tabelas ---------------------------------------------------------------------
# Painel A: coeficientes (no probit e no logit, do índice: não são probabilidades).
# Painel B: efeitos de casa e do anfitrião (MQO do saldo: gols; demais: probabilidade).
estrelas <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
fmt <- function(e) sprintf("%.3f%s (%.3f)", e[1], estrelas(e[3]), e[2])
vars_tab <- c("casa", "anf_copa", "anf_continental", "anf_regional", "dif_elo", "dif_forma_pts",
              "dif_forma_sg", "dif_experiencia", "dif_descanso")
painel_b <- function(k) sapply(c("casa", names(sedes)), function(s) {
  e <- efeitos[[paste(k, s, sep = "_")]]; if (is.null(e)) "" else fmt(e) })
r2 <- function(m) if (inherits(m, "glm")) 1 - m$deviance / m$null.deviance else summary(m)$r.squared
rotulos_b <- c("Efeito de casa", paste("Vantagem total do anfitrião:", names(sedes)))

# Tabela 3: principal
coluna3 <- function(k) {
  m <- modelos[[k]]; ct <- coeftest(m, vcov. = vcovs[[k]])
  c("", sapply(vars_tab, function(v) if (v %in% rownames(ct)) fmt(ct[v, c(1, 2, 4)]) else ""),
    "", painel_b(k), if (grepl("bruto", k)) "Não" else "Sim", nobs(m), sprintf("%.3f", r2(m)))
}
tab3 <- data.frame(
  variavel = c("Painel A: coeficientes", vars_tab, "Painel B: efeitos", rotulos_b,
               "Efeitos fixos de torneio e década", "Jogos", "R² / pseudo-R²"),
  sapply(c("mqo_bruto", "mqo", "probit_bruto", "probit"), coluna3), row.names = NULL)

# Tabela 4: especificações alternativas (só os efeitos). Com efeito fixo de par,
# só os pares com jogos em casa e neutros identificam o efeito de casa.
par_varia <- ave(am$casa, par_dirigido, FUN = function(x) length(unique(x)) == 2) == 1
coluna4 <- function(k) {
  fe <- k %in% names(fe_par)
  m  <- if (fe) fe_par[[k]]$m else modelos[[k]]
  c(painel_b(k), if (fe) "Sim" else "Não", nobs(m),
    if (fe) c(length(unique(par_dirigido[par_varia])), sum(par_varia)) else c("–", "–"))
}
tab4 <- data.frame(variavel = c(rotulos_b, "Efeitos fixos de par", "Jogos",
                                "Pares que identificam o efeito", "Jogos nesses pares"),
                   sapply(c("mqo", "mqo_par", "probit", "mpl", "logit", "mpl_par"), coluna4),
                   row.names = NULL)

write.csv(tab3, file.path(pasta_est, "tab3_principal.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab4, file.path(pasta_est, "tab4_alternativas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(modelos = modelos, vcovs = vcovs, fe_par = fe_par, efeitos = efeitos),
        file.path(pasta_est, "modelos.rds"))             # usados nas seções 8 a 10

# 7.4 Conferências e resumo -----------------------------------------------------
stopifnot(all(sapply(modelos, nobs) == nrow(am)),
          !anyNA(unlist(lapply(modelos[-c(1, 3)], coef))), modelos$probit$converged,
          all(sapply(c(vcovs, lapply(fe_par, `[[`, "V")), function(V) all(diag(V) > 0))),  # sem variância negativa
          abs(efeitos$mqo_copa[1] - sum(coef(modelos$mqo)[c("casa", "anf_copa")])) < 1e-8)
print(tab3, row.names = FALSE)
print(tab4, row.names = FALSE)


# -----------------------------------------------------------------------------
# 8. Diagnóstico e testes
# -----------------------------------------------------------------------------
# Só testes que respondem a uma hipótese prévia ou mudam uma decisão. Os testes
# de casa e de cada anfitrião já estão na Tabela 3 (t ou z com erro agrupado).
pasta_diag <- "resultados/diagnostico"
dir.create(pasta_diag, recursive = TRUE, showWarnings = FALSE)
mqo <- modelos$mqo; prb <- modelos$probit

# Teste de Wald (H0: coeficientes = 0), qui-quadrado com a matriz agrupada
wald <- function(m, V, nomes) {
  b <- coef(m)[nomes]; W <- drop(t(b) %*% solve(V[nomes, nomes]) %*% b)
  data.frame(estatistica = round(W, 2), gl = length(nomes),
             p = signif(pchisq(W, length(nomes), lower.tail = FALSE), 3))
}

# 8.1 Testes ----------------------------------------------------------------------
# Forma: as duas medidas andam juntas, então o teste relevante é o conjunto.
# RESET: potências do valor ajustado indicam se falta não linearidade no MQO.
forma <- c("dif_forma_pts", "dif_forma_sg")
am$aj2 <- fitted(mqo)^2; am$aj3 <- fitted(mqo)^3
reset <- lm(update(formula(mqo), . ~ . + aj2 + aj3), data = am)
testes <- cbind(
  teste  = c("Forma (pontos e saldo) = 0", "Forma (pontos e saldo) = 0", "RESET (ajustado² e ³ = 0)"),
  modelo = c("MQO", "Probit", "MQO"),
  rbind(wald(mqo, vcovs$mqo, forma), wald(prb, vcovs$probit, forma),
        wald(reset, vcovCL(reset, cluster = ~ home_team + away_team), c("aj2", "aj3"))))
am$aj2 <- am$aj3 <- NULL                                  # colunas auxiliares do RESET

# 8.2 Sinais esperados ------------------------------------------------------------
esperado <- c(casa = "+", anf_copa = "+", anf_continental = "+", anf_regional = "≈ 0",
              dif_elo = "+", dif_forma_pts = "+", dif_forma_sg = "+", dif_experiencia = "+",
              dif_descanso = "ambíguo")
sinal <- function(m, V) {
  ct <- coeftest(m, vcov. = V)[names(esperado), ]
  ifelse(ct[, 4] >= 0.05, "n.s.", ifelse(ct[, 1] > 0, "+", "−"))
}
sinais <- data.frame(variavel = names(esperado), esperado,
                     mqo = sinal(mqo, vcovs$mqo), probit = sinal(prb, vcovs$probit), row.names = NULL)

# 8.3 Multicolinearidade (VIF) e previsões do MPL ------------------------------
# VIF = 1 / (1 - R²) de cada controle contínuo contra as demais colunas do
# modelo (incluindo torneio e década); acima de 10 indica problema.
X <- model.matrix(mqo)[, -1]
vif <- sapply(names(esperado)[5:9], function(v) {
  res <- lm.fit(cbind(1, X[, colnames(X) != v]), X[, v])$residuals
  sum((X[, v] - mean(X[, v]))^2) / sum(res^2)
})
fora_01 <- mean(fitted(modelos$mpl) < 0 | fitted(modelos$mpl) > 1)

# Ajuste do probit: pseudo-R² de McFadden e fração corretamente prevista (corte
# 0,5), no total e por resultado, contra a regra ingênua (sempre "não vitória").
pr <- fitted(prb); y <- am$vitoria
ajuste_probit <- round(c(pseudo_r2 = 1 - prb$deviance / prb$null.deviance,
                         acerto_total = mean((pr > 0.5) == y), acerto_vitorias = mean(pr[y == 1] > 0.5),
                         acerto_nao_vitorias = mean(pr[y == 0] <= 0.5), regra_ingenua = mean(y == 0)), 3)

# 8.4 Salvar e resumo -------------------------------------------------------------
write.csv(testes, file.path(pasta_diag, "testes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sinais, file.path(pasta_diag, "sinais.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(testes, row.names = FALSE)
print(sinais, row.names = FALSE)
cat("VIF:", paste(names(vif), round(vif, 2), sep = " = ", collapse = "; "), "\n")
cat("MPL: previsões fora de [0, 1]:", round(100 * fora_01, 1), "% dos jogos\n")


# -----------------------------------------------------------------------------
# 9. Robustez
# -----------------------------------------------------------------------------
# Muda uma coisa por vez e guarda só o efeito de casa: beta1 (saldo, MQO) e
# EMM (vitória, probit).
pasta_rob <- "resultados/robustez"
dir.create(pasta_rob, recursive = TRUE, showWarnings = FALSE)

# 9.1 Função de reestimação ------------------------------------------------------
# troca: antigo = novo ("" tira); extra: termos a mais; dados; cluster;
# so_saldo: não reestima o probit (quando a vitória não muda).
reestimar <- function(troca = NULL, extra = NULL, dados = am, cluster = ~ home_team + away_team,
                      so_saldo = FALSE) {
  ctrl <- controles_mod
  for (v in names(troca)) ctrl <- if (troca[[v]] == "") setdiff(ctrl, v) else replace(ctrl, ctrl == v, troca[[v]])
  dados <- droplevels(dados)
  f <- function(y) reformulate(c("casa", ctrl, extra), y)
  m_saldo <- lm(f("saldo"), data = dados)
  V  <- vcovCL(m_saldo, cluster = cluster)
  b  <- c(coef(m_saldo)["casa"], sqrt(V["casa", "casa"]))
  b  <- c(b, 2 * pt(-abs(b[1] / b[2]), df.residual(m_saldo)))
  vit <- "= principal"
  if (!so_saldo) {
    m_vit <- glm(f("vitoria"), family = binomial("probit"), data = dados)
    e <- avg_comparisons(m_vit, variables = list(casa = 0:1), vcov = vcovCL(m_vit, cluster = cluster),
                         newdata = transform(dados, anf_copa = 0, anf_continental = 0, anf_regional = 0))
    vit <- fmt(c(e$estimate, e$std.error, e$p.value))
  }
  c(saldo = fmt(b), vitoria = vit, jogos = nobs(m_saldo))
}

# 9.2 Exercícios ------------------------------------------------------------------
ano_am <- as.integer(format(am$data, "%Y"))
sem_forma <- c(dif_forma_pts = "", dif_forma_sg = "")
exercicios <- list(
  "Especificação principal"            = list(),
  # medida e forma funcional
  "Elo com bônus de 0"                 = list(troca = c(dif_elo = "dif_elo_b0")),
  "Elo com bônus de 150"               = list(troca = c(dif_elo = "dif_elo_b150")),
  "Elo com termo quadrático"           = list(extra = "I(dif_elo^2)"),
  "Forma: janela de 5 jogos"           = list(troca = c(dif_forma_pts = "dif_forma_pts_5", dif_forma_sg = "dif_forma_sg_5")),
  "Forma: janela de 20 jogos"          = list(troca = c(dif_forma_pts = "dif_forma_pts_20", dif_forma_sg = "dif_forma_sg_20")),
  "Sem forma"                          = list(troca = sem_forma),
  # amostra
  "Mínimo de 20 jogos anteriores"      = list(dados = am[am$n_ant_m >= 20 & am$n_ant_v >= 20, ]),
  "Mínimo de 30 jogos anteriores"      = list(dados = am[am$n_ant_m >= 30 & am$n_ant_v >= 30, ]),
  "A partir de 1950"                   = list(dados = am[ano_am >= 1950, ]),
  "A partir de 1990"                   = list(dados = am[ano_am >= 1990, ]),
  "Saldo limitado a ±5"                = list(dados = transform(am, saldo = pmax(pmin(saldo, 5), -5)),
                                              so_saldo = TRUE),
  # erros-padrão
  "Agrupamento por par (confronto)"    = list(cluster = ~ confronto))

tab5 <- do.call(rbind, lapply(exercicios, function(a) do.call(reestimar, a)))
tab5 <- data.frame(exercicio = names(exercicios), tab5, row.names = NULL)

# 9.3 Salvar e resumo -------------------------------------------------------------
stopifnot(tab5$saldo[1]   == fmt(coeftest(modelos$mqo, vcov. = vcovs$mqo)["casa", c(1, 2, 4)]),
          tab5$vitoria[1] == fmt(efeitos$probit_casa))     # reproduz a Tabela 3
write.csv(tab5, file.path(pasta_rob, "tab5_robustez.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(tab5, row.names = FALSE)


# -----------------------------------------------------------------------------
# 10. Extensões
# -----------------------------------------------------------------------------
# Mesmas regras da análise principal (controles da Tabela 3, agrupamento duplo,
# EMM no probit). Tabela 6: altitude e período; Tabela 7: pênaltis.
pasta_ext <- "resultados/extensoes"
dir.create(pasta_ext, recursive = TRUE, showWarnings = FALSE)

# 10.1 Funções --------------------------------------------------------------------
# p-valores e IC pela t com G - 1 graus de liberdade (G = clusters na menor
# dimensão), o que importa nas amostras pequenas das disputas. hc1 = TRUE usa o
# erro-padrão robusto simples, como verificação.
ajustar <- function(x, dados, y, probit = FALSE) {         # modelo + matriz agrupada
  m <- if (probit) glm(reformulate(x, y), family = binomial("probit"), data = dados)
  else lm(reformulate(x, y), data = dados)
  list(m = m, V = vcovCL(m, cluster = ~ home_team + away_team), Vhc = vcovHC(m, type = "HC1"),
       gl = min(lengths(lapply(dados[c("home_team", "away_team")], unique))) - 1, dados = dados)
}
comb <- function(a, w, hc1 = FALSE) {                      # combinação linear (est, ep, p, gl)
  x <- setNames(rep(0, length(coef(a$m))), names(coef(a$m))); x[names(w)] <- w
  V <- if (hc1) a$Vhc else a$V
  est <- sum(x * coef(a$m)); ep <- sqrt(drop(x %*% V %*% x))
  c(est, ep, 2 * pt(-abs(est / ep), a$gl), a$gl)
}
emm <- function(a, v, dados = a$dados, hc1 = FALSE) {      # EMM: dummy 0 -> 1 ou +1 km
  vars <- if (v %in% c("dif_subida", "dif_descida")) v else setNames(list(0:1), v)
  r <- avg_comparisons(a$m, variables = vars, newdata = dados, vcov = if (hc1) a$Vhc else a$V)
  c(r$estimate, r$std.error, 2 * pt(-abs(r$estimate / r$std.error), a$gl), a$gl)
}
p_ind  <- function(e) sprintf("p = %.3f (índice)", e[3])   # teste no índice do probit
sem_anf_d <- function(d) transform(d, anf_copa = 0, anf_continental = 0, anf_regional = 0)
rhs <- c("casa", controles_mod)

# 10.2 Altitude -------------------------------------------------------------------
# Com e sem altitude na mesma amostra; robustez sem os casos aproximados.
alt <- function(dados) {
  dados <- droplevels(dados); x1 <- c(rhs, "dif_subida", "dif_descida")
  s0 <- ajustar(rhs, dados, "saldo"); p0 <- ajustar(rhs, dados, "vitoria", TRUE)
  s1 <- ajustar(x1, dados, "saldo");  p1 <- ajustar(x1, dados, "vitoria", TRUE)
  rbind("Casa, sem altitude" = c(fmt(comb(s0, c(casa = 1))), fmt(emm(p0, "casa", sem_anf_d(dados)))),
        "Casa, com altitude" = c(fmt(comb(s1, c(casa = 1))), fmt(emm(p1, "casa", sem_anf_d(dados)))),
        "Subida (km)"        = c(fmt(comb(s1, c(dif_subida = 1))), fmt(emm(p1, "dif_subida"))),
        "Descida (km)"       = c(fmt(comb(s1, c(dif_descida = 1))), fmt(emm(p1, "dif_descida"))),
        "Simetria (subida + descida = 0)" = c(fmt(comb(s1, c(dif_subida = 1, dif_descida = 1))),
                                              p_ind(comb(p1, c(dif_subida = 1, dif_descida = 1)))),
        "Jogos" = nobs(s1$m))
}
aprox <- am$cidade_sem_alt | am$alt_hab_aprox_m | am$alt_hab_aprox_v
tab_alt2 <- alt(am[am$am_altitude & !aprox, ])[-c(1, 5), ]
rownames(tab_alt2) <- paste(rownames(tab_alt2), "(sem aproximados)")

# 10.3 Evolução no tempo ----------------------------------------------------------
# Casa dividido em uma dummy por período; as dummies de década continuam.
am_per <- am; cp <- paste0("casa_p", 1:4)
for (k in 1:4) am_per[[cp[k]]] <- am$casa * (as.integer(am$periodo) == k)
sp <- ajustar(c(cp, controles_mod), am_per, "saldo")
pp <- ajustar(c(cp, controles_mod), am_per, "vitoria", TRUE)
iguais <- function(a) {                                    # Wald: os quatro períodos iguais
  R <- sapply(2:4, function(k) as.numeric(names(coef(a$m)) == cp[1]) - (names(coef(a$m)) == cp[k]))
  d <- crossprod(R, coef(a$m)); W <- drop(t(d) %*% solve(t(R) %*% a$V %*% R) %*% d)
  sprintf("χ²(3) = %.2f; p = %.3f", W, pchisq(W, 3, lower.tail = FALSE))
}
tab_per <- rbind(
  t(sapply(1:4, function(k) c(fmt(comb(sp, setNames(1, cp[k]))),
                              fmt(emm(pp, cp[k], sem_anf_d(am_per[as.integer(am_per$periodo) == k, ])))))),
  c(fmt(comb(sp, c(casa_p1 = 1, casa_p4 = -1))), p_ind(comb(pp, c(casa_p1 = 1, casa_p4 = -1)))),
  c(iguais(sp), paste(iguais(pp), "(índice)")))
rownames(tab_per) <- c(paste("Casa,", levels(am$periodo)), "Primeiro período − último", "Quatro períodos iguais")

tab6 <- rbind(alt(am[am$am_altitude, ]), tab_alt2, tab_per)
tab6 <- data.frame(linha = rownames(tab6), saldo_mqo = tab6[, 1], vitoria_probit = tab6[, 2], row.names = NULL)

# 10.4 Pênaltis a favor e disputas ---------------------------------------------------
# Pênaltis a favor: referência do torneio = eliminatórias da Copa (não há
# amistosos com todos os gols registrados). Disputas: poucas explicativas;
# cobrar primeiro só onde a ordem foi sorteada (até jul/2003).
d_pen <- droplevels(am[am$am_penalti, ]); d_pen$torneio <- relevel(d_pen$torneio, "Eliminatórias Copa")
d_disp <- am[am$am_disputa, ]
d_sort <- d_disp[!is.na(d_disp$cobrou_primeiro) & d_disp$ordem_sorteada, ]
m7 <- list(pen  = ajustar(rhs, d_pen, "dif_penalti"),
           mpl  = ajustar(c("casa", "dif_elo"), d_disp, "venceu_disputa"),
           prb  = ajustar(c("casa", "dif_elo"), d_disp, "venceu_disputa", TRUE),
           sort = ajustar(c("cobrou_primeiro", "casa", "dif_elo"), d_sort, "venceu_disputa"))
efeito7 <- function(hc1 = FALSE) list(
  "Pênaltis a favor: casa (MQO)"             = comb(m7$pen, c(casa = 1), hc1),
  "Disputa: casa (MPL)"                      = comb(m7$mpl, c(casa = 1), hc1),
  "Disputa: casa (probit, EMM)"              = emm(m7$prb, "casa", hc1 = hc1),
  "Disputa: cobrou primeiro (MPL, sorteada)" = comb(m7$sort, c(cobrou_primeiro = 1), hc1))
e7 <- efeito7(); e7_hc1 <- efeito7(hc1 = TRUE)
tab7 <- data.frame(linha = names(e7), efeito = sapply(e7, fmt),
                   ic95 = sapply(e7, function(e) { h <- qt(0.975, e[4]) * e[2]
                   sprintf("[%.3f; %.3f]", e[1] - h, e[1] + h) }),
                   ep_hc1 = sprintf("%.3f", sapply(e7_hc1, `[`, 2)),
                   gl_t = sapply(e7, `[`, 4), jogos = sapply(m7, function(a) nobs(a$m)), row.names = NULL)

# 10.5 Salvar e resumo ------------------------------------------------------------
write.csv(tab6, file.path(pasta_ext, "tab6_altitude_periodo.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab7, file.path(pasta_ext, "tab7_penaltis.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(tab6, row.names = FALSE)
print(tab7, row.names = FALSE)


# -----------------------------------------------------------------------------
# 11. Encerramento
# -----------------------------------------------------------------------------
# Versões do R e dos pacotes (cálculos podem mudar entre versões). Sem
# set.seed(): nada no script é aleatório.
sessionInfo()
writeLines(capture.output(sessionInfo()), "resultados/sessao_R.txt")
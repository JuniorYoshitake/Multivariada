# códigos

library(knitr)
library(kableExtra)
library(dplyr)
library(psych)

# QUESTÃO 66 (4.3 Artes e barroso

# verificação da fatorabilidade

pizza <- read.csv("Pizza.csv") %>% select(-brand, -id)

pizza_ord <- pizza %>%
    mutate(across(
        everything(),
        \(x) {
            cortes <- quantile(
                x,
                probs = seq(0, 1, length.out = 12),
                na.rm = TRUE
            )
            cortes <- unique(cortes)
            
            cut(
                x,
                breaks = cortes,
                include.lowest = TRUE,
                labels = FALSE
            ) - 1
        }
    ))

matriz_poly <- polychoric(pizza_ord, max.cat = 11)$rho

kable(
    round(matriz_poly, 4),
    format = "latex",
    booktabs = TRUE) %>%
    kable_styling(latex_options = c("hold_position"))

# definição do número de fatores

autovalores <- eigen(matriz_poly)$values

plot(autovalores, type = "b", pch = 19, col = "blue",
     main = "Scree Plot - Nutrientes de Pizza",
     xlab = "Número do Componente",
     ylab = "Autovalor (Variância)")

# Extração de Fatores

fit <- principal(matriz_poly, nfactors = 3, rotate = "varimax")

# cargas fatoriais
cargas <- as.data.frame(unclass(fit$loadings))

# comunalidades
cargas$h2 <- fit$communality

# unicidades
cargas$ei <- fit$uniquenesses

# soma dos quadrados das cargas
ss <- data.frame(
    RC1 = sum(cargas$RC1^2),
    RC2 = sum(cargas$RC2^2),
    RC3 = sum(cargas$RC3^2),
    h2 = NA,
    ei = NA
)

cargas <- cargas %>% mutate(h2 = round(h2, 4),
                            ei = round(ei, 4))

# adiciona linha final
tabela_final <- bind_rows(cargas, SS = ss)

tabela_final$h2[is.na(tabela_final$h2)] <- ""
tabela_final$ei[is.na(tabela_final$ei)] <- ""

rownames(tabela_final)[8] <- "SS"

# tabela latex
kable(
    tabela_final,
    format = "latex",
    booktabs = TRUE,
    align = "c",
    escape = FALSE,
    digits = 2,
    col.names = c("F1", "F2", "F3", "$h^2$", "$\\psi$")
) %>%
    kable_styling(
        latex_options = c("hold_position", "striped")
    )

p <- nrow(cargas)

var_f1 <- ss$RC1 / p
var_f2 <- ss$RC2 / p
var_f3 <- ss$RC3 / p

var_total <- var_f1 + var_f2 + var_f3

# Interpretação dos fatores e Análise das comunalidades 
#e Avaliação do ajuste do modelo

Psi <- diag(fit$uniquenesses)

Lambda <- as.matrix(unclass(fit$loadings))

Sigma_hat <- Lambda %*% t(Lambda) + Psi

Residual <- matriz_poly - Sigma_hat

kable(
    round(Residual, 4),
    format = "latex",
    booktabs = TRUE) %>%
    kable_styling(latex_options = c("hold_position"))

## Questão 68 (4.8 Arte e Barroso) #####################################

library(readxl)

dados <- read_xlsx(path = "Stress\\Stress.xlsx") %>% select(-Código)

# Tratamento da esparsidade dos dados

dados$X1 <- ifelse(dados$X1 == 6, 5, dados$X1)

dados$X5[dados$X5 == 1] <- 0

dados$X5 <- as.numeric(as.factor(dados$X5)) - 1

dados$X6 <- ifelse(dados$X6 == 0, 1, dados$X1)

dados$X7[dados$X7 == 0] <- 2

dados$X7 <- as.numeric(as.factor(dados$X7)) - 1

dados$X8 <- ifelse(dados$X8 == 7, 6, dados$X8)

dados$X8 <- ifelse(dados$X8 == 6, 5, dados$X8)

dados$X10 <- ifelse(dados$X10 == 6, 5, dados$X10)

dados$X11 <- ifelse(dados$X11 == 6, 5, dados$X11)

dados$X12 <- ifelse(dados$X12 == 6, 5, dados$X12)

dados$X13[dados$X13 == 1] <- 0

dados$X13 <- as.numeric(as.factor(dados$X13)) - 1

matriz_poly <- polychoric(dados)$rho

kable(
    round(matriz_poly[, 1:7], 4),
    format = "latex",
    booktabs = TRUE
) %>%
    kable_styling(
        latex_options = c("HOLD_position")
    )

kable(
    round(matriz_poly[, 8:14], 4),
    format = "latex",
    booktabs = TRUE
) %>%
    kable_styling(
        latex_options = c("HOLD_position")
    )

# AF com m=2

fit <- principal(matriz_poly, nfactors = 2, rotate = "varimax")

# cargas fatoriais
cargas <- as.data.frame(unclass(fit$loadings))

# comunalidades
cargas$h2 <- fit$communality

# unicidades
cargas$ei <- fit$uniquenesses

ss <- data.frame(
    RC1 = sum(cargas$RC1^2),
    RC2 = sum(cargas$RC2^2),
    h2 = NA,
    ei = NA
)

# arredondamento
cargas <- cargas %>%
    mutate(
        h2 = round(h2, 4),
        ei = round(ei, 4)
    )

# adiciona linha SS
tabela_final <- bind_rows(cargas, SS = ss)

# substitui NAs por vazio
tabela_final$h2[is.na(tabela_final$h2)] <- ""
tabela_final$ei[is.na(tabela_final$ei)] <- ""

# nome da última linha
rownames(tabela_final)[nrow(tabela_final)] <- "SS"

# evita erro com "_"
rownames(tabela_final) <- gsub(
    "_",
    "\\\\_",
    rownames(tabela_final)
)

# variância explicada
p <- nrow(cargas)

var_f1 <- ss$RC1 / p
var_f2 <- ss$RC2 / p

var_total <- var_f1 + var_f2

# tabela LaTeX
kable(
    tabela_final,
    format = "latex",
    booktabs = TRUE,
    align = "c",
    escape = FALSE,
    digits = 2,
    col.names = c("F1", "F2", "$h^2$", "$\\psi$")
) %>%
    kable_styling(
        latex_options = c("HOLD_position", "striped")
    )

# AF com m = 6

fit <- principal(matriz_poly, nfactors = 6, rotate = "varimax")

# cargas fatoriais
cargas <- as.data.frame(unclass(fit$loadings))

# comunalidades
cargas$h2 <- fit$communality

# unicidades
cargas$ei <- fit$uniquenesses

ss <- data.frame(
    RC1 = sum(cargas$RC1^2),
    RC2 = sum(cargas$RC2^2),
    RC3 = sum(cargas$RC3^2),
    RC4 = sum(cargas$RC4^2),
    RC5 = sum(cargas$RC5^2),
    RC6 = sum(cargas$RC6^2),
    h2 = NA,
    ei = NA
)

# arredondamento
cargas <- cargas %>%
    mutate(
        h2 = round(h2, 4),
        ei = round(ei, 4)
    )

# adiciona linha SS
tabela_final <- bind_rows(cargas, SS = ss)

# substitui NAs por vazio
tabela_final$h2[is.na(tabela_final$h2)] <- ""
tabela_final$ei[is.na(tabela_final$ei)] <- ""

# nome da última linha
rownames(tabela_final)[nrow(tabela_final)] <- "SS"

# evita erro com "_"
rownames(tabela_final) <- gsub(
    "_",
    "\\\\_",
    rownames(tabela_final)
)

# variância explicada
p <- nrow(cargas)

var_f1 <- ss$RC1 / p
var_f2 <- ss$RC2 / p
var_f3 <- ss$RC3 / p
var_f4 <- ss$RC4 / p
var_f5 <- ss$RC5 / p
var_f6 <- ss$RC6 / p

var_total <- var_f1 + var_f2 + var_f3 + var_f4 + var_f5 + var_f6

# tabela LaTeX
kable(
    tabela_final,
    format = "latex",
    booktabs = TRUE,
    align = "c",
    escape = FALSE,
    digits = 2,
    col.names = c("F1", "F2", "F3", "F4","F5", "F6", "$h^2$", "$\\psi$")
) %>%
    kable_styling(
        latex_options = c("HOLD_position", "striped")
    )

# Scree-plot

autovalores <- eigen(matriz_poly)$values

plot(autovalores, type = "b", pch = 19, col = "blue",
     main = "Scree Plot - BemEstarFin",
     xlab = "Número do Componente",
     ylab = "Autovalor (Variância)")

# Questão 81 (6.8 Artes e Barroso) #####################################

tab645 <- data.frame(
    Moradores = c("1", "2", "3", "4", "5", "6 ou mais", "Total"),
    Norte = c(627, 1105, 1269, 1119, 605, 561, 5286),
    Nordeste = c(2592, 4475, 5015, 3699, 1681, 1054, 18516),
    Sudeste = c(4918, 8301, 8228, 6067, 2144, 1082, 30740),
    Sul = c(1680, 3089, 2975, 1933, 686, 305, 10668),
    Centro_Oeste = c(772, 1447, 1409, 1120, 463, 223, 5434),
    Total = c(10589, 18417, 18896, 13938, 5579, 3225, 70644)
)

tab645 |>
    kable(
        format = "latex",
        booktabs = TRUE,
        align = c("l", rep("c", 6)),
        format.args = list(big.mark = ".")
    ) |>
    kable_styling(
        latex_options = "HOLD_position"
    )

tab646 <- data.frame(
    Regiao = c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste", "Total"),
    Branca = c(3279, 14190, 45438, 22813, 5746, 91466),
    Preta = c(1388, 6572, 8580, 1321, 1437, 19298),
    Parda = c(13640, 36323, 34642, 6041, 9256, 99902),
    Total = c(18307, 57085, 88660, 30175, 16439, 210666)
)

tab646 |>
    kable(
        format = "latex",
        booktabs = TRUE,
        align = c("l", rep("c", 4)),
        format.args = list(big.mark = ".")
    ) |>
    kable_styling(
        latex_options = "HOLD_position"
    )

tab1 <- t(as.matrix(tab645[1:6, 2:6]))
tab2 <- as.matrix(tab646[1:5, 2:4])

tab_just <- cbind(tab1, tab2)

res_just <- ca(tab_just)

plot(res_just)

tab645 <- data.frame(
    Moradores = c("1", "2", "3", "4", "5", "6 ou mais", "Total"),
    Norte = c(627, 1105, 1269, 1119, 605, 561, 5286),
    Nordeste = c(2592, 4475, 5015, 3699, 1681, 1054, 18516),
    Sudeste = c(4918, 8301, 8228, 6067, 2144, 1082, 30740),
    Sul = c(1680, 3089, 2975, 1933, 686, 305, 10668),
    Centro_Oeste = c(772, 1447, 1409, 1120, 463, 223, 5434),
    Total = c(10589, 18417, 18896, 13938, 5579, 3225, 70644)
)

rownames(tab645) <- tab645$Moradores
tab645$Moradores <- NULL

tab645_ca <- tab645[-nrow(tab645), -ncol(tab645)]

res_ca <- ca(tab645_ca)

plot(res_ca)

tab646 <- data.frame(
    Regiao = c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste", "Total"),
    Branca = c(3279, 14190, 45438, 22813, 5746, 91466),
    Preta = c(1388, 6572, 8580, 1321, 1437, 19298),
    Parda = c(13640, 36323, 34642, 6041, 9256, 99902),
    Total = c(18307, 57085, 88660, 30175, 16439, 210666)
)

rownames(tab646) <- tab646$Regiao
tab646$Regiao <- NULL

tab646 <- tab646 %>% select(-Total)

res_646 <- ca(tab646)

plot(res_646)

# Questão 82 (6.9 Artes e Barroso) #####################################

dados <- data.frame(
    Grau_Instrucao = c(
        "Fundamental incompleto",
        "Médio incompleto",
        "Superior incompleto",
        "Superior completo",
        "Total"
    ),
    Norte = c(125, 56, 91, 16, 288),
    Nordeste = c(356, 108, 240, 55, 759),
    Sudeste = c(155, 124, 259, 76, 614),
    Sul = c(53, 60, 103, 28, 244),
    COeste = c(69, 59, 79, 18, 225),
    Total = c(758, 407, 772, 193, 2130)
)


dados |>
    kable(
        format = "latex",
        booktabs = TRUE,
        align = c("l", rep("c", 6))
    ) |>
    kable_styling(
        latex_options = c("HOLD_position")
    )

tab <- as.matrix(
    dados[dados$Grau_Instrucao != "Total",
          c("Norte", "Nordeste", "Sudeste", "Sul", "COeste")]
)

rownames(tab) <- dados$Grau_Instrucao[dados$Grau_Instrucao != "Total"]

resultado_ca <- ca(tab)

plot(resultado_ca)

summary(resultado_ca)
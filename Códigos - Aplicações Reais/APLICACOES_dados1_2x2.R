library('tidyverse')
library('combinat')
library('glue')


###########################################################
######################### FUNÇÕES #########################
###########################################################
path <- "xxxxxx"
setwd(path)


source(glue('{path}/00. TH_kx2_funcoes_dist.R'))
source(glue('{path}/00. TH_kx2_funcoes_hpd.R'))
source(glue('{path}/05.1. TH_kx2_resultados_funcoes_aplicacao.R'))

#####################################################################
################### DADOS 1 - Aspirina vs Infarto ###################
#####################################################################

# Livro Agresti - Página 72 - Caso 2x2
# População 1: Placebo
# População 2: Aspirina
# Variável Binomial: teve ou não teve infarto
# Dados: 
    # Placebo - Infarto = 1: 28; Infarto = 0: 656 (Total = 684)
    # Aspirina - Infarto = 1: 18; Infarto = 0: 658 (Total = 676)



name_dist = 'BP'
k = 2
n = c(10,10)
confianca = 0.8
prior_par_list = list(c(1, 1), c(1, 1))
dados = as.tibble(t(matrix(c(28, 18)))) # número de infartos
m = c(684, 676)                         # tamanhos amostrais
name_dist = name_dist
eps = 0.1
N = 100000

caminho <- glue("{path}/DadosAplicacao/Aplicacao1/Test_Homog_{k}x2/{name_dist}")
if (!dir.exists(caminho)) {
  dir.create(caminho, recursive = TRUE)
}


# start <- Sys.time()
# gerar_dados_prag(k = k,
#                  n = n,
#                  name_dist = name_dist,
#                  eps = eps,
#                  N = N,
#                  path = path)
# end <- Sys.time()
# end-start


start <- Sys.time()
result_hpd <- gerar_hpd(k=k,
                        n = n,
                        confianca = confianca,
                        prior_par_list = prior_par_list,
                        dados = dados,
                        m = m,
                        name_dist = name_dist,
                        eps = eps,
                        N = N,
                        path = caminho)

end <- Sys.time()
end-start
## Funcionou, mas demora 18 min pra rodar (N=500) / 3h (N=1000)


result_hpd$resultado
# Mesmo resultado para N=500 e N=1000
# # A tibble: 1 × 4
# V1    V2 teste result   
# <dbl> <dbl> <dbl> <chr>    
#   1    28    18     0 Aceita H0
# browseURL(result_hpd$caminho_graf)
result_hpd$grafico_hpd

nome_dados <- 'dados1'
saveRDS(result_hpd, glue("{caminho}/resultado_aplicacao_full_{nome_dados}_N{N}.rds"))


##########################
#caminho_graf = glue("{path}/DadosAplicacao1/Test_Homog_{k}x2/{name_dist}/")
dados_gerados_prag = read_rds(glue("{caminho}/dados_gerados_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))
start <- Sys.time()
estima_propb_prag_hpd(k
                      ,n
                      ,confianca
                      ,prior_par_list
                      ,dados_gerados_prag
                      ,m
                      ,dados_full = dados
                      ,N = 1000
                      ,eps = 0.1)

end <- Sys.time()
end-start
# # A tibble: 1 × 4
# V1    V2 p_intersec_prag_hpd prob_hat_prag
# <dbl> <dbl>               <dbl>         <dbl>
#   1    28    18                   0.808         0.993



############################
####### FREQUENTISTA #######
############################


dados_matriz <- cbind(
  Infarto = as.numeric(dados[1, ]),
  NaoInfarto = m - as.numeric(dados[1, ])
)
rownames(dados_matriz) <- c("Placebo", "Aspirina")
dados_matriz

#### Teste qui-quadrado de homogeneidade (equivalente ao de independência em 2x2)
chisq.test(dados_matriz)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  dados_matriz
# X-squared = 1.7146, df = 1, p-value = 0.1904
# ----- ACEITO H0


##### Teste exato de Fisher (importante para 2x2)
fisher.test(dados_matriz)
# Fisher's Exact Test for Count Data
# 
# data:  dados_matriz
# p-value = 0.1768
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.8236804 3.0256940
# sample estimates:
# odds ratio 
#   1.559785 
# ---------- ACEITO H0


#### Odds ratio com IC
if(!require(epitools)) install.packages("epitools")
library(epitools)
oddsratio(dados_matriz, method = "wald")
# $data
# Infarto NaoInfarto Total
# Placebo       28        656   684
# Aspirina      18        658   676
# Total         46       1314  1360
# 
# $measure
# NA
# odds ratio with 95% C.I. estimate     lower    upper
# Placebo  1.000000        NA       NA
# Aspirina 1.560298 0.8546703 2.848502
# 
# $p.value
# NA
# two-sided  midp.exact fisher.exact chi.square
# Placebo          NA           NA         NA
# Aspirina  0.1481369    0.1767816  0.1444434
# 
# $correction
# [1] FALSE
# 
# attr(,"method")
# [1] "Unconditional MLE & normal approximation (Wald) CI"
# -------------- ODD RATION IC CONTÉM O ZERO - ACEITO H0
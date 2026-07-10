library('tidyverse')
library('combinat')
library('glue')


path <- "xxxxxxxxxxx"
# setwd(path)

source(glue('{path}/00. TH_kx2_funcoes_dist.R'))

####################################################################################
########### SIMULA CASO 2X2 ########################################################
####################################################################################

###### Grid de valores

k <- 2 # número de subpops
n <- rep(50, k)
N <- 10000 # tamanho do grid
eps <- 0.1

thetas <- seq(0.00000001, 0.99999999, length.out = N) # grid de thetas
## criar grid com todas as combinações de k thetas
aux <- expand_grid(thetas, thetas)
if(k>2){
  for(i in 3:(k)){
    aux <- expand_grid(aux, thetas)
    names(aux) <- glue('theta{1:i}')
  }
}
grid_thetas <- aux
grid_thetas

# grid_thetas <- grid_thetas %>% 
#   mutate(theta3 = 0.5) %>% 
#   distinct()
# grid_thetas




############################ LER DADOS
theta0 <- 0.5

dados_gerados_BP_simplex <- readRDS(
  glue("{path}/DadosAplicacao/Simulacoes/dados_gerados_BP_simples_k_{k}_theta0_{theta0}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))

dados_gerados_BP_full <- readRDS(
  glue("{path}/DadosAplicacao/Simulacoes/dados_gerados_BP_full_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))



#######################################
###### PLOTAR GRÁFICOS NO R2 (projeção)
#######################################

## ESCOLHER QUAL DISTÂNCIA e TIPO
dist <- "BP"
tipo <- "simples"

if(dist == "BP"){
  if(tipo == 'simples'){
    dados_gerados <- dados_gerados_BP_simplex
    tipo1 <- paste(tipo, glue('_theta0_{theta0}'), sep = "")
  }
  if(tipo == 'full'){
    dados_gerados <- dados_gerados_BP_full  
    tipo1 <- tipo
  }
}

# if(dist == "KL"){
#   if(tipo == 'simples'){
#     dados_gerados <- dados_gerados_KL_simplex
#     tipo1 <- paste(tipo, glue('_theta0_{theta0}'), sep = "")
#   }
#   if(tipo == 'full'){
#     dados_gerados <- dados_gerados_KL_full 
#     tipo1 <- tipo
#   }
# }


caminho_graf <- glue("{path}/DadosAplicacao/Simulacoes/Graf_testeHomog_k_{k}_{tipo1}_n_{paste(n, collapse = '_')}_eps_{eps}")
fig2 <- plot_prag(k, dados_gerados, thetas, col_3d = 'Viridis') + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_text(aes(x = 0.99, y = 0.9, label = glue(" t1 = t2")))
fig2
# htmlwidgets::saveWidget(fig2, glue("{caminho_graf}.html"))
# browseURL(glue("{caminho_graf}.html"))  
ggsave(glue("{caminho_graf}.png"),
       width = 5,
       height = 5)







####################################################################################
########### SIMULA CASO 3X2 ########################################################
####################################################################################

###### Grid de valores

k <- 3 # número de subpops
n <- rep(10, k)
N <- 100 # tamanho do grid
eps <- 0.1

thetas <- seq(0.00000001, 0.99999999, length.out = N) # grid de thetas
## criar grid com todas as combinações de k thetas
aux <- expand_grid(thetas, thetas)
for(i in 3:(k)){
  aux <- expand_grid(aux, thetas)
  names(aux) <- glue('theta{1:i}')
}
grid_thetas <- aux
grid_thetas

# grid_thetas <- grid_thetas %>% 
#   mutate(theta3 = 0.5) %>% 
#   distinct()
# grid_thetas




############################ LER DADOS
theta0 <- 0.5

dados_gerados_BP_simplex <- readRDS(
  glue("{path}/DadosAplicacao/Simulacoes/dados_gerados_BP_simples_k_{k}_theta0_{theta0}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))

dados_gerados_BP_full <- readRDS(
  glue("{path}/DadosAplicacao/Simulacoes/dados_gerados_BP_full_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))



#######################################
###### PLOTAR GRÁFICOS NO R2 (projeção)
#######################################

## ESCOLHER QUAL DISTÂNCIA e TIPO
dist <- "BP"
tipo <- "full"

if(dist == "BP"){
  if(tipo == 'simples'){
    dados_gerados <- dados_gerados_BP_simplex
    tipo1 <- paste(tipo, glue('_theta0_{theta0}'), sep = "")
  }
  if(tipo == 'full'){
    dados_gerados <- dados_gerados_BP_full  
    tipo1 <- tipo
  }
}



caminho_graf <- glue("{path}/DadosAplicacao/Simulacoes/Graf_testeHomog_k_{k}_{tipo1}_n_{paste(n, collapse = '_')}_eps_{eps}")
fig2 <- plot_prag(k, dados_gerados, thetas, col_3d = 'Viridis') 
fig2
htmlwidgets::saveWidget(fig2, glue("{caminho_graf}.html"))
browseURL(glue("{caminho_graf}.html"))  


# ggsave(glue("{caminho_graf}"),
#        width = 5,
#        height = 5)



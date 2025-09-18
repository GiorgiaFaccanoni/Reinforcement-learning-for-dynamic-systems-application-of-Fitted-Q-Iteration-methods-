#librerie
library(splines)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(e1071)
library(mgcv)
library(gridExtra)

#caricamento dataset e creazione del dataset in un formato tabella
House<-readRDS("RECORD_30.RDS")  #dati ogni 30 minuti per 60 giorni 
House[["Hs"]]     #stato al tempo t, la temperatura al tempo t
House[["Xs"]]     #azione al tempo t 
House[["CSIs"]]   #temperature esterne al tempo t

s<-seq(from=0.5, to=24, by=0.5)
min_30<-rep(s, 60)
giorno<-c()
for (i in 1:60) {giorno<-c(giorno, rep(i, 48))}

house_heat<-data.frame(
  stato<-as.vector(t(House[["Hs"]][,-49])),
  stato_succ<-as.vector(t(House[["Hs"]][,-1])),
  azione<-as.vector(t(House[["Xs"]])),
  giorno<-giorno, 
  min_30<-min_30
)

colnames(house_heat)<-c("stato"  ,
                        "stato_succ",
                        "azione"  ,
                        "giorno", "min_30"  )
#calcolo dei costi 
tariffa_ora<-list(time=c(0,14,20, 34, 42), 
                  cost=c(1,10,5,10,1))

compute.cost <- function( xs, tariff ) {
  COST<-c()
  COST[1] <- 0
  for( i in 2:48) {
    min_30 <- i
    ix <- which( tariffa_ora$time<=min_30)
    ix <- rev(ix)[1]
    COST[i] <-tariffa_ora$cost[ix] * xs[i] 
  }
  return(COST)
}

hour.cost <- matrix(0,nrow=60, ncol=48)
for( i in 1:nrow(hour.cost)){  
  hour.cost[i,] <- compute.cost(House$Xs[i,], tariffa_ora)
}
rowSums(hour.cost)  #costi orari della policy safe
costo<- -as.vector(t(hour.cost))
done<-ifelse(house_heat$min_30==24, 1, 0)
house_heat<-cbind(house_heat,costo, done)
colnames(house_heat)<-c("stato"  ,
                        "stato_succ",
                        "azione"  ,
                        "giorno"   ,
                        "min_30"   ,
                        "costo" ,
                        "done")

write.csv(house_heat, "house_heat_30.csv")


#Analisi esplorativa del dataset "house heating"
dati<-read.csv("house_heat_30.csv")[,-1]
House<- RECORD_30 #dataset senza pre processing
dati$CSIs<-as.vector(t(House[["CSIs"]])) 
#definizione nuova reward
lambda<-600
dati$reward_new<-dati$costo - lambda*abs(dati$stato-18) 

ggplot(dati, aes(x=stato, y=reward_new))+
  geom_point(col="salmon", pch=19)+
  #coord_cartesian(xlim = c(17,19))+ 
  theme_minimal()+
  xlab("Stato corrente") +
  ylab("Reward") + 
  labs(title="Reward in funzione dello stato")

ggplot(dati, aes(x=azione, y=reward_new))+
  geom_point(col="salmon", pch=19)+
  #coord_cartesian(xlim = c(17,19))+ 
  theme_minimal()+
  xlab("Azione") +
  ylab("Reward") + 
  labs(title="Reward in funzione dell'azione")

#divisione del dataset in training e test
index<-1:1920 
training<-dati[index, ] #training formato dai primi 40 giorni
stato_succ<-matrix(training$stato_succ, byrow=T, ncol=48, nrow=nrow(training)/48)
stato<- matrix(training$stato, byrow=T, ncol=48, nrow=nrow(training)/48)
test<- dati[-index, ] #test formato dagli ultimi 20 giorni

#inizializzazione per algoritmo FQI
gamma <- 0.99 
azioni_possibili <- seq(0, 70, by = 0.5)
n_iter <-50
q_mod <- NULL
media_reward<-c()
std_reward<-c()
gruppi <- rep(1:(nrow(training)/48), each = 48)
costo_temp<-matrix(0, ncol=48, nrow=nrow(training)/48)
reward_temp<-matrix(0, ncol=48, nrow=nrow(training)/48)
mse<-c()
q_next_mat<-matrix(0, nrow=n_iter, ncol=length(unique(training$stato_succ)))

#Algoritmo Fitted Q-iteration
set.seed(123)
for (k in 1:n_iter) {
  print(k)
  
  target_q <- training$reward_new
  
  if (!is.null(q_mod)) {
    # Pre-calcolo tutti i Q-values per stati unici 
    stati_unici <- unique(training$stato_succ)
    q_next <-  sapply(stati_unici, function(s) {
      new_data <- data.frame(stato = s, azione = azioni_possibili)
      p<-predict(q_mod, new_data)
      out<-c(max(p),azioni_possibili[which.max(p)])
      return (out)
    }) 
    
    idx <- match(training$stato_succ, stati_unici)
    temp<-data.frame(stato=training$stato_succ, azione=q_next[2,][idx])
    
    azione_mat<-matrix(temp$azione, ncol=48, nrow=nrow(training)/48, byrow=T)
    for( i in 1:nrow(azione_mat)){  
      costo_temp[i,] <- compute.cost(azione_mat[i,], tariffa_ora)
    }
    
    for ( i in 1:nrow(azione_mat)){
      reward_temp[i,]<- -costo_temp[i,] - lambda*abs(18 - stato[i,])}
    reward_daily <- rowSums(reward_temp)
    media_reward[k]<-mean(reward_daily)
    std_reward[k]<-sd(reward_daily)
    
    target_q[training$done != 1] <- target_q[training$done != 1] + gamma * q_next[1,][idx][training$done != 1]
    q_next_mat[k,]<-q_next[1,]
    }
  
  #Modelli utilizzati: 
  
  #LM
  #q_mod<-lm(target_q~stato + azione + stato*azione, data=training)
  
  #SPLINES CUBICHE
  #q_mod<-mgcv::gam(target_q~s(stato, bs="cr", k=7)+ s(azione, bs="cr", k=7)+ azione*stato,  data=training)
  
  #SPLINES SMOOTHING (thin plate)
  q_mod<-mgcv::gam(target_q~s(stato)+ s(azione)+ azione*stato,  data=training, method="REML")
  
  #B-SPLINE
  #q_mod<-mgcv::gam(target_q~s(stato, bs="bs", k=7)+ s(azione, bs="bs", k=7)+ azione*stato,  data=training)
  
  #NATURAL CUBIC SPLINE
  #q_mod<-lm(target_q~ns(stato,  df=3,  knots=quantile(training$stato, c(0.25,0.5,0.75)))+ ns(azione, df=3,   knots=quantile(training$azione, c(0.25,0.5,0.75)))+ azione*stato,  data=training)
  
  #RANDOM FOREST
  #q_mod <- randomForest::randomForest(target_q ~ stato + azione, data = training, ntree = 200, nodesize = 5)
  
  mse[k] <- mean((predict(q_mod) - target_q)^2)
}

#valori Q
q_next_mat[50, ]

# Funzione per calcolo dell'azione ottimale
optimal_action <- function(state) {
  q_vals <- predict(q_mod, data.frame(stato = state, azione = azioni_possibili))
  azioni_possibili[which.max(q_vals)]
}

#Aggiornamento del nuovo stato, azioni ottime sul test set
#parametri pre-impostati
m = 1470
c = 1005.4
M = 3600/20 
R = 4.329 * 10^-5 

h.ref <- 18
h0 <- 16
dt <- 1/2
Tf <- 24
ts <- seq(0,Tf-dt,by=dt)

Kp=17.5; Ki=5.9; Kd=1.4; max.heat=70

stato_nuovo<-numeric(nrow(test))
stato_nuovo[1]<- test$stato[1]
azione_ottima<-c()
stato_oss<- test$stato

for (i in 1:nrow(test)){
  if (test$done[i]==0){
    azione_ottima[i]<-optimal_action(stato_nuovo[i])
    dQg = M * c * (azione_ottima[i]-stato_nuovo[i])
    dQl = (stato_nuovo[i]-test$CSIs[i]) / R  
    stato_nuovo[i+1] = stato_nuovo[i] + (dt/(m*c)) * (dQg-dQl)}
  else {
    azione_ottima[i]<-optimal_action(stato_nuovo[i])
    stato_nuovo[i+1]<-stato_oss[i+1]
  }
}

stato_nuovo<-stato_nuovo[-961]
round(quantile(stato_nuovo,c(0.01, 0.25,0.5,0.75, 1)),2) #stato
table(azione_ottima)  #azioni ottime
sum(stato_nuovo<17.5) #limite inferiore
sum(stato_nuovo>18.5) #limite superiore

#Calcolo dei costi con la nuova policy
azione_mat<-matrix(azione_ottima, ncol=48, nrow=20, byrow=T)
costo<-matrix(0, ncol=48, nrow=20)

for( i in 1:nrow(azione_mat)){  
  costo[i,] <- compute.cost( azione_mat[i,], tariffa_ora)
}
costo<-as.vector(t(costo))
sum(costo)      #costi toali con nuova policy
sum(test$costo) #costi totali con policy safe-by-design
sum(costo- - test$costo) #differenza tra i costi 

#--------------------------------------------------
#Aggiornamento stato e azioni nuove sul dataset completo
stato_nuovo<-numeric(nrow(dati))
stato_nuovo[1]<- dati$stato[1]
azione_ottima<-c()
stato_oss<- dati$stato
for (i in 1:nrow(dati)){
  if (dati$done[i]==0){
    azione_ottima[i]<-optimal_action(stato_nuovo[i])
    dQg = M * c * (azione_ottima[i]-stato_nuovo[i])   
    dQl = (stato_nuovo[i]-dati$CSIs[i]) / R   
    stato_nuovo[i+1] = stato_nuovo[i] + (dt/(m*c)) * (dQg-dQl)}
  else {
    azione_ottima[i]<-optimal_action(stato_nuovo[i])
    stato_nuovo[i+1]<-stato_oss[i+1]
  }
}
stato_nuovo<-stato_nuovo[-2881]
round(quantile(stato_nuovo,c(0.01, 0.25,0.5,0.75, 1)),2)
table(azione_ottima)
sum(stato_nuovo<17.5)
sum(stato_nuovo>18.5)
#Costi
azione_mat<-matrix(azione_ottima, ncol=48, nrow=60, byrow=T)
costo<-matrix(0, ncol=48, nrow=60)
for( i in 1:nrow(azione_mat)){  
  costo[i,] <- compute.cost( azione_mat[i,], tariffa_ora)
}
costo<-as.vector(t(costo))
sum(costo)
sum(dati$costo)
sum(costo--dati$costo)


#---------------------------------------------------
#Grafico media della reward cumulata giornaliera
ggplot_with_ci <- function(media, se, 
            xlab = "Iterazioni", 
            ylab = "Valore", 
            title = "Andamento con intervalli di confidenza",
            level = 0.95) {
  n <- length(media)
  z_value <- qnorm(1 - (1 - level)/2)
  df <- data.frame(
    iterazione = 1:n,
    media = media,
    lower = media - z_value * se,
    upper = media + z_value * se
  )
  
  ggplot(df, aes(x = iterazione, y = media)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "skyblue", alpha = 0.3) +  # IC
    geom_line(color = "dodgerblue4", size = 1) +  # media
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal()
}

ggplot_with_ci(media_reward, std_reward,
    ylab = "Media delle rewards cumulate giornaliere",
    title = "Andamento delle rewards durante l'addestramento")


#grafico dell'andamento della temperatura nel tempo
dati_grafico <- data.frame(
tempo = 1:length(stato_nuovo),
stato = stato_nuovo
)

ggplot(dati_grafico, aes(x = tempo, y = stato)) +
  geom_line(color = "dodgerblue4", linewidth = 0.8) +
  geom_hline(yintercept = 17.5, color = "salmon", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 18.5, color = "salmon", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Andamento della temperatura nel tempo",
    subtitle = "Con vincoli di sicurezza (17.5°C - 18.5°C)",
    x = "Tempo (ogni 30 minuti)",
    y = "Next state (Temperatura in °C)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(dati_grafico$tempo), by = 48),
    minor_breaks = NULL  
  ) +
  scale_y_continuous(
    breaks = sort(unique(c(
      seq(floor(min(dati_grafico$stato)),
          ceiling(max(dati_grafico$stato)),
          by = 0.5),
      17.5, 18.5
    )))
  )


#grafico della densità dei costi giornalieri
gruppi <- rep(1:(length(test$costo)/48), each = 48)
costi_oss <- tapply(-test$costo, gruppi, sum) #con policy safe-by-design
costi_new <- tapply(costo, gruppi, sum) #con nuova policy

df <- data.frame(
  valore = c(costi_oss, costi_new),
  gruppo = factor(rep(c("Costi policy safe-by design", "Costi nuova policy"),
                      times = c(length(costi_oss), length(costi_new))))
)

p <- ggplot(df, aes(x = valore, fill = gruppo)) + 
  geom_density(alpha = 0.5, linewidth = 0.8, color = NA) + 
  scale_fill_manual(values = 
                  c("Costi policy safe-by design" = "salmon",
                  "Costi nuova policy" = "dodgerblue4"),
                  name = "Gruppi") +  
  
  labs(title = "Confronto distribuzioni dei costi con policy diverse",
       x = "Costi giornalieri (€/per giorno)",
       y = "Densità") +
  
  theme_bw() +
  theme(panel.border = element_blank(),  
   panel.grid.major = element_line(color = "grey85", size = 0.3),
   panel.grid.minor = element_line(color = "grey90", size = 0.2))
# Aggiunta mediane 
medie <- aggregate(valore ~ gruppo, df, median) 
p <- p + geom_vline(data = medie,
      aes(xintercept = valore, color = gruppo),
      linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  scale_color_manual(values = 
              c("Costi policy safe-by design" = "salmon",
              "Costi nuova policy" = "dodgerblue4"))
print(p)

#CONFORMAL 
dati<-read.csv("house_heat_30.csv")[,-1]
lambda<-600
dati$reward_new<-dati$costo - lambda*abs(dati$stato-18)
House<- RECORD_30 
dati$CSIs<-as.vector(t(House[["CSIs"]]))

#DIVIDO IN TRAIN CALIBRATIONE TEST 
index<-1:1440 #30 giorni
training<- dati[index, ]
cal<- dati[1441:1920, ] #10 giorni
test<- dati[1921:2880, ] #20 giorni
stato_succ<-matrix(training$stato_succ, byrow=T, ncol=48, nrow=nrow(training)/48)
stato<- matrix(training$stato, byrow=T, ncol=48, nrow=nrow(training)/48)
gamma <- 0.99
azioni_possibili <- seq(0, 70, by = 0.5)
n_iter <-50
q_mod <- NULL
media_reward<-c()
std_reward<-c()
gruppi <- rep(1:(nrow(training)/48), each = 48)
tariffa_ora<-list(time=c(0, 14,20, 34, 42), 
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
costo_temp<-matrix(0, ncol=48, nrow=nrow(training)/48)
reward_temp<-matrix(0, ncol=48, nrow=nrow(training)/48)
mse<-c()
q_next_mat<-matrix(0, nrow=n_iter, ncol=length(unique(training$stato_succ)))
set.seed(123)
#alleno il modello con lambda 600 e spline smoothing
for (k in 1:n_iter) {
  print(k)
  
  target_q <- training$reward_new
  
  if (!is.null(q_mod)) {
    # Pre-calcolo tutti i Q-values per stati unici 
    stati_unici <- training$stato #unique(training$stato_succ)
    tempo <- training$min_30
    q_next <-  matrix(NA, ncol=2, nrow=1440)
    for (i in 1:nrow(training)){
      new_data <- data.frame(stato = stati_unici[i], azione = azioni_possibili, min_30=tempo[i])
      p<-predict(q_mod, new_data)
      out<-c(max(p),azioni_possibili[which.max(p)])
      q_next[i,]<- out
      i<- i+1
    }
    #   sapply(c(stati_unici, training$min_30), function(s, tempo) {
    #   new_data <- data.frame(stato = s, azione = azioni_possibili, min_30=tempo)
    #   p<-predict(q_mod, new_data)
    #   out<-c(max(p),azioni_possibili[which.max(p)])
    #   return (out)
    # }) #mi crea una matrice con 2 righe e 6996 colonne
    
    idx <- match(training$stato_succ, stati_unici)
    temp<-data.frame(stato=training$stato_succ, azione=q_next[2,][idx])
    
    #che per la reward
    
    azione_mat<-matrix(temp$azione, ncol=48, nrow=nrow(training)/48, byrow=T)
    for( i in 1:nrow(azione_mat)){  
      costo_temp[i,] <- compute.cost(azione_mat[i,], tariffa_ora)
    }
    
    for ( i in 1:nrow(azione_mat)){
      reward_temp[i,]<- -costo_temp[i,] - lambda*abs(18 - stato[i,])}
    reward_daily <- rowSums(reward_temp)
    media_reward[k]<-mean(reward_daily)
    std_reward[k]<-sd(reward_daily)
    
    # Mappatura veloce
    target_q[training$done != 1] <- target_q[training$done != 1] + gamma * q_next[,1][idx][training$done != 1]
    q_next_mat[k,]<-q_next[1,]
  }
  
  #LM
  #q_mod<-lm(target_q~stato + azione + stato*azione, data=training)
  #SPLINES CUBICHE
  #q_mod<-mgcv::gam(target_q~s(stato, bs="cr", k=16) +s(azione, bs="cr", k=12)+ stato*azione, data=training, method="REML")
  q_mod<-mgcv::gam(target_q~s(stato)+ s(azione)+ azione*stato + min_30,  data=training, method="REML")
  
  #MARs
  #q_mod<-earth::earth(target_q~stato+ azione, data = training, degree = 2)
  
  #SPLINES SMOOTHING
  #q_mod<-lm(target_q~ +stato*azione, data=training)
  
  #RANDOM FOREST
  #q_mod <- randomForest::randomForest(target_q ~ stato + azione, data = training, ntree = 200, nodesize = 5)
  
  #SVM 
  #q_mod <- e1071::svm(target_q ~ stato + azione, data = training, kernel="radial",  type="eps-regression", gamma=1)
  
  mse[k] <- mean((predict(q_mod) - target_q)^2)
}


optimal_action <- function(state, tempo) {
  q_vals <- predict(q_mod, data.frame(stato = state, azione = azioni_possibili, min_30=tempo))
  azioni_possibili[which.max(q_vals)]
}

#campiono un sottoinsieme del training, considero 10 giorni 480 osservazioni
oss<- seq(1:480)
azioni<-training$azione[oss]
stati<- training$stato[oss]
done<- training$done[oss]
csis<- training$CSIs[oss]
#Aggiornamento 
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

stati_nuovo<-numeric(length(stati))
stati_nuovo[1]<- stati[1]
azioni_ottime<-c()
stati_oss<- stati

for (i in 1:length(stati)){
  if (done[i]==0){
    azioni_ottime[i]<-optimal_action(stati_nuovo[i], tempo[i])
    dQg = M * c * (azioni_ottime[i]-stati_nuovo[i])     # gain
    dQl = (stati_nuovo[i]-csis[i]) / R   # loss
    stati_nuovo[i+1] = stati_nuovo[i] + (dt/(m*c)) * (dQg-dQl)}
  else {
    azioni_ottime[i]<-optimal_action(stati_nuovo[i], tempo[i])
    stati_nuovo[i+1]<-stati_oss[i+1]
  }
}


azioni_tot<- c(azioni, azioni_ottime) #azioni osservate e azioni calcolate tramite modello
y<- c(rep(1, times=480),rep(0, times=480)) #1 se azione osservata, 0 se azione stimata
dd<-data.frame(azioni=azioni_tot,  y=y)
dd$azioni<- round(dd$azioni*2)/2
dd$azioni
#modello di regressione logistica
logistica<- glm(y~azioni, data=dd, family = binomial)


# --------------------- 
#modifico per calibration, salvo anche i q_val
optimal_action_cal <- function(state, tempo) {
  q_vals <- predict(q_mod, data.frame(stato = state, azione = azioni_possibili, min_30=tempo))
  out<- c(max(q_vals), azioni_possibili[which.max(q_vals)])
  out
}

tempo_cal=cal$min_30
stati_cal<-numeric(length(cal$stato))
stati_cal[1]<- cal$stato[1]
azioni_cal<-c()
stati_oss_cal<- cal$stato
q_val_cal<- c()
for (i in 1:length(cal$stato)){
  if (cal$done[i]==0){
    azioni_cal[i]<-optimal_action_cal(stati_cal[i], tempo_cal[i])[2]
    q_val_cal[i]<-optimal_action_cal(stati_cal[i], tempo_cal[i])[1]
    dQg = M * c * (azioni_cal[i]-stati_cal[i])     # gain
    dQl = (stati_cal[i]-cal$CSIs[i]) / R   # loss
    stati_cal[i+1] = stati_cal[i] + (dt/(m*c)) * (dQg-dQl)}
  else {
    azioni_cal[i]<-optimal_action_cal(stati_cal[i], tempo_cal[i])[2]
    q_val_cal[i]<-optimal_action_cal(stati_cal[i], tempo_cal[i])[1]
    stati_cal[i+1]<-stati_oss_cal[i+1]
  }
}


#calcolo i pesi per ogni azione del calibration
previsti<-predict.glm(logistica, newdata=data.frame(azioni=azioni_cal), type="response")
pesi<- previsti/(1-previsti) #pesi 

#calcolo della score function
giorni<- c(31,32, 33, 34, 35, 36, 37, 38, 39, 40)
score<-c()
score1<-c()
for(j in giorni ){
  reward_cal<- cal$reward_new[cal$giorno==j]
  for(i in 1:length(reward_cal)){
    score[i]<-abs(gamma*sum(reward_cal[i:length(reward_cal)])-q_val_cal[i])
  }
  score1<-c(score1, score)
  score<-c()
}

#peso la score function con i pesi
score1<- score1*pesi
n<-nrow(cal)
q<- ceiling((n+1)*(1-0.1)) #calcolo il quantile con alpha=0.1
score_ord<-sort(score1)
quantile<-score_ord[q]

#test: calcolo le azioni e i q-values ottimi
optimal_action_cal <- function(state, tempo) {
  q_vals <- predict(q_mod, data.frame(stato = state, azione = azioni_possibili, min_30=tempo))
  out<- c(max(q_vals), azioni_possibili[which.max(q_vals)])
  out
}

tempo_test=test$min_30
stati_test<-numeric(length(test$stato))
stati_test[1]<- test$stato[1]
azioni_test<-c()
stati_oss_test<- test$stato
q_val_test<- c()


for(i in 1:length(test$stato)){
  azioni_test[i]<-optimal_action_cal(stati_test[i], tempo_test[i])[2]
  q_val_test[i]<-optimal_action_cal(stati_test[i], tempo_test[i])[1]
  stati_test[i+1]<-stati_oss_test[i+1]
}

#intervallo
PI<- data.frame(lo=q_val_test-quantile, up=q_val_test+quantile)
quant<- rep(quantile, 960)

intervalli<- data.frame(q_values_test=q_val_test, quantile=quant, PI_lo=PI$lo, PI_up=PI$up )

#coverage con reward cumulata
giorni<- unique(test$giorno)
reward_cum<-c()
reward_cum1<-c()
for(j in giorni ){
  reward_test<- test$reward_new[test$giorno==j]
  for(i in 1:length(reward_test)){
    reward_cum[i]<- sum(reward_test[i:length(reward_test)])
  }
  reward_cum1<-c(reward_cum1, reward_cum)
  reward_cum<-c()
}


coverage <- mean((reward_cum1 > intervalli$PI_lo) & (reward_cum1 < intervalli$PI_up))
int_len <- round(mean( intervalli$PI_up - intervalli$PI_lo), 5)

#write.csv(intervalli, "PI_spline600.csv")
#----------------------------------------------------------------------------------------
#sub-sample del calibration
cal<- dati[1441:1920,]
set.seed(123)
ind<- sample(1:480, 360)
ind<-sort(ind)
cal<-cal[ind,]


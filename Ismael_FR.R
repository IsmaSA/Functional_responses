suppressMessages({
  library(dplyr, quiet = TRUE, warn.conflicts = FALSE)
  library(reshape, quiet = TRUE, warn.conflicts = FALSE)
  library(ggplot2)
  library(stringr)
  library(tidyr)  
  library(stringr)
  library(readxl)
  library(frair)
  library(patchwork)
  library(glmmTMB)
})
library(openxlsx)


cray=read.xlsx("EXP 1.xlsx")



Cray_0=cray%>%filter(Claw=="0")
Cray_1=cray%>%filter(Claw=="1")
Cray_2=cray%>%filter(Claw=="2")


Fair_0_res = frair_test(formula = Tot  ~ Density , data = Cray_0) 
Fair_1_res = frair_test(formula = Tot  ~ Density , data = Cray_1) 
Fair_2_res = frair_test(formula = Tot  ~ Density , data = Cray_2) 


# All the All the crayfish groups follows a type-II response


# PLOT 0 claw
g_fit_0_claw <- frair_fit(formula = Tot  ~ Density, data = Cray_0,
                    response = "rogersII",
                    start = list(a = 1, h = 0.1),
                    fixed = list(T = 4/24)) 
with(Cray_0, plot(Density, Eat, xlab="Prey Density",
               ylab="No. Prey killed"))
lines(g_fit_0_claw, lty = 1, col = "grey25")

print(g_fit_0_claw)

summary(g_fit_0_claw$fit)

# Coefficients:
#   a     h     T 
# 2.101 0.039 0.167


# PLOT 1 claw
g_fit_1_claw <- frair_fit(formula = Tot  ~ Density, data = Cray_1,
                          response = "rogersII",
                          start = list(a = 1, h = 0.1),
                          fixed = list(T = 4/24)) 
with(Cray_1, plot(Density, Eat, xlab="Prey Density",
                  ylab="No. Prey killed"))
lines(g_fit_1_claw, lty = 1, col = "grey25")

print(g_fit_1_claw)
summary(g_fit_1_claw$fit)

# Coefficients:
#   a     h     T 
# 7.286 0.030 0.167  

# PLOT 2 claw
g_fit_2_claw <- frair_fit(formula = Tot  ~ Density, data = Cray_2,
                          response = "rogersII",
                          start = list(a = 1, h = 0.1),
                          fixed = list(T = 4/24)) 
with(Cray_2, plot(Density, Eat, xlab="Prey Density",
                  ylab="No. Prey killed"))
lines(g_fit_2_claw, lty = 1, col = "grey25")


print(g_fit_2_claw)

# Coefficients:
#   a     h     T 
# 10.094  0.017  0.167  


fitp_fit1 <- frair_boot(g_fit_0_claw)
confint(fitp_fit1 , citypes ='bca')


# Coefficient  CI Type        Lower   Upper   
# a            BCa            1.081   4.758   
# h            BCa            0.019   0.069


fitp_fit2 <- frair_boot(g_fit_1_claw)
confint(fitp_fit2 , citypes ='bca')

# Coefficient  CI Type        Lower   Upper   
# a            BCa            2.986   15.906  
#h            BCa            0.015   0.043  

fitp_fit3 <- frair_boot(g_fit_2_claw)
confint(fitp_fit3 , citypes ='bca')

# Coefficient  CI Type        Lower   Upper   
# a            BCa            5.279   18.697  
# h            BCa            0.012   0.022 


# PLOT

plot(x = 1,
     type = "n",
     xlim = c(0, 36), 
     ylim = c(0, 15),
     xlab="Prey Density", ylab="No. of prey killed")

drawpoly(fitp_fit1, col=("royalblue1"), border=NA, tozero=TRUE)
drawpoly(fitp_fit2, col=("palegreen2"), border=NA, tozero=TRUE)
drawpoly(fitp_fit3, col=("orangered1"), border=NA, tozero=TRUE)


points(x = Cray_1$Density, y = Cray_1$Tot,   cex=0.9, col= "palegreen2",alpha = 0.5)
points(x = Cray_2$Density, y = Cray_2$Tot,  cex=0.9, col= "orangered1",alpha = 0.5)
points(x = Cray_0$Density, y = Cray_0$Tot,  cex=0.9, col= "royalblue1",alpha = 0.5)
#?points

lines(g_fit_0_claw, col = "royalblue4", lty=2, lwd=3) ##Group Non-reproductive
lines(g_fit_1_claw, col = "palegreen4",lty=2,  lwd=3 ) #Glands
lines(g_fit_2_claw, col = "orangered4",lty=2,  lwd=3 ) #Eggs

legend(1, 15, legend=c("0 claw", "1 claw","2 claw"),
       col=c("royalblue1", "palegreen2","orangered1"), lty=1:2, cex=0.8)




##Rate attack

a<- fitp_fit1$bootcoefs %>% as.data.frame() %>% mutate(Lower="1.081",
                                                       Upper="4.758",
                                                       point="2.101") 
a$Claw <- "0 Claw"
b<- fitp_fit2$bootcoefs %>% as.data.frame()%>% mutate(Lower="2.986",
                                                      Upper="15.906",
                                                      point="7.286")
b$Claw <- "1 Claw"
c<- fitp_fit3$bootcoefs %>% as.data.frame()%>% mutate(Lower="5.279",
                                                      Upper="18.697",
                                                      point="10.094")
c$Claw <- "2 Claw"

df<- bind_rows(a,b,c)
str(df)


df$Lower<- as.numeric(df$Lower)
df$Upper<- as.numeric(df$Upper)
df$point<- as.numeric(df$point)

df$Clawless <- factor(df$Claw, 
                      levels = c("0 claw","1 claw","2 claw"))

df <- df[!duplicated(df$Claw),]

p<- ggplot(df, aes(Claw,a,colour = Claw))
p1<- p + geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +theme_bw() + ylab("Attack rate (± SE)")+
  ylim(0,20)+ geom_point(data=df,aes(Claw,point)) +xlab("")
p1


###Handling time

d<- fitp_fit1$bootcoefs %>% as.data.frame() %>% mutate(Lower="0.019",
                                                       Upper="0.069",
                                                       point="0.039") 
d$Claw <- "0 Claw"
e<- fitp_fit2$bootcoefs %>% as.data.frame()%>% mutate(Lower="0.015",
                                                      Upper="0.043",
                                                      point="0.030")
e$Claw <- "1 Claw"
f<- fitp_fit3$bootcoefs %>% as.data.frame()%>% mutate(Lower="0.012",
                                                      Upper="0.022",
                                                      point="0.017")
f$Claw <- "2 Claw"

df1<- bind_rows(d,e,f)
str(df1)
?log

df1$Lower<- as.numeric(df1$Lower)
df1$Upper<- as.numeric(df1$Upper)
df1$point<- as.numeric(df1$point)

df1$Claw <- factor(df1$Claw, 
                       levels = c("0 Claw","1 Claw","2 Claw"))

df1 <- df1[!duplicated(df1$Claw),]



p3<- ggplot(df1, aes(Claw,h,colour = Claw)) 
p2<-p3 + geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +theme_bw() + 
  ylab("Handling time (± SE)") + geom_point(data=df1,aes(Claw,point))+
  xlab("")

p2

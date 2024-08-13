#Multiple soil conditioner addition - Incubation experiment -- Special issue NHESS
#Juan F. Dueñas Ph.D. https://github.com/juan-duenas

# This code reproduces the analysis and the main figure presented in the paper.
# ATTENTION: there are some lines of code within functions that need to be (un)comment by the user before running the code
# The lines to be (un)comment as well as the objects where the results are stored depend on the model object and variable analyzed. 
# every instance in which this is the case is signaled by a message next to the line of code
# I am sure there is a more efficient way to do this. You are welcome to improve the code.

#load packages
pkgs <- c("tidyverse", "ggpubr", "car", "MASS", "emmeans", "betareg", "boot", "ggrepel")

vapply(pkgs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)

#path Github repo
URL <- "https://raw.githubusercontent.com/juan-duenas/NHESS/main/MADdb.csv" 
'%notin%' <- Negate('%in%') #custom function

#load dataset
MAD <- read_csv2(URL) %>%
         mutate(whc1=whc/100)%>%
         mutate(wsa1=wsa/100)%>%
         mutate(Nfactors = as.factor(Nfactors))%>%
         mutate(Amendment = fct_relevel(Amendment, c("Control","Biochar","Compost", "Microbial wash",
                                              "Silica","Straw","Control Biochar", "3 Factors","5 Factors")))

MAD1 = MAD%>%filter(Amendment %notin% c("Control Biochar"))%>%droplevels()
MAD2 = MAD%>%filter(Amendment %notin% c("Biochar","Compost", "Microbial wash",
                                        "Silica","Straw"))%>%droplevels()


# multiple amendments against control 1####
#Models against control 1 
m1 <- betareg(whc1~Nfactors, MAD1)
print(m1)

m2 <- betareg(wsa1~Nfactors, MAD1)
summary(m2)

m3 <- lm(pH~Nfactors, MAD1)
print(m3)

m4 <- lm(ratio~Nfactors, MAD1)
print(m4)

# boot function - !!ATENTION!! replace models and (un)comment lines when necessary
bstat <- function (data, i)
{
  d <- data [i,]
  fit <- update(m3, data=d) # replace model object here
  #ko <- predict(fit, newdata = with(d, expand.grid(Nfactors = fit$levels$mean)), type = "response") # (un)comment when passing betareg objects - comment when other models are run
  ko <- predict(fit, newdata = with(d, expand.grid(Nfactors = fit$xlevels$Nfactors)), type = "response") # (un)comment when passing lm objects - comment when other models are run
  #store updated model mean predictions
  rtn <- c("1"=NA, "2"=NA,  "3"=NA, "4"=NA) #dummy logical vector
  rtn[names(ko)] <- ko #replace missing values in predict() vector, leaving NAs when necessary
  rtn 
}

set.seed(1234) # run to make results reproducible
br <- boot(statistic = bstat, data = MAD1, R =5000, parallel = "multicore", ncpus = 4) # option multicore does not work on windows
plot(br)

bres <- list() # create list to store CIs from boot
for (i in 1:ncol(br$t[,1:4])){
  bres[[i]] <- boot.ci(br, level = .95, type = c("bca"),index = i) # extract CI of the BCA type
}

# loop to extract CI from boot object
media <- br$t0[1:4] # 'media' means mean in Spanish
cdown <- c() 
cup <- c()
Nfactors <- c(0,1,3,5)

for (i in 1:length(bres[[i]])){
  cdown[i] <- bres[[i]]$bca[,4]
  cup[i] <- bres[[i]]$bca[,5]
}

#get results together in a data frame in order to plot - Attention - run each chunk of code for the appropriate model results signaled by the name of the object
CIs.m1 <- cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(media=round(media*100,2))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cdown=round(cdown*100,2))%>%
  mutate(cup=as.numeric(cup))%>%
  mutate(cup=round(cup*100,2))

CIs.m2 <- cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(media=round(media*100,2))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cdown=round(cdown*100,2))%>%
  mutate(cup=as.numeric(cup))%>%
  mutate(cup=round(cup*100,2))

CIs.m3 <- cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cup=as.numeric(cup))

CIs.m4 <- cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cup=as.numeric(cup))

#Plot first part of Fig. 1 - Attention - use the appropiate CIs.xx object for each panel
colorp <- c('#000000','#E69F00','#56B4E9','#009E73') # colorblind palette
labelsx_legend <- c("None","One (SD)","Three","Five")

#WHC
pwhc1  <- ggplot() +
  geom_point(CIs.m1, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m1, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("A")+
  ylab(expression("WHC (% gg"^"-1"*")")) + xlab("")+
  ggdist::stat_halfeye(data=MAD1,
                       density = "bounded",
                       scale=0.65, # (un)comment when using density="histogram"
                       adjust = 1,
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=whc, fill=Nfactors),
                       alpha=0.5
  ) +
  scale_fill_manual(values=colorp) +
  geom_hline(yintercept=CIs.m1[1,1], linetype="dashed")+
  scale_y_continuous(breaks=c(45,55,65), limits = c(40,67.5))+
  scale_x_discrete(labels= labelsx_legend)+
  guides(color="none", fill="none")+
  theme_bw()
pwhc1

#WSA
pwsa1  <- ggplot() +
          geom_point(CIs.m2, mapping=aes(Nfactors, media, color=Nfactors))+
          geom_errorbar(CIs.m2, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
          scale_color_manual(values=colorp) + 
          ggtitle("C")+
          ylab("WSA (%)") + xlab("")+
          ggdist::stat_halfeye(data=MAD1,
                       density = "bounded",
                       scale=0.65, # (un)comment when using density="histogram"
                       adjust = 1,
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=wsa, fill=Nfactors),
                       alpha=0.5
                      ) +
          scale_fill_manual(values=colorp) +
          geom_hline(yintercept=CIs.m2[1,1], linetype="dashed")+
          scale_y_continuous(breaks=c(40,60,80), limits = c(28,82))+
          scale_x_discrete(labels= labelsx_legend)+
          guides(color="none", fill="none")+
        theme_bw()
pwsa1

#pH
pph1  <- ggplot() +
  geom_point(CIs.m3, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m3, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("E")+
  ylab(expression("pH")) + xlab("")+
  ggdist::stat_halfeye(data=MAD1,
                       density = "bounded",
                       scale=0.65, # (un)comment when using density="histogram"
                       adjust = 1,
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=pH, fill=Nfactors),
                       alpha=0.5
  ) +
  scale_y_continuous(breaks=c(6.0, 6.5, 7.0), limits = c(5.65, 7.1))+
  scale_x_discrete(labels= labelsx_legend)+
  scale_fill_manual(values=colorp) +
  geom_hline(yintercept=CIs.m3[1,1], linetype="dashed")+
  guides(color="none", fill="none")+
  theme_bw()
pph1

#B:F ratio
prat1  <- ggplot() +
  geom_point(CIs.m4, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m4, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("G")+
  ylab(expression("B:F ratio")) + xlab("")+
  ggdist::stat_halfeye(data=MAD1,
                       density = "bounded",
                       scale=0.65, # (un)comment when using density="histogram"
                       adjust = 1,
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=ratio, fill=Nfactors),
                       alpha=0.5
  ) +
  scale_fill_manual(values=colorp) +
  scale_x_discrete(labels= labelsx_legend)+
  geom_hline(yintercept=CIs.m4[1,1], linetype="dashed")+
  guides(color="none", fill="none")+
  theme_bw()
prat1

#emmeans
con.m1 <- emmeans(m1, specs = trt.vs.ctrl~Nfactors, type="response")
con.m1$contrasts %>% summary(infer=TRUE)

con.m2 <- emmeans(m2, specs = trt.vs.ctrl~Nfactors, type="response")
con.m2$contrasts %>% summary(infer=TRUE)

con.m3 <- emmeans(m3, specs = trt.vs.ctrl~Nfactors, type="response")
con.m3$contrasts %>% summary(infer=TRUE)

con.m4 <- emmeans(m4, specs = trt.vs.ctrl~Nfactors, type="response")
con.m4$contrasts %>% summary(infer=TRUE)

#Table with broom
ts1<- rbind(broom::tidy(m1)[,2:6],broom::tidy(m2)[,2:6],broom::tidy(m3), broom::tidy(m4))
write.table(ts1, paste(path,"/ts1.txt",sep=""))
# multiple amendments against control 2 ####

#Models against control 2 
m1.1 <- betareg(whc1~Nfactors, MAD2)
print(m1.1)

m2.1 <- betareg(wsa1~Nfactors, MAD2)
summary(m2.1)

m3.1 <- lm(pH~Nfactors, MAD2)
summary(m3.1)

m4.1 <- lm(ratio~Nfactors, MAD2)
summary(m4.1)

## boot function - !!ATENTION!! replace models and (un)comment lines when necessary
bstat2 <- function (data, i)
{
  d <- data [i,]
  fit <- update(m3.1, data=d)
  #ko <- predict(fit, newdata = with(d, expand.grid(Nfactors = fit$levels$mean)), type = "response") # (un)comment when passing betareg objects - comment when other models are run
  ko <- predict(fit, newdata = with(d, expand.grid(Nfactors = fit$xlevels$Nfactors)), type = "response") # (un)comment when passing lm objects - comment when other models are run
  #store updated model mean predictions
  rtn <- c("1"=NA, "2"=NA, "3"=NA, "4"=NA) #dummy logical vector
  rtn[names(ko)] <- ko #replace missing values in predict() vector, leaving NAs when necessary
  rtn 
}

set.seed(123) # run to make results reproducible
br1.2 <- boot(statistic = bstat2, data = MAD2, R =5000 , parallel = "multicore", ncpus = 4) # option multicore does not work on windows
plot(br1.2)

bres2 <- list() # create list to store CIs from boot
for (i in 1:ncol(br1.2$t[,1:4])){
  bres2[[i]] <- boot.ci(br1.2, level = .95, type = c("bca"),index = i) # extract CI of the BCA type
}

# loop to extract CI from boot object
media <- br1.2$t0[1:4] # 'media' means mean in Spanish
cdown <- c() 
cup <- c()
Nfactors <- c(0,1,3,5)

for (i in 1:4){
  cdown[i] <- bres2[[i]]$bca[,4]
  cup[i] <- bres2[[i]]$bca[,5]
}

#get everything together in a data frame in order to plot - Attention - run each chunk of code for the appropriate model results signaled by the name of the object
CIs.m1.1 <- cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(media=round(media*100,2))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cdown=round(cdown*100,2))%>%
  mutate(cup=as.numeric(cup))%>%
  mutate(cup=round(cup*100,2))

CIs.m2.1 <- cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(media=round(media*100,2))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cdown=round(cdown*100,2))%>%
  mutate(cup=as.numeric(cup))%>%
  mutate(cup=round(cup*100,2))

CIs.m3.1<-  cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cup=as.numeric(cup))

CIs.m4.1<-  cbind(media, cdown, cup, Nfactors)%>% data.frame()%>%
  mutate(Nfactors=as.factor(Nfactors))%>%
  mutate(media=as.numeric(media))%>%
  mutate(cdown=as.numeric(cdown))%>%
  mutate(cup=as.numeric(cup))

#Plot first part of Fig. 1 - Attention - use the appropiate CIs.xx object for each panel
labelsx_legend <- c("None", "One (TD)", "Three", "Five") # run this code for the second batch of panels

#WHC
pwhc2  <- ggplot() +
  geom_point(CIs.m1.1, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m1.1, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("B")+
  xlab("") +
  ylab("") +
  ggdist::stat_halfeye(data=MAD2,
                       density = "bounded",
                       scale=0.65, # uncomment when using density="histogram"
                       ## custom bandwidth
                       adjust = 1,
                       ## move geom to the right
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=whc, fill=Nfactors),
                       alpha=0.5
  )+
  scale_y_continuous(breaks=c(45,55,65), limits = c(40,67.5))+
  scale_x_discrete(labels= labelsx_legend)+
  scale_fill_manual(values=colorp) +
  geom_hline(yintercept=CIs.m1.1[1,1], linetype="dashed")+
  guides(color="none", fill="none")+
  theme_bw()+
  theme()
pwhc2

#WSA
pwsa2  <- ggplot() +
  geom_point(CIs.m2.1, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m2.1, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("D")+
  xlab("")+
  ylab("")+
  ggdist::stat_halfeye(data=MAD2,
                       density = "bounded",
                       scale=0.65, # uncomment when using density="histogram"
                       ## custom bandwidth
                       adjust = 1,
                       ## move geom to the right
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=wsa, fill=Nfactors),
                       alpha=0.5
  )+
  scale_fill_manual(values=colorp) +
  geom_hline(yintercept=CIs.m2.1[1,1], linetype="dashed")+
  scale_y_continuous(breaks=c(40,60,80), limits = c(28,82))+
  scale_x_discrete(labels= labelsx_legend)+
  guides(color="none", fill="none")+
  theme_bw()+
  theme()
pwsa2

#pH
pph2  <- ggplot() +
  geom_point(CIs.m3.1, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m3.1, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("F")+
  xlab("")+
  ylab("")+
  ggdist::stat_halfeye(data=MAD2,
                       density = "bounded",
                       scale=0.65, # uncomment when using density="histogram"
                       ## custom bandwidth
                       adjust = 1,
                       ## move geom to the right
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=pH, fill=Nfactors),
                       alpha=0.5
  )+
  scale_fill_manual(values=colorp) +
  geom_hline(yintercept=CIs.m3.1[1,1], linetype="dashed")+
  scale_y_continuous(breaks=c(6.0, 6.5,7.0), limits = c(5.65, 7.1))+
  scale_x_discrete(labels= labelsx_legend)+
  guides(color="none", fill="none")+
  theme_bw()+
  theme()
pph2

#B:F ratio
prat2  <- ggplot() +
  geom_point(CIs.m4.1, mapping=aes(Nfactors, media, color=Nfactors))+
  geom_errorbar(CIs.m4.1, mapping=aes(x=Nfactors, y=media, ymin=cdown, ymax=cup, color=Nfactors), width=0.1, linewidth=1, alpha=0.9)+
  scale_color_manual(values=colorp) + 
  ggtitle("H")+
  ylab("")+
  xlab("")+
  ggdist::stat_halfeye(data=MAD2,
                       density = "bounded",
                       scale=0.65, # (un)comment when using density="histogram"
                       ## custom bandwidth
                       adjust = 1,
                       ## move geom to the right
                       justification = -.1,
                       ## remove slab interval
                       .width = 0,
                       point_colour = NA,
                       aes(x=Nfactors, y=ratio, fill=Nfactors),
                       alpha=0.5
  )+
  scale_fill_manual(values=colorp) +
  scale_x_discrete(labels= labelsx_legend)+
  geom_hline(yintercept=CIs.m4.1[1,1], linetype="dashed")+
  guides(color="none", fill="none")+
  theme_bw()+
  theme()
prat2

#emmeans

con.m1.1 <- emmeans(m1.1, specs = trt.vs.ctrl~Nfactors, type="response")
con.m1.1$contrasts %>% summary(infer=TRUE)

con.m2.1 <- emmeans(m2.1, specs = trt.vs.ctrl~Nfactors, type="response")
con.m2.1$contrasts %>% summary(infer=TRUE)

con.m3.1 <- emmeans(m3.1, specs = trt.vs.ctrl~Nfactors, type="response")
con.m3.1$contrasts %>% summary(infer=TRUE)

con.m4.1 <- emmeans(m4.1, specs = trt.vs.ctrl~Nfactors, type="response")
con.m4.1$contrasts %>% summary(infer=TRUE)


#Tables with broom
ts2<- rbind(broom::tidy(m1.1)[,2:6],broom::tidy(m2.1)[,2:6],broom::tidy(m3.1), broom::tidy(m4.1))
write.table(ts3, paste(path,"/ts2.txt",sep=""))

#Combination of pannels in one plot -- Figure 1
fig1 <- ggarrange(pwhc1, pwhc2, pwsa1, pwsa2, pph1, pph2, prat1, prat2,
          ncol = 2, nrow = 4)
annotate_figure(fig1, bottom = text_grob("Number of conditioners in the mix"))
ggsave(filename = "Fig.1_300dpi.png", width = 5.8, height = 8.3, dpi = 300, bg="white")
ggsave(filename = "Fig.1_150dpi.png", width = 5.8, height = 8.3, dpi = 150, bg="white")

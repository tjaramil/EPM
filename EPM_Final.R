#Load Library
library(dplyr)

#Load Data
#Click Session > Set working directory(the file your data is in)>ok

EPM <- read.csv("EPMtest.csv", skip = 4, header = T, stringsAsFactors = F) #EPM is the name used here but replace it with the name of your file, also we skiped the first 4 rows

#Checkout the data
summary(EPM)
glimpse(EPM)
names(EPM)

#Determine WT, HET, KO
EPM$Geno[EPM$Geno== "s"] <- "WT"  #Change "s" to whatever your letter is for WT
EPM$Geno[EPM$Geno== "r"] <- "HET" #Change "r" to whatever your letter is for HET
EPM$Geno[EPM$Geno== "t"] <- "KO"  #Change "t" to whatever your letter is for KO

#Create your dataframes
#View your data frame and find the row #'s specific for open and closed arm data, should be the bottom of dataset.  Replace the numbers with the corresponding numbers   
EPM.Open <- EPM[c(331:384),  ]
EPM.Closed <- EPM[c(386:439),]

#take a peak at your data 
summary(EPM.Open)
summary(EPM.Closed)


#1st set the order you want data arranged to x
x <- c("WT", "HET", "KO")


# Determine stats for Arms
names(EPM.Open)

Open.Arm.Stats <- EPM.Open %>% 
        group_by(Geno) %>% 
        summarise(Count = n(),
                  Open.Bouts.mean = mean(Bouts),
                  Open.Bouts.sd = sd(Bouts),
                  Open.Bouts.sem = Open.Bouts.sd/sqrt(Count),
                  Open.Dur.mean = mean(Duration.Second.),  # R converts () imported as .
                  Open.Dur.sd = sd(Duration.Second.),
                  Open.Dur.sem = Open.Dur.sd/sqrt(Count)) %>% 
        ungroup() %>%
        mutate(Geno = factor(Geno, levels = x)) %>% 
        arrange(Geno)
names(EPM.Closed)                        

Closed.Arm.Stats <- EPM.Closed %>% 
        group_by(Geno) %>% 
        summarise(Count = n(),
                  Closed.Bouts.mean = mean(Bouts),
                  Closed.Bouts.sd = sd(Bouts),
                  Closed.Bouts.sem = Closed.Bouts.sd/sqrt(Count),
                  Closed.Dur.mean = mean(Duration.Second.),
                  Closed.Dur.sd = sd(Duration.Second.),
                  Closed.Dur.sem = Closed.Dur.sd/sqrt(Count)) %>% 
        ungroup() %>% 
        mutate(Geno = factor(Geno, levels = x)) %>% 
        arrange(Geno)


Arm.stats <- left_join(Open.Arm.Stats, Closed.Arm.Stats, by = c("Geno", "Count"))

# If you want to save the data as is
write.csv(Arm.stats, "Arm.Stats.csv") #1st arguement is name of dataframe you created, 2nd is what you want to call it +.csv

##################################################################################################################
# Stats between Genotypes 

# If you have one independent variable (eg.Genotype) do a t-test. If you have 2 (Genotype and Sex) skip to the ANOVA 

## Prepare the data
## Bouts and Duration: Open and Closed

#Open

WT.Open <- subset(EPM.Open, Geno == "WT")
WT.Open.Bouts <- WT.Open$Bouts       
WT.Open.Dur <- WT.Open$Duration.Second.

HET.Open <- subset(EPM.Open, Geno == "HET")
HET.Open.Bouts <- HET.Open$Bouts
HET.Open.Dur <- HET.Open$Duration.Second.

KO.Open <- subset(EPM.Open, Geno == "KO")
KO.Open.Bouts <- KO.Open$Bouts
KO.Open.Dur <- KO.Open$Duration.Second.

#Closed

WT.Closed <- subset(EPM.Closed, Geno == "WT")
WT.Closed.Bouts <- WT.Closed$Bouts       
WT.Closed.Dur <- WT.Closed$Duration.Second.

HET.Closed <- subset(EPM.Closed, Geno == "HET")
HET.Closed.Bouts <- HET.Closed$Bouts
HET.Closed.Dur <- HET.Closed$Duration.Second.

KO.Closed <- subset(EPM.Closed, Geno == "KO")
KO.Closed.Bouts <- KO.Closed$Bouts
KO.Closed.Dur <- KO.Closed$Duration.Second.


###########################################################################################################
## Perform the t-Test 
library(stats)

# H0: mean WT Bouts in the open arm is = mean Bouts in HET or KO
# two-sided ttest
# assume non-equal variance
# mu(Ho hypothesis)
# alternative is 2sided
# confidence interval is 95%
# variance is not equal
# It is not a paired tTest

# Bouts tTest: Open
t.test(WT.Open.Bouts,HET.Open.Bouts, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs HET                
t.test(WT.Open.Bouts,KO.Open.Bouts, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs KO
t.test(HET.Open.Bouts,KO.Open.Bouts, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #HET vs KO

# Bouts tTest: Closed
t.test(WT.Closed.Bouts,HET.Closed.Bouts, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs HET                
t.test(WT.Closed.Bouts,KO.Closed.Bouts, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs KO
t.test(HET.Closed.Bouts,KO.Closed.Bouts, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #HET vs KO

#Duration tTest: Open
t.test(WT.Open.Dur,HET.Open.Dur, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs HET                
t.test(WT.Open.Dur,KO.Open.Dur, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs KO
t.test(HET.Open.Dur,KO.Open.Dur, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #HET vs KO

#Duration tTest: Closed
t.test(WT.Closed.Dur,HET.Closed.Dur, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs HET                
t.test(WT.Closed.Dur,KO.Closed.Dur, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #WT vs KO
t.test(HET.Closed.Dur,KO.Closed.Dur, mu=0, alternative = "two.sided", conf=0.95, var.equal = F, paired = F) #HET vs KO

# Welch's t-test performs better than Student t-test whenever sample sizes and variances are unequal between groups, 
# and gives the same result when sample sizes and variances are equal.

#########################################################################################

## ANOVA

# If you have 2 indep variables (eg. Genotype and Sex) conducted a 2-way ANOVA
#indep variable are Genotype and Sex
#Dep variables is arm frequency
names(EPM.Open)
Anova.1 <- aov(Bouts ~ Sex * Geno, data = EPM.Open) #Change out bouts for other dependent variables
summary(Anova.1)

#Post hoc test 
library(agricolae)

tukey <- TukeyHSD(Anova.1, trt='group')
tukey

scheffe.Geno <- scheffe.test(Anova.1, "Geno", group = F, console = F)
scheffe.Geno 

scheffe.Sex <- scheffe.test(Anova.1, c("Sex", "Geno"), group = F, console = F)
scheffe.Sex 

LSD.test <- LSD.test(Anova.1, "Geno", p.adj = "hochberg", group = F, console = F)
LSD.test


# If you find that sex is significant then explore further why.  If not continue to plotting.  

Arm.stats.Closed.bySex <- EPM.Closed %>% 
        group_by(Geno,Sex) %>% 
        summarise(Closed.bouts = n(),
                  Closed.bouts.mean = mean(Bouts),
                  Closed.bouts.sd = sd(Bouts),
                  Closed.bouts.sem = Closed.bouts.sd/sqrt(Closed.bouts),
                  Closed.dur.mean = mean(Duration.Second.),
                  Closed.dur.sd = sd(Duration.Second.),
                  Closed.dur.sem = Closed.dur.sd/sqrt(Closed.bouts)) %>% 
        ungroup() %>% 
        mutate(Geno = factor(Geno, levels = x)) %>% 
        arrange(Geno,desc(Sex))  #Descending will put Males 1st, remove it and Females will be 1st. 

Arm.stats.Open.bySex <- EPM.Open %>% 
        group_by(Geno,Sex) %>% 
        summarise(Open.bouts = n(),                  
                  Open.bouts.mean = mean(Bouts),
                  Open.bouts.sd = sd(Bouts),
                  Open.bouts.sem = Open.bouts.sd/sqrt(Open.bouts),
                  Open.dur.mean = mean(Duration.Second.),
                  Open.dur.sd = sd(Duration.Second.),
                  Open.dur.sem = Open.dur.sd/sqrt(Open.bouts)) %>% 
        ungroup() %>% 
        mutate(Geno = factor(Geno, levels = x)) %>% 
        arrange(Geno,desc(Sex))  #Descending will put Males 1st, remove it and Females will be 1st. 

#######################################################################################################################

# Plotting 
names(Arm.stats)
str(Arm.stats)
library(ggplot2)

Open.bouts <- ggplot(Arm.stats, aes(Geno, Open.Bouts.mean, fill = Geno))+
        geom_bar(stat = "identity", width = 0.7, colour = "black")+
        scale_fill_manual(values = c("white", "grey", "black"))+
        geom_errorbar(aes(ymin = Open.Bouts.mean, ymax = Open.Bouts.mean + Open.Bouts.sem), width = .4)+
        xlab(NULL)+
        ylab("Number of Bouts in Open Arms")+
        theme(panel.border = element_blank())+
        theme(axis.line = element_line(color = "black", size = 0.5)) +
        scale_y_continuous(expand = c(0,0))+
        theme(aspect.ratio = 9/16)+
        theme(legend.position = "none")
        

Open.dur <- ggplot(Arm.stats, aes(Geno, Open.Dur.mean, fill=Geno))+
        geom_bar(stat = "identity", width = 0.7, colour = "black")+
        scale_fill_manual(values = c("white", "grey", "black"))+
        geom_errorbar(aes(ymin = Open.Dur.mean, ymax = Open.Dur.mean + Open.Dur.sem), width = .4)+
        xlab(NULL)+
        ylab("Duration in Open Arms (sec)")+
        theme(panel.border = element_blank())+
        theme(axis.line = element_line(color = "black", size = 0.5))+
        scale_y_continuous(expand = c(0,0))+
        theme(aspect.ratio = 9/16)+
        theme(legend.position = "none")
        
Closed.bouts <- ggplot(Arm.stats, aes(Geno, Closed.Bouts.mean, fill = Geno))+
        geom_bar(stat = "identity", width = 0.7, colour = "black")+
        scale_fill_manual(values = c("white", "grey", "black"))+
        geom_errorbar(aes(ymin = Closed.Bouts.mean, ymax = Closed.Bouts.mean + Closed.Bouts.sem), width = .4)+
        xlab(NULL)+
        ylab("Number of Bouts in Closed Arms")+
        theme(panel.border = element_blank())+
        theme(axis.line = element_line(color = "black", size = 0.5))+
        scale_y_continuous(expand = c(0,0))+
        theme(aspect.ratio = 9/16)+
        theme(legend.position = "none")

Closed.dur <- ggplot(Arm.stats, aes(Geno, Closed.Dur.mean, fill = Geno))+
        geom_bar(stat = "identity", width = 0.7, colour = "black")+
        scale_fill_manual(values = c("white", "grey", "black"))+
        geom_errorbar(aes(ymin = Closed.Dur.mean, ymax = Closed.Dur.mean + Closed.Dur.sem), width = .4)+
        xlab(NULL)+
        ylab("Duration in Closed Arms (sec)")+
        theme(panel.border = element_blank())+
        theme(axis.line = element_line(color = "black", size = 0.5))+
        scale_y_continuous(expand = c(0,0))+
        theme(aspect.ratio = 9/16)+
        theme(legend.position = "none")

# ALL Plots on one sheet
#load mulitplot function

multiplot(Open.bouts, Open.dur, Closed.bouts, Closed.dur, cols=2)


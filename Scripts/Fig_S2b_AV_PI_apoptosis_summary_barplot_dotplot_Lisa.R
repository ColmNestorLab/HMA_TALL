################################################################################
######################## ANNEXIN/PI APOPTOSIS SAUMMARY ########################
################################################################################

# Lisa Haglund (lisa.haglund@liu.se)
# 2024-07-15


# plot barplot of total cell death + dotplot representing datapoints of replicates.
# Errorbars representing +/- SD of biological triplicates calculated in Exel

# Figures: Extended Fig. 2b


##### load packages #####


library(ggplot2)
library(readxl)

##### set working directory and check package versions #####

# set working directory
setwd("C:/Users/lizah/OneDrive/Documents")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()


##### define color palettes #####

# No palettes used


#### import data #####

### cell cycle analysis results
# mean from biological triplicates
# column names: "Cell_line","drug","concentration","cell_cycle","percent" 


apto2 <- read_excel("Degree project/Flowcytometry/apto2.xlsx")
dotplot <- read_excel("Degree project/dotplot.xlsx")

#reformating data
order.t <- c("Control", "5AZA300", "DEC300", "GSK300", "GSK3000")

colnames(apto2)[1] ="cell.line"
apto2[is.na(apto2)] <- 0

data_alt = pivot_longer(apto2, 
                        cols = c(-cell.line, -Treatment))
data_alt = na.omit(data_alt)

colnames(dotplot)[1] ="cell.line"

long <- left_join(data_alt, apto2[,c(1,2,5)], by = c("cell.line", "Treatment"))

##### plot cell total cell death: barplot #####

### Extended Figure 2b
# bar plot
# mean of biological triplicates
# dotplot showing datapoints for replicates


ggplot(data = long [which(long$name == "total"),], aes(x = as.factor(Treatment), y = value)) +
  geom_bar(stat = "identity", position = position_dodge2(), width = 0.60, fill = "grey70") +
  facet_wrap(~cell.line) +
  scale_x_discrete(limits = order.t) +
  geom_errorbar(data = apto2[,c(1,2,5,6)], aes(x = as.factor(Treatment), y = total, ymin = total - SD.total, ymax = total + SD.total), width = 0.3)+
  geom_dotplot(data = dotplot, aes(x= as.factor(Treatment), y = total, group = interaction(factor(Treatment), factor(cell.line))),alpha = .8, position_dodge(width = .8), binaxis = "y", 
               stackdir = "center", stackratio = 1)+
  theme_bw()




##### statistics #####

#Errorars represent standard diviation of total cell death of 3 biological replicates
# calculated by "STDAV.S" in Exel
#No further statictical analysis in this figure






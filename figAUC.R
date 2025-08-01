### figure for dyslipidemia
library(stringr)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggplot2)


library(ggprism)
library(patchwork)
library(magrittr)

setwd("C:/Users/jfrenc49/OneDrive - St. Jude Children's Research Hospital/Dyslipidemia")


fig2 <- data.frame(group=c("St. Jude Lifetime Cohort","St. Jude Lifetime Cohort",
                           "Childhood Cancer Survivor Study","Childhood Cancer Survivor Study",
                           "Dutch Childhood Cancer Survivor (DCCSS)- LATER Study",
                           "Dutch Childhood Cancer Survivor (DCCSS)- LATER Study"),
                   AUC = c(0.856,0.865,0.753,0.769,0.63,0.67),
                   AUC_l = c(0.810,0.821,0.712,0.731,0.59,0.65),
                   AUC_U = c(0.901,0.908,0.794,0.808,0.65,0.69),
                   model = c("Clinical","Clinical + PRS","Clinical","Clinical + PRS",
                             "Clinical","Clinical + PRS"),
                   P = c(NA,0.048,NA,0.0010,NA,0.000035))
    
fig2$model <- as.factor(fig2$model)
fig2$group = relevel(as.factor(fig2$group), ref = 'St. Jude Lifetime Cohort')
levels(fig2$group) = c('St. Jude Lifetime Cohort', 'Childhood Cancer Survivor Study',
                       "Dutch Childhood Cancer Survivor (DCCSS)- LATER Study")


test = fig2 %>% group_by(model)
test$group1 = rep(c(NA, as.character(test$model[1])),3)
test$group2 = rep(c(NA, as.character(test$model[2])),3)
test$y.position = test$AUC_U + 0.01
#test$y.position[c(4:6, 8:12)] = test$y.position[c(4:6, 8:12)] + c(0.02, 0.026, 0.031, -0.001, 0.013, 0.023, 0.038, 0.04)

library(ggplot2)





# dot plot
test_clean <- test %>% filter(!is.na(P))
  
  
pd = position_dodge(0.8)

out <- 
  ggplot(fig2, aes(x = model,
                 y = AUC, color = model)) +
  facet_wrap(~factor(group, levels = c('St. Jude Lifetime Cohort', 'Childhood Cancer Survivor Study',
                                       "Dutch Childhood Cancer Survivor (DCCSS)- LATER Study")),
             labeller = label_wrap_gen(width=20)) +
  geom_point(position = pd, size = 7) +
  geom_errorbar(aes(ymin=AUC_l, ymax=AUC_U), width=0, size=3, position = pd) +
  scale_color_manual(values = c("Clinical" = "red", "Clinical + PRS" = "blue")) +
  geom_text(aes(label=sprintf("%0.3f", AUC)), position=pd,
            hjust=-0.4, vjust=-0.55, angle=0, size=4, show.legend = F) +
  labs(x = NULL,
       y = "Area Under the ROC Curve (95% CI)") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        strip.text = element_text(size = 14, color= 'black'),
        legend.box = "vertical",
        legend.key.width = unit(1.1, 'cm'), 
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black", size = 14),
        text = element_text(size=14,colour="black"),
        axis.title.y = element_text(vjust = +1.5, size = 15)) + 
  stat_pvalue_manual(test_clean,label = "P", label.size = 4, parse = TRUE)
ggsave("CAD_AUC.png", plot = out, width = 10, height = 6, dpi = 300)

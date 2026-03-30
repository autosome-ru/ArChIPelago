
# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggrepel)
library(gapminder)
library(ggExtra)
library(extrafont)
font_import()


theme_set(theme_get() + theme(text = element_text(family = 'Arial')))


Figure 1
##### CHS ####

Model_key_list=c("RandomForestClassifier", "LogisticRegression", "XGBClassifier", "BaggingClassifier_XGBClassifier", "BaggingClassifier_LogisticRegression")




my_pal <- function(range = c(1, 6)) {
  force(range)
  function(x) scales::rescale(x, to = range, from = c(0, 1))
}

# 
# data <- read.csv("/Users/pavel/Desktop/Archipelago_trees/HUMAN_MOUSE_total_100k.csv", sep="\t")
# 
# data1 = data %>% 
#   dplyr::filter(Model=="RandomForestClassifier" & PWM=="mono+di") %>% 
#   select(c("TF_name", 
#            "roc_auc_test_H_PWM", "roc_auc_test_H",
#            "pr_auc_test_H_PWM", "pr_auc_test_H", "Count", "Seq_count"))
# 
# data = data1 %>% 
#   mutate(ROC_delta_H = roc_auc_test_H-roc_auc_test_H_PWM, PR_delta_H = pr_auc_test_H-pr_auc_test_H_PWM) %>% 
#   group_by(TF_name) %>% 
#   summarise(ROC_delta_H=mean(ROC_delta_H), PR_delta_H=mean(PR_delta_H), Count_PWM=sum(unique(Count)), Seq_count_min=min(Seq_count)) %>% 
#   ungroup()



data <- read.csv("./HUMAN_MOUSE_total_100k.csv", sep="\t")
data1 <- data %>%
  dplyr::filter(Model == "RandomForestClassifier") %>%
  select(
    TF_name, PWM,
    roc_auc_test_H_PWM_mono, roc_auc_test_H,
    pr_auc_test_H_PWM_mono, pr_auc_test_H,
    Count, Seq_count
  )
data2 <- data1 %>%
  group_by(TF_name) %>%
  mutate(Count_PWM = sum(unique(Count))) %>%
  ungroup()
data3 <- data2 %>%
  filter(PWM=="mono+di")
data <- data3 %>%
  mutate(
    ROC_delta_H = roc_auc_test_H - roc_auc_test_H_PWM_mono,  # compare mono+di model vs mono baseline
    PR_delta_H  = pr_auc_test_H  - pr_auc_test_H_PWM_mono    # compare mono+di model vs mono baseline
  ) %>% 
  group_by(TF_name) %>% 
  summarise(
    ROC_delta_H = mean(ROC_delta_H, na.rm = TRUE),
    PR_delta_H  = mean(PR_delta_H,  na.rm = TRUE),
    Count_PWM  = first(Count_PWM),
    Seq_count_min = min(Seq_count, na.rm = TRUE)
  ) %>% 
  ungroup()


c = data %>%
  arrange(desc(Seq_count_min)) %>%
  mutate(TF_name = factor(TF_name)) %>%
  ggplot(aes(x=PR_delta_H, y=ROC_delta_H, size=Seq_count_min, fill=Count_PWM)) + # 
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(range = c(1, 15), breaks= c(500, 1000, 5000, 10000, 50000),
             #labels = c(5*10^2, 10^3, 5*10^3, 10^4, 5*10^4),
             name="Size of \nthe positive \nset") +
  scale_fill_viridis(discrete=F, option="C", name="Number of PWMs") +
  theme_minimal(base_size = 21) +
  theme(legend.position="right",
        legend.box="vertical", 
        legend.margin=margin())+
  ylab("\u0394auROC") +
  xlab("\u0394auPRC") + 
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  geom_text_repel(data = subset(data, (ROC_delta_H < 0)|(ROC_delta_H>0.05)|(PR_delta_H<0)|(PR_delta_H>0.2)),
                  aes(label = TF_name), 
                  point.padding = 1,
                  min.segment.length = 3,
                  max.time = 1, max.iter = 1e5,
                  #box.padding = 0.3, 
                  segment.curvature = -0.1,
                  #segment.ncp = 3,
                  segment.angle = 20,
                  size=5)+
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2),
         shape = "none")



data <- read.csv("./HUMAN_MOUSE_total_100k.csv", sep="\t")

data = data %>%
  dplyr::filter(Model=="RandomForestClassifier") %>%
  select(c("TF_name",
           "roc_auc_test_H_PWM_mono", "roc_auc_test_H",
           "pr_auc_test_H_PWM_mono", "pr_auc_test_H", "Count", "Seq_count", "PWM"))

a = data %>%
  arrange(desc(Seq_count)) %>%
  ggplot(aes(x=roc_auc_test_H_PWM_mono, y=roc_auc_test_H, color=PWM)) + # 
  geom_point(alpha=0.5, shape=20, size=3)+
  theme_minimal(base_size = 21) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, vjust = 0.5))+
  ylab("auROC") +
  xlab("auROC best PWM") + 
  expand_limits(x=0, y=0) +
  geom_abline(intercept = 0, linetype="dotted") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

a = ggMarginal(a, type = "histogram", 
               #margins = "x",
               #color = "gray",
               fill = "white", 
               alpha=0.5,
               #bins = 60, 
               size=8, groupColour = TRUE)

b = data %>%
  arrange(desc(Seq_count))  %>%
  ggplot(aes(x=pr_auc_test_H_PWM_mono, y=pr_auc_test_H, color=PWM)) + # 
  geom_point(alpha=0.5, shape=20, size=3)+
  theme_minimal(base_size = 21) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, vjust = 0.5))+
  ylab("auPRC") +
  xlab("auPRC best PWM") + 
  expand_limits(x=0, y=0) +
  geom_abline(intercept = 0, linetype="dotted") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

b = ggMarginal(b, type = "histogram", 
               #margins = "x",
               #color = "gray",
               fill = "white", 
               alpha=0.5,
               size=8, groupColour = TRUE)

plot = ggarrange(
  
  ggarrange(a, b, nrow = 2, ncol = 1, labels = c("A", "B"),
            font.label=list(size=25), vjust = -0.3 
            #heights = c(1, 1.5)
  ),
  
  ggarrange(c, 
            labels = c("C"),
            ncol = 1, nrow = 1, 
            font.label=list(size=25), vjust = -0.3),  
  widths = c(1, 2.5), ncol = 2, nrow = 1
)

image = annotate_figure(plot, top = text_grob("", 
                                      color = "black", face = "bold", size = 25))

library(Cairo)

Cairo(file='Figure_2_20012026_mono_di_HUMAN_HUMAN_RandomForestClassifier.pdf', type="pdf", width=310, height=185, units="mm")
image
dev.off()


data4 <- data3 %>%
  mutate(
    ROC_delta_H = roc_auc_test_H - roc_auc_test_H_PWM_mono,  # compare vs mono baseline
    PR_delta_H  = pr_auc_test_H  - pr_auc_test_H_PWM_mono    # compare vs mono baseline
  )
median_comparison <- data4 %>%
  summarise(
    median_roc_H      = median(roc_auc_test_H, na.rm = TRUE),
    median_roc_H_PWM  = median(roc_auc_test_H_PWM_mono, na.rm = TRUE),
    median_pr_H       = median(pr_auc_test_H, na.rm = TRUE),
    median_pr_H_PWM   = median(pr_auc_test_H_PWM_mono, na.rm = TRUE)
  )
median_gain <- data4 %>% 
  summarise(
    median_ROC_delta_H = median(ROC_delta_H, na.rm = TRUE),
    median_PR_delta_H  = median(PR_delta_H,  na.rm = TRUE)
  )
median_summary <- bind_cols(median_comparison, median_gain)





PhyloObj <- pseq
index = "Observed"
strata = "LP"
y_label = "Observed species"
add_legend = FALSE
gg_title = NULL
m_Kolors = Kolors
  #group.colors <- c(ctrl = "dodgerblue3", Gem = "firebrick2")
condition_names <- sample_data(PhyloObj)[[strata]] %>% table() %>% names()
group.colors <- m_Kolors[1:length(condition_names)]
names(group.colors) <- condition_names
#my.labels <- c("Ctrl", "Gem")
my.labels <- condition_names
## calculate alpha diversity

tmp2 <- PhyloObj %>% estimate_richness()
##------------------------------
new_names <- rownames(tmp2) %>% 
  gsub("^X", "", .) %>% 
  gsub("\\.", "-", .)
rownames(tmp2) <- new_names
rich_meta <- merge(PhyloObj %>% sample_data(), tmp2, by = "row.names")
##-- only using statistics test
ob <- rich_meta %>% 
  t_test(as.formula(paste0(index, " ~ ", strata))) %>% 
  adjust_pvalue(method = "BH") %>%  
  add_significance("p.adj") %>% 
  add_xy_position() %>% 
  mutate(p_show = case_when(
    p.adj < 0.001 ~ "<0.001***",
    p.adj < 0.01 ~ sprintf("%.3f**", round(p.adj,3)),
    p.adj < 0.05 ~ sprintf("%.3f*", round(p.adj,3)),
    TRUE ~ sprintf("%.3f(ns)", round(p.adj,3))
  ))

# Logistic regression
model <- glm(
  Observed ~ LP,
  data = rich_meta
)

summary(model)
#mutate(p_show = sprintf("%.3f(%s)", round(p.adj,3), p.adj.signif))

p1 <- ggplot (rich_meta, aes(x = !!sym(strata), y = !!sym(index), fill=!!sym(strata)))+ 
  geom_violin() + #geom_boxplot()
  geom_boxplot(width=0.3) +
  geom_point (position=position_jitterdodge( jitter.width = 0.05), size = 0.3, col = alpha('grey', 0.2))+
  theme_minimal() + # theme_gray(), theme_bw()
  scale_x_discrete(labels= my.labels)+
  #xlab(strata)+
  ylab(y_label)+
  # facet_grid(.~tp)+
  #scale_x_discrete(labels= my.labels)+
  scale_fill_manual(values=group.colors, labels = my.labels)+
  ggtitle(gg_title) +
  #theme(panel.grid =  element_blank()) + 
  #ggtitle("KPC tumor vs. Healthy pancreas - Observed species") +
  theme(axis.text.y = element_text (size=12),
        axis.title = element_text(size=12, face="bold"))

if(!add_legend){    
  p1 <- p1 + theme(legend.position = "none",
                   # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                   legend.key.size = unit(4,"mm"),
                   axis.text.x = element_blank(),
                   plot.title = element_text(size = 12))
  
} else {
  p1 <- p1 + theme(legend.text = element_text(size = 12),
                   #legend.title = element_blank(),
                   # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                   legend.key.size = unit(4,"mm"),
                   axis.text.x = element_blank(),
                   plot.title = element_text(size = 12))
}
p1 <- p1 + stat_pvalue_manual(ob, label = "p_show", inherit.aes = FALSE, tip.length = 0.01, step.increase = 0.1)
return(p1)

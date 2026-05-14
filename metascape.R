# GO和KEGG富集分析双向柱状图

rm(list = ls())
library(stringr)
options(stringsAsFactors = F)

dt = readxl::read_xlsx("KEGG.xlsx")
dt <- dt[c(1:15),]
#View(dt)
# 创建示例数据
colnames(dt)
dt <- dt[,c("Category","Description","Log(q-value)","Enrichment")]
colnames(dt) = c("Classification", "Pathways", "logP.value", "Enrichment")
dt$Classification[str_detect(dt$Classification,"Bio")] <- "GO_BP"
dt$Classification[str_detect(dt$Classification,"Mol")] <- "GO_MF"
dt$Classification[str_detect(dt$Classification,"Cel")] <- "GO_CC"
dt$Classification[str_detect(dt$Classification,"KEGG")] <- "KEGG"


#dt$Classification[str_detect(dt$Classification,"Hallmark")] <- "Hallmark_Gene"

dt = dt[order(dt$Classification,dt$logP.value),]
dt = dt[!(duplicated(dt$Pathways)),]
dt$Pathways = factor(dt$Pathways, levels =  dt$Pathways)
library(ggplot2)
library(ggpubr)
# 左图'#ee4c58','#56c1ab','#80c5d8','#437eb8',

p1 <- ggplot(dt, aes(x = Pathways,
                     y = logP.value,
                     fill = Classification)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#c4a2c5')) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  ylab('logP.value') +
  xlab('') +
  #labs(title = "Mature neurons") +  # 添加标题
  theme_pubr() +
  theme(
    # 统一字体设置
    text = element_text( size = 16),
    axis.text = element_text(size = 14.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    
    # 标题样式设置（可选）
    plot.title = element_text(
      size = 12,           # 标题字体稍大
      face = "bold",       # 加粗
      hjust = 0.5,         # 水平居中
      margin = margin(b = 10)  # 标题与图表的间距
    ),
    
    # Nature风格的其他调整
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    
    # 调整图例和边距
    legend.position = "right",
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# 显示图形
print(p1)
# 右图
p2 <- ggplot(dt, aes( x = Pathways,
                      y = Enrichment),fill = '#595758')+
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(stat = 'identity')+
  coord_flip()+
  ylab('Enrichment Score')+
  xlab('') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p2
library(patchwork)
p1+ p2 + plot_layout(guides = 'collect')+ plot_layout(widths = c(4, 3))
0# 保存
ggsave('./NY071_Result/Day7_KEGG&GO.pdf',width = 8,height = 8)
ggsave('./NY071_Result/Day7_KEGG&GO.png',width = 8,height = 8)

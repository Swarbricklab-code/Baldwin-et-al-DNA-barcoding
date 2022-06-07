# script for the identificatoin of common DEGs between groups

# Louise Baldwin
# 20220511

#input: weighted list of DEGs in csv format
#output: csv, venn diagrams

# Session set up -------------------------------------------

library(VennDiagram)
library(tidyverse)

# ----------------------------------------------------------

IE1<-read.csv("../logFCs_FC_pvalue_IE1vsBulkBatchcorr.csv", header = TRUE)
IE2<-read.csv("../logFCs_FC_pvalue_IE2vsBulkBatchcorr.csv", header = TRUE)
NT1<-read.csv("../logFCs_FC_pvalue_NT1vsBulkBatchcorr.csv", header = TRUE)
NT2<-read.csv("../logFCs_FC_pvalue_NT2vsBulkBatchcorr.csv", header = TRUE)

IE1_sig <- (IE1 %>% filter(FDR<0.05))
IE2_sig <- (IE2 %>% filter(FDR<0.05))
NT1_sig <- (NT1 %>% filter(FDR<0.05))
NT2_sig <- (NT2 %>% filter(FDR<0.05))

IE1_up <- (IE1_sig %>% filter(logFC>0))
IE1_down <- (IE1_sig %>% filter(logFC<0))
IE2_up <- (IE2_sig %>% filter(logFC>0))
IE2_down <- (IE2_sig %>% filter(logFC<0))
NT1_up <- (NT1_sig %>% filter(logFC>0))
NT1_down <- (NT1_sig %>% filter(logFC<0))
NT2_up <- (NT2_sig %>% filter(logFC>0))
NT2_down <- (NT2_sig %>% filter(logFC<0))

IE1_up_noNT <- IE1_up[!(IE1_up$Column1 %in% NT1_up$X),]
IE1_up_noNT <- IE1_up_noNT[!(IE1_up_noNT$Column1 %in% NT2_up$X),]
IE1_down_noNT <- IE1_down[!(IE1_down$Column1 %in% NT1_down$X),]
IE1_down_noNT <- IE1_down_noNT[!(IE1_down_noNT$Column1 %in% NT2_down$X),]


IE2_down_noNT <- IE2_up[!(IE2_up$Column1 %in% NT1_up$X),]
IE2_up_noNT <- IE2_up_noNT[!(IE2_up_noNT$Column1 %in% NT2_up$X),]
IE2_down_noNT <- IE2_down[!(IE2_down$Column1 %in% NT1_down$X),]
IE2_down_noNT <- IE2_down_noNT[!(IE2_down_noNT$Column1 %in% NT2_down$X),]

write.csv(IE1_up_noNT, "./IE1up_not_NT.csv", row.names = FALSE)
write.csv(IE1_down_noNT,"./IE1down_not_NT.csv", row.names = FALSE)
write.csv(IE2_up_noNT,"./IE2up_not_NT.csv", row.names = FALSE)
write.csv(IE2_down_noNT,"./IE2down_not_NT.csv", row.names = FALSE)

#IE1 and IE2 unique genes only
IE1up_unique <- IE1_up_noNT[!(IE1_up_noNT$Column1 %in% IE2_up_noNT$Column1),]
IE1down_unique <- IE1_down_noNT[!(IE1_down_noNT$Column1 %in% IE2_down_noNT$Column1),]

IE2up_unique <- IE2_up_noNT[!(IE2_up_noNT$Column1 %in% IE1_up_noNT$Column1),]
IE2down_unique <- IE2_down_noNT[!(IE2_down_noNT$Column1 %in% IE1_down_noNT$Column1),]

write.csv(IE1up_unique, "./IE1_up_unique.csv", row.names = FALSE)librar
write.csv(IE1down_unique, "./IE1_down_unique.csv", row.names = FALSE)
write.csv(IE2up_unique, "./IE2_up_unique.csv", row.names = FALSE)
write.csv(IE2down_unique, "./IE2_down_unique.csv", row.names = FALSE)


#find genes common to IE1 and IE2, but not NT1 and NT2

IE1_up <- read.csv("./IE1up_not_NT.csv")
IE2_up <- read.csv("./IE2up_not_NT.csv")

common <- merge(IE1_up, IE2_up, by.x = "Column1", by.y = "Column1")
write.csv(common, "./IE1IE2_up_common.csv", row.names = FALSE)

IE1_down <- read.csv("./IE1down_not_NT.csv")
IE2_down <- read.csv("./IE2down_not_NT.csv")
common_down <- merge(IE1_down, IE2_down, by.x = "Column1", by.y = "Column1")
write.csv(common_down, "./IE1IE2_down_common.csv", row.names = FALSE)



##redo venns
pink <- rgb (244, 185, 185, max = 255, alpha = 200, names = "pink")
yellow <- rgb(249, 235, 174, max = 255, alpha = 200, names = "yellow")
blue <- rgb (174, 206, 249, max = 255, alpha = 200, names = "blue")
purple <- rgb(174, 183, 246, max = 255, alpha = 200, names = "purple")

two.colours <- c(blue, purple)
three.colours <- c(blue, purple, yellow)
colours <- c(pink, yellow, blue, purple)

IE1_up_list <- as.vector(IE1_up$Column1)
IE1_down_list <- as.vector(IE1_down$Column1)
IE2_up_list <- as.vector(IE2_up$Column1)
IE2_down_list <- as.vector(IE2_down$Column1)
NT1_up_list <- as.vector(NT1_up$X)
NT1_down_list <- as.vector(NT1_down$X)
NT2_up_list<- as.vector(NT2_up$X)
NT2_down_list <- as.vector(NT2_down$X)

#gene overlap between IE1up and IE2up

v<-venn.diagram(
  x=list(IE1_up_list, IE2_up_list),
  category.names = c("IE1_up", "IE2_up"),
  filename ="IE1IE2_up_venn.png", 
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  lty = "blank",
  fill = two.colours,
  fontface = "bold",
  cat.default.pos = "outer"
)

v<-venn.diagram(
  x=list(IE1_down_list, IE2_down_list),
  category.names = c("IE1_down", "IE2_down"),
  filename ="IE1IE2_down_venn.png", 
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  lty = "blank",
  fill = two.colours,
  fontface = "bold",
  cat.default.pos = "outer"
)

q <- venn.diagram(
  x = list(IE1_up_list, IE2_up_list, NT1_up_list, NT2_up_list),
  category.names = c("IE1_up","IE2_up", "NT1_up", "NT2_up"),
  filename = "quad_up_venn.png",
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  fill = colours, 
  lty = "blank",
  fontface = "bold",
  cat.default.pos = "outer"
)

q <- venn.diagram(
  x = list(IE1_down_list, IE2_down_list, NT1_down_list, NT2_down_list),
  category.names = c("IE1_down","IE2_down", "NT1_down", "NT2_down"),
  filename = "quad_down_venn.png",
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  fill = colours, 
  lty = "blank",
  fontface = "bold",
  cat.default.pos = "outer"
)

q <- venn.diagram(
  x = list(IE1_down_list, NT1_down_list, NT2_down_list),
  category.names = c("IE1_down", "NT1_down", "NT2_down"),
  filename = "IE1_NT1_NT2_down_venn.png",
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  fill = three.colours, 
  lty = "blank",
  fontface = "bold",
  cat.default.pos = "outer"
)

q <- venn.diagram(
  x = list(IE1_up_list, NT1_up_list, NT2_up_list),
  category.names = c("IE1_up", "NT1_up", "NT2_up"),
  filename = "IE1_NT1_NT2_up_venn.png",
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  fill = three.colours, 
  lty = "blank",
  fontface = "bold",
  cat.default.pos = "outer"
)



q <- venn.diagram(
  x = list(IE2_down_list, NT1_down_list, NT2_down_list),
  category.names = c("IE2_down", "NT1_down", "NT2_down"),
  filename = "IE2_NT1_NT2_down_venn.png",
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  fill = three.colours, 
  lty = "blank",
  fontface = "bold",
  cat.default.pos = "outer"
)

q <- venn.diagram(
  x = list(IE2_up_list, NT1_up_list, NT2_up_list),
  category.names = c("IE2_up", "NT1_up", "NT2_up"),
  filename = "IE2_NT1_NT2_up_venn.png",
  height = 6000, 
  width = 6000,
  output = TRUE,
  imagetype = "png",
  resolution = 1200,
  compression = "lzw",
  fill = three.colours, 
  lty = "blank",
  fontface = "bold",
  cat.default.pos = "outer"
)




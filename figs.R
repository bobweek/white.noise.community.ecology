#
#
# creating figures for community dynamics.
#
#

require(mvtnorm)
require(grid)
require(gridExtra)
require(ggplot2)
require(ggthemes)
require(ggpubr)

#
# community dynamics for strong competition
#

df.mc <- read.csv("/home/bb/Gits/white.noise.community.ecology/sample_path_mc.csv")

endT <- df.mc$time[2]

spp = unique(df.mc$spp)
extinct.mc = 0
for(i in 1:length(spp)){
  if( df.mc$ext[which(df.mc$spp==spp[i])[1]]=="Extinct" ) extinct.mc = extinct.mc+1
}

px.mc = ggplot(df.mc,aes(x=time,y=x,group=interaction(spp,ext),color=ext)) + geom_line(size=0.2) + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  scale_color_manual(values=c('black','red')) + xlab("") + ylab("Mean Trait") + labs(color = "") + theme_minimal() + theme(axis.title = element_text(size=15), legend.position = "none") #+ xlim(c(0,endT))

pG.mc = ggplot(df.mc,aes(x=time,y=G,group=interaction(spp,ext),color=ext)) + geom_line(size=0.2) + scale_y_continuous(trans='log10') + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  scale_color_manual(values=c('black','red')) + xlab("") + ylab("Trait Variance") + labs(color = "") + theme_minimal() + theme(axis.title = element_text(size=15), legend.position = "none") #+ xlim(c(0,endT))

pN.mc = ggplot(df.mc,aes(x=time,y=N,group=interaction(spp,ext),color=ext)) + geom_line(size=0.2) + scale_y_continuous(trans='log10') +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  scale_color_manual(values=c('black','red')) + xlab("Time") + ylab("Abundance") + labs(color = "") + theme_minimal() + theme(axis.title = element_text(size=15), legend.position = c(0.9,0.5)) #+ xlim(c(0,endT))

p.mc = grid.arrange(px.mc,pG.mc,pN.mc,nrow=3,top = textGrob("Moderate Competition",gp=gpar(fontsize=20,font=1)))

ggsave("/home/bb/Gits/white.noise.community.ecology/mcp.png",p.mc,width=5,height=8)

#
# community dynamics for weak competition
#

df.wc <- read.csv("/home/bb/Gits/white.noise.community.ecology/sample_path_wc.csv")

spp = unique(df.wc$spp)
extinct.wc = 0
for(i in 1:length(spp)){
  if( df.wc$ext[which(df.wc$spp==spp[i])[1]]=="Extinct" ) extinct.wc = extinct.wc+1
}

endT <- df.wc$time[2]

px.wc = ggplot(df.wc,aes(x=time,y=x,group=interaction(spp,ext),color=ext)) + geom_line(size=0.2) + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  scale_color_manual(values=c('black','red')) + xlab("") + ylab("") + labs(color = "") + theme_minimal() + theme(legend.position = "none") #+ xlim(c(0,endT))

pG.wc = ggplot(df.wc,aes(x=time,y=G,group=interaction(spp,ext),color=ext)) + geom_line(size=0.2) + scale_y_continuous(trans='log10') + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  scale_color_manual(values=c('black','red')) + xlab("") + ylab("") + labs(color = "") + theme_minimal() + theme(legend.position = "none") #+ xlim(c(0,endT))

pN.wc = ggplot(df.wc,aes(x=time,y=N,group=interaction(spp,ext),color=ext)) + geom_line(size=0.2) + scale_y_continuous(trans='log10') +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  scale_color_manual(values=c('black','red')) + xlab("Time") + ylab("") + labs(color = "") + theme_minimal() + theme(axis.title = element_text(size=15), legend.position = c(0.9,0.5)) #+ xlim(c(0,endT))

p.wc = grid.arrange(px.wc,pG.wc,pN.wc,nrow=3,top = textGrob("Weak Competition",gp=gpar(fontsize=20,font=1)))

ggsave("/home/bb/Gits/white.noise.community.ecology/wcp.png",p.wc,width=5,height=8)

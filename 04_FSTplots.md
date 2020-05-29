# Plot Space vs time

MOD core vs exp
```
p2 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761674.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("Mod Core vs Mus LR761674.1 - PX02") +  theme(axis.title.x = element_blank())
p3 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761666.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761666.1 - PX03") +  theme(axis.title.x = element_blank())
p4 <-ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761647.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761647.1 - PX04") +  theme(axis.title.x = element_blank())
p5 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761648.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761648.1 - PX05") +  theme(axis.title.x = element_blank())
p6 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761657.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761657.1 - PX06") +  theme(axis.title.x = element_blank())
p7 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761668.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761668.1 + PX07") +  theme(axis.title.x = element_blank())
p8 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761662.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761662.1 - PX08") +  theme(axis.title.x = element_blank())
p9 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761659.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761659.1 - PX09") +  theme(axis.title.x = element_blank())
p10 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761649.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761649.1 - PX10") +  theme(axis.title.x = element_blank())
p11 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761675.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
p12 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761655.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761655.1 - PX12") +  theme(axis.title.x = element_blank())
p13 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761653.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761653.1 - PX13") +  theme(axis.title.x = element_blank())
p14 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761669.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761669.1 - PX14") +  theme(axis.title.x = element_blank())
p15 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761651.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761651.1 - PX15") +  theme(axis.title.x = element_blank())
p16 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761656.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761656.1 - PX30/16") +  theme(axis.title.x = element_blank())
p17 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761660.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761660.1 - PX17") +  theme(axis.title.x = element_blank())
p18 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761663.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761663.1 - PX63") +  theme(axis.title.x = element_blank())
p19 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761667.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761667.1 - PX19") +  theme(axis.title.x = element_blank())
p20 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761654.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761654.1 - PX20/28") +  theme(axis.title.x = element_blank())
p21 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761664.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761664.1 - PX21") +  theme(axis.title.x = element_blank())
p22 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761652.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761652.1 - PX22") +  theme(axis.title.x = element_blank())
p23 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761658.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761658.1 - PX23") +  theme(axis.title.x = element_blank())
p24 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761673.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761673.1 - PX24") +  theme(axis.title.x = element_blank())
p25 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761665.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761665.1 - PX25") +  theme(axis.title.x = element_blank())
p26 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761671.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761671.1 - PX05/26") +  theme(axis.title.x = element_blank())
p27 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761670.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761670.1 - PX27") +  theme(axis.title.x = element_blank())

p29 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761661.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761661.1 - PX60") +  theme(axis.title.x = element_blank())
p31 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761672.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761672.1 - PX31") +  theme(axis.title.x = element_blank())

p32 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761650.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761650.1 - PX1Z") +  theme(axis.title.x = element_blank())


```




Mod vs Mus
```
p2mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761674.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("Mod Core vs Mus LR761674.1 - PX02") +  theme(axis.title.x = element_blank())
p3mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761666.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761666.1 - PX03") +  theme(axis.title.x = element_blank())
p4mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761647.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761647.1 - PX04") +  theme(axis.title.x = element_blank())
p5mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761648.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761648.1 - PX05") +  theme(axis.title.x = element_blank())
p6mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761657.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761657.1 - PX06") +  theme(axis.title.x = element_blank())
p7mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761668.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761668.1 + PX07") +  theme(axis.title.x = element_blank())
p8mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761662.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761662.1 - PX08") +  theme(axis.title.x = element_blank())
p9mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761659.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761659.1 - PX09") +  theme(axis.title.x = element_blank())
p10mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761649.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761649.1 - PX10") +  theme(axis.title.x = element_blank())
p11mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761675.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
p12mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761655.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761655.1 - PX12") +  theme(axis.title.x = element_blank())
p13mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761653.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761653.1 - PX13") +  theme(axis.title.x = element_blank())
p14mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761669.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761669.1 - PX14") +  theme(axis.title.x = element_blank())
p15mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761651.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761651.1 - PX15") +  theme(axis.title.x = element_blank())
p16mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761656.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761656.1 - PX30/16") +  theme(axis.title.x = element_blank())
p17mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761660.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761660.1 - PX17") +  theme(axis.title.x = element_blank())
p18mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761663.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761663.1 - PX63") +  theme(axis.title.x = element_blank())
p19mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761667.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761667.1 - PX19") +  theme(axis.title.x = element_blank())
p20mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761654.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761654.1 - PX20/28") +  theme(axis.title.x = element_blank())
p21mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761664.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761664.1 - PX21") +  theme(axis.title.x = element_blank())
p22mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761652.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761652.1 - PX22") +  theme(axis.title.x = element_blank())
p23mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761658.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761658.1 - PX23") +  theme(axis.title.x = element_blank())
p24mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761673.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761673.1 - PX24") +  theme(axis.title.x = element_blank())
p25mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761665.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761665.1 - PX25") +  theme(axis.title.x = element_blank())
p26mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761671.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761671.1 - PX05/26") +  theme(axis.title.x = element_blank())
p27mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761670.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761670.1 - PX27") +  theme(axis.title.x = element_blank())

p29mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761661.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761661.1 - PX60") +  theme(axis.title.x = element_blank())
p31mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761672.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761672.1 - PX31") +  theme(axis.title.x = element_blank())

p32mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761650.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761650.1 - PX1Z") +  theme(axis.title.x = element_blank())
```


Write all plots to pdf and combine using an online pdf combiner. 
```
pdf("FST1.pdf")

multiplot(p2,p2mus, p3, p3mus, cols=1)

dev.off()



pdf("FST2.pdf")

multiplot(p4,p4mus, p5, p5mus, cols=1)

dev.off()



pdf("FST3.pdf")

multiplot(p6,p6mus, p7, p7mus, cols=1)

dev.off()



pdf("FST4.pdf")

multiplot(p8,p8mus, p9, p9mus, cols=1)

dev.off()



pdf("FST5.pdf")

multiplot(p10,p10mus, p11, p11mus, cols=1)

dev.off()



pdf("FST6.pdf")

multiplot(p12,p12mus, p13, p13mus, cols=1)

dev.off()



pdf("FST7.pdf")

multiplot(p14,p14mus, p15, p15mus, cols=1)

dev.off()


pdf("FST8.pdf")

multiplot(p17,p17mus, p18, p18mus, cols=1)

dev.off()



pdf("FST9.pdf")

multiplot(p19,p19mus, p20, p20mus, cols=1)

dev.off()



pdf("FST10.pdf")

multiplot(p21,p21mus, p22, p22mus, cols=1)

dev.off()



pdf("FST11.pdf")

multiplot(p23,p23mus, p24, p24mus, cols=1)

dev.off()



pdf("FST12.pdf")

multiplot(p25,p25mus, p26, p26mus, cols=1)

dev.off()



pdf("FST13.pdf")

multiplot(p27,p27mus, p29, p29mus, cols=1)

dev.off()



pdf("FST14.pdf")

multiplot(p30,p30mus, p31, p31mus, cols=1)

dev.off()

pdf("FST15.pdf")
multiplot(p32,p32mus, cols=1)

```
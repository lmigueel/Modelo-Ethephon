########## ETHEPHON ##########

vtx12_eth <- c(0,0,0,0.366666667,0.85,1.672222222,	2.616666667,	6.716666667,
               7.605555556,	9.288888889,	10.71666667,	13.71666667,	14.72777778,
               15.86111111,	16.97777778,	19.37777778,	21.14444444,	25.16666667,
               25.83333333)

vtx12_ct <- c(0,0.027777778,	0.061111111,	0.622222222,	1.044444444,	2.105555556,
              3.533333333,	8.083333333,	8.905555556,	10.63333333,	12.62222222,	
              16.72222222,	18.17222222,	19.15,	20.78888889,	24.74444444,	27.25555556,
              32.35,	34.12777778)

ctc4_ct <- c(0,	0,	0,	0,	0.05,	0.13125,	0.25625,	0.36875,	0.493333333,	0.625,
             0.7875,	1.46875,	1.76875,	2.04375,	2.33125,	3.24375,	3.9125,	6.2125,
             7.4375)

ctc4_eth <- c(0,	0,0,	0.129411765,	0.217647059,	0.229411765,	0.417647059,	0.8,
              0.894117647,	1.3, 1.617647059,	3.476470588,	3.452941176,	4.252941176,
              5.317647059,	7.847058824,	9.464705882,	14.09411765,	16.45294118)

dias <- c(1,2,3,4,5,6,7,8,12,13,14,15,18,20,22,24,26,28,30)

windows()

setwd("../Ethefon/")
pdf("dots.pdf")
par(mfrow=c(2,2))
plot(dias,vtx12_ct,las=1,pch=16,main="VTX12 - Controle")
plot(dias,vtx12_eth,las=1,pch=16,main="VTX12 - Ethephon")
plot(dias,ctc4_ct,las=1,pch=16,main="CTC4 - Controle")
plot(dias,ctc4_eth,las=1,pch=16,main="CTC4 - Ethephon")
dev.off()

#vtx12_0
logisticModel <- nls(vtx12_ct~M/(1+exp(Po+k*dias)), start=list(Po=0, k=-0.21, M=17))
summary(logisticModel)
coef(logisticModel)

pdf("curvas.pdf")
par(mfrow=c(2,2))

x1 <- seq(min(dias), max(dias), length=100)
y1 <- predict(logisticModel, list(dias=x1))
plot(dias,vtx12_ct,las=1,pch=16,main="VTX12 - Controle")
points(x1, y1, type='l', col='blue')

#vtx12_0
logisticModel1 <- nls(vtx12_eth~M/(1+exp(Po+k*dias)), 
                      start=list(Po=0, k=-0.21, M=17))
summary(logisticModel)
coef(logisticModel)

x2 <- seq(min(dias), max(dias), length=100)
y2 <- predict(logisticModel1, list(dias=x2))
plot(dias,vtx12_eth,las=1,pch=16,main="VTX12 - Ethephon")
points(x2, y2, type='l', col='blue')

#ctc4_0
logisticModel2 <- nls(ctc4_ct~M/(1+exp(Po+k*dias)), 
                      start=list(Po=5, k=-0.5, M=10),
                      control = list(maxiter = 100,
                                     minFactor = 1/4096),trace = TRUE)
summary(logisticModel)
coef(logisticModel)

x3 <- seq(min(dias), max(dias), length=100)
y3 <- predict(logisticModel2, list(dias=x3))
plot(dias,ctc4_ct,las=1,pch=16,main="CTC4 - Controle")
points(x3, y3, type='l', col='blue')

#ctc4_100
logisticModel3 <- nls(ctc4_eth~M/(1+exp(Po+k*dias)), 
                      start=list(Po=5, k=-0.5, M=30),
                      control = list(maxiter = 50,minFactor = 1/4096),
                      trace=TRUE)
summary(logisticModel)
coef(logisticModel)

x4 <- seq(min(dias), max(dias), length=100)
y4 <- predict(logisticModel3, list(dias=x4))
plot(dias,ctc4_eth,las=1,pch=16,main="CTC4 - Ethephon")
points(x4, y4, type='l', col='blue')

dev.off()

############ Teste estatistico #########

# eh normal?
shapiro.test(ctc4_ct)
shapiro.test(ctc4_eth)
shapiro.test(vtx12_ct)
shapiro.test(vtx12_eth)

hist(ctc4_ct)
hist(vtx12_eth)

res <- wilcox.test(ctc4_ct,vtx12_ct,exact = FALSE)
res

res <- wilcox.test(ctc4_eth,vtx12_eth,exact = FALSE)
res

res <- wilcox.test(ctc4_ct,ctc4_eth,exact = FALSE)
res

res <- t.test(vtx12_ct,vtx12_eth,exact = FALSE,paired = FALSE)
res

boxplot(ctc4_ct, ctc4_eth, vtx12_ct, vtx12_eth,
        main = "Comparação entre variedades de cana",
        at = c(1,2,4,5),
        names = c("ctc4_ct", "ctc4_eth", "vtx12_ct", "vtx12_eth"),
        las = 1,
        col = c("blue","red"),
        border = "black",
        horizontal = FALSE,
        notch = TRUE,
        xlab="Variedades",
        ylab="Altura (cm)",
        cex.lab=1.5,
        cex.axis=1.5
)

library(tidyverse)
library(foreign)
library(qgraph)
library(psychonetrics)
install.packages("devtools")
devtools::install_github("SachaEpskamp/bootnet")
download.file('https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.2/bootnet_1.5.3.tgz')
install.packages("bootnet")
library(bootnet)

#Jordan et al.(2017)のデータをダウンロードして，読み込む
download.file("https://doi.org/10.1371/journal.pone.0182162.s004","pone.0182162.s004.sav")
data <- read.spss("pone.0182162.s004.sav", to.data.frame=TRUE)

#GAD-7（不安に関する７項目）のデータを使う。renameで変数名を整理する
data_gad <- data %>% 
  rename(gad7a = S_GAD7_a, gad7b = S_GAD7_b, gad7c = S_GAD7_c, gad7d = S_GAD7_d,gad7e = S_GAD7_e, gad7f = S_GAD7_f, gad7g = S_GAD7_g) %>% 
  select(gad7a, gad7b, gad7c, gad7d, gad7e, gad7f, gad7g)

#bootnetのestimateNetwork()関数で，thresholdとalphaを指定する
results_gad <- estimateNetwork(data_gad,default = "pcor", threshold = "sig", alpha = 0.05, corMethod = "cor_auto")
#GAD-7は4件法で順序変数なので，「corMethod = "cor_auto"」と指定して，自動的にポリコリック相関で推定
plot(results_gad, theme = "colorblind", cut = 0,layout = "spring",labels = TRUE)

install.packages('qgraph')
library(qgraph)
# dplyr パッケージをインストール（すでにインストール済みの場合は不要）
install.packages("dplyr")
# dplyr パッケージをロード
library(dplyr)
install.packages("magrittr")
library(magrittr)
install.packages("getmatrix")
library(getmatrix)

#モデル選択：刈り込み
ggm(data_gad) %>% runmodel %>%
  modelsearch(criterion = "bic") %>%
  getmatrix("omega") %>% qgraph(theme = "colorblind", cut = 0,layout = "spring",labels = TRUE)

#モデル選択：正規化
results_gad <- estimateNetwork(data_gad,default = "EBICglasso", corMethod = "cor_auto")
plot(results_gad, theme = "colorblind", cut = 0,layout = "spring",labels = TRUE)

#エッジの重みの正確度の推定結果の算出
accuracy_edge <- bootnet(results_gad, nBoots = 2500, nCores =4, statistics =   c("edge", "strength", "closeness", "betweenness"))
plot(accuracy_edge, labels = FALSE, order = "sample")

#中心性指標の算出
centralityPlot(results_gad, include = c("Strength", "Betweenness", "Closeness"))
#中心性指標の安定性の算出
stability_centrality <- bootnet(results_gad, nBoots = 2500, type = "case", nCores =4, statistics =  c("strength", "closeness", "betweenness"))
corStability(stability_centrality)

#有意差検定：gad7dとgad7fで差があるか検討する
differenceTest(accuracy_edge, 2, 6, "strength")

#エッジ間の差の検定結果をプロット
plot(accuracy_edge, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")

#Strengthにおけるノード間の差の検定結果をプロット
plot(accuracy_edge,"strength")

library(tidyverse)
library(foreign)
library(kableExtra)
library(bootnet)
library(qgraph)
library(psychonetrics)
library(NetworkComparisonTest)
install.packages('kableExtra')
install.packages('psychonetrics')


download.file("https://doi.org/10.1371/journal.pone.0182162.s004",
              "pone.0182162.s004.sav")
data <- read.spss("pone.0182162.s004.sav",to.data.frame=TRUE)

#PHQ9とGAD7のみ使用変数として抽出
data %>%
  select_if(grepl("PHQ9|GAD7", names(.))) ->usedata
cns1<-colnames(usedata)
v_names<-chartr("1-9","a-i",c(1:9))
cns2<-c(paste0("gad7",v_names[1:7]),
        paste0("phq9",v_names))
cols<-setNames(cns1, cns2)
usedata %>% 
  rename_(.dots=cols) -> usedata

#EBICglassoで正則化モデルで推定
library(bootnet)
results <- estimateNetwork(usedata,
                           default = "EBICglasso",
                           corMethod = "cor_auto")

#プロット
results_bridge<-plot(results, 
                     theme = "colorblind", 
                     cut = 0,layout = "spring",labels = TRUE, 
                     groups=list("GAD7"=c(1:7),"PHQ9"=c(8:16)))

#ブリッジ中心性指標の算出
library(networktools)
bridge_centrality <- bridge(results_bridge, 
                            communities=c(rep(1,7),rep(2,9)))
plot(bridge_centrality)

#ブリッジ中心性指標の安定性
bridge_stability <- bootnet(results, nBoots = 2500, 
                            type = "case", nCores =4, 
                            statistics =  c("bridgeStrength", 
                                            "bridgeCloseness", 
                                            "bridgeBetweenness"), 
                            communities = list("GAD7"=c(1:7), 
                                               "PHQ9"=c(8:16)))

corStability(bridge_stability)

#ブリッジ中心性指標の安定性のプロット
plot(bridge_stability,
     c("bridgeStrength","bridgeCloseness",
       "bridgeBetweenness"))

#性別変数をusedataに追加し、性別ごとのデータセットを抽出
usedata$gender<-data$gender
usedata %>% 
  filter(gender=="Weiblich")->women
usedata %>% 
  filter(gender=="Männlich")->men

#性別ごとのデータセットでネットワーク分析を実施
network_male <- estimateNetwork(women, 
                                default = "EBICglasso", 
                                corMethod = "cor_auto")
network_female <- estimateNetwork(men, 
                                  default = "EBICglasso", 
                                  corMethod = "cor_auto")

#群間でレイアウトが一致するように調整
#2群のネットワークの平均をとって比較するプロットのノードの位置を固定 (averageLayout)
#両群のエッジの中から最大値をとって、比較するプロットのエッジの最大値を固定 (maximum=Max) 　

L <- averageLayout(network_male, network_female)
Max <- max(abs(c(getWmat(network_male), getWmat(network_female)))) 
layout(t(1:2))
plot(network_male, layout = L, title = "Males", maximum = Max) 
plot(network_female, layout = L, title = "Females", maximum = Max)

#ネットワーク比較検定
install.packages('NetworkComparisonTest')
library(NetworkComparisonTest)
results_NCT<-NCT(na.omit(women[,-17]),na.omit(men[,-17]), 
                 test.edges=TRUE, test.centrality=TRUE, 
                 centrality=c("expectedInfluence","bridgeStrength"), 
                 communities=c(rep(1,7),rep(2,9)),
                 progressbar=F)

summary(results_NCT)

plot(results_NCT,
     what="network")

plot(results_NCT,
     what="strength")

results_NCT$einv.pvals%>% 
  kbl(.) %>%
  kable_styling()

results_NCT$diffcen.pval %>% 
  kbl(.) %>%
  kable_styling()



# network_meta
# We can use stata to make the network plot, contribution plot, forest plot, funnel plot, and cluster rank plot, while we can use R to make the Beyesian-based network meta analysis, plot the forest, get the relative effect for the efficacy-safety matrix, ranking probability distribution, SUCRA for both the efficacy and safety, inconsistancy and heterogeneity. 


#Stata coding
----------package one networkplot-------------
*install the package of networkplot
help networkplot   then install on website
*import the data with the variables of study, t1, t2, diff, stderr

*network plot
gen bias = .
replace bias = 1 if study == "Burger 2015"
replace bias = 2 if study == "Burger 2019"
replace bias = 3 if study == "Huang 2018"
*ROB with Low, unclear, high coded as 1,2,3 respectively
networkplot t1 t2, edgecolor(by bias)

*contribution plot
gen lnHR = ln(diff)
gen selnHR = ln(stderr)

netweight lnHR selnHR t1 t2

*inconsistency and funnel plot
ifplot lnHR selnHR t1 t2 study
netfunnel lnHR selnHR t1 t2, bycomparison

ifplot lnRR selnRR t1 t2 study
netfunnel lnRR selnRR t1 t2, bycomparison

*clusterank plot
clusterank safety efficacy treatment
label variable efficacy "SUCRAS for efficacy(Progression-free survival)"
label variable safety "SUCRAS for safety(Adverse events)"
db clusterank


---------package two mvmeta------------
*install the package of mvmeta
ssc install mvmeta
*set the network data
network setup r n, studyvar(study) trtvar(treatment) format(augment) rr
network map,improve
network meta i
network meta c
network forest
network sidesplit all, tau
network rank min, all zero reps(5000) gen(prob)

netleague, lab() sort() export("") eform

network convert pairs

netfunnel _y _stderr _t1 _t2 , random bycomp add(lfit _stderr _ES_CEN)noalpha ylabel(0 0.1 0.2 0.3)

netweight _y _stderr _t1 _t2
ifplot _y _stderr _t1 _t2 study, tau2(loop)




#R coding
```{r install the package of gemtc}
install.packages("gemtc")
library("gemtc")
```
```{r network analysis for count data}
network_count <- read.table(textConnection('
study responders sampleSize treatment
"Burger 2015" 68 132 Chlorambucil
"Burger 2015" 71 135 Ibrutinib
"Burger 2019" 67 104 Ibrutinib
"Burger 2019" 68 104 Ibrutinib_Rituximab
"Woyach 2018" 107 176 Bendamustine_Rituximab
"Woyach 2018" 74 180 Ibrutinib
"Sharman 2021" 48 58 Ibrutinib
"Sharman 2021" 45 59 Ibrutinib_Ublituximab
"Huang 2018" 31 52 Rituximab
"Huang 2018" 86 104 Ibrutinib
"Zucca 2017" 13 138 Rituximab
"Zucca 2017" 15 131 Chlorambucil'), header = TRUE)

network_count_data <- mtc.network(data.ab = network_count, description = "Network_count", treatments = NULL)

count_model <-mtc.model(network_count_data, type="consistency", n.chain=4,likelihood="binom",link="log",linearModel="random")

count_results <- mtc.run(count_model, n.adapt = 20000, n.iter = 50000, thin = 1)

summary(count_results)

summary(relative.effect(count_results, "Rituximab"))

forest(relative.effect(count_results, "Rituximab"))

gelman.plot(count_results)
count_ranks <- rank.probability(count_results)
print(count_ranks)
plot(count_ranks)

sucra_safety <- sucra(count_ranks)

write.csv(count_ranks, "count_ranks.csv")

rownames(count_ranks) <- count_ranks$intervention
count_ranks <- subset (count_ranks, select = -intervention)
count_ranks1 <- stack(count_ranks)
View(count_ranks1)
write.csv(count_ranks1, "ranks1.csv")

ggplot(data=ranks1,aes(x=intervention, y=probability, fill=rank))+geom_bar(stat="identity", position=position_dodge())+scale_fill_brewer(palette="Paired") + scale_x_discrete(limits=c("Ibrutinib", "Chlorambucil", "Rituximab","Ibrutinib_Rituximab","Ibrutinib_Ublituximab","Bendamustine_Rituximab")) + coord_flip()+labs(title="Safety(Adverse events)")

safety_network <- relative.effect.table(count_results)
write.csv(safety_network, "safety_network.csv")

#inconsistancy
resultnodesplit_count <- mtc.nodesplit(network_count_data,likelihood="binom",link="log",linearModel="random")
print(summary(resultnodesplit_count))
plot(summary(resultnodesplit_count))

#heterogeneity
resultanohe_count <- mtc.anohe(network_count_data,likelihood="binom",link="log",linearModel="random")
print(summary(resultanohe_count))
plot(summary(resultanohe_count))
```

```{r network analysis for survival data}
network_sample <- read.table(textConnection('
study diff std.err treatment
"Burger 2015" 6.25 0.048469387755102 Chlorambucil
"Burger 2015" NA NA Ibrutinib
"Burger 2020" 6.84931506849315 0.0306122448979592 Chlorambucil
"Burger 2020" NA NA Ibrutinib
"Burger 2019" 1.162 0.482142857142857 Ibrutinib
"Burger 2019" NA NA Ibrutinib_Rituximab
"Woyach 2018" 2.7027027027027 0.0790816326530612 Bendamustine_Rituximab
"Woyach 2018" NA NA Ibrutinib
"Sharman 2021" 2.17391304347826 0.160714285714286 Ibrutinib
"Sharman 2021" NA NA Ibrutinib_Ublituximab
"Huang 2018" 5.55555555555556 0.0517857142857143 Rituximab
"Huang 2018" NA NA Ibrutinib
"Zucca 2017" 1.1 0.211734693877551 Rituximab
"Zucca 2017" NA NA Chlorambucil'), header = TRUE)

network_sample$diff <- as.numeric(network_sample$diff)
network_sample_data <- mtc.network(data.re = network_sample, description = "Network", treatments = NULL)

plot(network_sample_data)

sample_model <-mtc.model(network_sample_data, type="consistency", n.chain=4,likelihood="binom",link="cloglog",linearModel="random")
sample_results <- mtc.run(sample_model, n.adapt = 20000, n.iter = 50000, thin = 1)

summary(sample_results)
summary(relative.effect(sample_results, "Rituximab"))

forest(relative.effect(sample_results, "Rituximab"))
forest(relative.effect(sample_results, "Rituximab",col.square="green",col.diamond="blue"))

forest(relative.effect(sample_results, "Rituximab"), xlim=c(0.01,70))


gelman.plot(sample_results)
ranks <- rank.probability(sample_results)
print(ranks)
plot(ranks)

write.csv(ranks, "ranks_pfs_original.csv")

rownames(ranks) <- ranks$intervention
ranks <- subset (ranks, select = -intervention)
ranks_pfs <- stack(ranks)
View(ranks_pfs)
write.csv(ranks_pfs, "ranks_pfs.csv")

ggplot(data=ranks_pfs,aes(x=intervention, y=probability, fill=rank))+geom_bar(stat="identity", position=position_dodge())+scale_fill_brewer(palette="Paired") + scale_x_discrete(limits=c("Ibrutinib", "Chlorambucil", "Rituximab","Ibrutinib_Rituximab","Ibrutinib_Ublituximab","Bendamustine_Rituximab")) + coord_flip()+labs(title="Efficacy(Progression-free survival)")

sucra_pfs <- sucra(ranks)
sucra_pfs <- sucra(ranks,lower.is.better=TRUE)
sucra_pfs


a <- relative.effect.table(sample_results)
write.csv(a, "network_meta_pfs.csv")

#inconsistancy
resultnodesplit <- mtc.nodesplit(network_sample_data,likelihood="binom",link="cloglog",linearModel="random")
print(summary(resultnodesplit))
plot(summary(resultnodesplit))

#heterogeneity
resultanohe <- mtc.anohe(network_sample_data,likelihood="binom",link="cloglog",linearModel="random")
print(summary(resultanohe))
plot(summary(resultanohe))
```





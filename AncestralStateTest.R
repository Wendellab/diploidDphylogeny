# inital code tests for reconstructing ancestral states (of genome size and TE total Mb) using both ace (pic) and contMap
# exemplary clusters were selected for the test
# results are used for GS size presentation at PAG

####################################################################################################################################
# R version 3.2.3 (2015-12-10)                                                                                                     #
# Platform: x86_64-w64-mingw32/x64 (64-bit)                                                                                        #
# Running under: Windows 8.1 x64 (build 9600)                                                                                      #
#                                                                                                                                  #
# locale:                                                                                                                          #
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252         #
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252                                                    #
#                                                                                                                                  #
# attached base packages:                                                                                                          #
# [1] stats     graphics  grDevices utils     datasets  methods   base                                                             #
#                                                                                                                                  #
# other attached packages:                                                                                                         #
# [1] phytools_0.5-10 maps_3.0.1      ape_3.4                                                                                      #
#                                                                                                                                  #
# loaded via a namespace (and not attached):                                                                                       #
#  [1] igraph_1.0.1            XVector_0.8.0           magrittr_1.5            BiocGenerics_0.14.0     splines_3.2.3               #
#  [6] zlibbioc_1.14.0         nnls_1.4                MASS_7.3-45             IRanges_2.2.7           scatterplot3d_0.3-36        #
# [11] mnormt_1.5-3            lattice_0.20-33         quadprog_1.5-5          tools_3.2.3             parallel_3.2.3              #
# [16] grid_3.2.3              nlme_3.1-122            msm_1.6                 clusterGeneration_1.3.4 phangorn_2.0.1              #
# [21] plotrix_3.6-1           survival_2.38-3         numDeriv_2014.2-1       Matrix_1.2-3            S4Vectors_0.6.6             #
# [26] animation_2.4           Biostrings_2.36.4       stats4_3.2.3            expm_0.999-0            mvtnorm_1.0-3               #
#                                                                                                                                  #
####################################################################################################################################


library(ape)
library(phytools)
GStree <- read.nexus("GSandCluster.nex")
GStable <- read.table("GSandCluster.table")
GSv <- GStable$genome_size
C3v <- GStable$cluster3_gypsy_Mb
C5v <- GStable$cluster5_gypsy_Mb
C24v <- GStable$cluster24_gypsy_Mb
C55v <- GStable$cluster55_copia_Mb
C62v <- GStable$cluster62_copia_Mb
names(GSv) <-row.names(GStable)
names(C3v) <-row.names(GStable)
names(C5v) <-row.names(GStable)
names(C24v) <-row.names(GStable)
names(C55v) <-row.names(GStable)
names(C62v) <-row.names(GStable)

GSgradient <- contMap(GStree, GSv, res=1000, plot=TRUE)
picGS <- ace(GSv, GStree, type = "continuous", method = "pic")
ancGS <- format(picGS$ace, digits = 2)
plot(GSgradient)
nodelabels(ancGS, adj=c(1.1,-0.6), frame="none", cex=0.8)

C3g <- contMap(GStree, C3v, res=1000, plot=TRUE)
picC3 <- ace(C3v, GStree, type = "continuous", method = "pic")
ancC3 <- format(picC3$ace, digits = 2)
plot(C3g)
nodelabels(ancC3, adj=c(1.5,-0.6), frame="none", cex=0.8)

C5g <- contMap(GStree, C5v, res=1000, plot=TRUE)
picC5 <- ace(C5v, GStree, type = "continuous", method = "pic")
ancC5 <- format(picC5$ace, digits = 2)
plot(C5g)
nodelabels(ancC5, adj=c(1.5,-0.6), frame="none", cex=0.8)

C24g <- contMap(GStree, C24v, res=1000, plot=TRUE)
picC24 <- ace(C24v, GStree, type = "continuous", method = "pic")
ancC24 <- format(picC24$ace, digits = 2)
plot(C24g)
nodelabels(ancC24, adj=c(1.5,-0.6), frame="none", cex=0.8)

C55g <- contMap(GStree, C55v, res=1000, plot=TRUE)
picC55 <- ace(C55v, GStree, type = "continuous", method = "pic")
ancC55 <- format(picC55$ace, digits = 2)
plot(C55g)
nodelabels(ancC55, adj=c(1.5,-0.6), frame="none", cex=0.8)

C62g <- contMap(GStree, C62v, res=1000, plot=TRUE)
picC62 <- ace(C62v, GStree, type = "continuous", method = "pic")
ancC62 <- format(picC62$ace, digits = 2)
plot(C62g)
nodelabels(ancC62, adj=c(1.5,-0.6), frame="none", cex=0.8)

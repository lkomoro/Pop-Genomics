#This script is for QA/QC and basic population structure analyses using the package StrataG (see E. Archer's github for details on this package)
#This is meant to be a template/example that can be amended for other datasets-currently code matches a hawksbill turtle microsatellite dataset
#Plan to mutate this into a Markdown that is easier to follow


# make sure you have Rtools installed
if (!require('devtools')) install.packages('devtools')#can comment out if/once you have installed
# install from GitHub
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)#need internet connection. don't need to run
#every time, just if you think Eric's made changes and you want to update

#intro####
setwd('~/Dropbox/home-NOAAtransitionfiles/Ei_Msat_StrataG')#change this accordingly
getwd()#this just double checks that you are where you wanted to be
library(strataG)
rm(list = ls())
options(stringsAsFactors = F)

description <- "EiMSat_April_20"
strata.file <- "Ei_Msat_strata.csv"
genotypes.file <- "EiMsatData020915.csv"
strata <- "Population"  #name of stratum to use in analysis, not really using currently below, but could

#read in data ####
Msat.Ei.geno <- readGenData(genotypes.file) #analysis set export with population column removed
    #missing genotype must be "NA" (not "0")
Msat.Ei.strata <- readGenData(strata.file) #LABID and population from analysis set.
strata.schemes <- Msat.Ei.strata[, c("Population","LABID")]#note-here you can have additional columns, and then this makes it easy to switch between strata levels below
rownames(strata.schemes) <- Msat.Ei.strata$LABID

head(Msat.Ei.geno)
#change 0s to NAs (if needed)
Msat.Ei.geno[Msat.Ei.geno==0]=NA
#check that they are NAs
View(Msat.Ei.geno)

#Run these lines if you want to take out the samples that have missing data over a certain percent: #####
numNAs<-apply(Msat.Ei.geno, 1, function(z) sum(is.na(z)))
Msat.Ei.geno.a<-Msat.Ei.geno[!(numNAs>((ncol(Msat.Ei.geno)-(NCOL(Msat.Ei.strata)-1))*0.75)),]
#this takes out rows that are missing more data for than 75% of the markers
#you can change the percent in the formula as you like (right now is 0.75 above)
#note that we are counting the NAs, so the percent is the inversion, 
#i.e. if more than 75% of the loci have NAs-we have data for less than 25% of the markers

#run if want to remove samples where have less than 3 in a strata or missing data for an entire locus:
library(plyr)
#Msat.Ei.geno.short<-merge (Msat.Ei.geno.a,Msat.Ei.strata)
#table(Msat.Ei.geno.short$Population)
#Msat.Ei.geno.short<-subset(Msat.Ei.geno.short,Population!="MEXICO")
#Msat.Ei.geno.short<-subset(Msat.Ei.geno.short,Population!="PALAU") 
#Msat.Ei.geno.short<-subset(Msat.Ei.geno.short,Population!="PANAMA") 
#Msat.Ei.geno.short<-Msat.Ei.geno.short[,1:51]

#--- create a diploid gtypes object####
Msat.g <- new("gtypes", gen.data = Msat.Ei.geno.a[, -1], ploidy = 2,
                     ind.names = Msat.Ei.geno.a[, 1], schemes = strata.schemes, strata="Population")
summary(Msat.g) #this gives a quick summary
save(Msat.g, file = paste(description, "_gtypes.rdata", sep=""))

#can do the same for ones with small groups removed...copy/amend code above as needed

# Test for all QA/QC tests ####
#(locus summaries, check for duplicates, etc.)
qaqc(Msat.g, label = paste(description, "_QC_results", sep=""), num.shared = 0.9) 

#Note-in the locus summary files:
#1 the propUniqueAlleles is more useful for haplotype (mitochondrial) bc it is the percentage of alleles at that locus that only occur once (so for diploid markers, this=one individual is heterozygous)
#2 the allelic richness of each locus calculated as the number of alleles divided by the number of samples without missing data at that locus.

#check for private alleles:
pA<-privateAlleles(Msat.g)
write.csv(pA, file = paste(description, "_private_alleles.csv", sep = ""))

#check for LD:
LD<-LDgenepop(Msat.g)#for all together
write.csv(LD, file = paste(description, "_all_LD.csv", sep = ""))
#x <- strataSplit(Msat.g)#split up by strata, and then can call it to do something
#by strata-
#test first with something easy like just printing num Alleles
for(g in strataSplit(Msat.g)) print(numAlleles(g))#this is how you loop it through all the strata, use print if you just want to see it
#to write them to file:
for(g in strataSplit(Msat.g)) {
  LD <- LDgenepop(g, use.genepop = T)
  fname <- paste(description,"_LD_", strataNames(g)[1], ".csv", sep = "")
  write.csv(LD, fname)
}
# HWE-
#for all together-
HWE<-hweTest(Msat.g, use.genepop = TRUE, label = "HWE.genepop")#note the default is genepop is false bc eric didnt want to assume people had that, but is actually more accurate so change to true
write.csv(HWE, file = paste(description, "_all_HWE.csv", sep = ""))

#looping through each strata...
for(g in strataSplit(Msat.g)) {
  hwe <- hweTest(g, use.genepop = T)
  na <- numAlleles(g)
  result <- cbind(num.alleles = na, hwe)
  fname <- paste(description,"_HWE_ ", strataNames(g)[1], ".csv", sep = "")
  write.csv(result, fname)
}# or use this if want to run and bind multiple things together that have the same dimensions
  

#Allele Frequencies:
z<-alleleFreqs(Msat.g, by.strata = FALSE)# list of allele frequencies for each locus. Each element is a matrix or array with 
#frequencies by count (freq) and proportion (prop) of each allele.
#can change to strata=TRUE if want by to do for each strata
z#just look at on screen for now to check-different loci have different lengths (# of alleles) 
z$ERIM04t# to look at one-I need to come back and fix code to write to file by locus
write.csv(z$ERIM04t, file = paste(description, "_ERIM04t_alleleFreq.csv", sep = ""))#example-can do one by one but I will put the loci names in a list and loop them

#this jackknifes n samples at a time (set to 1) to see how affects HWE 
Msat.Ei.JackHWE<-jackHWE(Msat.g, exclude.num = 1, min.hwe.samples = 3, show.progress = TRUE,
        use.genepop = TRUE)#performs a HWE jackknife where all combinations of exclude.num samples are left out and HWE is recalculated
Msat.Ei.JackHWE.influential<-jackInfluential(Msat.Ei.JackHWE, alpha = 0.05)
Msat.Ei.JackHWE.influential #HWE and identifies "influential" samples. Samples are "influential" if the observed HWE p-value is < alpha, but is > alpha when the samples are not present.
write.csv(Msat.Ei.JackHWE.influential$influential, file = paste(description, "_JackHWE.influential.csv", sep = ""))
write.csv(Msat.Ei.JackHWE.influential$allele.freqs, file = paste(description, "_JackHWE.allele.freqs.csv", sep = ""))
write.csv(Msat.Ei.JackHWE.influential$odds.ratio, file = paste(description, "_JackHWE.odds.ratio.csv", sep = ""))

pdf(paste(description, "_JackHWE.influential", ".pdf", sep = ""), width = 15, height = 8)
#creates a cumulative frequency plot of all odds-ratios from jack.influential. 
#A vertical dashed line marks the smallest influential exclusion
plot(Msat.Ei.JackHWE.influential, main = "Msat.Ei.JackHWE.influential")#need to check in with Eric to fix dimensions so can read #
dev.off()

# Run just overall, pairwise and/or both tests for all or specified metrics####
#All:
#popStructTest(Msat.Ei.g.Pop, nrep = 10, stats = "all", type = c("both", "overall","pairwise"), keep.null = FALSE, quietly = FALSE, num.cores = 1,write.output = FALSE)
#Overall:
overall.1 <- overallTest(Msat.g, stats = "all", nrep = 10000, write.output=FALSE)
write.csv(overall.1$result, file = paste(description, "_overall_test.csv", sep = ""))
#test with nrep <100; use nrep=1000 or 10000
#write.output=TRUE prints pairwise matrix for each test, with test value in the upper right
#quadrant and p-values in the lower left quadrant (only makes sense for pairwise)

# Run just pairwise tests for specified metrics####
pairwise.1 <- pairwiseTest(Msat.g, stats = c("fst", "fst.prime","gst.dbl.prime"),
                            nrep = 10000, write.output=TRUE) #set nrep = 10000 after testing at <100
write.csv(pairwise.1$result, file = paste(description, "_pairwise_results.csv", sep = ""))
warnings()
#Can save to file Warning messages about dropped loci to remember if applicable

#Calculate observed heterozygosity for each locus ()####
oh <- obsvdHet(Msat.g)
write.csv(oh, paste(description, "_Ho.csv", sep = ""))
#oh here is the same as that found in the qaqc output files, this just makes it available for easy graphing manipulation

#quick exploratory plots of some metrics:####
#read in heterozygosity data
library(ggplot2)
data<-read.csv("EiMSat_April_20_Ho.csv")
data<-rename(data, c("X"="locus", "x"="mean.heterozygosity"))
head(data)
#plot the mean.heterozygosity####
pdf(paste(description, "_mean_Het_by_locus", ".pdf", sep = ""), width = 15, height = 8)

ggplot(data = data, aes(x = locus, y = mean.heterozygosity)) + geom_point(size=4, shape = 18, fill = "black") + 
  labs(x = "Locus", y = "Observed Heterozygosity\n") + theme_classic()
dev.off()

#add more plots as desired for other metrics later-See Peter's request of by allele drawn on board
#Run STRUCTURE:#### see other script, add in as needed here.
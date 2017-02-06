rm(list = ls())

data(dolph.msats)
data(dolph.strata)
dolph.strata$Offshore <- dolph.strata$broad #copy the values into a new column
dolph.strata$Offshore[dolph.strata$Offshore == "Coastal"] <- NA #make all the 'coastal' values NA, so now this column is a subset denoting just the offshore animals

#First way to make a gtypes object:
d.msats <- new("gtypes", gen.data = dolph.msats[, -1], 
               ploidy = 2, ind.names = dolph.msats[, 1],
               schemes = dolph.strata, description = "my stuff"
)

#Make a gtypes object without any strata:
d.msats2 <- df2gtypes(dolph.msats, ploidy = 2, id.col = 1, 
                      strata.col = NULL, loc.col = 2)

#Merge 'fine' strata column into genotype file, then make gtypes (not sure of rationale here?)
strat.msats <- merge(dolph.strata[, c("id", "fine")], 
                     dolph.msats, by = "id", all.y = T)

d.msats3 <- df2gtypes(strat.msats, ploidy = 2, schemes = dolph.strata)


# haploid data

#create from df values only/directly:
d.dloop <- df2gtypes(dolph.strata[, c(1, 4, 2)], ploidy = 1)

data(dolph.haps)
data(dolph.seqs)

#create with hap IDs from df, but then link haplotypes to sequences in as well from other file:
#(depends on types of analyses wish to run)
dloop.haps <- df2gtypes(dolph.strata[, c(1, 4, 2)], ploidy = 1,
                        sequences = dolph.haps, description = "haplotypes")

#link IDs to sequences
dloop.seqs <- df2gtypes(dolph.strata[, c(1, 4, 1)], ploidy = 1,
                        sequences = dolph.seqs, description = "ids")

# label id sequences to haplotypes
x <- labelHaplotypes(as.DNAbin(dolph.seqs), prefix = "Ttru.")

new.haps <- labelHaplotypes(dloop.seqs, prefix = "Ttru.")$gtypes


# gtypes from sequences only
seq.g <- sequence2gtypes(dolph.seqs)

# add schemes
schemes(seq.g) <- dolph.strata

# stratify to fine
fine.seq <- stratify(seq.g, "fine")
off.seq <- stratify(seq.g, "Offshore")
off.seq.no.drop <- stratify(seq.g, "Offshore", drop = FALSE)

# indexing
msats.g[c(T, F), , ]

sex <- sample(c("M", "F"), nInd(msats.g), rep = T) #creating new variable (in real data would have this information in input file)

f.msats <- msats.g[sex == "F", , ]#parse out

id.msats <- msats.g[c("26318", "41538", "78062"), , ]

id.msats <- msats.g[c("26318", "41538", "78062"), c("EV37", "D11t") , ]

st.msats <- msats.g[, , "Coastal"]

# export to matrix
msats.mat <- as.matrix(msats.g)
msats.mat <- cbind(id = rownames(msats.mat), strata = strata(msats.g), msats.mat)
write.csv(msats.mat, file = "msats.g.csv", row.names = F)

# reading data in
gen.data <- readGenData("msats.g.csv")


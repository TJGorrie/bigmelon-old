### R code from vignette source 'bigmelon.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: UnevaluatedCode (eval = FALSE)
###################################################
## source('http://bioconductor.org/biocLite.R')
## biocLite('wateRmelon', 'gdsfmt')


###################################################
### code chunk number 3: UnevaluatedCode (eval = FALSE)
###################################################
## install.packages('bigmelon_0.99.4.tar.gz', repos = NULL, type = 'source')


###################################################
### code chunk number 4: code-block
###################################################
library(bigmelon)


###################################################
### code chunk number 5: UnevaluatedCode (eval = FALSE)
###################################################
## # read in an IDAT file with barcode 'sentrixid_rnncnn'
## gfile <- iadd('sentrixid_rnncnn', gds = 'melon.gds')
## gfile <- iadd2('Data/IDATLocations/dataset', gds = 'melon.gds', chunksize = 100)


###################################################
### code chunk number 6: code-block
###################################################
data(melon)
gfile <- es2gds(melon, 'melon.gds')


###################################################
### code chunk number 7: UnevaluatedCode (eval = FALSE)
###################################################
## library(methylumi)
## # read Illumina methylation data into a MethyLumiSet object
## melon <- methyLumiR('finalreport.txt')
## # read Illumina methylation final report into a gds.class object.
## gfile <- finalreport2gds('finalreport.txt', gds='melon.gds')


###################################################
### code chunk number 8: UnevaluatedCode (eval = FALSE)
###################################################
## # convert a MethyLumiSet object to a gds.class object
## gfile <- es2gds(melon, 'melon.gds')


###################################################
### code chunk number 9: code-block
###################################################
print(gfile)


###################################################
### code chunk number 10: code-block
###################################################
index.gdsn(gfile, 'betas')
class(index.gdsn(gfile, 'betas'))
# Access nodes with additional nodes inside
index.gdsn(gfile, 'fData/TargetID')


###################################################
### code chunk number 11: code-block
###################################################
betas(gfile)
class(betas(gfile))


###################################################
### code chunk number 12: code-block
###################################################
ls.gdsn(gfile)
# Look into nodes with additional nodes
ls.gdsn(index.gdsn(gfile, 'fData'))


###################################################
### code chunk number 13: code-block
###################################################
# Call a gdsn.class node
anode <- betas(gfile)
anode
class(anode)
dat <- read.gdsn(anode)
dim(dat)
head(dat)
datsub <- readex.gdsn(anode, sel = list(1:5, 1:3))
dim(datsub)
datsub


###################################################
### code chunk number 14: code-block
###################################################
# Re-using node from previous example
anode
datsub <- anode[1:5,1:3]
dim(datsub)
datsub
# Additionally, the row and col names can be turned
# for faster indexing.
anode[1:5, 1:3, name = FALSE]


###################################################
### code chunk number 15: code-block
###################################################
# Logical Indexing
anode[1:5,c(T,F,F)]
# Ordering calls
anode[c(5,9,1,500,345), c(8,4,1,3)]
# Indexing by characters (and drop)
anode[c('cg00000029', 'cg00000236'), '6057825008_R02C01', drop=F]
# Loading entire data (no indexing)
dat <- anode[ , ]
dim(dat)


###################################################
### code chunk number 16: code-block
###################################################
gfile[1:5, 1:3, node = 'betas', name = TRUE]
gfile[1:5, 1:3, node = 'methylated', name = TRUE]


###################################################
### code chunk number 17: code-block
###################################################
read.gdsn(index.gdsn(gfile, "paths"))
head(read.gdsn(index.gdsn(gfile, "fData/TargetID")))
head(read.gdsn(index.gdsn(gfile, "pData/sampleID")))


###################################################
### code chunk number 18: code-block
###################################################
rawmet <- methylated(gfile)[,]
rawume <- unmethylated(gfile)[,]


###################################################
### code chunk number 19: IncludeGraphic
###################################################
boxplot(log(rawmet), las=2, cex.axis=0.8)


###################################################
### code chunk number 20: IncludeGraphic1
###################################################
boxplot(log(rawume), las=2, cex.axis=0.8)


###################################################
### code chunk number 21: code-block
###################################################
rawbet <- betas(gfile)[,]


###################################################
### code chunk number 22: IncludeGraphic2
###################################################
outlyx(rawbet, plot = TRUE)


###################################################
### code chunk number 23: IncludeGraphic3
###################################################
outlyx(gfile, plot = TRUE, perc = 0.01)


###################################################
### code chunk number 24: code-block
###################################################
pfilter(gfile)


###################################################
### code chunk number 25: code-block
###################################################
backup.gdsn(gds = NULL, node = index.gdsn(gfile, 'betas'))
ls.gdsn(index.gdsn(gfile, 'backup'))


###################################################
### code chunk number 26: code-block
###################################################
f <- createfn.gds('melon2.gds')
backup.gdsn(gds = f, node = index.gdsn(gfile, 'betas'))
f
copyto.gdsn(node = f, source = index.gdsn(gfile, 'betas'), name = 'betacopy') 
f
copyto.gdsn(node = gfile, source = index.gdsn(gfile, 'betas'), name = 'betacopy')
# Close File
closefn.gds(f)


###################################################
### code chunk number 27: UnseenCodeAndOutput
###################################################
unlink('melon2.gds')


###################################################
### code chunk number 28: code-block
###################################################
 dasen(gfile)
 # Alternatively it is possible to store normalized betas to a separate node
 dasen(gfile, node="normbeta")
 index.gdsn(gfile, "normbeta")


###################################################
### code chunk number 29: code-block
###################################################
# Example of apply.gdsn
apply.gdsn(betas(gfile), margin = 2, as.is='double', FUN = function(x,y){
  mean(x, na.rm=y)
  }, y = TRUE)


###################################################
### code chunk number 30: code-block
###################################################
gds2mlumi(gfile)
gds2mset(gfile, anno="450k")


###################################################
### code chunk number 31: code-block
###################################################
# Closing the connection
closefn.gds(gfile)


###################################################
### code chunk number 32: UnseenCodeAndOutput
###################################################
unlink('melon.gds', force = TRUE)


###################################################
### code chunk number 33: code-block
###################################################
sessionInfo()



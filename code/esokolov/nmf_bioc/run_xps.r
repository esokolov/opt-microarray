library(xps)

### define directories:
# directory containing Affymetrix annotation files
annDir <- "/home/affy/annotation"
# directory to store ROOT scheme files
scmDir <- "/home/affy/schemes"

annotationVersion<-"32"
libraryVersion<-"4" # 4th version is for treating HuGene 1.0 as exon arrays

if( !file.exists(paste(scmDir,"/Scheme_HuGene10stv1r",libraryVersion,"_na",annotationVersion,".root",sep=""))){
	scheme.hugene10stv1 <- import.exon.scheme(paste("Scheme_HuGene10stv1r",libraryVersion,"_na",annotationVersion,sep=""),filedir=scmDir,
	layoutfile=paste(annDir,"/HuGene-1_0-st-v1.r",libraryVersion,".clf",sep=""),
	schemefile=paste(annDir,"/HuGene-1_0-st-v1.r",libraryVersion,".pgf",sep=""),
	probeset=paste(annDir,"/HuGene-1_0-st-v1.na",annotationVersion,".hg19.probeset.csv",sep=""),
	transcript=paste(annDir,"/HuGene-1_0-st-v1.na",annotationVersion,".hg19.transcript.csv",sep="")
	)
}else{
	scheme.hugene10stv1<-root.scheme(paste(scmDir,"/Scheme_HuGene10stv1r",libraryVersion,"_na",annotationVersion,".root",sep=""))
}

root.data(scheme.hugene10stv1, rootfile = "/home/sokolov/tmp/sample1000_cel.root", celnames = "*")

files <- list.files("/mnt/teradisk/affy/HuGene10_ST_CELS/sample/", "*.CEL", ignore.case = F)
sampleCels1 <- import.data(xps.scheme = scheme.hugene10stv1, filedir = "/mnt/teradisk/affy/HuGene10_ST_CELS/sample",
							filename = "sample1000_part1", celdir = "/mnt/teradisk/affy/HuGene10_ST_CELS/sample", celfiles = files[1:500])
sampleCels2 <- import.data(xps.scheme = scheme.hugene10stv1, filedir = "/mnt/teradisk/affy/HuGene10_ST_CELS/sample",
							filename = "sample1000_part2", celdir = "/mnt/teradisk/affy/HuGene10_ST_CELS/sample", celfiles = files[501:1000])
							
sampleExpr1 <- xps::rma(xps.data = sampleCels1, filename = "sample1000_part1_expr", filedir = "/home/affy/HuGene10_ST_CELS/sample/",
					option = "transcript", background = "antigenomic", normalize = TRUE, exonlevel = "core+metacore")
sampleExpr2 <- xps::rma(xps.data = sampleCels1, filename = "sample1000_part2_expr", filedir = "/home/affy/HuGene10_ST_CELS/sample/",
					option = "transcript", background = "antigenomic", normalize = TRUE, exonlevel = "core+metacore")

sampleExpr1 <- attachExpr(sampleExpr1)
sampleExpr2 <- attachExpr(sampleExpr2)


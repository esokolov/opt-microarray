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

files = c("Alisov_1.1_62_(HuGene-1_0-st-v1).CEL",
			"Alisov_1.2_63_(HuGene-1_0-st-v1).CEL",
			"Antonov_1.1_59_(HuGene-1_0-st-v1).CEL",
			"Antonov_1.2_60_(HuGene-1_0-st-v1).CEL",
			"Blohin_1.1_26_(HuGene-1_0-st-v1).CEL",
			"Blohin_1.2_27_(HuGene-1_0-st-v1).CEL",
			"Chekalenko_1.1_65_(HuGene-1_0-st-v1).CEL",
			"Chekalenko_1.2_66_(HuGene-1_0-st-v1).CEL",
			"Doronin_1.1_68_(HuGene-1_0-st-v1).CEL",
			"Doronin_1.2_69_(HuGene-1_0-st-v1).CEL",
			"Kolpakov_1.1_53_(HuGene-1_0-st-v1).CEL",
			"Kolpakov_1.2_54_(HuGene-1_0-st-v1).CEL",
			"Korolkov_1.1_32_(HuGene-1_0-st-v1).CEL",
			"Korolkov_1.2_33_(HuGene-1_0-st-v1).CEL",
			"Slavutskij_1.1_1_(HuGene-1_0-st-v1).CEL",
			"Slavutskij_1.2_2_(HuGene-1_0-st-v1).CEL",
			"Stepanov_1.1_35_(HuGene-1_0-st-v1).CEL",
			"Stepanov_1.2_36_(HuGene-1_0-st-v1).CEL",
			"Sveridov_1.1_56_(HuGene-1_0-st-v1).CEL",
			"Sveridov_1.2_57_(HuGene-1_0-st-v1).CEL")

cels <- import.data(xps.scheme = scheme.hugene10stv1, filedir = "/home/sokolov/tmp",
							filename = "skiers", celdir = "/home/affy/CEL/skiers/", celfiles = files)
cels <- attachInten(cels)
inten <- intensity(cels)
inten$Probe.ID <- inten$Y * 1050 + inten$X + 1
probes<-read.table("probes_worthy.tab",sep="\t")
row.names(probes)<- probes$Probe.ID

inten <- subset(inten, Probe.ID %in% probes$Probe.ID)[,-c(1:2)]
inten$Transcript.Cluster.ID <- probes$Transcript.Cluster.ID[order(probes$Probe.ID)]
inten <- inten[order(inten$Transcript.Cluster.ID, inten$Probe.ID),]
write.table(inten, file="raw_inten_skiers.tab", sep="\t", row.names=F, col.names = F)
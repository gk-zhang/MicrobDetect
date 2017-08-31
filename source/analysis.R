load("gc.RData")


########################################
#
# do the analysis for each level of taxonomy
#
########################################


################################################################
#
#create the abundance files of tp,nt,nb, for each taxonomy level
#
################################################################

# tissue.group: one of the three groups: tp, nt, nb
# tax.level: one of the six taxonomy levels "species", "genus", "family", "order", "class", "phylum", "ggenome_file_name"
# return a data frame containing the abundance at given taxonomy level for the given tissue group

# need package "data.table", if not installed, do:
# install.packages("data.table")

library(data.table)

tax.abundance <- function(tissue.group, tax.level){

		#dt <- data.table(taxonomy_micro)
		#dt[, grp := .GRP, by = paste(tax.level)]
		#setkey(dt, grp)
		#level <- dt[, list(list(.SD)), by = grp]$V1
	
		level <- split(taxonomy_micro,taxonomy_micro[,paste(tax.level)])
			
		data <- data.frame(tissue.group[1])
		for (i in 1:length(level)) {
			tax <- subset(tissue.group, select=c(match(level[[i]][,1],list_micro)+1))
			dat <- data.frame(rowMeans(tax))
			colnames(dat) <- names(level)[i]
			data <- cbind(data,dat)
		}
		return(data) 
}

###########################################################
#
# package "multtest" is required, if not installed yet, do:
#
## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("multtest")
#
###########################################################
library(multtest)

#do the paired t-test at different taxonomy level
#parameter:
#tissue.group.1, one of the three groups: tp, nt, nb
#tissue.group.2, one of the three groups: tp, nt, nb
#tax.level: one of the six taxonomy levels "species", "genus", "family", "order", "class", "phylum"
#mt.test: set the number (1 to 9) for different methods for p-values adjustment, we used "BY Adjusted p-values for the Benjamini & Yekutieli (2001)" (7)
### 1: Bonferroni Bonferroni single-step adjusted p-values for strong control of the FWER.
### 2: Holm Holm (1979) step-down adjusted p-values for strong control of the FWER.
### 3: Hochberg Hochberg (1988) step-up adjusted p-values for strong control of the FWER (for raw (unadjusted) p-values satisfying the Simes inequality).
### 4: SidakSS Sidak single-step adjusted p-values for strong control of the FWER (for positive orthant dependent test statistics).
### 5: SidakSD Sidak step-down adjusted p-values for strong control of the FWER (for positive orthant dependent test statistics).
### 6: BH Adjusted p-values for the Benjamini & Hochberg (1995) step-up FDR-controlling procedure (independent and positive regression dependent test statistics).
### 7: BY Adjusted p-values for the Benjamini & Yekutieli (2001) step-up FDR-controlling procedure (general dependency structures).
### 8: ABH Adjusted p-values for the adaptive Benjamini & Hochberg (2000) step-up FDR-controlling procedure.  This method ammends the original step-up procedure using an estimate of the number of true null hypotheses obtained from p-values.
### 9: TSBH Adjusted p-values for the two-stage Benjamini & Hochberg (2006) step-up FDR-controlling procedure.
#
#return a list with the following objects:
#"spe.sig":a list of microbes names with the significant different between the given groups at the given taxonomy level
#"spe.all": a list of microbes names at the given taxonomy level
#"gp1":a data frame with the abundance for the tissue group 1 at the given taxonomy level 
#"gp2":a data frame with the abundance for the tissue group 2 at the given taxonomy level 
#"gp1.sig":a data frame with the abundance for the tissue group 1 at the given taxonomy level, only for the significant species
#"gp2.sig":a data frame with the abundance for the tissue group 2 at the given taxonomy level, only for the significant species
#"gp1.gp2.p.ordered": a ordered list of p values between group1 and group2, NA removed
paired.testing.rfc.sig<-function (tissue.group.1,tissue.group.2,tax.level, mt.test) 
{
		gp1 <- tax.abundance(tissue.group.1,tax.level)
		gp2 <- tax.abundance(tissue.group.2,tax.level)
		level <- split(taxonomy_micro,taxonomy_micro[,paste(tax.level)])
	
		gp1.gp2.res <- c()
		#gp1.ga2.conf.int <- c()
	
		for (i in 1:length(level)) {
			xx.gp1<-gp1[order(gp1$pat_id),c(1,i+1)]
			yy.gp2<-gp2[order(gp2$pat_id),c(1,i+1)]
			iii.gp1<-(xx.gp1$pat_id%in%yy.gp2$pat_id)
			iii.gp2<-(yy.gp2$pat_id%in%xx.gp1$pat_id)
			xx.gp1<-xx.gp1[iii.gp1,]
			yy.gp2<-yy.gp2[iii.gp2,]
			#res<-wilcox.test(xx.gp1[,2],yy.gp2[,2],paired=T)
			res<-wilcox.test(xx.gp1[,2],yy.gp2[,2],paired=T,conf.level=0.95,conf.int=T)
			gp1.gp2.res<-append(gp1.gp2.res,res$p.value)
			#gp1.ga2.conf.int <- append(gp1.ga2.conf.int,res$conf.int[1:2])
					
		}
		
		#gp1.ga2.conf.int <- matrix(gp1.ga2.conf.int,nrow = 2)
		#colnames(gp1.ga2.conf.int) <- names(gp1[-1])
				
		names(gp1.gp2.res)<-names(gp1[-1])
	
		gp1.gp2.p.ordered<-gp1.gp2.res[order(gp1.gp2.res)]
		gp1.gp2.p.ordered<-gp1.gp2.p.ordered[!is.na(gp1.gp2.p.ordered)]
	
		yy<-mt.rawp2adjp(gp1.gp2.p.ordered)
	
		ll<-length(yy$index[yy$adjp[,mt.test]<0.05])
		if (ll == 0) {
			spe.sig = ""
			
			newList <- list("spe.sig"=spe.sig, "spe.all"=names(gp1[-1]), "gp1"=gp1, "gp2"=gp2, "gp1.gp2.p.ordered"=gp1.gp2.p.ordered, "gp1.sig"="", "gp2.sig"="", "gp1.pg2.conf.int"="")
		}
		else {
			gp1.gp2.p.ordered[1:ll]
			nn<-names(gp1.gp2.p.ordered[1:ll])
			spe.sig<-nn[order(nn)]
	
			newList <- list("spe.sig"=spe.sig, "spe.all"=names(gp1[-1]), "gp1"=gp1, "gp2"=gp2, "gp1.gp2.p.ordered"=gp1.gp2.p.ordered, "gp1.sig"=gp1[spe.sig], "gp2.sig"=gp2[spe.sig], "gp1.pg2.conf.int"=)
		}
		
		return(newList)
}

#for "tp" and "nt" groups, do the paired t-test at "genome_file_name", "genome_name", "species", "genus", "family", "order", "class", "phylum" levels
tp.nt.genome_file_name <- paired.testing.rfc.sig(tp, nt, "genome_file_name", 7)
tp.nt.genome_name <- paired.testing.rfc.sig(tp, nt, "genome_name", 7)
tp.nt.species <- paired.testing.rfc.sig(tp, nt, "species", 7)
tp.nt.genus<- paired.testing.rfc.sig(tp, nt, "genus", 7)
tp.nt.family <- paired.testing.rfc.sig(tp, nt, "family", 7)
tp.nt.order <- paired.testing.rfc.sig(tp, nt, "order", 7)
tp.nt.class <- paired.testing.rfc.sig(tp, nt, "class", 7)
tp.nt.phylum <- paired.testing.rfc.sig(tp, nt, "phylum", 7)


# do some statistic of 
clinical_tp_nt <-  clinical[clinical$bcr_patient_barcode%in%nt$pat_id,]
drug_tp_nt <- drug[drug$bcr_patient_barcode%in%nt$pat_id,]
radiation_tp_nt <- radiation[radiation$bcr_patient_barcode%in%nt$pat_id,] ### 0 rows returned 
consent_tp_nt <- consent[consent$bcr_patient_barcode%in%nt$pat_id,] 
follow_tp_nt <- follow[follow$bcr_patient_barcode%in%nt$pat_id,] 

clinical_tp <-  clinical[clinical$bcr_patient_barcode%in%tp$pat_id,]
drug_tp <- drug[drug$bcr_patient_barcode%in%tp$pat_id,]
radiation_tp <- radiation[radiation$bcr_patient_barcode%in%tp$pat_id,] 
consent_tp <- consent[consent$bcr_patient_barcode%in%tp$pat_id,] 
follow_tp <- follow[follow$bcr_patient_barcode%in%tp$pat_id,] 

write.table(clinical_tp_nt, file="clinical_tp_nt.csv", quote=FALSE, row.names=FALSE, sep=";")
write.table(drug_tp_nt, file="drug_tp_nt.csv", quote=FALSE, row.names=FALSE, sep=";")
write.table(consent_tp_nt, file="consent_tp_nt.csv", quote=FALSE, row.names=FALSE, sep=";")
write.table(follow_tp_nt, file="follow_tp_nt.csv", quote=FALSE, row.names=FALSE, sep=";")

summary(droplevels(clinical_tp_nt$age_at_initial_pathologic_diagnosis))
summary(droplevels(clinical_tp_nt$anatomic_neoplasm_subdivision))


#e.g. to have a list of the significant phylum:
#tp.nt.phylum$spe.sig

#do taxa distribution (bar charts)
#gp1: the abundance data frame of the tissue group 1
#gp2: the abundance data frame of the tissue group 2
#main1: the title of the left bar plot
#main2: the title of the right bar plot
#mainall: the main title of the graphic
#outpdf: the outpdf file name
tax.distribution.plot <- function (gp1, gp2, main1, main2, mainall, outpdf)
{
		gp1.match <- gp1[gp1$pat_id%in%gp2$pat_id,]
		gp2.match <- gp2[gp2$pat_id%in%gp1$pat_id,]
	
		gp1.match <- gp1.match[order(gp1.match$pat_id),]
		gp2.match <- gp2.match[order(gp2.match$pat_id),]
	
		gp1.matrix <- data.matrix(gp1.match[-1])
		rownames(gp1.matrix)=unname(unlist(gp1.match[1]))

		prop1 <- prop.table(gp1.matrix,margin=1)
		prop1 <- t(prop1)

		gp2.matrix <- data.matrix(gp2.match[-1])
		rownames(gp2.matrix)=unname(unlist(gp2.match[1]))

		prop2 <- prop.table(gp2.matrix,margin=1)
		prop2 <- t(prop2)
	
		pdf(outpdf,width=18,height=8)

		layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
		par(mai=c(0.5, 1, 0.5, 1), las=1)
		barplot(prop1, xlab="proportional abundance", main=main1, horiz=TRUE, col=heat.colors(length(rownames(prop1))), width=2, cex.names=0.5)
    barplot(prop2, xlab="proportional abundance", main=main2, horiz=TRUE, col=heat.colors(length(rownames(prop2))), width=2, cex.names=0.5)
	
		par(mai=c(0,0,0,0))
		plot.new()
		list.spe = rownames(prop1)
		#list.spe[list.spe=="Var.2"] <- "unknown"
		legend("center",fill=heat.colors(length(rownames(prop1))), ncol = 5, legend=list.spe, cex=1, pt.cex = 1)
	
	  mtext(mainall, cex=1)
		dev.off()
}

# do it at the phylum level
tax.distribution.plot(tp.nt.phylum$gp1,tp.nt.phylum$gp2,"tp","nt", "Composition of gastric microbiota at the phylum level.", "gc_phylum.pdf")

# do it at the class level
tax.distribution.plot(tp.nt.class$gp1,tp.nt.class$gp2,"tp","nt", "Composition of gastric microbiota at the class level.", "gc_class.pdf")

# do it at the order level
tax.distribution.plot(tp.nt.order$gp1,tp.nt.order$gp2,"tp","nt", "Composition of gastric microbiota at the order level.", "gc_order.pdf")

# do it at the family level
tax.distribution.plot(tp.nt.family$gp1,tp.nt.family$gp2,"tp","nt", "Composition of gastric microbiota at the family level.", "gc_family.pdf")

# do it at the genus level
tax.distribution.plot(tp.nt.genus$gp1,tp.nt.genus$gp2,"tp","nt", "Composition of gastric microbiota at the genus level.", "gc_genus.pdf")

# do it at the species level
tax.distribution.plot(tp.nt.species$gp1,tp.nt.species$gp2,"tp","nt", "Composition of gastric microbiota at the species level.", "gc_species.pdf")


#given a matrix, normalize each row to [-0.5,0.5]
#x: a given matrix
normalization <- function(x) {

  (x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x, na.rm=TRUE)) -0.5

}

# do box plot for the abundance difference between two groups
#gp1: the abundance data frame of the tissue group 1
#gp2: the abundance data frame of the tissue group 2
#ifNor: if the normalization of one vector is required, "T" or "F"
#ifOut: if leave the outliers in the plots, "T" or "F"
#outpdf: the outpdf file name
abundance.diff.plot <- function(gp1, gp2, ifNor, ifOut, outpdf){

		gp1.match <- gp1[gp1$pat_id%in%gp2$pat_id,]
		gp2.match <- gp2[gp2$pat_id%in%gp1$pat_id,]
	
		gp1.match <- gp1.match[order(gp1.match$pat_id),]
		gp2.match <- gp2.match[order(gp2.match$pat_id),]
		
		gp1.matrix <- data.matrix(gp1.match[-1])
		rownames(gp1.matrix)=unname(unlist(gp1.match[1]))
		#colnames(gp1.matrix)[colnames(gp1.matrix)=="Var.2"] <- "unknown"

		gp2.matrix <- data.matrix(gp2.match[-1])
		rownames(gp2.matrix)=unname(unlist(gp2.match[1]))
    #colnames(gp2.matrix)[colnames(gp2.matrix)=="Var.2"] <- "unknown"
		
		gp1.gp2.dif.matrix <- gp1.matrix - gp2.matrix
		
		#normalization 
		if (ifNor == "T") {
			gp1.gp2.dif.matrix <- apply(gp1.gp2.dif.matrix, 2, normalization)
		}
		
		pdf(outpdf, width=18,height=8)
 	  par(mar=c(12, 4, 2, 2),las=1)
		if (ifOut == "T" ) {
			#boxplot(gp1.gp2.dif.matrix, ylab ="abundance difference", outline=TRUE, cex.axis = 1, las=2)
			boxplot(asinh(gp1.gp2.dif.matrix), ylab ="abundance difference (arc-sine)", outline=TRUE, cex.axis = 1, las=2)
		}
		else if (ifOut == "F") {
			#boxplot(gp1.gp2.dif.matrix, ylab ="abundance difference", outline=FALSE, cex.axis = 1, las=2)
			boxplot(asinh(gp1.gp2.dif.matrix), ylab ="abundance difference (arc-sine)", outline=FALSE, cex.axis = 1, las=2)
		}
		
		#mtext("Phylums", side = 1, line = 8)
		
		dev.off()

}

#print the boxplot at phylum level, without normalization, with outliers
abundance.diff.plot(tp.nt.phylum$gp1, tp.nt.phylum$gp2, "F", "T", "boxplot_no_normalization_with_outliers.pdf")

#print the boxplot at phylum level, without normalization, removing outliers
abundance.diff.plot(tp.nt.phylum$gp1, tp.nt.phylum$gp2, "F", "F", "boxplot_no_normalization_no_outliers.pdf")

#print the boxplot at phylum level, with normalization, with outliers
abundance.diff.plot(tp.nt.phylum$gp1, tp.nt.phylum$gp2, "T", "T", "boxplot_with_normalization_with_outliers.pdf")

#print the boxplot at phylum level, with normalization, removing outliers
abundance.diff.plot(tp.nt.phylum$gp1, tp.nt.phylum$gp2, "T", "F", "boxplot_with_normalization_no_outliers.pdf")

# return the matrix of abundance difference
#gp1: the abundance data frame of the tissue group 1
#gp2: the abundance data frame of the tissue group 2
abundance.diff <- function(gp1, gp2) {

		gp1.match <- gp1[gp1$pat_id%in%gp2$pat_id,]
		gp2.match <- gp2[gp2$pat_id%in%gp1$pat_id,]
	
		gp1.match <- gp1.match[order(gp1.match$pat_id),]
		gp2.match <- gp2.match[order(gp2.match$pat_id),]
		
		gp1.matrix <- data.matrix(gp1.match[-1])
		rownames(gp1.matrix)=unname(unlist(gp1.match[1]))
		#colnames(gp1.matrix)[colnames(gp1.matrix)=="Var.2"] <- "unknown"

		gp2.matrix <- data.matrix(gp2.match[-1])
		rownames(gp2.matrix)=unname(unlist(gp2.match[1]))
    #colnames(gp2.matrix)[colnames(gp2.matrix)=="Var.2"] <- "unknown"
		
		gp1.gp2.dif.matrix <- gp1.matrix - gp2.matrix
		
		return(gp1.gp2.dif.matrix)
}

#get the matrix with the abundance difference at "genome_file_name", "genome_name", "species", "genus", "family", "order", "class", "phylum" level
genome_file_name.diff.matrix  <- abundance.diff(tp.nt.genome_file_name$gp1, tp.nt.genome_file_name$gp2)
genome_name.diff.matrix <- abundance.diff(tp.nt.genome_name$gp1, tp.nt.genome_name$gp2)
species.diff.matrix <- abundance.diff(tp.nt.species$gp1, tp.nt.species$gp2)
genus.diff.matrix <- abundance.diff(tp.nt.genus$gp1, tp.nt.genus$gp2)
family.diff.matrix <- abundance.diff(tp.nt.family$gp1, tp.nt.family$gp2)
order.diff.matrix <- abundance.diff(tp.nt.order$gp1, tp.nt.order$gp2)
class.diff.matrix <- abundance.diff(tp.nt.class$gp1, tp.nt.class$gp2)
phylum.diff.matrix <- abundance.diff(tp.nt.phylum$gp1, tp.nt.phylum$gp2)


#create the data frame with columsn: pid, NT abundance, TP abundance, abundance difference, clinical data (each have a column)
#gp1: the abundance data frame of the tissue group 1
#gp2: the abundance data frame of the tissue group 2
#spe_name: the name of the species at corresponding level
gp1=tp.nt.species$gp1
gp2=tp.nt.species$gp2
spe_name="Helicobacter pylori"
clinical.prep <- function(gp1, gp2, spe_name) {

		gp1 <- gp1[,c("pat_id", spe_name)]
		gp2 <- gp2[,c("pat_id", spe_name)]

		#for the samples pairs with both from tp and nt
		gp1.match <- gp1[gp1$pat_id%in%gp2$pat_id,]
		gp2.match <- gp2[gp2$pat_id%in%gp1$pat_id,]
	
		gp1.match <- gp1.match[order(gp1.match$pat_id),]
		gp2.match <- gp2.match[order(gp2.match$pat_id),]

		#for the samples with only from tp
		gp1.not.match <- gp1[!(gp1$pat_id%in%gp2$pat_id),]
		gp1.not.match <- gp1.not.match[order(gp1.not.match$pat_id),]
		
		#for the samples with only from nt
		gp2.not.match <- gp2[!(gp2$pat_id%in%gp1$pat_id),]
		gp2.not.match <- gp2.not.match[order(gp2.not.match$pat_id),]
		
		#create the dataframe
		df.tmp.match <- data.frame("pat_id"=gp1.match$pat_id, "nt"=gp2.match[[spe_name]], "tp"=gp1.match[[spe_name]],"diff"=(gp1.match[[spe_name]] - gp2.match[[spe_name]]))
		df.tmp.not.match.1 <- data.frame("pat_id"=gp1.not.match$pat_id, "nt"=rep(NA,length(gp1.not.match$pat_id)), "tp"=gp1.not.match[[spe_name]],"diff"=rep(NA,length(gp1.not.match$pat_id)))
		df.tmp.not.match.2 <- data.frame("pat_id"=gp2.not.match$pat_id, "nt"=gp2.not.match[[spe_name]], "tp"=rep(NA,length(gp2.not.match$pat_id)),"diff"=rep(NA,length(gp2.not.match$pat_id)))
		
		df.tmp <- rbind(df.tmp.match, df.tmp.not.match.1, df.tmp.not.match.2)
			
		#match the clinical data 
		df.tmp.match <- df.tmp[df.tmp$pat_id%in%clinical$bcr_patient_barcode,]		
		clinical.match <- clinical[clinical$bcr_patient_barcode%in%df.tmp$pat_id,]
		
		clinical.match[clinical.match=="[Not Available]"] <- NA
		clinical.match[clinical.match=="[Not Applicable]"] <- NA
		clinical.match[clinical.match=="[Not Available]|[Not Available]"] <- NA

		clinical.match <- clinical.match[match(df.tmp.match$pat_id,clinical.match$bcr_patient_barcode),]

		df.final.tmp <- data.frame(df.tmp.match, clinical.match[c(6,7,9,12,14,24,26,31,32,34,35,37,38,39,41,42,47,49,50,51,52,53,56,57,60,63,65,67,69,70,73,74,75,77,78)])
		
		write.table(df.final.tmp, file='abundance_clinical.tsv', quote=FALSE, sep='\t', row.names = FALSE)
		
		#match the consent data
		df.tmp.match <- df.tmp[df.tmp$pat_id%in%consent$bcr_patient_barcode,]
		consent.match <- consent[consent$bcr_patient_barcode%in%df.tmp$pat_id,]
		
		consent.match[consent.match=="[Not Available]"] <- NA
		consent.match[consent.match=="[Not Applicable]"] <- NA
		consent.match[consent.match=="[Not Available]|[Not Available]"] <- NA
		
		df.final.consent <- data.frame(df.tmp.match, consent.match[c(3,4,6,7,8,9,10,11,13,15)])
		
		write.table(df.final.consent, file='abundance_consent.tsv', quote=FALSE, sep='\t', row.names = FALSE)
		
		#match the radiation data
		df.tmp.match <- df.tmp[df.tmp$pat_id%in%radiation$bcr_patient_barcode,]
		radiation.match <- radiation[radiation$bcr_patient_barcode%in%df.tmp$pat_id,]
		
		radiation.match[radiation.match=="[Not Available]"] <- NA
		radiation.match[radiation.match=="[Not Applicable]"] <- NA
		radiation.match[radiation.match=="[Not Available]|[Not Available]"] <- NA
		
		df.final.consent <- data.frame(df.tmp.match, radiation.match[c(3,6,7,8,9,10,11,12,13,15,16,17)])
		
		write.table(df.final.consent, file='abundance_radiation.tsv', quote=FALSE, sep='\t', row.names = FALSE)		
		
		#match the drug data
		drug.match <- drug[drug$bcr_patient_barcode%in%df.tmp$pat_id,]
		
		drug.match[drug.match=="[Not Available]"] <- NA
		drug.match[drug.match=="[Not Applicable]"] <- NA
		drug.match[drug.match=="[Not Available]|[Not Available]"] <- NA
		drug.match[drug.match=="[Unknown]"] <- NA
		
		df.final.drug <- data.frame("pat_id"=drug.match$bcr_patient_barcode, drug.match[c(2,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22)])
		
		#df.final.drug <- merge(df.tmp, df.final.drug, all = TRUE)
		df.final.drug <- merge(df.tmp, df.final.drug)
		df.final.drug <- df.final.drug[with(df.final.drug, order(pat_id)), ]
		
		write.table(df.final.drug, file='abundance_drug.tsv', quote=FALSE, sep='\t', row.names = FALSE)
		
		#match the followup data
		follow.match <- follow[follow$bcr_patient_barcode%in%df.tmp$pat_id,]
		
		follow.match[follow.match=="[Not Available]"] <- NA
		follow.match[follow.match=="[Not Applicable]"] <- NA
		follow.match[follow.match=="[Not Available]|[Not Available]"] <- NA
		follow.match[follow.match=="[Unknown]"] <- NA
		
		df.final.follow <- data.frame("pat_id"=follow.match$bcr_patient_barcode, follow.match[c(2,3,4,5,6,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,27,28)])
		
		df.final.follow <- merge(df.tmp, df.final.follow)
		df.final.follow <- df.final.follow[with(df.final.follow, order(pat_id)), ]
		
		write.table(df.final.follow, file='abundance_follow.tsv', quote=FALSE, sep='\t', row.names = FALSE)
		
}

df_clinical <- read.csv("abundance_clinical.tsv", header=T, sep="\t")
df_consent <- read.csv("abundance_consent.tsv", header=T, sep="\t")
df_radiation <- read.csv("abundance_radiation.tsv", header=T, sep="\t")
df_drug <- read.csv("abundance_drug.tsv", header=T, sep="\t")
df_follow <- read.csv("abundance_follow.tsv", header=T, sep="\t")

#in clinical table, include the following columns (using factor() to remove levels):
age_at_initial_pathologic_diagnosis 
anatomic_neoplasm_subdivision
antireflux_treatment
barretts_esophagus
city_of_procurement
days_to_death
days_to_last_followup
family_history_of_stomach_cancer
gender
h_pylori_infection
histological_type
icd_10
icd_o_3_histology
icd_o_3_site
lymph_node_examined_count
neoplasm_histologic_grade
number_of_lymphnodes_positive_by_he
number_of_relatives_with_stomach_cancer
pathologic_M
pathologic_N
pathologic_N
pathologic_stage
neoplasm_cancer_status
primary_lymph_node_presentation_assessment
prior_dx
race
reflux_history
residual_tumor
state_province_country_of_procurement
system_version
tissue_prospective_collection_indicator
tissue_retrospective_collection_indicator
tissue_source_site
tumor_tissue_site
vital_status
year_of_initial_pathologic_diagnosis


# a function to create dataframe with the first column "pat_id" and the other columns=log(abundance(nt)/abundance(tp)), for the samples with tp and nt available.
#gp1: the abundance data frame of the tissue group 1
#gp2: the abundance data frame of the tissue group 2
abundance.diff.log <- function(gp1, gp2) {

		#for the samples pairs with both from tp and nt
		gp1.match <- gp1[gp1$pat_id%in%gp2$pat_id,]
		gp2.match <- gp2[gp2$pat_id%in%gp1$pat_id,]
	
		gp1.match <- gp1.match[order(gp1.match$pat_id),]
		gp2.match <- gp2.match[order(gp2.match$pat_id),]
		
		#for the samples with only from tp
		gp1.not.match <- gp1[!(gp1$pat_id%in%gp2$pat_id),]
		gp1.not.match <- gp1.not.match[order(gp1.not.match$pat_id),]
		
		#for the samples with only from nt
		gp2.not.match <- gp2[!(gp2$pat_id%in%gp1$pat_id),]
		gp2.not.match <- gp2.not.match[order(gp2.not.match$pat_id),]
		
		#create the dataframe
		df.tmp.match <- data.frame("pat_id"=gp1.match$pat_id, log(gp2.match[,-1]/gp1.match[,-1]))
		df.tmp.not.match.1 <- data.frame("pat_id"=gp1.not.match$pat_id, log(NA/gp1.not.match[,-1]))
		df.tmp.not.match.2 <- data.frame("pat_id"=gp2.not.match$pat_id, log(gp2.not.match[,-1]/NA))
		
		df.tmp <- rbind(df.tmp.match, df.tmp.not.match.1, df.tmp.not.match.2)
		
    return(df.tmp)
}

# get the dataframe with log(abundance(nt)/abundance(tp)) for each species level:
genome_file_name.diff.log  <- abundance.diff.log(tp.nt.genome_file_name$gp1, tp.nt.genome_file_name$gp2)
genome_name.diff.log <- abundance.diff.log(tp.nt.genome_name$gp1, tp.nt.genome_name$gp2)
species.diff.log <- abundance.diff.log(tp.nt.species$gp1, tp.nt.species$gp2)
genus.diff.log <- abundance.diff.log(tp.nt.genus$gp1, tp.nt.genus$gp2)
family.diff.log <- abundance.diff.log(tp.nt.family$gp1, tp.nt.family$gp2)
order.diff.log <- abundance.diff.log(tp.nt.order$gp1, tp.nt.order$gp2)
class.diff.log <- abundance.diff.log(tp.nt.class$gp1, tp.nt.class$gp2)
phylum.diff.log <- abundance.diff.log(tp.nt.phylum$gp1, tp.nt.phylum$gp2)

# write to file for the files above
write.table(genome_file_name.diff.log, file='genome_file_name.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(genome_name.diff.log, file='genome_name.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(species.diff.log, file='species.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(genus.diff.log, file='genus.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(family.diff.log, file='family.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(order.diff.log, file='order.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(class.diff.log, file='class.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(phylum.diff.log, file='phylum.diff.log.tsv', quote=FALSE, sep='\t', row.names = FALSE)

# get the matched clinical data for each pat_id, output the dataframe with the same structure as the log files above, for those without data, leave NA
 
clinical.list = list()
consent.list = list()
radiation.list = list()
drug.list = list()
follow.list = list()
ind = 1
for(pat_id in phylum.diff.log$pat_id){

	#clinical data
	if(is.element(pat_id, clinical$bcr_patient_barcode)) {
		tmp.df <- data.frame("pat_id"=pat_id, clinical[clinical$bcr_patient_barcode==pat_id,c(6,7,9,12,14,24,26,31,32,34,35,37,38,39,41,42,47,49,50,51,52,53,56,57,60,63,65,67,69,70,73,74,75,77,78)])
		tmp.df[tmp.df=="[Not Available]"] <- NA
		tmp.df[tmp.df=="[Not Applicable]"] <- NA
		tmp.df[tmp.df=="[Not Available]|[Not Available]"] <- NA
		tmp.df[tmp.df=="[Unknown]"] <- NA
		clinical.list[[ind]] <- tmp.df
	}
	else {
		tmp.df.rep.na <- as.data.frame(matrix(0, ncol = 36, nrow = 1))
		colnames(tmp.df.rep.na) <- c("pat_id", colnames(clinical[,c(6,7,9,12,14,24,26,31,32,34,35,37,38,39,41,42,47,49,50,51,52,53,56,57,60,63,65,67,69,70,73,74,75,77,78)]))
		tmp.df.rep.na[1,] <- c(pat_id, rep(NA, 35))
		clinical.list[[ind]] <- tmp.df.rep.na
	}
	
	#consent data
	if(is.element(pat_id, consent$bcr_patient_barcode)) {
		tmp.df <- data.frame("pat_id"=pat_id, consent[consent$bcr_patient_barcode==pat_id,c(3,4,6,7,8,9,10,11,13,15)])
		tmp.df[tmp.df=="[Not Available]"] <- NA
		tmp.df[tmp.df=="[Not Applicable]"] <- NA
		tmp.df[tmp.df=="[Not Available]|[Not Available]"] <- NA
		tmp.df[tmp.df=="[Unknown]"] <- NA
		consent.list[[ind]] <- tmp.df
	}
	else {
		tmp.df.rep.na <- as.data.frame(matrix(0, ncol = 11, nrow = 1))
		colnames(tmp.df.rep.na) <- c("pat_id", colnames(consent[,c(3,4,6,7,8,9,10,11,13,15)]))
		tmp.df.rep.na[1,] <- c(pat_id, rep(NA, 10))
		consent.list[[ind]] <- tmp.df.rep.na
	}
	
	radiation data
	 if(is.element(pat_id, radiation$bcr_patient_barcode)) {
		 tmp.df <- data.frame("pat_id"=pat_id, radiation[radiation$bcr_patient_barcode==pat_id,c(3,6,7,8,9,10,11,12,13,15,16,17)])
		 tmp.df[tmp.df=="[Not Available]"] <- NA
		 tmp.df[tmp.df=="[Not Applicable]"] <- NA
		 tmp.df[tmp.df=="[Not Available]|[Not Available]"] <- NA
		 tmp.df[tmp.df=="[Unknown]"] <- NA
		 radiation.list[[ind]] <- tmp.df
	 }
	 else {
		 tmp.df.rep.na <- as.data.frame(matrix(0, ncol = 13, nrow = 1))
		 colnames(tmp.df.rep.na) <- c("pat_id", colnames(radiation[,c(3,6,7,8,9,10,11,12,13,15,16,17)]))
		 tmp.df.rep.na[1,] <- c(pat_id, rep(NA, 12))
		 radiation.list[[ind]] <- tmp.df.rep.na
	 }
	
	#drug data, one patient can have more than one rows
	if(is.element(pat_id, drug$bcr_patient_barcode)) {
		tmp.df <- data.frame("pat_id"=pat_id, drug[drug$bcr_patient_barcode==pat_id,c(2,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22)])
		tmp.df[tmp.df=="[Not Available]"] <- NA
		tmp.df[tmp.df=="[Not Applicable]"] <- NA
		tmp.df[tmp.df=="[Not Available]|[Not Available]"] <- NA
		tmp.df[tmp.df=="[Unknown]"] <- NA
		drug.list[[ind]] <- tmp.df
	}
	else {
		tmp.df.rep.na <- as.data.frame(matrix(0, ncol = 19, nrow = 1))
		colnames(tmp.df.rep.na) <- c("pat_id", colnames(drug[,c(2,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22)]))
		tmp.df.rep.na[1,] <- c(pat_id, rep(NA, 18))
		drug.list[[ind]] <- tmp.df.rep.na
	}
	
	#followup data, one patient can have more than one rows
	if(is.element(pat_id, follow$bcr_patient_barcode)) {
		tmp.df <- data.frame("pat_id"=pat_id, follow[follow$bcr_patient_barcode==pat_id,c(2,3,4,5,6,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,27,28)])
		tmp.df[tmp.df=="[Not Available]"] <- NA
		tmp.df[tmp.df=="[Not Applicable]"] <- NA
		tmp.df[tmp.df=="[Not Available]|[Not Available]"] <- NA
		tmp.df[tmp.df=="[Unknown]"] <- NA
		follow.list[[ind]] <- tmp.df
	}
	else {
		tmp.df.rep.na <- as.data.frame(matrix(0, ncol = 24, nrow = 1))
		colnames(tmp.df.rep.na) <- c("pat_id", colnames(follow[,c(2,3,4,5,6,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,27,28)]))
		tmp.df.rep.na[1,] <- c(pat_id, rep(NA, 23))
		follow.list[[ind]] <- tmp.df.rep.na
	}

	ind = ind + 1
}

clinical.dataframe <- do.call(rbind, clinical.list)
consent.dataframe <- do.call(rbind, consent.list)
radiation.dataframe <- do.call(rbind, radiation.list)
follow.dataframe <- do.call(rbind, follow.list)
drug.dataframe <- do.call(rbind, drug.list)
 
# write to file for the files above
write.table(clinical.dataframe, file='clinical.df.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(consent.dataframe, file='consent.df.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(radiation.dataframe, file='radiation.df.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(follow.dataframe, file='follow.df.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(drug.dataframe, file='drug.df.tsv', quote=FALSE, sep='\t', row.names = FALSE)

# output the file for circus histogram plotting 
output_circus_file <- function(outfile){

	diff.matrix <- abundance.diff(tp.nt.genome_file_name$gp1, tp.nt.genome_file_name$gp2)
	tax.level="genome_file_name" 
	
	df <- data.frame(diff.matrix)
	rownames(df) <- rownames(diff.matrix)
	colnames(df) <- colnames(diff.matrix)
	 
	df <- df[taxonomy_micro[,paste(tax.level)]]
	diff.matrix <- data.matrix(df)
	
	v.mean <- apply(diff.matrix, 2, mean)

	df.new <- cbind(taxonomy_micro,v.mean)
	
	# go through the new dataframe with abundance difference, calculate for each taxonomy level, save in a new matrix "out.matrix"
	out.matrix <- matrix(,nrow=16045,ncol=7)
	
	# simpy do 7 loops, i.e. go through 7 columns, read line by line, once taxonomy ssignment changed, average abundance was given
  for(i in 3:9) {
	
		# initialize with the value from 1st line
		cur_tax = df.new[1,i]
		abu_tax = df.new[1,10]
		n_tax = 1
		
		for(j in 2:(length(df.new[,1])-1)) {
			
			if(df.new[j,i] == cur_tax) {
			
				abu_tax = abu_tax + df.new[j,10]
				n_tax = n_tax + 1
			}
			else {
				
				out.matrix[(j-n_tax):(j-1), (i-2)] = abu_tax/n_tax
				cur_tax = df.new[j,i]
				abu_tax = df.new[j,10]
				n_tax = 1
				
			}
		}
		
		# do for the last line
		j = length(df.new[,1])
		if(df.new[j,i] == cur_tax) {
			
			abu_tax = abu_tax + df.new[j,10]
			n_tax = n_tax + 1
			out.matrix[(j-n_tax+1):j, (i-2)] = abu_tax/n_tax
			
		}
		else {
			
			out.matrix[(j-n_tax):(j-1), (i-2)] = abu_tax/n_tax
			cur_tax = df.new[j,i]
			abu_tax = df.new[j,10]
			n_tax = 1
				
			out.matrix[j, (i-2)] = abu_tax/n_tax
		}
	}
	
	# do it in another way, let the abundance equal for the same species
	out.eq.matrix <- matrix(,nrow=16045,ncol=7)
		
	# go through df.new, for each level tax. assign the corresponding abundance in out.eq.matrix
	for (i in 3:9) {
		
		cur_tax_level = colnames(df.new)[i]
		
    # get the corresponding variable name		
		gp1 = eval(parse(text=paste("tp.nt.", cur_tax_level, "$gp1", sep="")))
		
		gp2 = eval(parse(text=paste("tp.nt.", cur_tax_level, "$gp2", sep="")))
		
		cur_diff.matrix = abundance.diff(gp1, gp2)
	
		cur_v.mean = apply(cur_diff.matrix, 2, mean)
		
		for (j in 1:length(df.new[,1])) {
			
			out.eq.matrix[j,i-2] = 	cur_v.mean[df.new[j,i]]
		}
	}
	
	# write df.new, out.matrix and out.eq.matrix to file
	write.table(df.new, file="for_circos/df.new", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	write.table(out.matrix, file="for_circos/out.matrix", sep="\t", col.names=FALSE, row.names=FALSE)
	write.table(out.eq.matrix, file="for_circos/out.eq.matrix", sep="\t", col.names=FALSE, row.names=FALSE)
	
}




# output the file for the matrix, no use
output_matrix_file <- function(diff.matrix, tax.level, outfile){
 
  
	v.mean <- apply(diff.matrix, 2, mean)
	
	ma.mean <- cbind(as.vector(taxonomy_micro[match(names(v.mean), taxonomy_micro[, paste(tax.level)]),]$phylum),names(v.mean),v.mean)
	colnames(ma.mean) <- c("phylum", tax.level, "mean")
	
	# reorder the column at given taxonomy level, and 
	#ma.mean <- ma.mean[match(taxonomy_micro[, paste(tax.level)], ma.mean[, paste(tax.level)]),]
	
	write.table(ma.mean, file = outfile, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	
}

#output the matrix at in a file at "genome_file_name", "genome_name", "species", "genus", "family", "order", "class", "phylum" level
output_matrix_file(genome_file_name.diff.matrix, "genome_file_name", "test_genome_file_name.txt")
output_matrix_file(genome_name.diff.matrix, "genome_name", "test_genome_name.txt")
output_matrix_file(species.diff.matrix, "species", "test_species.txt")
output_matrix_file(genus.diff.matrix, "genus", "test_genus.txt")
output_matrix_file(family.diff.matrix, "family", "test_family.txt")
output_matrix_file(order.diff.matrix, "order", "test_order.txt")
output_matrix_file(class.diff.matrix, "class", "test_class.txt")
output_matrix_file(phylum.diff.matrix, "phylum", "test_phylum.txt")


# given a vector, claculate the 95% confidence interval
# x: a given vector
confidencen.interval <- function(x) {

	error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
	
	left <- mean(x) - error
	
	right <- mean(x) + error
	
	left <- round(left, 3)
	
	right <- round(right, 3)

	return(paste("(", left, ", ", right, ")", sep=""))
}

# given a vector, claculate the 95% confidence interval, uspng wilcox testing
# x: a given vector
confidencen.interval.wil <- function(x) {

  res<-wilcox.test(x,conf.level=0.95,conf.int=T)

	#left <- round(res$conf.int[1], 5)
	#right <- round(res$conf.int[2], 5)
	
	#use asinh() to transform
	left <- round(asinh(res$conf.int[1]), 3)
	right <- round(asinh(res$conf.int[2]), 3)
	
	return(paste("(", left, ", ", right, ")", sep=""))

}

# given a vector, claculate median, uspng wilcox testing
# x: a given vector
median.wil <- function(x) {

	res<-wilcox.test(x,conf.level=0.95,conf.int=T)
	
	#median.round <- round(res$estimate, 5)
	
	#use asinh() to transform
	median.round <- round(asinh(res$estimate), 3)
	
	return(median.round)
}




# output the statistic results at "phylum" level in latex format
# use package "xtable" to output latex format, more info. about package "xtable" see: https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
# install "xtable" package if not done
# install.packages("xtable")
library(xtable)

# build a new matrix containing the columns we need for xtable
#gp1: the abundance data frame of the tissue group 1
#gp2: the abundance data frame of the tissue group 2
print.table <- function(gp1, gp2){
	gp1 <- tp.nt.phylum$gp1
	gp2 <- tp.nt.phylum$gp2

	gp1.match <- gp1[gp1$pat_id%in%gp2$pat_id,]
	gp2.match <- gp2[gp2$pat_id%in%gp1$pat_id,]

	gp1.match <- gp1.match[order(gp1.match$pat_id),]
	gp2.match <- gp2.match[order(gp2.match$pat_id),]

	#cm1 <- colMeans(gp1.match[-1])
	#cm2 <- colMeans(gp2.match[-1])
	#cm1 <- round(cm1[names(tp.nt.phylum$gp1.gp2.p.ordered)], 3)
	#cm2 <- round(cm2[names(tp.nt.phylum$gp1.gp2.p.ordered)], 3)

	gp1.matrix <- data.matrix(gp1.match[-1])
	rownames(gp1.matrix)=unname(unlist(gp1.match[1]))

	gp2.matrix <- data.matrix(gp2.match[-1])
	rownames(gp2.matrix)=unname(unlist(gp2.match[1]))
		
	gp1.gp2.dif.matrix <- gp1.matrix - gp2.matrix

	mean.gp1.gp2.dif <- colMeans(gp1.gp2.dif.matrix)
	mean.gp1.gp2.dif <- round(mean.gp1.gp2.dif[names(tp.nt.phylum$gp1.gp2.p.ordered)], 3)
	
	median.gp1.gp2.dif <- apply(gp1.gp2.dif.matrix, 2, median.wil)
	median.gp1.gp2.dif <- median.gp1.gp2.dif[names(tp.nt.phylum$gp1.gp2.p.ordered)]

	#conf.inter <- apply(gp1.gp2.dif.matrix, 2, confidencen.interval)
  conf.inter <- apply(gp1.gp2.dif.matrix, 2, confidencen.interval.wil)
	conf.inter <- conf.inter[names(tp.nt.phylum$gp1.gp2.p.ordered)]

	yy<-mt.rawp2adjp(tp.nt.phylum$gp1.gp2.p.ordered)

	paired.t.p <- round(tp.nt.phylum$gp1.gp2.p.ordered,3) 

	adjust.fdr.bh <- round(yy$adjp[,6], 3)
	names(adjust.fdr.bh) <- names(tp.nt.phylum$gp1.gp2.p.ordered)

	adjust.fdr.by <- round(yy$adjp[,7], 3)
	names(adjust.fdr.by) <- names(tp.nt.phylum$gp1.gp2.p.ordered)

	names.micro <- names(tp.nt.phylum$gp1.gp2.p.ordered)
	#names.micro[names.micro=="Var.2"] <- "unknown"

	#new_mat <- cbind("phylum name"=names.micro, "mean of tp"=cm1, "mean of nt"=cm2, "paired t-test p value"=paired.t.p, "adjust FDR using BH"=adjust.fdr.bh, "adjust FDR using BY"=adjust.fdr.by)
	#new_mat <- cbind("phylum name"=names.micro, "mean"=mean.gp1.gp2.dif, "confidence interval"=conf.inter, "paired t-test p"=paired.t.p, "adj-FDR BH"=adjust.fdr.bh, "adj-FDR BY"=adjust.fdr.by)
	new_mat <- cbind("phylum name"=names.micro, "median"=median.gp1.gp2.dif, "confidence interval"=conf.inter, "paired t-test p"=paired.t.p, "adj-FDR BH"=adjust.fdr.bh, "adj-FDR BY"=adjust.fdr.by)

	#reorder them alphabetic
	new_mat <- cbind("phylum name"=sort(names.micro), "median"=median.gp1.gp2.dif[sort(names(median.gp1.gp2.dif))] , "confidence interval"=conf.inter[sort(names(conf.inter))], "paired t-test p"=paired.t.p[sort(names(paired.t.p))], "adj-FDR BH"=adjust.fdr.bh[sort(names(adjust.fdr.bh))], "adj-FDR BY"=adjust.fdr.by[sort(names(adjust.fdr.by))])
	
	print(xtable(new_mat, caption = 'The microbes abundance comparison between the group of tumor tissue and non adjacent non tumor tissue at phylum level. The comparison was done by using paired t-test. The median and the confidence interval were recalculated using arc-sine.'), include.rownames=FALSE)
}

#use this to print out without science notation
options(scipen=999)

#print out the table in latex format
print.table(tp.nt.phylum$gp1, tp.nt.phylum$gp2)



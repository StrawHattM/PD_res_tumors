#modified lollipopPLot to plot lollipop plots with lowercase/titlecase names

lollipopPlot_fixed<- function (maf, gene = NULL, AACol = NULL, labelPos = NULL, labPosSize = 0.9, 
          showMutationRate = TRUE, showDomainLabel = TRUE, cBioPortal = FALSE, 
          refSeqID = NULL, proteinID = NULL, roundedRect = TRUE, repel = FALSE, 
          collapsePosLabel = TRUE, showLegend = TRUE, legendTxtSize = 0.8, 
          labPosAngle = 0, domainLabelSize = 0.8, axisTextSize = c(1, 
                                                                   1), printCount = FALSE, colors = NULL, domainAlpha = 1, 
          domainBorderCol = "black", bgBorderCol = "black", 
          labelOnlyUniqueDoamins = TRUE, defaultYaxis = FALSE, titleSize = c(1.2, 
                                                                             1), pointSize = 1.5) {
  if (is.null(gene)) {
    stop("Please provide a gene name.")
  }
  geneID = gene
  gff = system.file("extdata", "protein_domains.RDs", 
                    package = "maftools")
  gff = readRDS(file = gff)
  data.table::setDT(x = gff)
  mut = subsetMaf(maf = maf, includeSyn = FALSE, genes = gene, 
                  query = "Variant_Type != 'CNV'", mafObj = FALSE)
  if (is.null(AACol)) {
    pchange = c("HGVSp_Short", "Protein_Change", 
                "AAChange")
    if (length(pchange[pchange %in% colnames(mut)]) > 0) {
      pchange = suppressWarnings(pchange[pchange %in% colnames(mut)][1])
      message(paste0("Assuming protein change information are stored under column ", 
                     pchange, ". Use argument AACol to override if necessary."))
      colnames(mut)[which(colnames(mut) == pchange)] = "AAChange_"
    }
    else {
      message("Available fields:")
      print(colnames(mut))
      stop("AAChange field not found in MAF. Use argument AACol to manually specifiy field name containing protein changes.")
    }
  }
  else {
    if (length(which(colnames(mut) == AACol)) == 0) {
      message("Available fields:")
      print(colnames(mut))
      stop(paste0("Column ", AACol, " not found."))
    }
    else {
      colnames(mut)[which(colnames(mut) == AACol)] = "AAChange_"
    }
  }
  prot.dat = mut[Hugo_Symbol %in% geneID, .(Variant_Type, Variant_Classification, 
                                            AAChange_)]
  if (nrow(prot.dat) == 0) {
    stop(paste(geneID, "does not seem to have any mutations!", 
               sep = " "))
  }
  prot = gff[HGNC %in% str_to_upper(geneID)]
  if (nrow(prot) == 0) {
    stop(paste("Structure for protein", geneID, "not found.", 
               sep = " "))
  }
  if (!is.null(refSeqID)) {
    prot = prot[refseq.ID == refSeqID]
    if (nrow(prot) == 0) {
      stop(paste0(refSeqID, " not found!"))
    }
  }
  else if (!is.null(proteinID)) {
    prot = prot[protein.ID == proteinID]
    if (nrow(prot) == 0) {
      stop(paste0(refSeqID, " not found!"))
    }
  }
  else {
    txs = unique(prot$refseq.ID)
    if (length(txs) > 1) {
      message(paste(length(txs), " transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.", 
                    sep = ""))
      print(prot[!duplicated(protein.ID), .(HGNC, refseq.ID, 
                                            protein.ID, aa.length)])
      prot = prot[which(prot$aa.length == max(prot$aa.length)), 
      ]
      if (length(unique(prot$refseq.ID)) > 1) {
        prot = prot[which(prot$refseq.ID == unique(prot[, 
                                                        refseq.ID])[1]), ]
        message(paste("Using longer transcript", 
                      unique(prot[, refseq.ID])[1], "for now.", 
                      sep = " "))
      }
      else {
        message(paste("Using longer transcript", 
                      unique(prot[, refseq.ID])[1], "for now.", 
                      sep = " "))
      }
    }
  }
  len = as.numeric(max(prot$aa.length, na.rm = TRUE))
  prot = prot[!is.na(Label)]
  prot = prot[, `:=`(domain_lenght, End - Start)][order(domain_lenght, 
                                                        decreasing = TRUE)][, `:=`(domain_lenght, NULL)]
  sampleSize = as.numeric(maf@summary[ID %in% "Samples", 
                                      summary])
  mutRate = round(getGeneSummary(x = maf)[Hugo_Symbol %in% 
                                            geneID, MutatedSamples]/sampleSize * 100, digits = 2)
  cbioSubTitle = geneID
  if (showMutationRate) {
    cbioSubTitle = substitute(paste(italic(cbioSubTitle), 
                                    " : [Somatic Mutation Rate: ", mutRate, "%]"))
  }
  if (cBioPortal) {
    vc = c("Nonstop_Mutation", "Frame_Shift_Del", 
           "Missense_Mutation", "Nonsense_Mutation", 
           "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del", 
           "In_Frame_Ins")
    vc.cbio = c("Truncating", "Truncating", "Missense", 
                "Truncating", "Truncating", "Truncating", 
                "In-frame", "In-frame")
    names(vc.cbio) = vc
    col = grDevices::adjustcolor(col = c("black", "#33A02C", 
                                         "brown"), alpha.f = 0.7)
    col = c(Truncating = col[1], Missense = col[2], `In-frame` = col[3])
  }
  else {
    if (is.null(colors)) {
      col = get_vcColors(alpha = 0.7, named = TRUE)
    }
    else {
      col = colors
    }
  }
  prot.spl = strsplit(x = as.character(prot.dat$AAChange_), 
                      split = ".", fixed = TRUE)
  prot.conv = sapply(sapply(prot.spl, function(x) x[length(x)]), 
                     "[", 1)
  prot.dat[, `:=`(conv, prot.conv)]
  pos = gsub(pattern = "Ter.*", replacement = "", 
             x = prot.dat$conv)
  pos = gsub(pattern = "[[:alpha:]]", replacement = "", 
             x = pos)
  pos = gsub(pattern = "\\*$", replacement = "", 
             x = pos)
  pos = gsub(pattern = "^\\*", replacement = "", 
             x = pos)
  pos = gsub(pattern = "\\*.*", replacement = "", 
             x = pos)
  pos = as.numeric(sapply(X = strsplit(x = pos, split = "_", 
                                       fixed = TRUE), FUN = function(x) x[1]))
  prot.dat[, `:=`(pos, abs(pos))]
  if (nrow(prot.dat[is.na(pos)]) > 0) {
    message(paste("Removed", nrow(prot.dat[is.na(prot.dat$pos), 
    ]), "mutations for which AA position was not available", 
    sep = " "))
    prot.dat = prot.dat[!is.na(pos)]
  }
  prot.snp.sumamry = prot.dat[, .N, .(Variant_Classification, 
                                      conv, pos)]
  colnames(prot.snp.sumamry)[ncol(prot.snp.sumamry)] = "count"
  maxCount = max(prot.snp.sumamry$count, na.rm = TRUE)
  prot.snp.sumamry = prot.snp.sumamry[order(pos), ]
  if (cBioPortal) {
    prot.snp.sumamry$Variant_Classification = vc.cbio[as.character(prot.snp.sumamry$Variant_Classification)]
  }
  if (maxCount <= 5) {
    prot.snp.sumamry$count2 = 1 + prot.snp.sumamry$count
    lim.pos = 2:6
    lim.lab = 1:5
  }
  else {
    prot.snp.sumamry$count2 = 1 + (prot.snp.sumamry$count * 
                                     (5/max(prot.snp.sumamry$count)))
    lim.pos = prot.snp.sumamry[!duplicated(count2), count2]
    lim.lab = prot.snp.sumamry[!duplicated(count2), count]
  }
  if (length(lim.pos) > 6) {
    lim.dat = data.table::data.table(pos = lim.pos, lab = lim.lab)
    lim.dat[, `:=`(posRounded, round(pos))]
    lim.dat = lim.dat[!duplicated(posRounded)]
    lim.pos = lim.dat[, pos]
    lim.lab = lim.dat[, lab]
  }
  if (!defaultYaxis) {
    lim.pos = c(min(lim.pos), max(lim.pos))
    lim.lab = c(min(lim.lab), max(lim.lab))
  }
  clusterSize = 10
  if (repel) {
    prot.snp.sumamry = repelPoints(dat = prot.snp.sumamry, 
                                   protLen = len, clustSize = clusterSize)
  }
  else {
    prot.snp.sumamry$pos2 = prot.snp.sumamry$pos
  }
  xlimPos = pretty(0:max(prot$aa.length, na.rm = TRUE))
  xlimPos[length(xlimPos)] = max(prot$aa.length)
  if (!is.null(labelPos)) {
    prot.snp.sumamry = data.table::data.table(prot.snp.sumamry)
    if (length(labelPos) == 1) {
      if (labelPos != "all") {
        prot.snp.sumamry$labThis = ifelse(test = prot.snp.sumamry$pos %in% 
                                            labelPos, yes = "yes", no = "no")
        labDat = prot.snp.sumamry[labThis %in% "yes"]
      }
      else {
        labDat = prot.snp.sumamry
      }
    }
    else {
      prot.snp.sumamry$labThis = ifelse(test = prot.snp.sumamry$pos %in% 
                                          labelPos, yes = "yes", no = "no")
      labDat = prot.snp.sumamry[labThis %in% "yes"]
    }
    if (nrow(labDat) == 0) {
      message(paste0("Position ", labelPos, " doesn't seem to be mutated. Here are the mutated foci."))
      print(prot.snp.sumamry[, .(pos, conv, count, Variant_Classification)][order(pos)])
      stop()
    }
    if (collapsePosLabel) {
      uniquePos = unique(labDat[, pos2])
      labDatCollapsed = data.table::data.table()
      for (i in 1:length(uniquePos)) {
        uniqueDat = labDat[pos2 %in% uniquePos[i]]
        if (nrow(uniqueDat) > 1) {
          maxDat = max(uniqueDat[, count2])
          maxPos = unique(uniqueDat[, pos2])
          toLabel = uniqueDat[, conv]
          toLabel = paste(toLabel[1], paste(gsub(pattern = "^[A-z]*[[:digit:]]*", 
                                                 replacement = "", x = toLabel[2:length(toLabel)]), 
                                            collapse = "/"), sep = "/")
          labDatCollapsed = rbind(labDatCollapsed, data.table::data.table(pos2 = maxPos, 
                                                                          count2 = maxDat, conv = toLabel))
        }
        else {
          labDatCollapsed = rbind(labDatCollapsed, data.table::data.table(pos2 = uniqueDat[, 
                                                                                           pos2], count2 = uniqueDat[, count2], conv = uniqueDat[, 
                                                                                                                                                 conv]))
        }
      }
      labDat = labDatCollapsed
    }
  }
  domains = unique(prot[, Label])
  domain_cols = get_domain_cols()
  if (length(domains) > length(domain_cols)) {
    domain_cols = sample(colours(), size = length(domains), 
                         replace = FALSE)
  }
  domain_cols = domain_cols[1:length(domains)]
  domain_cols = grDevices::adjustcolor(col = domain_cols, alpha.f = domainAlpha)
  names(domain_cols) = domains
  col = col[unique(as.character(prot.snp.sumamry[, Variant_Classification]))]
  if (showLegend) {
    lo = matrix(data = c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
    graphics::layout(mat = lo, heights = c(4, 1.25))
    par(mar = c(1, 2.5, 2, 1))
  }
  else {
    par(mar = c(2.5, 2.5, 2, 1))
  }
  plot(0, 0, pch = NA, ylim = c(0, 6.5), xlim = c(0, len), 
       axes = FALSE, xlab = NA, ylab = NA)
  rect(xleft = 0, ybottom = 0.2, xright = len, ytop = 0.8, 
       col = "#95a5a6", border = bgBorderCol)
  axis(side = 1, at = xlimPos, labels = xlimPos, lwd = 1.2, 
       font = 1, cex.axis = axisTextSize[1], line = -0.4)
  axis(side = 2, at = lim.pos, labels = lim.lab, lwd = 1.2, 
       font = 1, las = 2, cex.axis = axisTextSize[2])
  segments(x0 = prot.snp.sumamry[, pos2], y0 = 0.8, x1 = prot.snp.sumamry[, 
                                                                          pos2], y1 = prot.snp.sumamry[, count2 - 0.03], lwd = 1.2, 
           col = "gray70")
  point_cols = col[as.character(prot.snp.sumamry$Variant_Classification)]
  points(x = prot.snp.sumamry[, pos2], y = prot.snp.sumamry[, 
                                                            count2], col = point_cols, pch = 16, cex = pointSize)
  prot[, `:=`(domainCol, domain_cols[prot[, Label]])]
  if (roundedRect) {
    if (requireNamespace("berryFunctions", quietly = TRUE)) {
      for (i in 1:nrow(prot)) {
        berryFunctions::roundedRect(xleft = prot[i, Start], 
                                    ybottom = 0.1, xright = prot[i, End], ytop = 0.9, 
                                    col = prot[i, domainCol], border = domainBorderCol, 
                                    rounding = 0.08)
      }
    }
    else {
      rect(xleft = prot[, Start], ybottom = 0.1, xright = prot[, 
                                                               End], ytop = 0.9, col = prot[, domainCol], border = domainBorderCol)
    }
  }
  else {
    rect(xleft = prot[, Start], ybottom = 0.1, xright = prot[, 
                                                             End], ytop = 0.9, col = prot[, domainCol], border = domainBorderCol)
  }
  title(main = cbioSubTitle, adj = 0, font.main = 2, cex.main = titleSize[1], 
        line = 0.8)
  title(main = unique(prot[, refseq.ID]), adj = 0, font.main = 1, 
        line = -0.5, cex.main = titleSize[2])
  if (showDomainLabel) {
    if (labelOnlyUniqueDoamins) {
      prot = prot[!duplicated(Label)]
    }
    prot$pos = rowMeans(x = prot[, .(Start, End)], na.rm = FALSE)
    text(y = 0.5, x = prot$pos, labels = prot$Label, font = 3, 
         cex = domainLabelSize)
  }
  if (!is.null(labelPos)) {
    text(x = labDat[, pos2], y = labDat[, count2 + 0.45], 
         labels = labDat[, conv], font = 1, srt = labPosAngle, 
         cex = labPosSize)
  }
  if (showLegend) {
    par(mar = c(0, 0.5, 1, 0), xpd = TRUE)
    plot(NULL, ylab = "", xlab = "", xlim = 0:1, 
         ylim = 0:1, axes = FALSE)
    lep = legend("topleft", legend = names(col), col = col, 
                 bty = "n", border = NA, xpd = TRUE, text.font = 1, 
                 pch = 16, xjust = 0, yjust = 0, cex = legendTxtSize, 
                 y.intersp = 1.5, x.intersp = 1, pt.cex = 1.2 * legendTxtSize, 
                 ncol = ceiling(length(col)/4))
    x_axp = 0 + lep$rect$w
    if (!showDomainLabel) {
      if (length(domain_cols) <= 4) {
        n_col = 1
      }
      else {
        n_col = (length(domain_cols)%/%4) + 1
      }
      lep = legend(x = x_axp, y = 1, legend = names(domain_cols), 
                   col = domain_cols, border = NA, ncol = n_col, 
                   pch = 15, xpd = TRUE, xjust = 0, bty = "n", 
                   cex = legendTxtSize, title = "Domains", 
                   title.adj = 0, pt.cex = 1.2 * legendTxtSize)
    }
  }
  if (printCount) {
    print(prot.snp.sumamry[, .(pos, conv, count, Variant_Classification)][order(pos)])
  }
}

#needed for lollipopPlot_fixed()
get_vcColors = function(alpha = 1, websafe = FALSE, named = TRUE){
  if(websafe){
    col = c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
            "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
            "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E",
            "#607D8B")
  }else{
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060', '#535c68')
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }
  
  if(named){
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                           'RNA','Splice_Site','Intron','Frame_Shift_Ins','In_Frame_Del','ITD','In_Frame_Ins',
                           'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event', 'pathway')
  }
  
  col
}
get_domain_cols = function(){
  c("#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#f19066",
    "#f5cd79", "#546de5", "#e15f41", "#c44569", "#786fa6", "#f8a5c2",
    "#63cdda", "#ea8685", "#596275", "#574b90", "#f78fb3", "#3dc1d3",
    "#e66767", "#303952")
}


# Changes Hugo Symbols to uppercase.

capHS<-function(maf) {
  
  maf %>%
    subsetMaf(tsb=levels(maf@data$Tumor_Sample_Barcode), 
              mafObj = F) %>%
    mutate(Hugo_Symbol=str_to_upper(Hugo_Symbol)) %>%
    read.maf()
}

#function to pipe MAF into tibble extraction 

extract.maf <- function(maf) { 
  subsetMaf(maf, tsb=levels(maf@data$Tumor_Sample_Barcode), mafObj=FALSE)
}

# Extracts a table with number of genes/variants above a certain VAF threshold 
vaf.table<-function(maflist, thresholds=seq(0, 0.9, by=0.1), unique="genes", basename=NULL){
  
  if (!unique == "genes" & !unique == "variants"){
    stop("Please specify a correct value for \"unique\" argument. default=\"genes\", also admits \"variants\".")
  }
  else {
    if (unique == "genes") {
      
      #Wrapper to create the temp vector with results for each MAF in maflist
      for (maf in maflist) {  
        
        #create empty vector with name for temp list. 1/maf
        assign(paste("genes", maf, sep="_"), c()) 
        
        #applies to the maf the same function with different VAF thresholds
        for(vaf in thresholds){ 
          
          # "a" is a temporary variable that serves to store the length output for each VAF threshold in each cycle
          a<- get(maf)@data %>%
            filter(VAF>=vaf) %>% 
            pull(Hugo_Symbol) %>% 
            unique() %>% 
            length()
          
          # glues "a" at the end of the specific maf vector
          assign(x=paste("genes", maf, sep="_"), value = c(get(paste("genes", maf, sep="_")), a))
          
        }
      }
      
      #creating an identifiable vector of the temp vectors names 
      vafgenestablenames<- paste("genes", maflist, sep="_") 
      
      #skeleton (empty df) of the final table
      vaf.gene.table<- data.frame()
      
      #binds one by one the temp vectors into a data frame, rowwise in the order of maflist
      for(i in 1:length(vafgenestablenames)){
        vaf.gene.table<-rbind(vaf.gene.table, get(vafgenestablenames[i]))
      }
      
      #setting dimnames for the final table
      rownames(vaf.gene.table)<-maflist
      colnames(vaf.gene.table)<- paste0("VAF>", thresholds)
      
      if(!is.null(basename)){
        
        vaf.gene.table %>% write.csv(file=paste0(basename, "_VAFtable.csv"))
      }
      
      #return
      return(vaf.gene.table)
      
    }  
    else{
      if (unique == "variants") {
        
        for(maf in maflist){
          
          assign(paste("variants", maf, sep="_"), c())
          for(vaf in thresholds){
            
            a<- get(maf)@data %>%
              filter(VAF>=vaf) %>% 
              pull(Hugo_Symbol) %>% 
              length()
            
            assign(x=paste("variants", maf, sep="_"), value = c(get(paste("variants", maf, sep="_")), a))
          }
        }
        
        vafvariantstablenames<- paste("variants", maflist, sep="_") 
        vaf.variant.table<- data.frame()
        
        for(i in 1:length(vafvariantstablenames)){
          vaf.variant.table<-rbind(vaf.variant.table, get(vafvariantstablenames[i]))
        }
        
        rownames(vaf.variant.table)<-maflist
        colnames(vaf.variant.table)<- paste0("VAF>", thresholds)
        
        if(!is.null(basename)){
          
          vaf.variant.table %>% write.csv(file=paste0(basename, "_VAFtable.csv"))
        }
        
        vaf.variant.table
        
      }
    }
  }
}


# Adds Protein_Change column with standardized HGVS.short protein mutation annotations
add.protchange <- function(maf){
  
  if (maf@maf.silent %>% rownames() %>% is_empty() == FALSE){
    
    message("MAF file contains silent mutations that may not be displayed. Careful!!")
    
  } 
  
  maf %>% 
    subsetMaf(tsb=levels(maf@data$Tumor_Sample_Barcode), mafObj = FALSE) %>%
    separate(Amino_acids, 
             sep="/", 
             into=c("ref.aa", "alt.aa"),
             remove = FALSE) %>%
    mutate(Protein_position=str_replace_all(Protein_position, "-", "_")) %>%
    separate(Protein_position,
             sep="_",
             into=c("pos1", "pos2"),
             remove= FALSE) %>%
    mutate(Protein_Change = case_when(
      Variant_Classification == "In_Frame_Ins" & !ref.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
      Variant_Classification == "In_Frame_Ins" & ref.aa == "-" ~ paste0("p.", Protein_position, "ins", alt.aa),
      Variant_Classification == "In_Frame_Del" & !alt.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
      Variant_Classification == "In_Frame_Del" & alt.aa == "-" ~ paste0("p.", ref.aa, Protein_position, "del"),
      Variant_Classification == "Frame_Shift_Ins" ~ paste0("p.", ref.aa, pos1, "fs"),
      Variant_Classification == "Frame_Shift_Del" ~ paste0("p.", ref.aa, pos1, "fs"),
      Variant_Classification == "Missense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
      Variant_Classification == "Nonsense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
      Variant_Classification == "Nonstop_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
      Variant_Classification == "Splice_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
      Variant_Classification == "Translation_Start_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa))) %>%
    select(!ref.aa & !alt.aa & !pos1 & !pos2) %>%
    read.maf()
}

# Modified mafCompare to work with total number of mutations rather than samples altered in

mutCompare<-function (m1, m2, m1Name = NULL, m2Name = NULL, minMut = 5, useCNV = TRUE) {
  
  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)
  if (is.null(m1Name)) {
    m1Name = "M1"
  }
  if (is.null(m2Name)) {
    m2Name = "M2"
  }
  if (useCNV) {
    m1.genes = as.character(m1.gs[total >= minMut, 
                                  Hugo_Symbol])
    m2.genes = as.character(m2.gs[total >= minMut, 
                                  Hugo_Symbol])
    uniqueGenes = unique(c(m1.genes, m2.genes))
  }
  else {
    m1.genes = as.character(m1.gs[total >= minMut, 
                                  Hugo_Symbol])
    m2.genes = as.character(m2.gs[total >= minMut, 
                                  Hugo_Symbol])
    uniqueGenes = unique(c(m1.genes, m2.genes))
  }
  m1.sampleSize = as.numeric(m1@summary[4, summary])
  m2.sampleSize = as.numeric(m2@summary[4, summary])
  m1.gs.comGenes = m1.gs[Hugo_Symbol %in% uniqueGenes]
  m2.gs.comGenes = m2.gs[Hugo_Symbol %in% uniqueGenes]
  sampleSummary = data.table::data.table(Cohort = c(m1Name, 
                                                    m2Name), SampleSize = c(m1.sampleSize, m2.sampleSize))
  if (useCNV) {
    m.gs.meged = merge(m1.gs.comGenes[, .(Hugo_Symbol, total)], 
                       m2.gs.comGenes[, .(Hugo_Symbol, total)], 
                       by = "Hugo_Symbol", all = TRUE)
  }
  else {
    m.gs.meged = merge(m1.gs.comGenes[, .(Hugo_Symbol, total)], 
                       m2.gs.comGenes[, .(Hugo_Symbol, total)], 
                       by = "Hugo_Symbol", all = TRUE)
  }
  m.gs.meged[is.na(m.gs.meged)] = 0
  m.gs.meged = as.data.frame(m.gs.meged)
  
  fisherTable = lapply(seq_len(nrow(m.gs.meged)), function(i) {
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i, 2]
    m2Mut = m.gs.meged[i, 3]
    xf = fisher.test(matrix(c(m1Mut, m1.sampleSize - m1Mut, 
                              m2Mut, m2.sampleSize - m2Mut), byrow = TRUE, nrow = 2), 
                     conf.int = TRUE, conf.level = 0.95)
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut, 
                                  m2Mut, pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, 
                                      fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  fisherTable[, `:=`(adjPval, p.adjust(p = pval, method = "fdr"))]
  colnames(fisherTable)[2:3] = c(m1Name, m2Name)
  return(list(results = fisherTable, SampleSummary = sampleSummary))
}


#multiplot takes objects from ggplot and puts them together

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


`%nin%` = Negate(`%in%`)

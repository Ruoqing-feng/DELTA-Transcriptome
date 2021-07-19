#@param genelist # input dysregulated genes, should be entrez ID
#@universe #background gene list
#@dd #a dataframe with a column containing expression level information and gene symbol as rownames.

#@author Yan Zhao
#@version 2019-12-18
# require(clusterProfiler)
# require(ReactomePA)

GOKEGGPlot = function(genelist, dd = NULL, universe = NULL, ont = 'ALL', top = 10,
                      pvalueCutoff = 0.05, qvalueCutoff = 0.05, outdir = NULL, main = NULL,
                      width = 7, height =5, charLength = 60, kegg_organism = 'hsa', GO_orgdb = org.Hs.eg.db){
  ####################################### Go term &KEGG for up grouped genes ####################
  message(Sys.time(), " # Enrichment analysis : GO term ...")
  go <- enrichGO(   gene   = genelist,
                     OrgDb  = GO_orgdb,
                     ont   = ont,
                     universe = universe,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = pvalueCutoff,
                     readable = T,
                     qvalueCutoff  = qvalueCutoff,
                     minGSSize = 5)
  go = data.frame(go)
  if(is.null(go) || nrow(go)==0){
      p1 = noEnrichPlot("No enriched terms")
      p2 = noEnrichPlot("No enriched terms")
      p3 = noEnrichPlot("No enriched terms")
  }
  else{
  go$termRatioStr =paste0(go$Count,'/', unlist(strsplit(go$BgRatio, '/'))[seq(1,nrow(go)*2,2)])
  go$termRatio = round(go$Count/as.numeric(unlist(strsplit(go$BgRatio, '/'))[seq(1,nrow(go)*2,2)]),2)
  go = go[order(go$pvalue, decreasing=T),]
  write.table(go, file.path(outdir, paste0(main, 'group_GO_enrich.csv')), sep = ',',row.names = F)
  print(paste0('GO terms: ', nrow(go)))

  # else{
  #   go$termRatioStr =paste0(go$Count,'/', unlist(strsplit(go$BgRatio, '/'))[1])
  #   go$termRatio = go$termRatioStr
  # }
    if (nrow(go)>top){
    go = go[(nrow(go)-top+1):nrow(go),]
  }
  ## Normalize term description ##
  terms = as.character(go$Description)
  terms = lapply(terms, function(x,k){
    x = as.character(x)
    if(nchar(x)>k){x=substr(x,start=1,stop=k)}
    return(x)}, charLength)
  go$Description = do.call(rbind, terms)
  go = go[!duplicated(go$Description),]
  
  p1 <- ggplot(data=go) +  geom_bar(aes(x=Description,y=-log10(pvalue), fill=ONTOLOGY), stat='identity') + coord_flip() +
    scale_x_discrete(limits=go$Description) +
    theme(axis.text.x=element_text(angle=0,size=6, vjust=0.7), axis.text.y=element_text(angle=0,size=8, vjust=0.7))
  
  p2 <- ggplot(go, aes(x=-log10( pvalue ), y=Description)) +
    geom_point(aes( size= termRatio , colour = Count)  ) + scale_y_discrete(limits=go$Description)+
    ggtitle("GO enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
    theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
          plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
          panel.background = element_rect(fill="white", colour='gray'),
          panel.grid.major = element_line(size = 0.05, colour = "gray"),
          panel.grid.minor.y = element_line(size=0.05, colour="gray"),
          panel.grid.minor.x = element_line(size=0.05, colour="gray")
    )
  go = go[order(go$p.adjust, decreasing=T),]
  p3 <- ggplot(go, aes(x=-log10( p.adjust ), y=Description)) +
    geom_point(aes( size= termRatio, colour = Count)  ) + scale_y_discrete(limits=go$Description)+
    ggtitle("GO enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
    theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
          plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
          panel.background = element_rect(fill="white", colour='gray'),
          panel.grid.major = element_line(size = 0.05, colour = "gray"),
          panel.grid.minor.y = element_line(size=0.05, colour="gray"),
          panel.grid.minor.x = element_line(size=0.05, colour="gray")
    )
  }
  ggsave(file.path(outdir, paste0(main,"group_GO_bar.pdf")), p1, width = width, height=height)
  ggsave(file.path(outdir,paste0(main,"group_Go_Bubble.pdf")), p2, width = width, height=height)
  ggsave(file.path(outdir,paste0(main,"group_Go_Bubble_ty2_y-padj.pdf")), p3, width = width, height=height)
  
  ###########KEGG
  message(Sys.time(), " # Enrichment analysis : KEGG ...")
  
  kegg = enrichKEGG(genelist, universe = universe, organism = kegg_organism)
  kegg = data.frame(kegg)
  if(is.null(kegg) || nrow(kegg)==0){
      q1 = noEnrichPlot("No enriched terms")
      q2 = noEnrichPlot("No enriched terms")
      q3 = noEnrichPlot("No enriched terms")
  }
  else{
  kegg$termRatioStr =paste0(kegg$Count,'/', unlist(strsplit(kegg$BgRatio, '/'))[seq(1,nrow(kegg)*2,2)])
  kegg$termRatio =round(kegg$Count/as.numeric(unlist(strsplit(kegg$BgRatio, '/'))[seq(1,nrow(kegg)*2,2)]),2)
  
  # else{
  #   kegg$termRatioStr =paste0(kegg$Count,'/', unlist(strsplit(kegg$BgRatio, '/'))[1])
  #   kegg$termRatio = kegg$termRatioStr
  # }
  
  kegg = kegg[order(kegg$p.adjust, decreasing = T),]
  write.table(kegg, file.path(outdir, paste0(main, 'group_KEGG_enrich.csv')), sep = ',', row.names = F)
  print(paste0('KEGG terms: ', nrow(kegg)))
  
  if (nrow(kegg)>top){
    kegg = kegg[(nrow(kegg)-top + 1):nrow(kegg),]
  }
  terms = as.character(kegg$Description)
  terms = lapply(terms, function(x,k){
    x = as.character(x)
    if(nchar(x)>k){x=substr(x,start=1,stop=k)}
    return(x)}, charLength)
  kegg$Description = do.call(rbind, terms)
  kegg = kegg[!duplicated(kegg$Description),]
  
  q1 <- ggplot(data=kegg) +  geom_bar(aes(x=Description,y=-log10(pvalue)), stat='identity') + coord_flip() +
    scale_x_discrete(limits=kegg$Description) +
    theme(axis.text.x=element_text(angle=0,size=6, vjust=0.7), axis.text.y=element_text(angle=0,size=8, vjust=0.7))
  
  q2 <- ggplot(kegg, aes(x=-log10( pvalue ), y=Description)) +
    geom_point(aes( size= termRatio, colour = Count)  ) + scale_y_discrete(limits=kegg$Description)+
    ggtitle("KEGG enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
    theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
          plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
          panel.background = element_rect(fill="white", colour='gray'),
          panel.grid.major = element_line(size = 0.05, colour = "gray"),
          panel.grid.minor.y = element_line(size=0.05, colour="gray"),
          panel.grid.minor.x = element_line(size=0.05, colour="gray")
    )
  kegg = kegg[order(kegg$p.adjust, decreasing = T),]
  q3 <- ggplot(kegg, aes(x=-log10( p.adjust ), y=Description)) +
    geom_point(aes( size= termRatio, colour = Count)  ) + scale_y_discrete(limits=kegg$Description)+
    ggtitle("KEGG enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
    theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
          plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
          panel.background = element_rect(fill="white", colour='gray'),
          panel.grid.major = element_line(size = 0.05, colour = "gray"),
          panel.grid.minor.y = element_line(size=0.05, colour="gray"),
          panel.grid.minor.x = element_line(size=0.05, colour="gray")
    )
    }
  ggsave(file.path(outdir, paste0(main, "group_KEGG_bar.pdf")), q1, width = width, height=height)
  ggsave(file.path(outdir, paste0(main, "group_KEGG_Bubble.pdf")), q2, width = width, height=height)
  ggsave(file.path(outdir, paste0(main, "group_KEGG_Bubble_ty2_y-padj.pdf")), q3, width = width, height=height)
  
  # message(Sys.time(), " # Enrichment analysis : Reactome ...")
  # reactome = enrichPathway(genelist, organism = kegg_organism, pvalueCutoff=pvalueCutoff, minGSSize = 5,readable=T,universe = universe)
  # reactome = as.data.frame(reactome)
  # message(Sys.time(), " # Enrichment analysis : Reactome ...")
  # if(is.null(reactome) || nrow(reactome)==0){
  #   q3 = noEnrichPlot("No enriched terms")
  # }
  # else{
  #   reactome$termRatioStr =paste0(kegg$Count,'/', unlist(strsplit(reactome$BgRatio, '/'))[seq(1,nrow(reactome)*2,2)])
  #   reactome$termRatio =round(reactome$Count/as.numeric(unlist(strsplit(reactome$BgRatio, '/'))[seq(1,nrow(reactome)*2,2)]),2)
  #   
  #   if (nrow(reactome)>top){
  #     reactome = reactome[(nrow(reactome)-top + 1):nrow(reactome),]
  #   }
  #   terms = as.character(reactome$Description)
  #   terms = lapply(terms, function(x,k){
  #     x = as.character(x)
  #     if(nchar(x)>k){x=substr(x,start=1,stop=k)}
  #     return(x)}, charLength)
  #   reactome$Description = do.call(rbind, terms)
  #   reactome = reactome[!duplicated(reactome$Description),]
  #   
  #   reactome = reactome[order(reactome$p.adjust, decreasing = T),]
  #   q3 <- ggplot(reactome, aes(x=-log10( p.adjust ), y=Description)) +
  #     geom_point(aes( size= termRatio, colour = Count)  ) + scale_y_discrete(limits=kegg$Description)+
  #     ggtitle("Reactome enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
  #     theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
  #           plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
  #           panel.background = element_rect(fill="white", colour='gray'),
  #           panel.grid.major = element_line(size = 0.05, colour = "gray"),
  #           panel.grid.minor.y = element_line(size=0.05, colour="gray"),
  #           panel.grid.minor.x = element_line(size=0.05, colour="gray")
  #     )
    
  # }
  # ggsave(file.path(outdir, paste0(main, "group_Reactome_Bubble_padj.pdf")), q3, width = width, height=height)
  # 
  
  return(kegg)
  


# if (nrow(kegg)>0){
#     if (length(grep('hsa05206', kegg$ID))){
#         kegg = kegg[-grep('hsa05206', kegg$ID),]
#     }
#     pathview(dd, pathway.id = kegg$ID, gene.idtype = 'symbol', out.suffix = 'KEGG_All' )
#     
# }
}

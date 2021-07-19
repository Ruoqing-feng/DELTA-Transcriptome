
noEnrichPlot = function(main='no enriched terms'){
  
test = matrix(c(0,0,0,0),2,2)
test = as.data.frame(test)
ggplot(test)+geom_text(aes(V1,V2, label = main))
}
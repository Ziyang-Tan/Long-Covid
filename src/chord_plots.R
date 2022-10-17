library(dplyr)
library(circlize)



clone_vb_chord_plot <- function(data_vb_sub, sample_name){
  df <- data_vb_sub %>% 
    filter(Sample_Name == sample_name) %>%
    group_by(clone_id, V) %>%
    tally()
  chord_order <- c(df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F),
                   df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F))
  chordDiagram(df, 
               annotationTrack = "grid", 
               order = chord_order,
               preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  title(sample_name)
  circos.clear()
}

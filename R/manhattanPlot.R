# function created from https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
manhattanPlot <- function(x, colors=c("grey", "skyblue"),
                          significanceLine = 8, 
                          snpsOfInterest=NULL, ...) {

  don <- x %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(x, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    mutate( is_annotate=ifelse(-log10(P)>8, "yes", "no")) 
  
  # Prepare X axis
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # Make the plot
  ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(colors, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    
    # Change X-legend 
    xlab("Chromosome") +
    
    # Add genome-wide line
    geom_hline(yintercept = significantLine, linetype="dashed") +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

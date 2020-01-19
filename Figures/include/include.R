################################################################################
#
# Wout Megchelenbrink
# Feb 25, 2017
#
#
# Include
# Includes custom functions
################################################################################


# From: http://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot
#' x      numeric vector for each slice
#' group  vector identifying the group for each slice
#' labels vector of labels for individual slices
#' col    colors for each group
#' radius radius for inner and outer pie (usually in [0,1])
donut_plot <- function(dt, main="maincat", sub="subcat", values="vals", labels = NA, col = RColorBrewer::brewer.pal(8, "Set1"), radius = c(.7, 1)) 
{
  cols <- data.table(main=unique(dt[[main]]), col_main=col[1:length(unique(dt[[main]]))])
  dt <- merge(dt, cols, by.x=main, by.y="main")
  
  # Each subcat must have a different shade!
  dt[, subid:=1:.N, by=main]
  n <- max(dt$subid)
  
  print(dt)
  
  dt[, col_sub:=alpha(col_main, (0.85-0.15)/n * ((n+1)-subid))]
  
  plot.new()
  
  par(new = TRUE)
  pie(dt[[values]], border = NA, radius = radius[2L],
      col = dt$col_sub, labels = labels, cex=1)
  
  par(new = TRUE)
  pie(dt[[values]], border = NA, radius = radius[1L],
      col = dt$col_main, labels = NA)
  
  dt[, col_main:=NULL]
  dt[, col_sub:=NULL]
}

format_scientific <- function(l) 
{
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}



correlation.plot <- function(mat, title="R")
{
    col.range <- c(min(mat), (min(mat)+max(mat))/2, max(mat))  
    col_fun = circlize::colorRamp2(col.range, c( bs.col.dark[2], "white", bs.col.dark[1]))
    
    hm <- Heatmap(mat, name = title, col = col_fun, rect_gp = gpar(type = "none"), 
                  cell_fun = function(j, i, x, y, width, height, fill) 
                  {
                    if(i >= j) 
                    {
                      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "#DDDDDD", fill = NA))
                      grid.circle(x = x, y = y, r = 0.06, gp = gpar(fill = col_fun(mat[i,j])))
                      grid.text(gp = gpar(fontsize=16), sprintf("%2.2f", mat[i,j]), x = x, y = y,  vjust = .5)
                    } 
                  }, cluster_rows = FALSE, cluster_columns = FALSE, 
                  show_row_names = T, show_column_names = T, row_names_side = "left", column_names_side = "bottom", column_names_gp = gpar(fontsize = 14, fontface="bold"),
                  row_names_gp = gpar(fontsize = 14, fontface="bold"),   heatmap_legend_param =  list(color_bar = "continuous", legend_height = unit(4, "cm")))
    return(hm)
}
  

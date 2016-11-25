geom_tallrect <- function
### ggplot2 geom with xmin and xmax aesthetics that covers the entire
### y range.
(mapping=NULL,
 data=NULL,
 stat="identity",
 position="identity",
 ...){
  require(proto)
  require(grid)
  GeomTallRect <- proto(ggplot2:::GeomRect,{
    objname <- "tallrect"
    required_aes <- c("xmin", "xmax")
    draw <- draw_groups <- function(.,data,scales,coordinates,
                                    ymin=0,ymax=1,...){
      ymin <- unit(ymin,"npc")
      ymax <- unit(ymax,"npc")
      with(coord_transform(coordinates, data, scales),ggname(.$my_name(), {
        rectGrob(xmin, ymin, xmax - xmin, ymax-ymin,
                 default.units = "native", just = c("left", "bottom"), 
                 gp=gpar(
                   col=colour, fill=alpha(fill, alpha), 
                   lwd=size * .pt, lty=linetype, lineend="butt"
                   )
                 )
      }))
    }
  })
  GeomTallRect$new(mapping = mapping, data = data, stat = stat,
                   position = position, ...)
}


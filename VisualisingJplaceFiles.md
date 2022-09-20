# Visualising placement (.jplace) files

In R: (note: do not treat this like a script...)

```r
#Loading libraries----
library(treeio)
library(ggtree)
library(phytools)
library(ape)

#Loading file(s)----
EPAtree<-read.jplace("/path/to/jplace/file")
print(EPAtree)

placetree <- ggtree(EPAtree) + geom_tiplab(size=3) + geom_nodepoint((aes(subset=nplace>0, size=nplace)), color="salmon", alpha=0.5) + geom_tippoint((aes(subset=nplace>0, size=nplace)), color="salmon", alpha=0.5)

#NOTE that the ggtree object loses placement data!!!

#function for grabbing all placements for a node----
grab_placements <- function(nodenumber,tree) {
  placementlist <- select(filter(tree@placements, node == nodenumber), name, like_weight_ratio) %>% arrange(desc(like_weight_ratio))
  return(placementlist)
}

#plotting----
quartz("MetVAMPS", width=8, height=12) #create new window (change height to fit tree)
plot.new()
Metplace
last_plot

ggtree:::identify.gg(x=Metplace, tolerance = 0.01, n=3, order = TRUE) %>% grab_placements(MetEPA) #has errors...

#saving PDF with placements----
quartz.save("MetVAMPS.pdf", type="pdf")
```

## How ggtree parses .jplace files

Part of series “Why provide adequate documentation when users can reverse-engineer how it works?”

ggtree version 2.4.1

```r
MetEPA <- read.jplace("RAxML_portableTree.MetSoilMetaT-EPA.jplace")
```

Creates a fairly strictly-defined S4 object consisting of:`@phylo` : strictly-defined S3 object describing the tree as a table; messing with this makes R scream. `phylo` is focused on nodes. It contains a `2xnnodes` table listing connected nodes in tips-base order, with associated `edge.lengths` and `tip.labels` for nodes that are tips. Edge numbering is hidden in attributes, to access use

```r
getNode <- function(S4tree, edge){
  attr(S4tree@phylo, 'edgeNum') %>% filter(edgeNum == edge) %>% select(node)
}
```

Or, more reliably, create a reference tibble:

```r
PrepareEdgeNodeS4 <- function(tree){
  EdgeNode<-tree@phylo$edge %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="EdgeName") %>% 
    rename("Node1"="V1","Node2"="V2")
  return(EdgeNode)
}
EdgeNode<-PrepareEdgeNodeS4(MetEPA)
```

And then refer to EdgeName, Node1, Node2 as needed in EdgeNode.

`@data` node, nplace:

can add column, associated with node, but tricky! Use `dplyr::mutate`, etc (to preserve data type or else `ggtree()` complains)

`@placements`

edge_num DOES NOT correspond to the edges in `@phylo`, but the original edge numberings in the jplace files. These are translated to nodes (see node column), which DO correspond to nodes in the `@phylo` tree. R phylo parsers number nodes by going through tips first, then deeper and deeper in the branching order; RAxML-EPA numbers edges in order of appearance in the output Newick file.

`@file`: name of source file, useful for auto title generation. Use filestrings package to create a shortname.

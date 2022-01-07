These scripts can be used to generate the Supplemental dotplots.

To generate the PDF images, you run the dotplot_maker.sh file which does the following:

The first step is to run minimap2 of the two assemblies, using -x asm5 to match regions of higher divergence, and -X to exclude self mappings

Then, the pafCoordsDotPlots_updated.R produces the images, using code adapted from https://github.com/tpoorten/dotPlotly

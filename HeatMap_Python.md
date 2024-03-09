# Heat Maps in python using seaborn

Assuming distance matrix already calculated

```python
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

heat =  pd.read_csv("SkolioBarth_identity_heatmap.csv", sep = ",")
heat.pop("Unnamed: 0")  # weird import issue
heat.index = heat.columns

# heatplot = sns.heatmap(heat, annot = True)  # default
# plt.show()

# removing top half
maskmatrix = np.triu(heat)  # upper triangle of array
# heatplot = sns.heatmap(heat, annot = True, mask = maskmatrix) # this excludes the diagonal
# but include diagonal
maskmatrix[np.diag_indices_from(maskmatrix)] = False
# heatplot = sns.heatmap(heat, annot = True, mask = maskmatrix)
heatplot = sns.heatmap(heat, annot = True, mask = maskmatrix, fmt = ".2f", cmap = sns.cm.mako_r)

plt.rcParams['svg.fonttype'] = 'none' #  render text as actual text
# plt.show()
plt.savefig("heatmap.svg", dpi=250)



```

```python
heatplot = sns.heatmap(heat, 
annot = True,  # include numbers inside each cell 
mask=maskmatrix,  # exclude upper triangle in this case
fmt=".2f", # keep 2 decimal places, e.g. 1.00 not 1
cmap=sns.cm.rocket_r,  # colour map
vmin=0, vmax=1)  # sets scale to 0-1; not always necessary


```

[python - How to plot only the lower triangle of a seaborn heatmap? - Stack Overflow](https://stackoverflow.com/questions/57414771/how-to-plot-only-the-lower-triangle-of-a-seaborn-heatmap)

[python 3.x - Round decimal places seaborn heatmap labels - Stack Overflow](https://stackoverflow.com/questions/63964006/round-decimal-places-seaborn-heatmap-labels)

Citation for Seaborn: `Waskom, M. L., (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021, https://doi.org/10.21105/joss.03021.`

[Citing and logo &#8212; seaborn 0.13.2 documentation](https://seaborn.pydata.org/citing.html)

[Choosing color palettes &#8212; seaborn 0.13.2 documentation](https://seaborn.pydata.org/tutorial/color_palettes.html)

-   notes on colour theory in diagram design: [Elegant Figures - Subtleties of Color (Part 1 of 6)](https://earthobservatory.nasa.gov/blogs/elegantfigures/2013/08/05/subtleties-of-color-part-1-of-6/)

Calculating (and displaying) pairwise distance matrices: [Making a pairwise distance matrix in pandas | Drawing from Data](https://drawingfromdata.com/pandas/clustering/making-a-pairwise-distance-matrix-in-pandas.html)



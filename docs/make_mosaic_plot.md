# Make mosaic plot
This tutorial is to use a python script to draw a mosaic plot for visualizing frequency distribution of two variables.

#### Python packages required

* [Pandas](https://pandas.pydata.org/)
* [SciPy](https://scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [statsmodels](https://www.statsmodels.org/stable/index.html)

##### Drawing a mosaic plot using `mosaic_plot.py`

You will use a python script [mosaic_plot.py](../scripts/mosaic_plot.py) in the path `path_to_the_package/KunDH-2023-CRM-MSM_metagenomics/scripts/`, and a table containing species associated with *MSM* and *Non-MSM* individuals which were identified as Gram-negative or not in [two_variable_mosaic.tsv](../example_data/two_variable_mosaic.tsv).

```{python}
usage: mosaic_plot.py [-h] [--input [INPUT]] [--facecolor_map [FACECOLOR_MAP]] [--font_style [FONT_STYLE]] [--output [OUTPUT]]

This program is to draw a mosaic plot.

optional arguments:
  -h, --help            show this help message and exit
  --input [INPUT]       Input a file containing two variable information regarding each individual subject.
  --facecolor_map [FACECOLOR_MAP]
                        Specify the the pathway to SCFA metabolisms database. default: /vol/projects/khuang/databases/SCFA/SCFA_pathways.tsv
  --font_style [FONT_STYLE]
                        Specify the font style, font family and font type is delimited by a comma. default: [sans-serif,Arial]
  --output [OUTPUT]     Specify the output figure name.

examples: mosaic_plot.py --input input_file.tsv --facecolor_map facecolor_mapfile.tsv --output mosaic_plot.png   
```

Example command:
```{bash}
mosaic_plot.py --input two_variable_mosaic.tsv --facecolor_map facecolor_map.tsv --output mosaic_plot.png
```
![Mosaic plot](../images/)

*Note*
The face color of mosaic plot should be specified as in the example [mapping file](../example_data/facecolor_map.tsv).


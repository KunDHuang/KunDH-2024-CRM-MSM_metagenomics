# Machine learning analysis
This tutorial is to use Random Forest model to evaluate the predictive capability of microbiome taxonomic composition.

## Evaluate the predictive power of microbiome taxonomic composition

#### Python package required

* [scikit-learn](https://scikit-learn.org/stable/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)

#### Random forest training and results evaluation with ROC-AUC curve

Here, we would like to introduce a python script [evaluation_kfold.py](../scripts/evaluation_kfold.py) which implements [random forest model](https://en.wikipedia.org/wiki/Random_forest) to evaluate the predictive capability of information encoded in microbiome taxonomic composition for classifying different individuals.

```{python}
usage: evaluation_kfold.py [-h] [--mpa_df [MPA_DF]] [--md_rows [MD_ROWS]] [--target_row [TARGET_ROW]] [--pos_feature [POS_FEATURE]] [--neg_feature [NEG_FEATURE]] [--fold_number [FOLD_NUMBER]] [--repeat_time [REPEAT_TIME]] [--output [OUTPUT]]
                           [--output_values [OUTPUT_VALUES]] [--nproc [NPROC]] [--transform [TRANSFORM]]

This script is to estimate ROC AUC based on metaphlan-style table with metadata being inserted.

optional arguments:
  -h, --help            show this help message and exit
  --mpa_df [MPA_DF]     Input a mpa-style table with metadata being inserted.
  --md_rows [MD_ROWS]   Input row numbers for specifying metadata without considering header row, zero-based, comma delimited. for example, 0,1,2,3,4.
  --target_row [TARGET_ROW]
                        Specify the row number for indicating target metadata to examine, zero-based without considering header row.
  --pos_feature [POS_FEATURE]
                        Specify the feature name to be labeled as posive, e.g. 1.
  --neg_feature [NEG_FEATURE]
                        Specify the feature name to be labeld as negative, e.g. 0.
  --fold_number [FOLD_NUMBER]
                        Specify the fold number you want split the whole dataset.
  --repeat_time [REPEAT_TIME]
                        Specify the repeat time you want to split the dataset.
  --output [OUTPUT]     Specify the output figure name.
  --output_values [OUTPUT_VALUES]
                        Specify the output file name for storing ROC-AUC values.
  --nproc [NPROC]       Specify the number of processors you want to use. 4 by default.
  --transform [TRANSFORM]
                        Transform values in the matrix, [arcsin_sqrt] or [binary] or [None]. [None] by default

examples: evaluation_kfold.py --mpa_df <mpa_df.tsv> --md_rows 0,1,2,3,4 --target_row 3 --pos_feature <CRC> --neg_feature <Healthy> --fold_number 10 --repeat_time 20 --output ROC_AUC.svg  
```

To demostrate the tutorial, we use a microbiome dataset [machine_learning_input.tsv](../example_data/machine_learning_input.tsv) containing 52 subjects (28 subjects have >3 sexual partners and 24 sibjects have 0-3 sexual partners) and corresponding relative abundances of gut microbiome species.

Example command:
```{bash}
evaluation_kfold.py --mpa_df machine_learning_input.tsv --md_rows 0 --target_row 0 --pos_feature ">3" --neg_feature "0_3" --fold_number 3 --repeat_time 50 --output roc_auc_npartners.png --output_values roc_auc_npartners_values.tsv --nproc 10
```

It generates a ROC-AUC curve to show the overall predictive capability of random forest model fitting to our input microbiome taxonomic data.
![roc_auc_npartners.png](../images/roc_auc_machine_learning.png)

Optionally, it can also generate the raw output [roc_auc_npartners_values.tsv](../example_data/roc_auc_npartners_values.tsv) used to generate the plot above. One can use it for other purposes.

**Note:** The figure displayed above had been edited using [inkscape](https://inkscape.org/) on the base of the crude output in order to enhance the readability and aesthetic sense.
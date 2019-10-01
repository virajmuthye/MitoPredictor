## Mitopredictor v1 
A novel machine-learning tool for prediction of mt-proteins in animals
#########################################################################

## How to run Mitopredictor

To run mitopredictor with default settings, follow the following steps.

## Step 0 Downloading software for mitopredictor pipeline
```markdown
1] Proteinortho v5.16b
2] TargetP v1.1
3] MitoFates
4] CD-HIT
5] pfam_scan.pl along with the Pfam-A libraries
6] R
```
## Step 1 Prep file for analysis

./step1_prep.bash

This step processes the input proteome file for the next three steps. 
It does three things:
1] Uses CD-HIT to remove redundancy. A default value of 98% similarity is used.
This can be changed in the script in the script "prep.bash" in the Prep folder.
<-------------------------
# removes redundancy by running CD-hit
# currently removes all proteins at 98% similarity, this can be changed by changing the -c option
cd-hit -i renamed.pep -o renamed.pep.98 -c 1 -n 5
-------------------------->

2] Removes all proteins below 100 amino-acids
3] Removes all proteins which do no begin with a methionine (potential fragments). This step is meant for MTS predictions (Step3) since MTS predictions are sensitive to proteins which do not have a complete N-terminus end.


## Step 2 Find orthologous groups using Proteinortho


! IMPORTANT: Steps 2, 3, 4 can be run in parallel. Once these three steps are complete, proceed with Step 5.


./step2_ortho.bash
In this step, Proteinortho is used to identify groups of orthologous proteins in the reference species and the inpur query proteome.
Default settings are used for Proteinortho. These can be changed in the "run_proteinortho.bash" script in the orthology folder.


## Step 3 Find proteins with MTS

./step3_mts.bash
In this step, TargetP and MitoFates predicts MTS in the query proteome.
Default settings are used for TargetP and MitoFates. These can be changed in the "runTargetP.bash" and "runMitoFates.bash" scripts in the mts folder.


## Step 4 Protein domain analysis for identification of mt-proteins


./step4_domain.bash
Protein domain information is used to identify mt-proteins.


## Step 5 Identify mt-proteins using information from steps 2,3,4

./step5_identify_and_concat.bash

Uses a Random Forest classifier to identify mt-proteins.
Training datasets for this model included human, mouse, C. elegans and D. melanogaster. If the query proteome is a mammal, user can change the training dataset to include only mammalian mt- and non-mt-proteins. This can be changed in the R script "rf.R".

<----------------------------------------------------------
#read in training dataset and the mito_score_matrix files
training <- read.table("all.training.csv", sep = ",", header = T, stringsAsFactors = F)

#if you want to use a mammal-only training dataset,uncomment the next line and comment the previous line
#training <- read.table("onlymammals.training.csv", sep = ",", header = T, stringsAsFactors = F)

--------------------------------------------------------->

Creates the FINAL.MATRIX file whicn includes information for all proteins in this study.
This FINAL.MATRIX file is the input for the "mito_app.R": an R Shiny applet for exploration of the inferred mt-proteome.

## Step 6 Make basic statistics of the inferred mt-proteome

./step6_make_stats.bash
The basic statistics of the inferred mt-proteome are stored in the "Stats" folder.

## Step 7 Cleans intermediate files and prepares mitopredictor for next analysis

./step7_clean.bash


---------------------------------------------------------------------
For any queries please contact Viraj Muthye at viraj.muthye@gmail.com
---------------------------------------------------------------------
























You can use the [editor on GitHub](https://github.com/virajmuthye/mitopredictor/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/virajmuthye/mitopredictor/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.

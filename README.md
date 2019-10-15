# MitoPredictor 
A novel machine-learning tool for prediction of mt-proteins in animals
Source code for MitoPredictor can be downloaded on GitHub https://github.com/virajmuthye/mitopredictor.github.io

## How to run MitoPredictor

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
```markdown
./step1_prep.bash
```
This step processes the input proteome file for the next three steps. 
It does three things:
1] Uses CD-HIT to remove redundancy. A default value of 98% similarity is used.
This can be changed in the script in the script "prep.bash" in the Prep folder.

```markdown
cd-hit -i renamed.pep -o renamed.pep.98 -c 1 -n 5
```


2] Removes all proteins below 100 amino-acids
3] Removes all proteins which do no begin with a methionine (potential fragments). This step is meant for MTS predictions (Step3) since MTS predictions are sensitive to proteins which do not have a complete N-terminus end.

## Step 2 Find orthologous groups using Proteinortho
### IMPORTANT: Steps 2, 3, 4 can be run in parallel. Once these three steps are complete, proceed with Step 5.

```markdown
./step2_ortho.bash
```

In this step, Proteinortho is used to identify groups of orthologous proteins in the reference species and the inpur query proteome.
Default settings are used for Proteinortho. These can be changed in the "run_proteinortho.bash" script in the orthology folder.

## Step 3 Find proteins with MTS
```markdown
./step3_mts.bash
```

In this step, TargetP and MitoFates predicts MTS in the query proteome.
Default settings are used for TargetP and MitoFates. These can be changed in the "runTargetP.bash" and "runMitoFates.bash" scripts in the mts folder.


## Step 4 Protein domain analysis for identification of mt-proteins
```markdown
./step4_domain.bash
```
Protein domain information is used to identify mt-proteins.


## Step 5 Identify mt-proteins using information from steps 2,3,4
```markdown
./step5_identify_and_concat.bash
```

Uses a Random Forest classifier to identify mt-proteins.

Creates the FINAL.MATRIX file whicn includes information for all proteins in this study.
This FINAL.MATRIX file is the input for the "mito_app.R": an R Shiny applet for exploration of the inferred mt-proteome.

## Step 6 Make basic statistics of the inferred mt-proteome
```markdown
./step6_make_stats.bash
```

The basic statistics of the inferred mt-proteome are stored in the "Stats" folder.

## Step 7 Cleans intermediate files and prepares mitopredictor for next analysis
```markdown
./step7_clean.bash
```

---------------------------------------------------------------------
For any queries please contact Viraj Muthye at viraj.muthye@gmail.com
---------------------------------------------------------------------





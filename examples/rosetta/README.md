This folder contains sample output from the Rosetta tutorial (https://docs.rosettacommons.org/demos/latest/tutorials/Protein-Protein-Docking/Protein-Protein-Docking) for the purpose of buildling Task #3 of PPInsight.
'Protein-Protein-Docking/output_files/expected_output' is a folder containg expected example output files from running the demo commands in 'demos/tutorials/Protein-Protein-Docking' after correctly installing Rosetta (https://docs.rosettacommons.org/docs/latest/getting_started/Getting-Started#local-installation-and-use-of-rosetta_installation-on-mac-linux):

**The values we want are in the score files** (the '.sc' files; https://docs.rosettacommons.org/demos/latest/tutorials/analysis/Analysis):
> "The score file contains the scores for all models generated in a run, broken down into the individual terms. You can extract certain columns (i.e. scoring terms) or sort them by any column.
> ```bash $ sort -n -k2 example_score_file.sc ```
> this will sort the entire file by column 2 (total score)
> ```bash $ sort -n -k2 example_score_file | awk '{print $2 "\t" $3}' ```
> this will do the same, but only print out columns 2 and 3

> Extract score and rmsd values for the best 1000 models!
> ```bash $ sort -n -k2 example_score_file.sc | head -n 1000 | awk '{print $2 "\t" $25 "\t" $NF}' > score_rmsd.dat ```
> this will sort by total score, take only the top 1000, extract columns 2 (score), 25 (rms) and the very last one (description,tag) and write this into a new file, called score_rmsd.dat."
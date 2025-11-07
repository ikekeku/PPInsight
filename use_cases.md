# This is the document for use cases that will be utilized by PPI 

## Describe what are the inputs and outputs, and how are the inputs transformed into outputs 


SETUP 
- Accessing a remote server (HYAK?)  

- Interpreting a protein name and pulling FASTA files from the internet:
The input of this function will be the protein name and the species of interest. The output of this function will be the according FASTA file for the protein of interest. The function itself will be able to input the two designations into the search bar of an according database in order to pull the file output. 

- Taking in a FASTA file and outputting a sequence (?) 
- Save the FASTA files of multiple proteins in one folder
- Taking the designated FASTA files and putting them into models to run
- Having chosen the proteins and saved the FASTA files, compare two models of the interaction between the selected proteins:
The inputs of this component are the models that the user wants to compare and the FASTA files of the selected proteins. This requires a menu of models from which the user is able to select at least (or exactly?) two. It also requires either a place where the user uploads the FASTA files, or where the program pulls the FASTA files from an earlier component. The output will be a CSV (or similar) file that contains the scores which will then serve as inputs for a visualization/plotting step.

- Output of the model result and keep as a value in your local server 
- Plot the scores that compare models
  

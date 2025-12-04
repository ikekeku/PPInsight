# This is the document for use cases that will be utilized by PPI 

## Describe what are the inputs and outputs, and how are the inputs transformed into outputs 




### TASK 1: Fetch proteins by search term and save files



### TASK 2: Run PPI predictors on file pairs and store scores



### Task 3: plot interaction scores (D3 visualizations)





##### EVERYONE RESTRUCTURE THE USE CASES CLEARLY UNDER THE ACCORDING TASK!!! 
SETUP 
- Accessing a remote server (HYAK?)  

- Interpreting a protein name and pulling FASTA files from the internet:
The input of this function will be the protein name and the species of interest. The output of this function will be the according FASTA file for the protein of interest. The function itself will be able to input the two designations into the search bar of an according database in order to pull the file output.
- Retrieve PDB structure from PDB ID or protein name
- If multiple FASTA files are returned (e.g., isoforms), display all candidates and let the user pick the one(s) of interest
- Option to upload a local FASTA file
-  Save the FASTA files of one or multiple proteins into one folder
- Taking the designated FASTA files and putting them into models to run
- Having chosen the proteins and saved the FASTA files, compare two models of the interaction between the selected proteins:
The inputs of this component are the models that the user wants to compare and the FASTA files of the selected proteins. This requires a menu of models from which the user is able to select at least (or exactly?) two. It also requires either a place where the user uploads the FASTA files, or where the program pulls the FASTA files from an earlier component. The output will be a CSV (or similar) file that contains the scores which will then serve as inputs for a visualization/plotting step.
- Output of the model result and keep as a value in your local server 
- Plot the scores that compare models
- Filter plots by score type or protein pair
- Export plot option: allow PNG or SVG export of graphs
- The system can compile all results into a text summary report, listing models and metrics, as well as key findings (input by user)
  

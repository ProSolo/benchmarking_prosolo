# Files in this folder

All files in this folders are *aggregational include files*. Using the snakemake *include* statement, they are meant to include all rules-files from the folder they are named after. I.e. they need to be updated whenever a new file is added in one of these folders.

# Subfolders

Each subfolder should contain rules for a distinct stage of any processing pipeline. I.e. input and output file format of a stage should be well defined and output filenames should contain a distinct tag that can be required by following stages of the pipeline. Pipeline parts are kept separate, to avoid excessive cluttering of filenames, as filename tags are reduced to one per stage. The *config.yaml* file should specify the rules to be applied in a stage.

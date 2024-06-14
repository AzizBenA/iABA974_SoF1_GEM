1-Follow the instructions provided in the carveme doc:
https://carveme.readthedocs.io/en/latest/installation.html

2-In my case im using a WSL so i need to request a new license key for Gurobi each time i run carveme.
Here are the instruction to do it :https://support.gurobi.com/hc/en-us/articles/7367019222929-How-do-I-set-up-Gurobi-in-WSL2-Windows-Subsystem-for-Linux

3-This is the command for creating my model:

carve GEM/files/FAA_merged.faa --output GEM/models/model_gramneg.xml -u gramneg --solver gurobi

This is the command to gapfill:

need to add he metabolites to the SF medium in the media_db.tsv :
home\aziz\anaconda3\envs\carveme\lib\python3.11\site-packages\carveme\data\input

gapfill GEM/models/updated_growth_gapfilled_autotrophic_1.sbml -m SF -o GEM/models/gapfilled_eth.xml


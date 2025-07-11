# pensa_analysis
my workflow for pensa


This will contain a workflow for doing PENSA analysis on simulations

## setup analysis directory
- If this is your first time doing this protocol, make a directory in your home (`cd ~`) called ‘pensa scripts’ (`mkdir ~/pensa_scripts` )  and also ‘pensa_setup_files’ (`mkdir ~/pensa_setup_files` )
- also set up pensa from here https://pensa.readthedocs.io/en/latest/installation.html
    - Copy files from the folders in my home (`/home/users/evejfine/pensa_scripts` and `/home/users/evejfine/pensa_setup_files`) into yours.
        - `cp /home/users/evejfine/pensa_scripts/* ~/pensa_scripts`
        - `cp /home/users/evejfine/pensa_setup_files/* ~/pensa_setup_files`
- make the parent analysis directory and then inside this, make a directory called scripts
    - `mkdir CB1_arr_Gi`  (replace CB1_arr_Gi with whatever you want to name it)
    - `mkdir CB1_arr_Gi/scripts`
- copy files from the pensa setup directory into the scripts directory:
    - in your parent directory, run `cp ~/pensa_setup_files/* scripts/`
 
## Preprocessing

- set up [1-preprocess.py](http://1-preprocess.py) (in the scripts folder)
    - Modify the root directories for each set of simulations, the paths to the pdb, psf, and trajectory files, the selection strings, and the number of replicates
    - then from your analysis directory, `bash ~/pensa_scripts/sub_1-preprocess.sh`
        - make sure that you activate the pensa environemnt `conda activate pensa` , which you presumably set up from https://pensa.readthedocs.io/en/latest/installation.html
- this will make a traj folder with the xtc files and gro files
    

  ## Compare conditions

- set up [2-compare.sh](http://2-compare.sh) (in the scripts folder) with the receptor name, and condition names (this will map to your gro and xtc files from preprocessing in the traj directory)
    - from your analysis directory, `bash ~/pensa_scripts/sub_2-compare.sh`
- this should create things in a results folder which has csvs, and a vispdb folder which you can visualize in pdb

## Visualize results from compare conditions

- in pymol, open any one of the torsions files
- in pymol: `cartoon putty`
- in pymol: `spectrum b, white_blue_green_yellow, minimum=0, maximum=1`

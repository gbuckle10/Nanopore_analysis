# Peter's Instructions :)

## To finish the analysis:
    - Basecalling:
        - If the basecalling isn't finished, you can continue from an existing bam file by using the --resume-from
        command.
        - In the 01_basecalling.sh file there is a basecalling_pod5() function. I have already fixed it so that it'll continue
        from the file that's already being done. It should be resumed from "data/basecalled_output/cyclomics2_continued.bam" and
        in the function call underneath it should be writing to "data/basecalled_output/cyclomics2_continued_forreal.bam"
        - The basecalling should also contain a demultiplex_bam() function. Make sure the local raw_bam variable assignment
        points to the right file. If you had to run it a second time it's already correct, if not you need to change the raw_bam
        variable to be cyclomics2_continued.bam
    - Merging bam files:
        - This is not part of the pipeline, so we need to do this manually. There should now be 2 bam files for each barcode, one in
        the set2 folder and one in the set1 folder. We need to do this before alignment because merging will destroy the sorting.
        - You can merge the bam files with the command:
            samtools merge -o analysis/combined/demultiplexed/barcode_X.bam path/to/file1.bam path/to/file2.bam
        - This might be quite annoying, but I think automating it at the moment would be more annoying.
    - Alignment:
        - Once the files are merged and put into the analysis/combined folder, you can run the alignment script. To just run the
        alignment script, go to config.yaml and comment out all of the run steps apart from align.
        - The 02_alignment.sh script should only run the align_and_index() method. In that method, the alignment_cmd will
        take the entire demultiplexed folder (demultiplexed), align and index, and send to the output directory (demux_sorted).
        - If there is an error here some of them are pretty self-explanatory, but if you call me I'll take a look.
        - After the alignment step, there should be a folder in nanopore_analysis/analysis/set2/demux_sorted containing
        a bam and a bam.bai file for each barcode.
        - There are more steps in our pipeline, but for the UXM deconvolution we don't need to do them. So from here we can do everything
        manually.
    - Now you have demultiplexed, sorted and indexed bam files for all of the barcodes. Now you can start deconvoluting.
##Deconvolution
    - You have a separate bam file for each barcode, you can run the steps we've already gone through.
    - Make sure wgbstools and uxm are in your PATH environment variable (you know how).
    - Make the pat files
        - bam2pat can take a whole directory of bam files. In our case you would want to do the following commands:
            mkdir -p analysis/combined/demux_pat/
            wgbstools bam2pat analysis/combined/demux_sorted/*.bam -output_dir analysis/combined/demux_pat/
        - Now there should be a load of pat files in the demux_pat folder.
    - Filter the pats using the atlas
        - You unfortunately need to do this one by one. First you can do mkdir -p analysis/combined/filtered/, then do
        the following for each barcode:
            wgbstools view analysis/combined/demux_pat/patname.pat.gz -L data/atlas/Atlas.U250.l4.hg38.full.tsv -o analysis/combined/filtered/filename.pat.
            You can give them any name you want, I would just name them after their barcode.
        - Now you have a load of pat files in the analysis/combined/filtered folder.
    - Index the pat files:
        - wgbstools index analysis/combined/filtered/*.pat
    - Use uxm to deconvolute the files one by one:
        uxm deconv analysis/combined/filtered/file.pat.gz --atlas data/atlas/Atlas.U250.l4.hg38.full.tsv -o analysis/combined/deconvoluted/barcodeXX.csv

To run uxm_deconv:
    - Get your pat.gz file and your atlas.tsv file.
    - Filter out the rows which aren't in the atlas using wgbstools view  path/to/pat.gz -L path/to/atlas -o path/to/output.pat
    - index the output
    - Make sure you have used the correct atlas (h g38/hg19) - it should be the same genome assembly used for the basecalling.
    - The atlas you used for the view function should be the same as the one you use for the deconv.
    - uxm deconv path/*pat.gz --atlas /path/to/atlas -o output.csv


The UXM and wgbs_tools files in the parent directory for the respective submodules should be symlinks, but sometimes they
end up as text files containing the path to the py file. To fix this:
    - rm path/to/file
    - ls -s src/uxm.py path/to/file

Process to install WGBS Tools:
    - At the moment this will be done seperately, but will be folded into the main project eventually.
    - Open conda and go to the directory you want to clone the repository into.
    - Run "git clone https://github.com/nloyfer/wgbs_tools.git" and do "cd wgbs_tools"
    - Make a new wgbs_tools conda environment:
        - conda create -n wgbs_build python=3.9
        - conda activate wgbs_build
        - conda install -c conda-forge gxx_linux-64
        - conda install -c conda-forge numpy pandas
    - You should be able to run python setup.py, but there is a good chance it won't work. If "which g++" gives you no output,
    you can fix it like so:
        - cd /home/{USERNAME}/miniconda3/envs/wgbs_build/bin/  (otherwise find out where your environments are).
        - ln -s x86_64-conda-linux-gnu-g++ g++
        - ln -s x86_64-conda-linux-gnu-gcc gcc
        - ln -s x86_64-conda-linux-gnu-cpp cpp
        - ln -s x86_64-conda-linux-gnu-c++ c++
    - The above commands point the g++, gcc, cpp and c++ commands to the relevant executables.
    - Now you can run setup.py and compile everything. Once setup.py is done you need to make wgbstools executable.


Before running the programme, activate the conda environment:
    conda activate nanopore_analysis
Sometimes, before running conda commands in a new WSL terminal you'll have to run:
    eval "$(path/to/conda shell.bash hook)"

Sample data can be gotten from s3://ont-open-data/. There might be a better way, but you can browse the data using the command:
    - aws s3 ls {folder_path} --no-sign-request
    - This will show the contents of the folder you provided the path for.

Use atlas from https://github.com/nloyfer/UXM_deconv/tree/main

Sometimes the meth_atlas submodule is an empty folder, and you'll get an error like "The directory <Project>/externals/meth_atlas
is registered as a Git root, but not Git repositories were found there". I'll find a more permanent fix, but to fix this:
    - git submodule sync
        - Reads the .gitmodules file and updates the local Git configuration to match
        - It should tell you "Synchronising submodule url for 'externals/meth_atlas"
    - git submodule update --init --recursive

When you're downloading a reference genome for alignment, you will need to download a reference genome with USCS-style
headers. e.g. chr1.

Download a full aligned "plate" using:
    aws s3 cp "s3://ont-open-data/gm24385_mod_2021.09/extra_analysis/alignment/20210510_1600_X1_FAQ32172_f02f2d1c.bam" "data/alignment_output/taken_alignment.bam" --no-sign-request

You need to install awscliv2:
    - Download the zip file using command: curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
    - Unzip using command: unzip awscliv2.zip
    - Install using command: sudo ./aws/install
This should be run on Linux/Ubuntu. On colab it's fine for testing purposes, but if you're running this on a Windows machine
Setup for WSL:
    Installation:
        run wsl --install in powershell
        Restart and open WSL. Create a username and password.
        Download and install miniconda. In the command line run:
         - wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
         - bash Miniconda3-latest-Linux-x86_64.sh
         - Accept the licence agreement. Say no when it asks you whether you want to always open in the conda environment
         - If you accidentally said yes to automatically initialise conda, you can undo that (after restarting WSL) with the command: conda config --set auto_activate_base false
        Close WSL to finish installation
    Setup the project to run in WSL:
        Open WSL again and go to the project folder
        Check that the environment.yml file is there, then run conda env create -f environment.yml to set up the Conda environment.
    If you have changed the environment.yml file or want to update the environment for some other reason, use the following command:
        - conda env update --file environment.yml --prune
    It should be now possible to run the pipeline by running main.py. There are some possible issues though:
        Errors like $'\r': command not found:
            I edited the bash code in Windows, but it's running from Ubuntu which means that the line endings are different.
            If you're lucky you can change the line endings in your code editor. For example, in Pycharm:
                On the bottom right status bar there is a button which has CRLF/LF/CR. Change this to LF and it'll have the right file endings.
            If not, you can do it manually:
                In the WSL terminal, install dos2unix so that you can convert from Windows to Unix-style LF line endings. Run the following in the terminal:
                    sudo apt-get update (connects to configured package repositories and gets the latest package lists)
                    sudo apt-get install dos2unix (installs dos2unix)
                When it's been installed, convert every bash file by running the following:
                    find . -type f -name "*.sh" -exec dos2unix {} \+;
                        find . - find something
                        -type f - find files
                        -name "*.sh" - find files which end with ".sh". These are bash files.
                        -exec dos2unix - execute the dos2unix command
                        {} - this is replaced with the list of filenames found
                        \+ - tells find to add all of the filenames
                        ; - end of the line.


With wgbs_tools:
    - To your aligned bam file, use bam2pat to generate pat and beta files
    - View any region you're interested in with:
        - wgbstools view -r chr1:910433-910476 path/to/file.pat.gz
        - wgbstools view -r chr1:910433-910476 --genome hg19 path/to/file.beta
    - Visualise with:
        - wgbstools vis *.beta -r chr1:22517933-22519650
        - wgbstools vis *.beta -r chr1:22517933-22519650 --heatmap
    - Visualise methylation patterns with:
        - wgbstools vis path/to/file.pat.gz -r chr1:22517933-22519650

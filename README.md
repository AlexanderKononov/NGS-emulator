# NGS-emulator
Tool to simulation of genome sequencing data with established properties to test and control different bioinformatics approaches.

## Idea
Main idea of the project is tool which can help biologist to create data semple with already known features such as clonality or chromosom abberations. Thereby scientist can estimate a works of some approaches, for instance compare prediction about clonality and originally clonality of sample given by this emulator.
So, Emulator should allowed users to create data with different DNA content (clonality, ploidy and crhromosome abberations), different levels of data (copy number data, B-allele fraction, single nucletide variation etc.) as well as differet size of sample data (lenghs and number of chromosomes). Opportunity to have small data size is especially important to test large computational approaches.

## User interfase
To realise noted flexibility, samle's features should be tune by wide spectrum of command and them attributes. To realise ease interaction between user and program, the process of tune was organized as "step by step" mode.

* First, user creates project, creates basic sample by command *'-normal'*.
    It can be imagined as the user take normal tissue sample on this step.
* Second, user changes basic sample by other commands such as *'-addclone', '-duplication', '-deletin'*.
    It is looked as the user molds or sculpts disered sample.
* Third, user get wishful data by command *'-seq'*.
    Given sample is used as matrix to simulation of sequencing process and after the simulation the user gets data of sample with known features.

For current version the execution script is *Cover_emulator.py* .
First step to create sample with 5 chromosomes 1000 dp lenght look such:

 `python Cover_emulator.py -normal 1000 5`

Second step can has next commands:

 `python Cover_emulator.py -addclone 0.5`

 `python Cover_emulator.py -duplication 1 2`

During this commands we add new clone in our sample with 0.5 fraction and then we add duplication of second chromosome in first clone added early.

 `python Cover_emulator.py -addclone 0.3`

 `python Cover_emulator.py -deletion 2 4`

During this commands we add new clone in our sample with 0.3 fraction. Now there are three cell fraction in our sample: clone 1 (f=0.5), clone 2 (f=0.3) and normal cells (f=0.2). Then we delete fourth cromosome in second clone added early

Third step look like next command:

 `python Cover_emulator.py -seq 15 4`

By this command we create CNA and BAF data with read depth about 15 per haploid (mean read cover of normal (diploid) sample is about 30) and varianse of read depth (nois) is 4.

Output data plases in program directori with extention *'.data.txt'*
For example : *file_CNA.data.txt* or *file_BAF.data.txt*

## List of available commands

  Command *'-normal'* is first command by which inner project files are created.
Originally this files have description for normal sample without any abberation. This script run with command *'-normal'* by 2-3 arguments: first is command *'-normal'* itself, second is length of DNA data per chromosome, third is optional argument and notes a number of chromosomes (default is 23).

  Command *'-seq'* is last command by which data from inner project files is transformed to output files.
This script run with command *'-seq'* by 3 arguments: first is command *'-seq'* itself, second is mean read depth of haploid sample and reflects number of read from one copy of DNA matrix, and third is level of noise. Then script creates data file by template created before.

  Command *'-addclone'* is command to adding clone.
This script run with command '-addclone' by 2 arguments: first command *'-addclone'* itself is noted and then new clone fraction (float value). Sum of all clone fraction must be less then 1.0 .

 Command *'-duplication'* is command to chromosome duplication.
This script run with command *'-duplication'* by 3 arguments: first is command *'-duplication'* itself then clone of changing (int value) and then chromosome of changing (int value). By default major chromosome is duplicated, but user can change it by additional argument *'minor'*.

  Command *'-deletion'* is command to chromosome deletion.
This script run with command *'-deletion'* by 3 arguments: first is command *'-deletion'* itself then clone of changing (int value) and then chromosome of changing (int value).By default minor chromosome is deleted, but user can change it by additional argument *'majore'*.

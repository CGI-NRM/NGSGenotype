# Fastrapper
- The FASTQ Bootstrapper

## Requirements
In order to run this program, erlang needs to be installed.

Arch based linux distribution installation:
```bash
$ sudo pacman -S erlang
```

## Execution
The program is run through a bash script. Navigate to this directory and run:
```bash
$ sh fastrapper.sh
```
The program will then look in the folder '../Raw_data' and load any FASTQ files
in there, and draw 10000 sequences (with repetition) from each forward and
reverse file and write the new files in the directory of the script. If your
FASTQ files are located somewhere else or if you want a different number of
bootstraps you can edit the file 'fastrapper.sh'.

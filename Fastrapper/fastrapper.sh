# Settings:
data_folder="../Raw_data" # Location of raw FASTQ files to be sampled
bootstraps=10000 # Number of lines (with repetition) that are to be sampled from each FASTQ file

# Compile and run erlang program:
erlc fastrapper.erl
erl -noshell -run fastrapper start $data_folder $bootstraps -s init stop

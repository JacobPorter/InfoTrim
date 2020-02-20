# InfoTrim

InfoTrim is a DNA read quality trimmer based on the Trimmomatic maximum information criterion model.  It adds an entropy term.
It works in both Python 2.7 and Python 3.X.  It is pure Python, so there is no pip installation.

If cited, use 

@inproceedings{porter2017infotrim,
  title={InfoTrim: A DNA read quality trimmer using entropy},
  author={Porter, Jacob and Zhang, Liqing},
  booktitle={Computational Advances in Bio and Medical Sciences (ICCABS), 2017 IEEE 7th International Conference on},
  pages={1--2},
  year={2017},
  organization={IEEE}
}

An extended paper can be found on the biorxiv at https://www.biorxiv.org/content/early/2017/10/11/201442

InfoTrim.py is the main program and can be run with Python.  Python must be installed.
InfoTrim_automate.py and the shell scripts in the shell directory were used to test the program with simulations, etc. as found in the papers.
The Utilities directory contains files necessary for the program to run such as a file with constants and a generic sequence iterator (SeqIterator.py)  that can iterate through FASTA and FASTQ files.

The following is the online help for the program.

Usage: InfoTrim.py [options] <FASTQ reads file location> 

Options:
  --version             show program's version number and exit
  
  -h, --help            show this help message and exit
  
  -p PROCESSES, --processes=PROCESSES
                        Number of processes to use. [default: 1]
                        
  -o OUTPUT, --output=OUTPUT
                        The file name for the output file. [default: stdout]
                        
  -s PHRED, --phred=PHRED
                        The tradeoff value for the phred score. [default: 0.1]
                        
  -r ENTROPY, --entropy=ENTROPY
                        The tradeoff value for the entropy feature. [default:
                        0.4]
                        
  -m MIN_LENGTH, --min_length=MIN_LENGTH
                        The minimum length target. [default: 25]
                        
  -t START_WEIGHT, --start_weight=START_WEIGHT
                        The weight for the exponential model for the starting
                        position of the read for the non-linear option.  If
                        set to 0, it will not be used. [default: 0]
                        
  -n, --non_linear      Use the non-linear approach to search for a substring
                        that maximizes the score.
                        
  -a ASCII_OFFSET, --ascii_offset=ASCII_OFFSET
                        The offset value for the ascii encoding of phred
                        scores. [default: 33]
                        
  -v, --verbose         Run in verbose.  The score for each read is printed to
                        stderr.


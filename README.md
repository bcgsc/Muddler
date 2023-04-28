# Muddler
Simulating Chromothripsis
# Dependencies and installation
Tested on python 3.9.12
Python packages: 
Biopython, Pybedtools
Versions tested:
Biopython=1.79
pybedtools=0.9.0

I used anaconda to install these into an environment with their dependencies.

`conda create -n muddler`

`conda activate muddler`

As biopython installs via pip rather than a conda recipe, you can install it to the conda environment

`conda install pip`

Then find the location of this installed pip. For example the path to pip could be:

`/anaconda/envs/muddler/bin/pip`

Installing biopython to the muddler environment

`/anaconda/envs/muddler/bin/pip install biopython`

For installing muddler

`git clone https://github.com/bcgsc/Muddler.git`

# Usage
Basic run command:

`python muddler_1.0.py --ref hg38.fa --out muddled_1 --chrome hg38chrom.bed --gaps hg38gaps.bed --frag hg38frag.bed --param param`

Paramater Info:
# Example usage with other simulators
-link to a wiki page with more details

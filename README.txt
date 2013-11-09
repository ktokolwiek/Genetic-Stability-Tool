## Genetic stability tool, v. 0.1 Alpha Americano.
## Gets you DNA sequences which are most unlike the original, but coding for
## the same protein.
## Go to http://2011.igem.org/Team:Edinburgh/Genetic_instability .
## Copyright (C) 2011 Lukasz Kopec, iGEM 2011 Team Edinburgh project
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## This tool includes bits of Biopython code, released under Biopython
## license. The license is available online: 
## http://www.biopython.org/DIST/LICENSE and in Bio folder.
## All Biopython code is (C) Copyright 2000 by Jeffrey Chang.
## All rights reserved.


Please enjoy this software responsibly, and report any feedback to:
e-mail:		L (dot) Kopec (at) sms (dot) ed (dot) ac (dot) uk
twitter:	@iGEMEdinburgh



usage: gen_stab.py [-h] [-f FILE] [-d DNA] [-n NOSEQ] [-g GEN] [-i IND]
                   [-c CROSS] [-m MUTATE] [-o CUTOFF] [-a {b,g,r}]

Genetic stability tool.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Take a DNA sequence from specified file in FASTA
                        format
  -d DNA, --dna DNA     Use the specified DNA sequence
  -n NOSEQ, --noseq NOSEQ
                        Number of sequences to generate
  -g GEN, --gen GEN     Number of generations of Genetic Algorithm
  -i IND, --ind IND     Number of individuals in Genetic Algorithm
  -c CROSS, --cross CROSS
                        Probability of crossing-over (0<c<1)
  -m MUTATE, --mutate MUTATE
                        Probability of mutation (0<m<1)
  -o CUTOFF, --cutoff CUTOFF
                        Cutoff point of inefficient E.coli codons (0<o<0.1)
  -a {b,g,r}, --action {b,g,r}
                        The mode of action - (b)est codon, (g)enetic algorthm,
                        (r)andom codon


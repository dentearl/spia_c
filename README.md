# 'SPIA IN C'
(c) October 2009 - 2012 The Author, see LICENSE.txt for details.

## Author
[Dent Earl](https://github.com/dentearl/)

## OVERVIEW
  SPIA is an acronym for Signalling Pathway Impact Analysis and is taken from
the 2009 paper by Tarca et al, A novel signaling pathway impact analysis. 
'SPIA in C', herein spia_c, is a C language implementation of the SPIA algorithm
based in part by the R code made available by Tarca et al. on their website,
http://bioinformaticsprb.med.wayne.edu/SPIA/

spia_c was written at UCSC by Dent Earl while working for Dr. Josh Stuart.
We wrote the version in C to facilitate as-fast-as-possible batch running of
SPIA for use with the UCSC Cancer Browser, or with the UCSC BioIntegrator.

IF your project is fairly simple, or doesn't require super-fast results, or
doesn't need to be in a pipeline, you may want to try out the original R version first.

## BACKGROUND
The SPIA method, described in the Tarca paper, involves a combination of a
traditional over-representation analysis (ORA) and a novel pathway analysis.
The pathway analysis is performed by grabbing pathways (the R version pulls
from Kegg, the C version requires you to provide pathways in the specific
Charlie-format), and then using the experiment results (both the test space,
meaning all genes that were on the array you used, and the genes that were
differentially expressed along with their log-fold change values) checks to
see if your experiment 'perturbed' a particular pathway.

## DETAILS
The ORA used in SPIA is the hypergeometeric. We get a p-value from the hypergeo
and from the pathway analysis and we combine them at the end.

<code>spia_c</code> loads all the tested genes and the differentially expressed genes into
hashes and then iteratively processes each pathway file in the pathway directory
(<code>--dir</code>), or in the UCSC Cancer Browser Analysis Team (UCSC CBAT)
format (<code>--pathFiles</code>).

First it is decomposed into 28 differentn by n matrices, where n is the 
number of genes in the pathway, and the number 28 comes from the number 
of different interaction types that Kegg uses (23) plus 5 custom 
types that Steve Benz and Charlie Vaske of (UCSC CBAT) use. The values in the matrices are
binary where a A_{ij} = 1 means that gene j effects gene i (the effect
being dependent on which of the 28 matrices. E.g. activation, inhibition, etc).
A_{ij} = 0 means that gene j does not effect gene i.

These matrices are compressed into a single matrix by normalizing the sum of each
column to 1 and then multiplying the matrix by the coeffecient of that relationship.
We use Tarca et al.'s coeffecient vector, but you could use anything. Additionally
they use binary values, but there is no reason why you couldn't use doubles should
you have cause to. HOWEVER, I never got around to writing the code for allowing
custom coefficent vectors.

The matrix B is then subtracted from the identity matrix, inverted, multiplied by the
B matrix and then multiplied by the differentially expressed vector:

Accumulation = B * (I-B)^{-1} * diffExpression

Further we know that

Accumulation = PerturbationFactor - diffExpression

And this allows us to solve for the perterbation factor. tA is the net accumulation,
and we sum the entire vector to calculate it.
tA = sum(Accumulation)

Next the bootstrap is applied to try to find an empirical distribution of tA scores.
This is where the PPert, or probability of perturbation, comes from.

## INSTALLATION
spia_c should be installed from the Makefile using

<code>$ make </code>

spia_c requires:

* uthash - (provided), a C-library for hashes.
* LAPACK and BLAS - Linear algebra libraries, used for inverting matrices.

## USE
spia_c requires a couple of things to run. It requires a directory
(<code>--dir</code>, or <code>--pathFiles</code> depending on format)
containing pathways in a tab-delimited, Charlie-formatted pathay files
(more information needed here). The files must end in '.tab'. It also requires a
tab-delimited file containing the differentially expressed genes from your 
experiment (geneName \t log-fold expression change). Lastly it requires a
tab-delimited file containing all genes that were tested in the experminent,
INCUDING the genes that were differentially expressed (order number 
{1,2,3,...} \t geneName). All geneNames may be strings. The net accumulation
data for each pathway may be silenced by using the <code>--quietNetAcc</code> command line
option. The number of bootstraps is specified with the <code>--nBoots</code> command line
option.

<code>bin/spia --dir testPathways/ --de testData/DE_Colorectal.tab  --array testData/ALL_Colorectal.tab --nBoots 2000</code>

<code>bin/spia --dir testPathways/ --de testData/DE_Colorectal.tab  --array testData/ALL_Colorectal.tab --nBoots 2000 --queitNetAcc</code>

There are also <code>--verbose</code> and <code>--debug</code> flags.

## CHARLIE FORMAT:
I don't really know what Charlie's format is exactly. In goes the Kegg xml, out
comes a tab delimited file, via kegg2tab.py (unknown location). One line looks like

<code>hsa:941	hsa:940	PPrel	activation	-->	path:hsa05330	Allograft rejection</code>

Here's my guess:
species:gene1 species:gene2 type(?) relationship relationshipPictogram path:speciesPathNumber commonName(??)
So from the example, Human gene 941 is activating gene 940 in pathway 05330.
We get the Charlie parse from
/projects/sysbio/Map/Data/Pathways/Human/KEGG/TabResults/

## UCSC CBAT FORMAT:
[fill this in!]

## REFERENCES
* Tarca et al. A novel signaling pathway impact analysis. Bioinformatics (2009) vol. 25 (1) pp. 75-82
* Draghici et al. A systems biology approach for pathway level analysis. Genome Res (2007) vol. 17 (10) pp. 1537-45
* Hassan et al. Signature pathways identified from gene expression profiles in the human uterine cervix before and after spontaneous term parturition. Am J Obstet Gynecol (2007) vol. 197 (3) pp. 250.e1-7
* Khatri et al. A System Biology Approach for the Steady-State Analysis of Gene Signaling Networks. 12th Iberoamerican Congress on Pattern Recognition, Valparasio, Chile (2007)

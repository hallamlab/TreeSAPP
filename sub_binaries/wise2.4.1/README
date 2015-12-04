

README for Wise2/Dynamite development


Wise2/Dynamite is written and maintained by Ewan Birney
<birney@ebi.ac.uk>.

The paper for Genewise is:

GeneWise and Genomewise. Birney E, Clamp M, Durbin R.
Genome Res. 2004 May;14(5):988-95.


The paper for Dynamite is:

Dynamite: a flexible code generating language for dynamic programming methods used in sequence comparison.
Birney E, Durbin R.  Proc Int Conf Intell Syst Mol Biol. 1997;5:56-64.


If you are using Genewise, almost certainly you want to
be also using Exonerate is protein2genome mode - for reasonably
close genomes (ie, inter-mammalian, or human/chicken), this
is as accurate as genewise and about 1,000 times faster.

Exonerate was written by Guy Slater and is available at:

http://www.ebi.ac.uk/~guy/exonerate/



INSTALLATION
-------------

cd into src and type

make all

possibly followed by

make test


there is not a make install. binaries are in src/bin after make.


The pthreads port no longer cleanly compiles. There was never
really that much point in using the pthreads port as (a) you
could trivially split databases and recombine the results 
for the same effect, and this was more sensible as it
would work on farm configurations and (b) it was excessively
long computation in anycase and probably you are better
off with Exonerate.


Development Notes:
------------------

Wise2 package is having a bit of a renaissance as my own (ie, Ewan's)
coding development package. You will see alot more experimental
code in these distributions and programs in development.


The old war-horse, genewise is being gently tweaked. There is far more
flexible splice site model. However currently this does not give
better results that plain old gt-ag for all cases. Because of that the
default is to use the gt-ag rule.


Other interesting programs that will come out from this include
promoterwise - a small region aligner designed for promoters
and scanwise - a new protein searching engine.


 





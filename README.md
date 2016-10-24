# Informed Proteomics

The Informed Proteomics project includes algorithms for analysis of both 
top down and bottom up MS/MS proteomics data.  Supported workflows include
both DDA and DIA.

## MSPathFinderT

MSPathFinder finds peptides in top-down LC-MS/MS datasets

Processing steps:

1. Run PbfGen.exe to convert the instrument file or .mzML file to an optimized binary file
	- Creates a .pbf file

2. Run ProMex on the .Pbf file to deisotope the data, including determine charge states
	- Creates a .ms1ft file

3. Run MSPathfinderT with the ? file and a fasta file to search for proteins
	- Creates _IcTda.tsv files

Example command lines:

PbfGen.exe -s MyDataset.raw

ProMex.exe -i MyDataset.pbf  -minCharge 2 -maxCharge 60 -minMass 3000 -maxMass 50000 -score n -csv n -maxThreads 0

MSPathFinderT.exe  -s MyDataset.pbf -feature MyDataset.ms1ft -d C:\FASTA\ProteinList.fasta -o C:\WorkFolder -t 10 -f 10 -m 1 -tda 1 -minLength 21 -maxLength 300 -minCharge 2 -maxCharge 30 -minFragCharge 1 -maxFragCharge 15 -minMass 3000 -maxMass 50000 -mod MSPathFinder_Mods.txt

## ProMex Syntax

ProMex version 1.0.5798 (November 16, 2015)
Usage: ProMex.exe
        [-i InputFolder or InputFile]
        [-o OutFolder (default : InputFolder)]
        [-minCharge MinCharge] (minimum charge state, default: 1)
        [-maxCharge MaxCharge] (maximum charge state, default: 60)
        [-minMass MinMassInDa] (minimum mass in Da, default: 600.0)
        [-maxMass MaxMassInDa] (maximum mass in Da, default: 50000.0)
        [-featureMap y or n (default: y)]
        [-score y or n (default: n)]
        [-maxThreads 0 (default: 0 (automatic set))]


## PbfGen Syntax

PbfGen version 1.0.5798 (November 16, 2015)
Usage: PbfGen.exe
        -s RawFilePath (*.raw or directory)
        [-o OutputDir]


## MSPathFinder Syntax ##

MSPathFinderT version 1.0.5798 (November 16, 2015)
Usage: MSPathFinderT.exe
        -s SpectrumFile (*.raw)
        -d DatabaseFile (*.fasta or *.fa)
        [-o OutputFolder]
        [-m SearchMode] (0: multiple internal cleavages, 1: single internal cleavage (default), 2: no internal cleavage)
        [-tagSearch 0/1] (1: include tag-based search (default), 0: skip tag-based search)
        [-mod ModificationFileName] (modification file, default: no modification)
        [-t PrecursorToleranceInPpm] (e.g. 10, Default: 10)
        [-f FragmentIonToleranceInPpm] (e.g. 10, Default: 10)
        [-tda 0/1] (0: don't search decoy database (default), 1: search shuffled decoy database)
        [-minLength MinSequenceLength] (minimum sequence length, default: 21)
        [-maxLength MaxSequenceLength] (maximum sequence length, default: 500)
        [-minCharge MinPrecursorCharge] (minimum precursor ion charge, default: 2)
        [-maxCharge MaxPrecursorCharge] (maximum precursor ion charge, default: 50)
        [-minFragCharge MinPrecursorCharge] (minimum fragment ion charge, default: 1)
        [-maxFragCharge MaxPrecursorCharge] (maximum fragment ion charge, default: 20)
        [-minMass MinSequenceMassInDa] (minimum sequence mass in Da, default: 3000.0)
        [-maxMass MaxSequenceMassInDa] (maximum sequence mass in Da, default: 50000.0)
        [-feature FeatureFile] (*.ms1ft, *_isos.csv, or *.msalign, default: Run ProMex)
        [-threads MaxNumThreads] (maximum number of threads, default: 0 automatic setting)

Enabling tag-based searching with "-tagSearch 1" can give 5% to 10% more matches, but can increase the runtime by 30% to 50%.

## MSPathFinder Mods File

See the Example_Files folder for sample modification definition files.

-------------------------------------------------------------------------------
Written by Sangtae Kim, Junkap Park, and Chris Wilkins for the Department of Energy (PNNL, Richland, WA)
Copyright 2015, Battelle Memorial Institute.  All Rights Reserved.

E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
Website: http://panomics.pnnl.gov/ or http://omics.pnl.gov
-------------------------------------------------------------------------------

Licensed under the Apache License, Version 2.0; you may not use this file except 
in compliance with the License.  You may obtain a copy of the License at 
http://www.apache.org/licenses/LICENSE-2.0

Notice: This computer software was prepared by Battelle Memorial Institute, 
hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
Department of Energy (DOE).  All rights in the computer software are reserved 
by DOE on behalf of the United States Government and the Contractor as 
provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
SOFTWARE.  This notice including this sentence must appear on any copies of 
this computer software.

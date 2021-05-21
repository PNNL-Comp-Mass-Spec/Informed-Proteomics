# Informed Proteomics

The Informed Proteomics project includes algorithms for proteomic mass spectrometry data analysis. Although the back-end data access and some of the scoring routines are general purpose, this repository is currently maintained for top down MS/MS datasets.

[![DOI](https://zenodo.org/badge/21950650.svg)](https://zenodo.org/badge/latestdoi/21950650)

## Install/Tutorials
#### [See the Informed Proteomics GitHub wiki for usage, tutorials, and more!](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/wiki)

## Downloads
https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/releases

### Continuous Integration

The latest versions of the Informed Proteomics tools are available on the [AppVeyor CI server](https://ci.appveyor.com/project/PNNLCompMassSpec/informed-proteomics/build/artifacts)

[![Build status](https://ci.appveyor.com/api/projects/status/j52ywc5d204gaxtp?svg=true)](https://ci.appveyor.com/project/PNNLCompMassSpec/informed-proteomics)

## MSPathFinderT

MSPathFinder finds peptides in top-down LC-MS/MS datasets. Similar to database search engines for bottom-up, it takes a FASTA file, a spectrum file, and a list of modifications as an input and reports proteoform spectrum matches (PsSMs) and their scores. These results are output in a tab-separated format and in a MzIdentML file.

Processing steps:

1. Run PbfGen.exe to convert the instrument file or .mzML file to an optimized binary file
  * Creates a .pbf file
    This file contains spectra information and full chromatograms for MS1 and MSn data, allowing fast access to extracted ion chromatograms during the search.
  * This step usually only needs to be performed once for a dataset, and by default will do nothing if it has previously been run on the dataset and can find the file.
  * ProMex and MSPathFinderT will perform this step automatically if they are given a spectrum file that is not a .pbf file.

2. Run ProMex on the .Pbf file to deisotope the data, including determine charge states
  * Creates a .ms1ft file
  * This file can be reused for multiple searches, as long as none of its parameters change
  * MSPathFinderT will perform this step automatically if not provided with a path to a feature file, using the respective parameters from the MSPathFinderT parameters.

3. Run MSPathfinderT with the .Pbf (or .raw) file and a FASTA file to search for proteins
  * Creates \_IcTda.tsv files and .mzid
  * Can also be run with the spectrum source file directly, with PbfGen and ProMex run as part of the process.

Example command lines:

`PbfGen.exe -s MyDataset.raw`

`ProMex.exe -i MyDataset.pbf  -minCharge 2 -maxCharge 60 -minMass 3000 -maxMass 50000 -score n -csv n -maxThreads 0`

`MSPathFinderT.exe  -s MyDataset.pbf -feature MyDataset.ms1ft -d C:\FASTA\ProteinList.fasta -o C:\WorkFolder -t 10 -f 10 -m 1 -tda 1 -minLength 21 -maxLength 300 -minCharge 2 -maxCharge 30 -minFragCharge 1 -maxFragCharge 15 -minMass 3000 -maxMass 50000 -mod MSPathFinder_Mods.txt`

### Results viewer / GUI
For viewing search results, you might want to consider [LCMS-Spectator](https://github.com/PNNL-Comp-Mass-Spec/LCMS-Spectator). It can also function as a GUI for running ProMex and MSPathFinder.

## Running on Linux

PbfGen, ProMex, and MSPathFinderT can be run on Linux using [Mono](https://www.mono-project.com/download/stable/)

Example command lines:
```
mono PbfGen.exe -s *.raw
mono PbfGen.exe -s *.mzML
mono PbfGen.exe -s *.mzML.gz

mono ProMex.exe -i *.pbf -minCharge 2 -maxCharge 60 -minMass 2000 -maxMass 50000 -score n -csv n -maxThreads 0

mono MSPathFinder/MSPathFinderT.exe -s *.pbf -d ID_006407_8F27399B.fasta -o . -ParamFile MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage.txt
```

## ProMex Syntax

```
ProMex version 1.0.6232 (January 23, 2017)
Usage: ProMex.exe
        [-i InputFolder or InputFile]
        [-o OutFolder (default : InputFolder)]
        [-minCharge MinCharge] (minimum charge state, default: 1)
        [-maxCharge MaxCharge] (maximum charge state, default: 60)
        [-minMass MinMassInDa] (minimum mass in Da, default: 2000.0)
        [-maxMass MaxMassInDa] (maximum mass in Da, default: 50000.0)
        [-featureMap y or n (default: y)]
        [-score y or n (default: n)]
        [-maxThreads 0 (default: 0 (automatic set))]

Syntax to create a PNG of the features in an existing ms1ft file
(requires both a .pbf file and a .ms1ft file)
ProMex.exe
        [-i PbfFile]
        [-o OutFolder (default : InputFolder)]
        [-minMass MinMassInDa] (minimum mass in Da, default: 2000.0)
        [-maxMass MaxMassInDa] (maximum mass in Da, default: 50000.0)
        [-ms1ft FeaturesFilePath (use a period to infer the name from the pbf file)]
```

## PbfGen Syntax

```
PbfGen version 1.0.6232 (January 23, 2017)
Usage: PbfGen.exe
        -s RawFilePath (*.raw or directory)
        [-o OutputDir]
```

## MSPathFinder Syntax

```
MSPathFinderT version 1.0.7569 (September 21, 2020)

Usage: MSPathFinderT.exe

  -i, -s, -specFile,     Required. Spectrum File (.raw or .pbf)
  arg#1

  -d, -database          Required. Database File (*.fasta or *.fa or *.faa)

  -o, -outputDir         Output Directory

  -ic                    Search Mode (Default: SingleInternalCleavage (or 1))
                         Possible values are:
                           0 or 'NoInternalCleavage': No Internal Cleavage
                           1 or 'SingleInternalCleavage': Single Internal Cleavage
                           2 or 'MultipleInternalCleavages': Multiple Internal
                           Cleavages

  -tagSearch             Include Tag-based Search (use true or false;
                         or use '0' for false or '1' for true) (Default: True)

  -memMatches            Number of matches to keep in memory; these matches are
                         used when computing spectral E-values (Default: 3)

  -n,                    Number of results to report for each mass spectrum
  -NumMatchesPerSpec     (Default: 1)

  -IncludeDecoy,         Include decoy results in the _IcTda.tsv file
  -IncludeDecoys         (Default: False)

  -mod                   Path to modification file that defines static and dynamic
                         modifications. Modifications can alternatively be defined
                         in a parameter file, as specified by /ParamFile or
                         -ParamFile
                         Modifications defined using the -mod switch take
                         precedence over modifications defined in a parameter file
                         (Default: empty string, meaning no modifications)

  -tda                   Database search mode:
                         0: don't search decoy database,
                         1: search shuffled decoy database
                         (Default: 0, Min: -1, Max: 1)

  -overwrite             Overwrite existing results. If false (default), looks for
                         files _IcTarget.tsv and _IcDecoy.tsv and uses the results
                         in the files if found (Default: False)

  -t, -precursorTol,     Precursor Tolerance (in PPM) (Default: 10)
  -PMTolerance

  -f, -fragmentTol,      Fragment Ion Tolerance (in PPM) (Default: 10)
  -FragTolerance

  -minLength             Minimum Sequence Length (Default: 21, Min: 0)

  -maxLength             Maximum Sequence Length (Default: 500, Min: 0)

  -minCharge             Minimum precursor ion charge (Default: 2, Min: 1)

  -maxCharge             Maximum precursor ion charge (Default: 50, Min: 1)

  -minFragCharge         Minimum fragment ion charge (Default: 1, Min: 1)

  -maxFragCharge         Maximum fragment ion charge (Default: 20, Min: 1)

  -minMass               Minimum sequence mass in Da (Default: 3000)

  -maxMass               Maximum sequence mass in Da (Default: 50000)

  -feature               .ms1ft, _isos.csv, or .msalign feature file (typically
                         the results from ProMex); leave blank/undefined if
                         processing multiple input files

  -threads               Maximum number of threads, or 0 to set automatically
                         (Default: 0, Min: 0)

  -act,                  Activation Method (Default: Unknown (or 6))
  -ActivationMethod      Possible values are:
                           0 or 'CID'
                           1 or 'ETD'
                           2 or 'HCD'
                           3 or 'ECD'
                           4 or 'PQD'
                           5 or 'UVPD'
                           6 or 'Unknown'

  -scansFile,            Optional text file with MS2 scans to process (tab, comma,
  -ScansFilePath         or space separated); any integer in the file is assumed
                         to be a scan number to process

  -flip                  If specified, FLIP scoring code will be used
                         (supports UVPD spectra) (Default: False)

  -ParamFile             Path to a file containing program parameters. Additional
                         arguments on the command line can supplement or override
                         the arguments in the param file. Lines starting with '#'
                         or ';' will be treated as comments; blank lines are
                         ignored. Lines that start with text that does not match a
                         parameter will also be ignored.

  -CreateParamFile       Create an example parameter file. Can supply a path; if
                         path is not supplied, the example parameter file content
                         will output to the console.

  NOTE:                  arg#1, arg#2, etc. refer to positional arguments, used
                         like "AppName.exe [arg#1] [arg#2] [other args]".
```

Enabling tag-based searching with `-tagSearch 1` can give 5% to 10% more matches, but can increase the runtime by 30% to 50%.

### Supported file formats

These programs, when used with no other software installed, only support the use of centroid mzML files as spectrum input.

Versions prior to February 1, 2019: If Thermo Finnigan MSFileReader is installed, it also supports reading from Thermo Finnigan .raw files ([Download here](https://thermo.flexnetoperations.com/control/thmo/download?element=6306677), requires registration to download).

Versions after February 1, 2019: If running as a 64-bit program, reading from Thermo Finnigan .raw files is supported via the included RawFileReader dlls (no additional software install required).

On Windows, several other formats are supported if an appropriate version of ProteoWizard is installed ([Download here](http://proteowizard.sourceforge.net/downloads.shtml), make sure the version downloaded matches system architecture)

On Linux, supported input files are .raw, .mzML, and .mzML.gz

## MSPathFinder Parameter Files

See the [Example_Files](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/tree/master/Example_Files) folder for sample parameter files
* Example command for invoking MSPathFinder with a parameter file:
```
MSPathFinderT.exe -s C:\WorkDir\Dataset.pbf -feature C:\WorkDir\Dataset.ms1ft -d C:\WorkDir\Proteins.fasta -o C:\WorkDir /ParamFile:C:\WorkDir\MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage_ReportTop2.txt
```

## MSPathFinder Mods File

See the [Example_Files](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/tree/master/Example_Files) folder for sample modification definition files.

## System Requirements

Minimum required:
* .NET 4.7.2

Minimum recommended:
* 2.4 GHz, quad-core CPU
* 16 GB RAM
* Windows 7 or newer
* 250 GB hard drive

-------------------------------------------------------------------------------
Written by Sangtae Kim, Junkap Park, and Chris Wilkins for the Department of Energy (PNNL, Richland, WA)
Copyright 2015, Battelle Memorial Institute.  All Rights Reserved.

E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/\

-------------------------------------------------------------------------------

Licensed under the Apache License, Version 2.0; you may not use this file except
in compliance with the License.  You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved.

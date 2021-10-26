# Informed Proteomics

The Informed Proteomics project includes algorithms for proteomic mass spectrometry data analysis. 
Although the back-end data access and some of the scoring routines are general purpose, 
this repository is currently maintained for top down MS/MS datasets.

[![DOI](https://zenodo.org/badge/21950650.svg)](https://zenodo.org/badge/latestdoi/21950650)

## Manuscript 

Implementation details are described in the manuscript "Informed-Proteomics: open-source software package for top-down proteomics", published in Nature Methods
* [PMID 28783154](https://pubmed.ncbi.nlm.nih.gov/28783154/)
* [Manuscript text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5578875/) in PMC
* [DOI 10.1038/nmeth.4388](https://doi.org/10.1038/nmeth.4388)

## Install/Tutorials
#### [See the Informed Proteomics GitHub wiki for usage, tutorials, and more!](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/wiki)

## Downloads
https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/releases

### Continuous Integration

The latest versions of the Informed Proteomics tools may be available on the [AppVeyor CI server](https://ci.appveyor.com/project/PNNLCompMassSpec/informed-proteomics/build/artifacts), though they get auto-deleted after 6 months.

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

4. Optionally use ProMexAlign to align MS1 features between datasets
  * Creates a .tsv file listing a consolidated list of MS1 features, plus the corresponding Feature IDs from the input files
  * A FeatureID of 0 means the feature was not present

Example command lines:

`PbfGen.exe -s MyDataset.raw`

`ProMex.exe -i MyDataset.pbf -minCharge 2 -maxCharge 60 -minMass 3000 -maxMass 50000 -score n -csv n -maxThreads 0`

`MSPathFinderT.exe -s MyDataset.pbf -feature MyDataset.ms1ft -d C:\FASTA\ProteinList.fasta -o C:\WorkDir -t 10 -f 10 -m 1 -tda 1 -minLength 21 -maxLength 300 -minCharge 2 -maxCharge 30 -minFragCharge 1 -maxFragCharge 15 -minMass 3000 -maxMass 50000 -mod MSPathFinder_Mods.txt`

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

## PbfGen Syntax

```
PbfGen.exe
    -s:RawFilePath (*.raw or directory)
    [-o:OutputDir]
    [-start:StartScan] [-end:EndScan]
```

`-i` or `-s` or `-InputFile`
* Input file or input directory
* Supports .pbf, .mzML, and several vendor formats (see documentation)

`-o` or `-OutputDirectory`
* Output directory.
* Default: directory containing input file

`-start`
* Start scan number
* Optionally use to limit scan range included in .pbf file

`-end`
* End scan number
* Optionally use to limit scan range included in .pbf file

## ProMex Syntax

```
ProMex.exe 
  -i:InputFile [-o:OutputDirectory]
  [-MinCharge:Value] [-MaxCharge:Value]
  [-MinMass:Value]   [-MaxMass:Value]
  [-FeatureMap:Flag] [-Score:Flag]
  [MaxThreads:Value] [-csv:Flag]
  [BinResPPM:Value]  [ScoreTh:Value]
```

`-i` or `-s` or `-InputFile`
* Input file or input directory
  * Supports .pbf, .mzML, and several vendor formats (see documentation)

`-o` or `-OutputDirectory`
* Output directory. 
* Default: directory containing input file

`-MinCharge`
* Minimum charge state. 
* Default: 1
* Min: 1, Max: 60

`-MaxCharge`
* Maximum charge state
* Default: 60
* Min: 1, Max: 60

`-MinMass`
* Minimum mass, in Daltons
* Default: 2000
* Min: 600, Max: 100000

`-MaxMass`
* Maximum mass, in Daltons
* Default: 50000
* Min: 600, Max: 100000

`-FeatureMap`
* Output the feature heatmap. Defaults to true.
* To disable, use `-FeatureMap:false` or include 'FeatureMap=False' in a parameter file

`-Score`
* Output extended scoring information
* Default: False

`-MaxThreads`
* Max number of threads to use
* Default: 0, meaning automatically determine the number of threads to use

`-csv`
* Also write feature data to a CSV file
* Default: False

`-BinResPPM` or `-BinningResolutionPPM`
* Binning resolution, in ppm
  * This is used when finding LC-MS features
  * Data points whose mass difference is less than the resolution are grouped together
* Allowed values are 1, 2, 4, 8, 16, 32, 64, or 128
* Default: 16

`-ScoreTh` or `-ScoreThreshold`
* Likelihood score threshold
* Default: -10

`-ms1ft`
* Name of a ms1ft feature file to use to create a feature plot (the corresponding .pbf file must be included in the same directory)
* Use `-ms1ft.` to infer the name from the pbf file

`-ParamFile`
* Path to a file containing program parameters. 
* Additional arguments on the command line can supplement or override the arguments in the param file. 
* Lines starting with '#' or ';' will be treated as comments; blank lines are ignored. 
* Lines that start with text that does not match a parameter will also be ignored.

`-CreateParamFile`
* Create an example parameter file.
* Can supply a path; if path is not supplied, the example parameter file content will output to the console.

### Create PNG using existing ms1ft file

Command to create a PNG of the features in an existing ms1ft file
(requires both a .pbf file and a .ms1ft file):

`ProMex.exe -i dataset.pbf -ms1ft dataset.ms1ft -featureMap`

## ProMexAlign syntax

```
ProMexAlign.exe DatasetInfoFile
```

* The dataset info file is a tab-delimited text file with either 3 or 4 columns of information
  * Expected columns: 
  * `Label`  `RawFilePath`  `Ms1FtFilePath`  `MsPathfinderIdFilePath`
* If data in the Label column is an empty string, dataset labels will be auto-assigned as `Dataset_1`, `Dataset_2`, etc.
* The raw files are either Thermo .raw files or .pbf files created by PbfGen
* The MsPathfinderIdFilePath column is optional
  * If \_IcTda.tsv files are listed in this column, additional columns will appear in the output file created by ProMexAlign

Example input file:

| Label | RawFilePath | Ms1FtFilePath | MsPathfinderIdFilePath |
|-------|-------------|---------------|------------------------|
| Intact\_Run8  | Intact\_100ng\_Run8.pbf  | Intact\_100ng\_Run8.ms1ft |
| Intact\_Run9  | Intact\_100ng\_Run9.pbf  | Intact\_100ng\_Run9.ms1ft |
| Intact\_Run10 | Intact\_100ng\_Run10.pbf | Intact\_100ng\_Run10.ms1ft |

## MSPathFinder Syntax

```
MSPathFinderT.exe
  -i:InputFile [-d:FastaFile] 
  [-o:OutputDirectory] [-overwrite:Flag]
  [-mod:ModificationFilePath]
  [-ic:InternalCleavageMode]
  [-tda:DatabaseSearchMode]
  [-tagSearch:Flag]
  [-t:PrecursorTolerance] 
  [-f:FragmentIonTolerance]
  [-MemMatches:Count] 
  [-n:NumMatchesPerSpec]
  [-IncludeDecoys:Flag] 
  [-MinLength:Value] [-MaxLength:Value]
  [-MinCharge:Value] [-MaxCharge:Value]
  [-MinFragCharge:Value] [-MaxFragCharge:Value]
  [-MinMass:Value] [-MaxMass:Value]
  [-Feature:FeatureFile]
  [-Threads:Value]
  [-act:ActivationMethod]
  [-ScansFile:MS2ScanFilterFile]
  [-flip:Flag]
  [-ParamFile]
  [-CreateParamFile]
```

`-i` or `-s` or `-specFile`
* Spectrum File (.raw or .pbf)

`-d` or `-database`
* Database File (*.fasta or *.fa or *.faa)

`-o` or `-outputDir`
* Output Directory

`-ic`
* Internal Cleavage Search Mode
* Default: SingleInternalCleavage (or 1)
* Possible values:
  * 0 or 'NoInternalCleavage': No Internal Cleavage
  * 1 or 'SingleInternalCleavage': Single Internal Cleavage
  * 2 or 'MultipleInternalCleavages': Multiple Internal Cleavages

`-TagSearch`
* Enable/disable Tag-based Search
* Use true or false
  * Or use '0' for false or '1' for true
  * Or use 'n' for false or 'y' for true
* Default: True

`-MemMatches`
* Number of matches to keep in memory
  * These matches are used when computing spectral E-values
* Default: 3

`-n` or `-NumMatchesPerSpec`
* Number of results to report for each mass spectrum
* Default: 1

`-IncludeDecoy` or `-IncludeDecoys`
* Include decoy results in the \_IcTda.tsv file
* Default: False

`-mod`
* Path to modification file that defines static and dynamic modifications.
* Modifications can alternatively be defined in a parameter file, as specified by `/ParamFile` or `-ParamFile`
* Modifications defined using the `-mod` switch take precedence over modifications defined in a parameter file
* Default: empty string, meaning no modifications

`-tda`
* Database search mode; can be 0, 1, or -1
  * 0: target search only
  * 1: target and shuffled decoy database
  * -1: only search shuffled decoy database
* Default: 0

`-Overwrite`
* Overwrite existing results
* If false, looks for files \_IcTarget.tsv and \_IcDecoy.tsv and uses the results 
* Default: False

`-t` or `-PrecursorTol` or `-PMTolerance`
* Precursor Tolerance (in PPM)
* Default: 10

`-f` or `-FragmentTol` or `-FragTolerance`
* Fragment Ion Tolerance (in PPM)
* Default: 10

`-MinLength`
* Minimum Sequence Length
* Default: 21

`-MaxLength`
* Maximum Sequence Length
* Default: 500

`-MinCharge`
* Minimum precursor ion charge
* Default: 2

`-MaxCharge`
* Maximum precursor ion charge
* Default: 50

`-MinFragCharge`
* Minimum fragment ion charge
* Default: 1

`-MaxFragCharge`
* Maximum fragment ion charge
* Default: 20

`-MinMass`
* Minimum sequence mass in Da
* Default: 3000

`-MaxMass`
* Maximum sequence mass in Da
* Default: 50000

`-Feature`
* Path to a .ms1ft, \_isos.csv, or .msalign feature file (typically the results from ProMex)
* Leave blank/undefined if processing multiple input files

`-threads`
* Maximum number of threads, or 0 to set automatically
* Default: 0

`-act` or `-ActivationMethod`
* Activation Method; possible values:
  * 0 or 'CID'
  * 1 or 'ETD'
  * 2 or 'HCD'
  * 3 or 'ECD'
  * 4 or 'PQD'
  * 5 or 'UVPD'
  * 6 or 'Unknown'
* Default: Unknown (or 6)


`-ScansFile` or `-ScansFilePath`
* Optional text file with MS2 scans to process (tab, comma, or space separated)
* Any integer in the file is assumed to be a scan number to process

`-Flip`
* If specified, FLIP scoring code will be used (experimental, supports UVPD spectra)
* Default: False

`-ParamFile`
* Path to a file containing program parameters. 
* Additional arguments on the command line can supplement or override the arguments in the param file. 
* Lines starting with '#' or ';' will be treated as comments; blank lines are ignored. 
* Lines that start with text that does not match a parameter will also be ignored.

`-CreateParamFile`
* Create an example parameter file.
* Can supply a path; if path is not supplied, the example parameter file content will output to the console.

Enabling tag-based searching with `-tagSearch 1` can give 5% to 10% more matches, but can increase the runtime by 30% to 50%.

### Supported file formats

These suggested input file format is centroided .mzML files.

For versions released after February 1, 2019:
* If running PbfGen and ProMex as 64-bit programs, reading from Thermo .raw files is supported via the included RawFileReader DLLs (no additional software is required).

For versions released prior to February 1, 2019: 
* If Thermo's MSFileReader is installed, PbfGen and ProMex also support reading from Thermo .raw files ([Download here](https://thermo.flexnetoperations.com/control/thmo/download?element=6306677), requires registration to download).

On Windows, several other formats are supported if an appropriate version of ProteoWizard is installed ([Download here](https://proteowizard.sourceforge.io/), make sure the version downloaded matches system architecture)

On Linux, supported input file formats are .raw, .mzML, and .mzML.gz

## MSPathFinder Parameter Files

See the [Example_Files](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/tree/master/Example_Files) directory on GitHub for sample parameter files
* Example command for invoking MSPathFinder with a parameter file:
```
MSPathFinderT.exe -s C:\WorkDir\Dataset.pbf -feature C:\WorkDir\Dataset.ms1ft -d C:\WorkDir\Proteins.fasta -o C:\WorkDir /ParamFile:C:\WorkDir\MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage_ReportTop2.txt
```

## MSPathFinder Mods File

See the [Example_Files](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics/tree/master/Example_Files) directory for sample modification definition files.

## System Requirements

Minimum required:
* .NET 4.7.2

Minimum recommended:
* 2.4 GHz, quad-core CPU
* 16 GB RAM
* Windows 7 or newer
* 250 GB hard drive

## Contacts

Written by Sangtae Kim, Junkap Park, and Chris Wilkins for the Department of Energy (PNNL, Richland, WA)\
Copyright 2015, Battelle Memorial Institute.  All Rights Reserved.\
E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov\
Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics/

## License

Licensed under the Apache License, Version 2.0; you may not use this program except
in compliance with the License.  You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved.

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using CsvHelper;
using CsvHelper.Configuration;
using InformedProteomics.Backend.Data.Sequence;
using PRISM;
using PSI_Interface.IdentData;

namespace InformedProteomics.Backend.SearchResults
{
    /// <summary>
    /// Container with utilities for reading/writing/storing database search results.
    /// </summary>
    public class DatabaseSearchResultData
    {
        // Ignore Spelling: Dehydro, desc, pre

        /// <summary>
        /// Scan number
        /// </summary>
        public int ScanNum { get; set; }

        /// <summary>
        /// Pre residue
        /// </summary>
        public string Pre { get; set; }

        /// <summary>
        /// Peptide sequence
        /// </summary>
        public string Sequence { get; set; }

        /// <summary>
        /// Post residue
        /// </summary>
        public string Post { get; set; }

        /// <summary>
        /// Name and location of modifications
        /// </summary>
        public string Modifications { get; set; }

        /// <summary>
        /// Match elemental composition (including modifications)
        /// </summary>
        public string Composition { get; set; }

        /// <summary>
        /// Name of Protein
        /// </summary>
        public string ProteinName { get; set; }

        /// <summary>
        /// Protein Description
        /// </summary>
        public string ProteinDescription { get; set; }

        /// <summary>
        /// Length of protein
        /// </summary>
        public int ProteinLength { get; set; }

        /// <summary>
        /// Start index of sequence in protein
        /// </summary>
        public int Start { get; set; }

        /// <summary>
        /// End index of sequence in protein
        /// </summary>
        public int End { get; set; }

        /// <summary>
        /// Charge
        /// </summary>
        public int Charge { get; set; }

        /// <summary>
        /// m/z of most abundant isotope
        /// </summary>
        public double MostAbundantIsotopeMz { get; set; }

        /// <summary>
        /// Calculated mass (monoisotopic m/z)
        /// </summary>
        public double Mass { get; set; }

        /// <summary>
        /// Number of matched fragments
        /// </summary>
        public int NumMatchedFragments { get; set; }

        /// <summary>
        /// Match Probability
        /// </summary>
        public double Probability { get; set; }

        /// <summary>
        /// SpecEValue
        /// </summary>
        public double SpecEValue { get; set; }

        /// <summary>
        /// EValue
        /// </summary>
        public double EValue { get; set; }

        /// <summary>
        /// QValue
        /// </summary>
        public double QValue
        {
            get => qValue;
            set
            {
                qValue = value;
                HasTdaScores = true;
            }
        }

        /// <summary>
        /// PepQValue
        /// </summary>
        public double PepQValue
        {
            get => pepQValue;
            set
            {
                pepQValue = value;
                HasTdaScores = true;
            }
        }

        /// <summary>
        /// If the FDR scores (QValue and PepQValue) have been set
        /// </summary>
        public bool HasTdaScores { get; private set; }

        /// <summary>
        /// Integer ID to associate with this result; initially 0
        /// </summary>
        public int ResultID { get; set; }

        /// <summary>
        /// The sequence, with the pre and post residues
        /// </summary>
        public string SequenceWithEnds
        {
            get
            {
                if (string.IsNullOrWhiteSpace(sequenceWithEnds))
                {
                    sequenceWithEnds = Pre + "." + Sequence + "." + Post;
                }
                return sequenceWithEnds;
            }
        }

        /// <summary>
        /// Ms1Feature id
        /// </summary>
        public int Ms1Feature { get; set; }

        private string sequenceWithEnds;
        private double qValue;
        private double pepQValue;

        /// <summary>
        /// Constructor
        /// </summary>
        public DatabaseSearchResultData()
        {
            Pre = "-";
            Post = "-";
            // Set to int.MaxValue instead of double.MaxValue to avoid bloat in the output if the value hasn't been set.
            SpecEValue = int.MaxValue;
            EValue = int.MaxValue;
            qValue = int.MaxValue;
            pepQValue = int.MaxValue;
            HasTdaScores = false;
        }

        private const double SmallestValueExcel = 9.99E-308;

        /// <summary>
        /// Write the resultData in TSV format to the specified path, possibly including FDR scores
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="resultData"></param>
        /// <param name="includeTdaScores">If FDR scores should be output also</param>
        public static void WriteResultsToFile(string filePath, IEnumerable<DatabaseSearchResultData> resultData, bool includeTdaScores = false)
        {
            var outputFile = new FileInfo(filePath);
            if (outputFile.Directory != null && !outputFile.Directory.Exists)
            {
                outputFile.Directory.Create();
            }

            using (var tsv = new CsvWriter(new StreamWriter(new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.Read)), CultureInfo.InvariantCulture))
            {
                SetCsvWriterConfig(tsv.Configuration);
                if (includeTdaScores)
                {
                    tsv.Configuration.RegisterClassMap<DatabaseSearchResultDataTdaMap>();
                }
                else
                {
                    tsv.Configuration.RegisterClassMap<DatabaseSearchResultDataMap>();
                }
                tsv.WriteRecords(resultData);
            }
        }

        /// <summary>
        /// Read in result data from the specified TSV format file.
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static List<DatabaseSearchResultData> ReadResultsFromFile(string filePath)
        {
            List<DatabaseSearchResultData> results;
            using (var stream = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                // check for QValue header, to determine the mapping used for reading
                var hasTda = false;
                var line = stream.ReadLine();
                if (!string.IsNullOrWhiteSpace(line) && line.IndexOf("QValue", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    hasTda = true;
                }
                stream.BaseStream.Seek(0, SeekOrigin.Begin);
                stream.DiscardBufferedData();
                using (var tsv = new CsvReader(stream, CultureInfo.InvariantCulture))
                {
                    SetCsvReaderConfig(tsv.Configuration);
                    if (hasTda)
                    {
                        tsv.Configuration.RegisterClassMap<DatabaseSearchResultDataTdaMap>();
                    }
                    else
                    {
                        tsv.Configuration.RegisterClassMap<DatabaseSearchResultDataMap>();
                    }
                    results = tsv.GetRecords<DatabaseSearchResultData>().ToList();
                }
            }
            if (results.Count == 0)
            {
                return null;
            }
            return results;
        }

        private static void SetCsvReaderConfig(IReaderConfiguration config)
        {
            config.Delimiter = "\t";

            config.PrepareHeaderForMatch = (header, _) => header?.Trim().ToLower();

            config.HeaderValidated = null;
            config.MissingFieldFound = null;
            //config.BadDataFound = null;
            config.Comment = '#';
            config.AllowComments = true;
        }

        private static void SetCsvWriterConfig(IWriterConfiguration config)
        {
            config.Delimiter = "\t";
            config.Comment = '#';
            config.AllowComments = true;
        }

        // ReSharper disable once ClassNeverInstantiated.Local
        private class DatabaseSearchResultDataMap : ClassMap<DatabaseSearchResultData>
        {
            protected int ColumnCount;

            // ReSharper disable once MemberCanBeProtected.Local
            public DatabaseSearchResultDataMap()
            {
                ColumnCount = 0;
                // ReSharper disable VirtualMemberCallInConstructor
                Map(x => x.ScanNum).Index(ColumnCount++).Name("Scan");
                Map(x => x.Pre).Index(ColumnCount++).Name("Pre");                                                                                                                          // Pre
                Map(x => x.Sequence).Index(ColumnCount++).Name("Sequence", "Peptide", "#Peptide");                                                                                         // Sequence
                Map(x => x.Post).Index(ColumnCount++).Name("Post");                                                                                                                        // Post
                Map(x => x.Modifications).Index(ColumnCount++).Name("Modifications");                                                                                                      // Modifications
                Map(x => x.Composition).Index(ColumnCount++).Name("Composition");                                                                                                          // Composition
                Map(x => x.ProteinName).Index(ColumnCount++).Name("ProteinName", "Protein");                                                                                               // ProteinName
                Map(x => x.ProteinDescription).Index(ColumnCount++).Name("ProteinDesc");                                                                                                   // ProteinDescription
                Map(x => x.ProteinLength).Index(ColumnCount++).Name("ProteinLength");                                                                                                      // ProteinLength
                Map(x => x.Start).Index(ColumnCount++).Name("Start");                                                                                                                      // Start position in protein
                Map(x => x.End).Index(ColumnCount++).Name("End");                                                                                                                          // End position in protein
                Map(x => x.Charge).Index(ColumnCount++).Name("Charge");                                                                                                                    // precursorCharge
                Map(x => x.MostAbundantIsotopeMz).Index(ColumnCount++).Name("MostAbundantIsotopeMz").ConvertUsing(x => StringUtilities.DblToString(x.MostAbundantIsotopeMz, 9, true));     // MostAbundantIsotopeMz
                Map(x => x.Mass).Index(ColumnCount++).Name("Mass").ConvertUsing(x => StringUtilities.DblToString(x.Mass, 9, true));                                                        // Mass
                Map(x => x.Ms1Feature).Index(ColumnCount++).Name("Ms1Features");
                Map(x => x.NumMatchedFragments).Index(ColumnCount++).Name("#MatchedFragments", "Score");                                                                                   // (Number of matched fragments)
                Map(x => x.Probability).Index(ColumnCount++).Name("Probability").ConvertUsing(x => StringUtilities.DblToString(x.Probability, 4));                                         // Probability
                Map(x => x.SpecEValue).Index(ColumnCount++).Name("SpecEValue").ConvertUsing(x => StringUtilities.DblToString(Math.Max(SmallestValueExcel, x.SpecEValue), 6, true, 0.001)); // SpecEValue; will be displayed using scientific notation if the value is less than 0.001
                Map(x => x.EValue).Index(ColumnCount++).Name("EValue").ConvertUsing(x => StringUtilities.DblToString(Math.Max(SmallestValueExcel, x.EValue), 6, true, 0.001));             // EValue; will be displayed using scientific notation if the value is less than 0.001
                // ReSharper restore VirtualMemberCallInConstructor
            }
        }

        // ReSharper disable once ClassNeverInstantiated.Local
        private sealed class DatabaseSearchResultDataTdaMap : DatabaseSearchResultDataMap
        {
            // Adds in the columns expected for TDA results
            public DatabaseSearchResultDataTdaMap()
            {
                Map(x => x.QValue).Index(ColumnCount++).Name("QValue").ConvertUsing(x => StringUtilities.DblToString(x.QValue, 7));
                Map(x => x.PepQValue).Index(ColumnCount++).Name("PepQValue").ConvertUsing(x => StringUtilities.DblToString(x.PepQValue, 7));
            }
        }

        /// <summary>
        /// Read results from TSV file into group of objects from PSI_Interface
        /// </summary>
        /// <param name="idFilePath"></param>
        /// <returns></returns>
        public static SimpleMZIdentMLReader.SimpleMZIdentMLData ReadResultsFromFileToMzIdData(string idFilePath)
        {
            var databaseSearchResultData = ReadResultsFromFile(idFilePath);

            var reader = new SimpleMZIdentMLReader();
            var simpleMzIdentMLData = reader.Read(idFilePath);

            foreach (var databaseSearchResult in databaseSearchResultData)
            {
                var peptide = new SimpleMZIdentMLReader.PeptideRef
                {
                    Sequence = databaseSearchResult.Sequence
                };

                var identification = new SimpleMZIdentMLReader.SpectrumIdItem
                {
                    Peptide = peptide,
                    Charge = databaseSearchResult.Charge,
                    ScanNum = databaseSearchResult.ScanNum,
                    SpecEv = databaseSearchResult.SpecEValue,
                };

                // Parse the modification list, e.g.
                // Dehydro 9,Dehydro 20,Dehydro 27,Dehydro 48

                var modParts = databaseSearchResult.Modifications.Split(',');
                if (modParts.Length > 0)
                {
                    foreach (var part in modParts)
                    {
                        var mod = part.Split(' ');
                        if (mod.Length < 2)
                        {
                            continue;
                        }

                        var modName = mod[0];
                        var modIndex = Convert.ToInt32(mod[1]);
                        var ipMod = Modification.Get(modName);

                        var modification = new SimpleMZIdentMLReader.Modification
                        {
                            Mass = ipMod.Mass,
                            Name = modName,
                        };

                        peptide.ModsAdd(modIndex, modification);
                    }
                }

                var proteinAccessions = databaseSearchResult.ProteinName.Split(',');
                foreach (var accession in proteinAccessions)
                {
                    var dbSequence = new SimpleMZIdentMLReader.DatabaseSequence
                    {
                        Accession = accession.Trim(),
                        ProteinDescription = databaseSearchResult.ProteinDescription
                    };
                    var peptideEvidence = new SimpleMZIdentMLReader.PeptideEvidence
                    {
                        DbSeq = dbSequence,
                        Pre = databaseSearchResult.Pre,
                        Post = databaseSearchResult.Post,
                        PeptideRef = identification.Peptide,
                        Start = databaseSearchResult.Start,
                        End = databaseSearchResult.End,
                    };

                    identification.PepEvidence.Add(peptideEvidence);
                }

                simpleMzIdentMLData.Identifications.Add(identification);
            }

            return simpleMzIdentMLData;
        }

        /// <summary>
        /// Show the scan number and the sequence
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return string.Format("Scan {0,4}, EValue {1:0.00000}, Probability {2:0.00000}: {3}", ScanNum, EValue, Probability, SequenceWithEnds);
        }
    }
}
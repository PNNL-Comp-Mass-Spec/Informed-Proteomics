using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;
using PNNLOmics.Utilities;
using PSI_Interface.IdentData;

namespace InformedProteomics.Backend.Results
{
    using System.Linq;

    /// <summary>
    /// Container with utilities for reading/writing/storing database search results.
    /// </summary>
    public class DatabaseSearchResultData
    {
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
            get { return qValue; }
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
            get { return pepQValue; }
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

        /// <summary>
        /// Construct using a string from a tsv file input
        /// </summary>
        /// <param name="line"></param>
        public DatabaseSearchResultData(string line)
            : this()
        {
            ParseTsvLine(line);
        }

        /// <summary>
        /// Header string for default TSV output, without FDR scores
        /// </summary>
        public const string TsvHeaderString = "Scan\tPre\tSequence\tPost\tModifications\tComposition\tProteinName\tProteinDesc" +
                                              "\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tMass\tMs1Features\t#MatchedFragments\tProbability\tSpecEValue\tEValue";

        /// <summary>
        /// Format string for default TSV output, without FDR scores
        /// </summary>
        public const string TsvFormatString = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}";

        /// <summary>
        /// Header string for default TSV output, with FDR scores
        /// </summary>
        public const string TdaTsvHeaderString = TsvHeaderString + "\tQValue\tPepQValue";

        /// <summary>
        /// Format string for default TSV output, with FDR scores
        /// </summary>
        public const string TdaTsvFormatString = TsvFormatString + "\t{18}\t{19}";

        /// <summary>
        /// Get the header string for default TSV output, with columns added for FDR scores if addTdaScores is true
        /// </summary>
        /// <param name="addTdaScores"></param>
        /// <returns></returns>
        public static string GetHeaderString(bool addTdaScores = false)
        {
            return addTdaScores ? TdaTsvHeaderString : TsvHeaderString;
        }

        /// <summary>
        /// Create a TSV format string from this object
        /// </summary>
        /// <param name="addTdaScores">True to add FDR results (QValue and PepQValue) to the output</param>
        /// <returns></returns>
        public string TsvFormattedString(bool addTdaScores = false)
        {
            var format = TsvFormatString;
            if (addTdaScores)
            {
                format = TdaTsvFormatString;
            }
            return String.Format(format,
                ScanNum,
                Pre,                  // Pre
                Sequence,             // Sequence
                Post,                 // Post
                Modifications,        // Modifications
                Composition,          // Composition
                ProteinName,          // ProteinName
                ProteinDescription,   // ProteinDescription
                ProteinLength,        // ProteinLength
                Start,                // Start position in protein
                End,                  // End position in protein
                Charge,               // precursorCharge
                StringUtilities.DblToString(MostAbundantIsotopeMz, 9, true), // MostAbundantIsotopeMz
                StringUtilities.DblToString(Mass, 9, true),                  // Mass
                Ms1Feature,
                NumMatchedFragments,                                         // (Number of matched fragments)
                StringUtilities.DblToString(Probability, 4),                 // Probability
                StringUtilities.DblToString(Math.Max(SmallestValueExcel, SpecEValue), 6, true, 0.001), // SpecEValue; will be displayed using scientific notation if the value is less than 0.001
                StringUtilities.DblToString(Math.Max(SmallestValueExcel, EValue), 6, true, 0.001),     // EValue; will be displayed using scientific notation if the value is less than 0.001
                StringUtilities.DblToString(QValue, 7),
                StringUtilities.DblToString(PepQValue, 7)
            );
        }

        private const double SmallestValueExcel = 9.99E-308;

        private static string lastSetHeaderLine = "";
        private static int scanNumIndex = 0;
        private static int preIndex = 1;
        private static int sequenceIndex = 2;
        private static int postIndex = 3;
        private static int modificationsIndex = 4;
        private static int compositionIndex = 5;
        private static int proteinNameIndex = 6;
        private static int proteinDescriptionIndex = 7;
        private static int proteinLengthIndex = 8;
        private static int startIndex = 9;
        private static int endIndex = 10;
        private static int chargeIndex = 11;
        private static int mostAbundantIsotopeMzIndex = 12;
        private static int massIndex = 13;
        private static int ms1FeaturesIndex = 14;
        private static int numMatchedFragmentsIndex = 15;
        private static int probabilityIndex = 16;
        private static int specEValueIndex = 17;
        private static int eValueIndex = 18;
        private static int qValueIndex = 19;
        private static int pepQValueIndex = 20;

        /// <summary>
        /// Set the file header string for a file that will be read in
        /// </summary>
        /// <param name="headerLine"></param>
        public static void SetInputFileHeader(string headerLine)
        {
            if (string.IsNullOrWhiteSpace(headerLine))
            {
                return;
            }
            if (!string.IsNullOrWhiteSpace(lastSetHeaderLine) && headerLine.EndsWith(lastSetHeaderLine))
            {
                return;
            }
            var tokens = headerLine.Split('\t');
            for (var i = 0; i < tokens.Length; i++)
            {
                switch (tokens[i])
                {
                    case "Scan":
                        scanNumIndex = i;
                        break;
                    case "Pre":
                        preIndex = i;
                        break;
                    case "Sequence":
                        sequenceIndex = i;
                        break;
                    case "Post":
                        postIndex = i;
                        break;
                    case "Modifications":
                        modificationsIndex = i;
                        break;
                    case "Composition":
                        compositionIndex = i;
                        break;
                    case "ProteinName":
                        proteinNameIndex = i;
                        break;
                    case "ProteinDesc":
                        proteinDescriptionIndex = i;
                        break;
                    case "ProteinLength":
                        proteinLengthIndex = i;
                        break;
                    case "Start":
                        startIndex = i;
                        break;
                    case "End":
                        endIndex = i;
                        break;
                    case "Charge":
                        chargeIndex = i;
                        break;
                    case "MostAbundantIsotopeMz":
                        mostAbundantIsotopeMzIndex = i;
                        break;
                    case "Mass":
                        massIndex = i;
                        break;
                    case "Ms1Features":
                        ms1FeaturesIndex = i;
                        break;
                    case "#MatchedFragments":
                        numMatchedFragmentsIndex = i;
                        break;
                    case "Probability":
                        probabilityIndex = i;
                        break;
                    case "SpecEValue":
                        specEValueIndex = i;
                        break;
                    case "EValue":
                        eValueIndex = i;
                        break;
                    case "QValue":
                        qValueIndex = i;
                        break;
                    case "PepQValue":
                        pepQValueIndex = i;
                        break;
                }
            }
        }

        /// <summary>
        /// Parse a string from a tsv file input
        /// </summary>
        /// <param name="line"></param>
        /// <param name="headerLine"></param>
        public void ParseTsvLine(string line, string headerLine = null)
        {
            SetInputFileHeader(headerLine);
            var tokens = line.Split('\t');

            if (tokens.Length > scanNumIndex)
            {
                ScanNum = int.Parse(tokens[scanNumIndex]);
            }
            if (tokens.Length > preIndex)
            {
                Pre = tokens[preIndex];
            }
            if (tokens.Length > sequenceIndex)
            {
                Sequence = tokens[sequenceIndex];
            }
            if (tokens.Length > postIndex)
            {
                Post = tokens[postIndex];
            }
            if (tokens.Length > modificationsIndex)
            {
                Modifications = tokens[modificationsIndex];
            }
            if (tokens.Length > compositionIndex)
            {
                Composition = tokens[compositionIndex];
            }
            if (tokens.Length > proteinNameIndex)
            {
                ProteinName = tokens[proteinNameIndex];
            }
            if (tokens.Length > proteinDescriptionIndex)
            {
                ProteinDescription = tokens[proteinDescriptionIndex];
            }
            if (tokens.Length > proteinLengthIndex)
            {
                ProteinLength = int.Parse(tokens[proteinLengthIndex]);
            }
            if (tokens.Length > startIndex)
            {
                Start = int.Parse(tokens[startIndex]);
            }
            if (tokens.Length > endIndex)
            {
                End = int.Parse(tokens[endIndex]);
            }
            if (tokens.Length > chargeIndex)
            {
                Charge = int.Parse(tokens[chargeIndex]);
            }
            if (tokens.Length > mostAbundantIsotopeMzIndex)
            {
                MostAbundantIsotopeMz = double.Parse(tokens[mostAbundantIsotopeMzIndex]);
            }
            if (tokens.Length > massIndex)
            {
                Mass = double.Parse(tokens[massIndex]);
            }
            if (tokens.Length > ms1FeaturesIndex)
            {
                Ms1Feature = Convert.ToInt32(tokens[ms1FeaturesIndex]);
            }
            if (tokens.Length > numMatchedFragmentsIndex)
            {
                NumMatchedFragments = int.Parse(tokens[numMatchedFragmentsIndex]);
            }
            if (tokens.Length > probabilityIndex)
            {
                Probability = double.Parse(tokens[probabilityIndex]);
            }
            if (tokens.Length > specEValueIndex)
            {
                SpecEValue = double.Parse(tokens[specEValueIndex]);
            }
            if (tokens.Length > eValueIndex)
            {
                EValue = double.Parse(tokens[eValueIndex]);
            }
            if (tokens.Length > qValueIndex)
            {
                QValue = double.Parse(tokens[qValueIndex]);
            }
            if (tokens.Length > pepQValueIndex)
            {
                PepQValue = double.Parse(tokens[pepQValueIndex]);
            }
        }

        /// <summary>
        /// Write the resultData in TSV format to the specified path, possibly including FDR scores
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="resultData"></param>
        /// <param name="includeTdaScores">If FDR scores should be output also</param>
        public static void WriteResultsToFile(string filePath, IEnumerable<DatabaseSearchResultData> resultData, bool includeTdaScores = false)
        {
            using (var stream = new StreamWriter(new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
            {
                stream.WriteLine(GetHeaderString(includeTdaScores));
                foreach (var result in resultData)
                {
                    stream.WriteLine(result.TsvFormattedString(includeTdaScores));
                }
            }
        }

        /// <summary>
        /// Read in result data from the specified TSV format file.
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static List<DatabaseSearchResultData> ReadResultsFromFile(string filePath)
        {
            var results = new List<DatabaseSearchResultData>();
            using (var stream = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                string line;
                var isFirstLine = true;
                while (!stream.EndOfStream && (line = stream.ReadLine()) != null)
                {
                    if (isFirstLine)
                    {
                        isFirstLine = false;
                        int test;
                        if (!int.TryParse(line.Substring(0, 1), out test))
                        {
                            SetInputFileHeader(line);
                            continue;
                        }
                    }

                    results.Add(new DatabaseSearchResultData(line));
                }
            }
            if (results.Count == 0)
            {
                return null;
            }
            return results;
        }

        /// <summary>
        /// Read results from tsv file into group of objects from PSI_Interface
        /// </summary>
        /// <param name="idFilePath"></param>
        /// <returns></returns>
        public static SimpleMZIdentMLReader.SimpleMZIdentMLData ReadResultsFromFileToMzIdData(string idFilePath)
        {
            var databaseSearchResultData = ReadResultsFromFile(idFilePath);

            var simpleMzIdentMLData = new SimpleMZIdentMLReader.SimpleMZIdentMLData(idFilePath);

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

                // Parse modification
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
                            Tag = modName,
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
    }
}
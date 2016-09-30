using System;
using System.Collections.Generic;
using System.IO;
using PNNLOmics.Utilities;

namespace InformedProteomics.Backend.Results
{
    public class DatabaseSearchResultData
    {
        public int ScanNum { get; set; }
        public string Pre { get; set; }
        public string Sequence { get; set; }
        public string Post { get; set; }
        public string Modifications { get; set; }
        public string Composition { get; set; }
        public string ProteinName { get; set; }
        public string ProteinDescription { get; set; }
        public int ProteinLength { get; set; }
        public int Start { get; set; }
        public int End { get; set; }
        public int Charge { get; set; }
        public double MostAbundantIsotopeMz { get; set; }
        public double Mass { get; set; }
        public int NumMatchedFragments { get; set; }
        public double Probability { get; set; }
        public double SpecEValue { get; set; }
        public double EValue { get; set; }

        public double QValue
        {
            get { return qValue; }
            set
            {
                qValue = value;
                HasTdaScores = true;
            }
        }

        public double PepQValue
        {
            get { return pepQValue; }
            set
            {
                pepQValue = value;
                HasTdaScores = true;
            }
        }
        public bool HasTdaScores { get; private set; }

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

        private string sequenceWithEnds;
        private double qValue;
        private double pepQValue;

        public DatabaseSearchResultData()
        {
            Pre = "-";
            Post = "-";
            SpecEValue = double.MaxValue;
            EValue = double.MaxValue;
            qValue = double.MaxValue;
            pepQValue = double.MaxValue;
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

        public const string TsvHeaderString = "Scan\tPre\tSequence\tPost\tModifications\tComposition\tProteinName\tProteinDesc" +
                                              "\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tMass\t#MatchedFragments\tProbability\tSpecEValue\tEValue";

        public const string TsvFormatString = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}";
        public const string TdaTsvHeaderString = TsvHeaderString + "\tQValue\tPepQValue";
        public const string TdaTsvFormatString = TsvFormatString + "\t{18}\t{19}";

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
        private static int numMatchedFragmentsIndex = 14;
        private static int probabilityIndex = 15;
        private static int specEValueIndex = 16;
        private static int eValueIndex = 17;
        private static int qValueIndex = 18;
        private static int pepQValueIndex = 19;

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
            if (headerLine.EndsWith(lastSetHeaderLine))
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

        public static List<DatabaseSearchResultData> ReadResultsFromFile(string filePath)
        {
            var results = new List<DatabaseSearchResultData>();
            using (var stream = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                var header = stream.ReadLine();
                string line;
                if (string.IsNullOrWhiteSpace(header))
                {
                    return null;
                }
                var headerTokens = header.Split('\t');
                int temp;
                if (int.TryParse(headerTokens[0], out temp))
                {
                    header = "";
                    line = header;
                }
                else
                {
                    SetInputFileHeader(header);
                    line = stream.ReadLine();
                }

                while (line != null)// && !stream.EndOfStream)
                {
                    results.Add(new DatabaseSearchResultData(line));
                    line = stream.ReadLine();
                }
            }
            if (results.Count == 0)
            {
                return null;
            }
            return results;
        }
    }
}
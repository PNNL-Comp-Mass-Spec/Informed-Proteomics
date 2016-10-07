using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.DIA.Search
{
    public class MsGfPostProcessor
    {
        const string DecoyPrefix = FastaDatabase.DecoyProteinPrefix;
        const double SpecEValueThreshold = 1E-8;

        //// Post-processing parameters
        //const double RelativeIsotopeIntensityThreshold = 0.7;
        private const int NumIsotopesToCheckForValidation = 3;

        public MsGfPostProcessor(string specFilePath, string msGfResultPath, Tolerance tolForBaseXic, Tolerance tolFromBaseXic)
            : this(new[] { specFilePath }, msGfResultPath, tolForBaseXic, tolFromBaseXic)
        {
        }

        public MsGfPostProcessor(IEnumerable<string> specFilePaths, string msGfResultPath, Tolerance tolForBaseXic, Tolerance tolFromBaseXic)
        {
            MsGfResultPath = msGfResultPath;
            ToleranceForBaseXic = tolForBaseXic;
            ToleranceFromBasicXic = tolFromBaseXic;

            Run = new Dictionary<string, LcMsRun>();
            foreach (var specFilePath in specFilePaths)
            {
                var run = InMemoryLcMsRun.GetLcMsRun(specFilePath);
                var specFileKey = Path.GetFileNameWithoutExtension(specFilePath);
                if(specFileKey != null) Run[specFileKey] = run;
            }
        }

        public Dictionary<string, LcMsRun> Run { get; private set; }
        public string MsGfResultPath { get; private set; }
        public Tolerance ToleranceForBaseXic { get; private set; }
        public Tolerance ToleranceFromBasicXic { get; private set; }

        public int PostProcessing(string outputFilePath)
        {
            // Parse MS-GF+ results
            var pepToResults = new Dictionary<string, MsGfMatch>();

            MsGfPlusHeaderInformation headerInfo = null;

            foreach (var line in File.ReadLines(MsGfResultPath))
            {
                if (line.StartsWith("#"))
                {
                    headerInfo = new MsGfPlusHeaderInformation(line);
                    continue;
                }

                var match = new MsGfMatch(line, headerInfo);
                if (!match.IsValid) continue;

                if (match.SpecEValue > SpecEValueThreshold) continue;

                if (!IsValid(match)) continue;

                MsGfMatch prevMatch;
                if (!pepToResults.TryGetValue(match.Peptide, out prevMatch))
                {
                    pepToResults[match.Peptide] = match;
                    match.NumMatches = 1;
                }
                else
                {
                    if (match.SpecEValue < prevMatch.SpecEValue)
                    {
                        pepToResults[match.Peptide] = match;
                        match.NumMatches += prevMatch.NumMatches;
                    }
                    else
                    {
                        ++prevMatch.NumMatches;
                    }
                }
            }

            //var filteredPsms = pepToResults.Select(entry => entry.Value).Where(IsValid).ToList();
            var filteredPsms = pepToResults.Select(entry => entry.Value).ToList();
            filteredPsms.Sort();

            // compute FDR
            var qValue = GetQValues(filteredPsms);

            var index = -1;
            var numId = 0;
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("#SpecFile\tPeptide\tScanNum\tPrecursorMz\tCharge\tProtein\tNumMatches\tDeNovoScore\tMSGFScore\tSpecEValue\tPepQValue");
                foreach (var match in filteredPsms)
                {
                    //if (match.Protein.StartsWith("DecoyPrefix")) continue;
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                        match.SpecFile
                        //, match.Peptide.Replace("C+57.021", "C")
                        , match.Peptide
                        , match.ScanNum
                        , new Ion(match.Formula, match.Charge).GetMonoIsotopicMz()
                        , match.Charge
                        , match.Protein
                        , match.NumMatches
                        , match.DeNovoScore
                        , match.MsgfScore
                        , match.SpecEValue
                        , qValue[++index]
                        );
                    if (!match.Protein.StartsWith(DecoyPrefix))
                    {
                        if (qValue[index] <= 0.01) ++numId;
                    }
                }
            }

            return numId;
        }

        private bool IsValid(MsGfMatch match)
        {
            var specFileKey = Path.GetFileNameWithoutExtension(match.SpecFile);
            if (specFileKey == null) return false;

            var precursorIon = new Ion(match.Formula, match.Charge);
            var prevMs1ScanNum = Run[specFileKey].GetPrecursorScanNum(match.ScanNum);
            var nextMs1ScanNum = Run[specFileKey].GetNextScanNum(prevMs1ScanNum, 1);

            //var isotopeMzs =
            //    precursorIon.GetIsotopes(NumIsotopesToCheckForValidation)
            //        .Select(isotope => precursorIon.GetIsotopeMz(isotope.Index)).ToArray();
            var isotopes = precursorIon.GetIsotopes(NumIsotopesToCheckForValidation).ToArray();
            Array.Sort(isotopes);   // sort by indices

            var prevMs1Spec = Run[specFileKey].GetSpectrum(prevMs1ScanNum);
            if (prevMs1Spec != null)
            {
                //if (prevMs1Spec.ContainsIon(precursorIon, ToleranceForBaseXic, 0.9)) return true;
                var isPrevIsotopeValid = false;
                foreach (var isotope in isotopes)
                {
                    var mz = precursorIon.GetIsotopeMz(isotope.Index);
                    if (prevMs1Spec.FindPeak(mz, ToleranceForBaseXic) != null)  // match
                    {
                        if (isPrevIsotopeValid) return true;
                        isPrevIsotopeValid = true;
                    }
                }
            }

            var nextMs1Spec = Run[specFileKey].GetSpectrum(nextMs1ScanNum);
            if (nextMs1Spec != null)
            {
                var isPrevIsotopeValid = false;
                foreach (var isotope in isotopes)
                {
                    var mz = precursorIon.GetIsotopeMz(isotope.Index);
                    if (nextMs1Spec.FindPeak(mz, ToleranceForBaseXic) != null)  // match
                    {
                        if (isPrevIsotopeValid) return true;
                        isPrevIsotopeValid = true;
                    }
                }
            }

            return false;
        }

        //private bool IsValidOld(MsGfMatch match)
        //{
        //    // Valid if spectral E-value is below 1E-15
        //    //if (match.SpecEValue < 1e-15) return true;

        //    var specFileKey = Path.GetFileNameWithoutExtension(match.SpecFile);
        //    if (specFileKey == null) return false;

        //    var precursorIon = new Ion(match.Formula, match.Charge);
        //    var basePeakIndex = precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex();
        //    var basePeakMz = precursorIon.GetIsotopeMz(basePeakIndex);
        //    var baseXic = Run[specFileKey].GetExtractedIonChromatogram(precursorIon.GetMostAbundantIsotopeMz(), ToleranceForBaseXic, match.ScanNum);

        //    // check whether the most abundant isotope peak exists
        //    var prevMs1ScanNum = Run[specFileKey].GetPrecursorScanNum(match.ScanNum);
        //    var nextMs1ScanNum = Run[specFileKey].GetNextScanNum(prevMs1ScanNum, 1);
        //    if (!baseXic.ContainsScanNum(prevMs1ScanNum) && !baseXic.ContainsScanNum(nextMs1ScanNum))
        //    {
        //        return false;
        //        //if (match.SpecEValue < 1e-15)
        //        //{
        //        //    Console.WriteLine("{0} {1} {2} ({3}) is rejected (no base peak)", match.Peptide, match.Charge, match.ScanNum, match.SpecEValue);
        //        //}
        //    }

        //    // Check whether all abundant isotope peaks exist
        //    //foreach (var isotope in precursorIon.GetIsotopes(RelativeIsotopeIntensityThreshold))
        //    foreach (var isotope in precursorIon.GetIsotopes(NumIsotopesToCheckForValidation))    // top 2 isotopes
        //    {
        //        Xic xic;
        //        if (isotope.Index == basePeakIndex)
        //        {
        //            xic = baseXic;
        //        }
        //        else
        //        {
        //            //if(match.Peptide.Equals("HGGEDYVFSLLTGYCEPPTGVSLR"))
        //            //    Console.WriteLine("Debug");
        //            var isotopeIndex = isotope.Index;
        //            var mzDifference = precursorIon.GetIsotopeMz(isotopeIndex) - basePeakMz;
        //            xic = Run[specFileKey].GetIsotopeExtractedIonChromatogram(baseXic, mzDifference, ToleranceFromBasicXic);
        //        }

        //        if (xic.Count == 0)
        //        {
        //            if (match.SpecEValue < 1e-15)
        //            {
        //                //Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", match.SpecFile, match.Peptide, match.Charge,
        //                //    match.ScanNum, new Ion(_aaSet.GetComposition(match.Peptide) + Composition.H2O, match.Charge).GetMonoIsotopicMz(), match.SpecEValue);
        //            }
        //            return false;
        //        }
        //    }

        //    //if (!hasValid && match.SpecEValue < 1e-15)
        //    //{
        //    //    Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", match.SpecFileKey, match.Peptide, match.Charge,
        //    //        match.ScanNum, new Ion(_aaSet.GetComposition(match.Peptide) + Composition.H2O, match.Charge).GetMonoIsotopicMz(), match.SpecEValue);
        //    //};

        //    //var monoIsotopeMz = precursorIon.GetMonoIsotopicMz();
        //    //var monoIsotopeXic = Run[specFileKey].GetExtractedIonChromatogram(monoIsotopeMz, Tolerance, match.ScanNum);

        //    //// check isotopes at index -1
        //    //var isotopeMzMinus1 = precursorIon.GetIsotopeMz(-1);
        //    //var xicAtIndexMinus1 = Run[specFileKey].GetExtractedIonChromatogram(isotopeMzMinus1, new Tolerance(5), match.ScanNum);
        //    //if (xicAtIndexMinus1.Count > 0 && xicAtIndexMinus1.GetCorrelation(monoIsotopeXic) > 0.7)
        //    //{
        //    //    // Check the ratio
        //    //    var abundanceMono = monoIsotopeXic.GetSumIntensities();
        //    //    var abundanceMinusOne = xicAtIndexMinus1.GetSumIntensities();
        //    //    var isoEnv = _aaSet.GetComposition(match.Peptide).GetIsotopomerEnvelopeFromNominalMass();
        //    //    var approxRato = isoEnv[1]/isoEnv[0];
        //    //    var ratioDiff = abundanceMono/abundanceMinusOne/(isoEnv[1]/isoEnv[0]);
        //    //    if (ratioDiff > 0.8 && ratioDiff < 1.2)
        //    //    {
        //    //        if (match.SpecEValue < 1e-15)
        //    //        {
        //    //            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", match.SpecFile, match.Peptide, match.Charge,
        //    //                match.ScanNum, new Ion(_aaSet.GetComposition(match.Peptide) + Composition.H2O, match.Charge).GetMonoIsotopicMz(), match.SpecEValue, ratioDiff);
        //    //        }
        //    //        return false;
        //    //    }
        //    //}

        //    //if (xic.Count > 0)
        //    //{
        //    //    var isotopeRatio = xic.GetSumIntensities() / baseIntensity / isotope.Item2;
        //    //    var correlation = xic.GetCorrelation(baseXic);

        //    //    if (isotopeRatio > 0.8 && isotopeRatio < 1.2
        //    //        && correlation > 0.8)
        //    //    {
        //    //        isValid = true;
        //    //    }
        //    //}

        //    // Check if isotope ratio is within tolerance
        //    //if (isotopeRatio > isotopeRatioTolerance || isotopeRatio < 1 / isotopeRatioTolerance)
        //    //{
        //    //    isValid = false;
        //    //    //Console.WriteLine("Off ratio\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", isDecoy, peptide, scanNum, charge, precursorIon.GetMonoIsotopicMz(), isotopeMz, isotopeRatio);
        //    //    break;
        //    //}

        //    // check isotopes at index 0.5
        //    //if (match.Charge == 4)
        //    //{
        //    //    var isotopeMzPointFive = precursorIon.GetIsotopeMz(0.5);
        //    //    var xicAtIndexPointFive = Run[specFileKey].GetExtractedIonChromatogram(isotopeMzPointFive, Tolerance, match.ScanNum);
        //    //    if (xicAtIndexPointFive.Count > 0 && xicAtIndexPointFive.GetCorrelation(monoIsotopeXic) > 0.7) return false;
        //    //}

        //    return true;
        //}

        // matchList must have been sorted
        private double[] GetQValues(IList<MsGfMatch> matchList)
        {
            var numPsms = matchList.Count;
            var numDecoy = new int[numPsms];
            var numTarget = new int[numPsms];
            var fdr = new double[numPsms];
            var qValue = new double[numPsms];

            for (var i = 0; i < matchList.Count; i++)
            {
                var match = matchList[i];
                var isDecoy = match.Protein.StartsWith(DecoyPrefix) ? 1 : 0;
                if (i == 0) numDecoy[i] = isDecoy;
                else numDecoy[i] = numDecoy[i - 1] + isDecoy;
                numTarget[i] = i + 1 - numDecoy[i];
                fdr[i] = numDecoy[i] / (double)numTarget[i];
            }

            qValue[numPsms - 1] = fdr[numPsms - 1];
            for (var i = numPsms - 2; i >= 0; i--)
            {
                qValue[i] = Math.Min(qValue[i + 1], fdr[i]);
            }
            return qValue;
        }
    }

    public class MsGfResults
    {
        public MsGfResults(string resultFilePath)
        {
            _msGfMatches = new List<MsGfMatch>();
            MsGfPlusHeaderInformation headerInfo = null;
            var prevScanNum = -1;
            foreach (var line in File.ReadLines(resultFilePath))
            {
                if (headerInfo == null && line.StartsWith("#"))
                {
                    headerInfo = new MsGfPlusHeaderInformation(line);
                    continue;
                }

                var match = new MsGfMatch(line, headerInfo);

                if (match.ScanNum == prevScanNum) continue;
                prevScanNum = match.ScanNum;

                if (!match.IsValid || match.Protein.StartsWith(FastaDatabase.DecoyProteinPrefix)) continue;

                _msGfMatches.Add(match);
            }
            _msGfMatches.Sort();
        }

        public List<MsGfMatch> MsGfMatches
        {
            get { return _msGfMatches; }
        }

        public List<MsGfMatch> GetMatchesAtPsmFdr(double psmFdrThreshold)
        {
            return _msGfMatches.TakeWhile(m => !(m.QValue > psmFdrThreshold)).ToList();
        }

        public List<MsGfMatch> GetMatchesAtPepFdr(double pepFdrThreshold)
        {
            var matches = new List<MsGfMatch>();
            var pepSet = new HashSet<string>();

            foreach (var m in _msGfMatches)
            {
                if (m.PepQValue > pepFdrThreshold) continue;
                if (pepSet.Contains(m.Peptide)) continue;
                pepSet.Add(m.Peptide);
                matches.Add(m);
            }
            return matches;
        }

        private readonly List<MsGfMatch> _msGfMatches;
    }

    public class MsGfMatch : IComparable<MsGfMatch>
    {
        public MsGfMatch(string specFile, string peptide, int scanNum, int charge, string protein, int deNovoScore, int msgfScore, double specEValue)
        {
            SpecFile = specFile;
            Peptide = peptide;
            ScanNum = scanNum;
            Charge = charge;
            Protein = protein;
            DeNovoScore = deNovoScore;
            MsgfScore = msgfScore;
            SpecEValue = specEValue;
        }

        public MsGfMatch(string resultStr, MsGfPlusHeaderInformation header)
        {
            var token = resultStr.Split('\t');
            if (token.Length != header.NumColumns)
            {
                IsValid = false;
            }
            else
            {
                SpecFile = token[header.SpecFileColNum];
                Peptide = token[header.PeptideColNum];
                if (header.FormulaColNum > 0)
                {
                    Formula = Composition.Parse(token[header.FormulaColNum]);
                }
                ScanNum = Convert.ToInt32(token[header.ScanNumColNum]);
                Charge = Convert.ToInt32(token[header.ChargeColNum]);
                Protein = token[header.ProteinColNum];
                DeNovoScore = Convert.ToInt32(token[header.DeNovoScoreColNum]);
                MsgfScore = Convert.ToInt32(token[header.MsgfScoreColNum]);
                SpecEValue = Convert.ToDouble(token[header.SpecEValueColNum]);
                if (header.QValueColNum > 0)
                {
                    QValue = Convert.ToDouble(token[header.QValueColNum]);
                }
                if (header.PepQValueColNum > 0)
                {
                    PepQValue = Convert.ToDouble(token[header.PepQValueColNum]);
                }
                IsValid = true;
            }
        }

        public bool IsValid { get; private set; }
        public string SpecFile { get; private set; }
        public string Peptide { get; private set; }
        public Composition Formula { get; private set; }
        public int ScanNum { get; private set; }
        public int Charge { get; private set; }
        public string Protein { get; private set; }
        public int DeNovoScore { get; private set; }
        public int MsgfScore { get; private set; }
        public double SpecEValue { get; private set; }
        public double QValue { get; private set; }
        public double PepQValue { get; private set; }
        public int NumMatches { get; set; }

        public int CompareTo(MsGfMatch other)
        {
            if (SpecEValue > other.SpecEValue) return 1;
            if (SpecEValue < other.SpecEValue) return -1;
            return 0;
        }
    }

    public class MsGfPlusHeaderInformation
    {
        public MsGfPlusHeaderInformation(string header)
        {
            var token = header.Split('\t');

            NumColumns = token.Length;

            SpecFileColNum = -1;
            PrecursorColNum = -1;
            PeptideColNum = -1;
            FormulaColNum = -1;
            ScanNumColNum = -1;
            ChargeColNum = -1;
            ProteinColNum = -1;
            DeNovoScoreColNum = -1;
            MsgfScoreColNum = -1;
            SpecEValueColNum = -1;
            QValueColNum = -1;
            PepQValueColNum = -1;

            for (var i = 0; i < token.Length; i++)
            {
                if (token[i].Equals("#SpecFile")) SpecFileColNum = i;
                else if (token[i].StartsWith("Precursor")) PrecursorColNum = i;
                else if (token[i].Equals("Peptide")) PeptideColNum = i;
                else if (token[i].Equals("Formula")) FormulaColNum = i;
                else if (token[i].Equals("ScanNum")) ScanNumColNum = i;
                else if (token[i].Equals("Charge")) ChargeColNum = i;
                else if (token[i].Equals("Protein")) ProteinColNum = i;
                else if (token[i].Equals("DeNovoScore")) DeNovoScoreColNum = i;
                else if (token[i].Equals("MSGFScore")) MsgfScoreColNum = i;
                else if (token[i].Equals("SpecEValue")) SpecEValueColNum = i;
                else if (token[i].Equals("QValue")) QValueColNum = i;
                else if (token[i].Equals("PepQValue")) PepQValueColNum = i;
            }
        }

        public int SpecFileColNum { get; private set; }
        public int ScanNumColNum { get; private set; }
        public int PrecursorColNum { get; private set; }
        public int PeptideColNum { get; private set; }
        public int FormulaColNum { get; private set; }
        public int ChargeColNum { get; private set; }
        public int ProteinColNum { get; private set; }
        public int DeNovoScoreColNum { get; private set; }
        public int MsgfScoreColNum { get; private set; }
        public int SpecEValueColNum { get; private set; }
        public int QValueColNum { get; private set; }
        public int PepQValueColNum { get; private set; }

        public int NumColumns { get; private set; }
    }
}

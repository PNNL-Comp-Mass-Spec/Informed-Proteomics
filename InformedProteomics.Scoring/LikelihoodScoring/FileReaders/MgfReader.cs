using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.FileReaders
{
    public class MgfReader : IDataFileReader
    {
        public MgfReader(string fileName, bool decoy)
        {
            _fileName = fileName;
            _decoy = decoy;
        }
        private enum MgfState { Label, Parameter, Peak };
        public IList<SpectrumMatch> Read()
        {
            var specMatches = new List<SpectrumMatch>();
            var file = File.ReadLines(_fileName);
            var mgfState = MgfState.Label;

            var sequence = string.Empty;
            int scanNum = 0, charge = 0;
            var peaks = new List<Peak>();

            var peptideSet = new HashSet<string>();
            var lineNumber = 0;

            foreach (var line in file)
            {
                lineNumber++;
                switch (mgfState)
                {
                    case MgfState.Label:
                        if (line == "BEGIN IONS")
                        {
                            mgfState = MgfState.Parameter;
                        }
                        else
                        {
                            throw new FormatException("Invalid MGF file, expected BEGIN IONS at line " + lineNumber);
                        }

                        break;
                    case MgfState.Parameter:
                        var parameter = line.Split('=');
                        if (parameter.Length < 2)
                        {
                            throw new FormatException(string.Format("Line {0} is invalid in the MGF file: {1}", lineNumber, line));
                        }

                        if (parameter[0] == "SEQ")
                        {
                            sequence = parameter[1];
                        }
                        else if (parameter[0] == "SCANS")
                        {
                            scanNum = Convert.ToInt32(parameter[1]);
                        }
                        else if (parameter[0] == "CHARGE")
                        {
                            var chargeStr = parameter[1].Substring(0, parameter[1].Length - 1);
                            charge = Convert.ToInt32(chargeStr);
                            mgfState = MgfState.Peak;
                            if (string.IsNullOrWhiteSpace(sequence) || scanNum == 0 || charge == 0)
                            {
                                throw new FormatException("Incomplete spectrum entry, line: " + lineNumber);
                            }
                        }
                        break;
                    case MgfState.Peak:
                        if (line == "END IONS")
                        {
                            if (peaks.Count == 0)
                            {
                                throw new FormatException("Empty peak list, line: " + lineNumber);
                            }

                            mgfState = MgfState.Label;
                            if (peptideSet.Contains(sequence))
                            {
                                sequence = string.Empty;
                                scanNum = 0;
                                charge = 0;
                                peaks.Clear();
                                continue;
                            }
                            peptideSet.Add(sequence);
                            var spectrum = new ProductSpectrum(peaks, scanNum) { MsLevel = 2 };
                            var sequenceReader = new MgfSequenceReader();
                            var seq = sequenceReader.GetSequence(sequence);
                            var specMatch = new SpectrumMatch(seq, spectrum, scanNum, charge, _decoy);
                            sequence = string.Empty;
                            scanNum = 0;
                            charge = 0;
                            specMatches.Add(specMatch);
                            peaks.Clear();
                        }
                        else
                        {
                            var parts = line.Split('\t');
                            if (parts.Length < 2)
                            {
                                throw new FormatException(string.Format("Line {0} is invalid in the MGF file: {1}", lineNumber, line));
                            }

                            var mz = Convert.ToDouble(parts[0]);
                            var intensity = Convert.ToDouble(parts[1]);
                            peaks.Add(new Peak(mz, intensity));
                        }
                        break;
                }
            }
            return specMatches;
        }

        private readonly string _fileName;
        private readonly bool _decoy;
    }
}

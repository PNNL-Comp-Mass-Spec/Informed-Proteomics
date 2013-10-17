using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.DIA.Search
{
    public class PeptideCentricAnalysis
    {
        public PeptideCentricAnalysis(
            LcMsRun run,
            IEnumerable<string> peptideEnumerator,
            AminoAcidSet aminoAcidSet): this(run, peptideEnumerator, aminoAcidSet, 1, 3)
        {
            
        }

        public PeptideCentricAnalysis(
            LcMsRun run, 
            IEnumerable<string> peptideEnumerator, 
            AminoAcidSet aminoAcidSet, 
            int minCharge, 
            int maxCharge)
        {
            Run = run;
            PeptideEnumerator = peptideEnumerator;
            AminoAcidSet = aminoAcidSet;
            MinCharge = minCharge;
            MaxCharge = maxCharge;
        }

        public LcMsRun Run { get; private set; }
        public IEnumerable<string> PeptideEnumerator { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }

        public void IntermediateSearch(string outputFilePath)
        {
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("Annotation\tCharge\tScanNum");
                foreach (var annotation in PeptideEnumerator)
                {
                    // annotation: pre + "." + peptide + "." + post (e.g. R.PEPTIDER.G)
                    var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, annotation);
                    foreach (var sequenceComposition in seqGraph.GetSequenceCompositions())
                    {
                        for (var precursorCharge = MinCharge; precursorCharge <= MaxCharge; precursorCharge++)
                        {
                            var precursorIon = new Ion(sequenceComposition + Composition.H2O, precursorCharge);
                            var matchingSpectra = Run.GetMatchingMs2Spectra(precursorIon);
                            if (matchingSpectra == null) continue;
                            foreach (var productSpectrum in Run.GetMatchingMs2Spectra(precursorIon))
                            {
                                writer.WriteLine("{0}\t{1}\t{2}", annotation, precursorCharge, productSpectrum.ScanNum);
                            }
                        }
                    }
                }
            }
        }
    }
}

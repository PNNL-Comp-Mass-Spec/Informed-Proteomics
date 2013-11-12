using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestDia
    {
        [Test]
        public void ComputeSpikedInPeptideMzHist()
        {
            const string pepListFile = @"C:\cygwin\home\kims336\Data\DIA\SpikedPeptides.txt";

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var charges = new[] {2};

            var hist = new int[4];

            var sum = 0;

            Console.WriteLine("Peptide\tCharge\tMz");
            foreach (var line in File.ReadLines(pepListFile))
            {
                if (line.Length == 0) continue;
                var peptide = line;
                var composition = aaSet.GetComposition(peptide) + Composition.H2O;

                foreach (var charge in charges)
                {
                    var precursorIon = new Ion(composition, charge);
                    var precursorIonMz = precursorIon.GetMz();

                    if (precursorIonMz < 400 || precursorIonMz >= 900) continue;
                    var histIndex = (int)((precursorIonMz - 400)/125);
                    hist[histIndex]++;

                    Console.WriteLine("{0}\t{1}\t{2}\t{3}", peptide, charge, precursorIonMz, histIndex);

                    sum++;
                }
            }

            Console.WriteLine("\nRange\tNum\tRatio");
            for (var i = 0; i < hist.Length; i++)
            {
                Console.WriteLine("{0}-{1}\t{2}\t{3}", 400+i*125, 525+i*125, hist[i], hist[i] / (float)sum);
            }

        }
    }
}

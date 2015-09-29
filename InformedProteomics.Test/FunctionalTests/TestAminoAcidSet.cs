﻿using System;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestAminoAcidSet
    {
        [Test]
        public void TestParsingManyMods()
        {
            const string modFilePath = @"\\protoapps\UserData\Jungkap\Lewy\db\Mods.txt";
            var aaSet = new AminoAcidSet(modFilePath);
            //aaSet.Display();


            //SequenceLocation.ProteinNTerm
            var residue = AminoAcid.ProteinNTerm.Residue;
            var location = SequenceLocation.ProteinNTerm;
            var aa = aaSet.GetAminoAcid(residue, location);
            Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
            foreach (var modIndex in aaSet.GetModificationIndices(residue, location))
            {
                var modification = aaSet.GetModificationParams().GetModification(modIndex);
                Console.WriteLine(modification.Mass);
                //Console.Write("\t" + _modificationParams.GetModification(modIndex));
            }
            Console.WriteLine();
            residue = AminoAcid.ProteinCTerm.Residue;
            location = SequenceLocation.ProteinCTerm;
            aa = aaSet.GetAminoAcid(residue, location);
            Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
            foreach (var modIndex in aaSet.GetModificationIndices(residue, location))
            {
                var modification = aaSet.GetModificationParams().GetModification(modIndex);
                Console.WriteLine(modification.Mass);
                //Console.Write("\t" + _modificationParams.GetModification(modIndex));
            }


            //foreach (var aa in AminoAcid.StandardAminoAcidArr)
            
                /*
                var keys = _locationSpecificResidueMap[location].Keys.ToArray();
                Array.Sort(keys);
                foreach (var residue in keys)
                {
                    var aa = GetAminoAcid(residue, location);
                    Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
                    foreach (var modIndex in GetModificationIndices(residue, location))
                    {
                        Console.Write("\t" + _modificationParams.GetModification(modIndex));
                    }
                    Console.WriteLine();
                }
            }     */            


        }

        [Ignore]
        public void TestParsingGlycoMods()
        {
            const string modFilePath = @"C:\cygwin\home\kims336\Data\Debug\MSPathFinder_Mods.txt";
            var aaSet = new AminoAcidSet(modFilePath);
            aaSet.Display();
        }
    }
}

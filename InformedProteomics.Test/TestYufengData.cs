using System;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestYufengData
    {
        public const string TestRawFilePath = @"H:\Research\Yufeng\TopDownYufeng\raw\yufeng_column_test2.raw";
        
        [Test]
        public void Test43KProtein()
        {
            //var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath);
            const string sequence = "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            const string annotation = "_." + sequence + "._";
            var aaSet = new AminoAcidSet();
            var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
            Console.WriteLine(seqGraph.GetUnmodifiedSequenceComposition().Mass);
        }
    }
}

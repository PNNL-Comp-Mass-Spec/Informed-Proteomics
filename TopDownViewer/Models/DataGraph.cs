using System;
using System.Linq;
using GraphX;
using GraphX.Logic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using QuickGraph;

namespace InformedProteomics.TopDownViewer.Models
{
    //Graph data class
    public class DataGraph : BidirectionalGraph<DataVertex, DataEdge> { }

    //Logic core class
    public class LogicCore : GXLogicCore<DataVertex, DataEdge, BidirectionalGraph<DataVertex, DataEdge>> { }

    //Vertex data object
    public class DataVertex : VertexBase
    {
        public DataVertex()
        {
            Text = "";
            Selected = false;
        }
        public string Text { get; set; }

        public bool Selected { get; set; }

        public Composition PrefixComposition { get; set; }
        public Composition SuffixComposition { get; set; }
        public ModificationCombination ModificationCombination { get; set; }

        public override string ToString()
        {
            return Text;
        }
    }

    //Edge data object
    public class DataEdge : EdgeBase<DataVertex>
    {
        public DataEdge(DataVertex source, DataVertex target, double weight = 1)
            : base(source, target, weight)
        {

            Modifications = ModificationCombination.NoModification;
        }

        public DataEdge()
            : base(null, null, 1)
        {
            Modifications = ModificationCombination.NoModification;
            SequenceIndex = 0;
        }

        public AminoAcid AminoAcid { get; set; }

        public ModificationCombination Modifications
        {
            get { return _modifications; }
            set
            {
                _modifications = value;
                if (_modifications.Modifications.Count > 0)
                    AminoAcid = new ModifiedAminoAcid(AminoAcid, _modifications.Modifications.LastOrDefault());
            }
        }

        public int SequenceIndex { get; set; }

        public string Text
        {
            get
            {
                var text =
                    (Modifications.GetNumModifications() == 0)
                        ? String.Format("{0}", AminoAcid.Residue)
                        : String.Format("{0}[{1}]", AminoAcid.Residue, Modifications);
                return text;
            }
        }

        public override string ToString()
        {
            return Text;
        }

        private ModificationCombination _modifications;
    }
}

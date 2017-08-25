using System.Collections.Generic;

namespace InformedProteomics.Backend.MathAndStats
{
    public class VennDiagram<T>
    {
        public VennDiagram(ISet<T> set1, ISet<T> set2)
        {
            Set1 = set1;
            Set2 = set2;
            ComputeVennDiagram();
        }

        public ISet<T> Set1 { get; }
        public ISet<T> Set2 { get; }

        public ISet<T> Intersection { get; private set; }
        public ISet<T> Set1Only { get; private set; }
        public ISet<T> Set2Only { get; private set; }

        private void ComputeVennDiagram()
        {
            Intersection = new HashSet<T>();
            Set1Only = new HashSet<T>();
            Set2Only = new HashSet<T>();

            foreach (var element in Set1)
            {
                if (Set2.Contains(element)) Intersection.Add(element);
                else Set1Only.Add(element);
            }

            foreach (var element in Set2)
            {
                if (!Intersection.Contains(element)) Set2Only.Add(element);
            }
        }
    }
}

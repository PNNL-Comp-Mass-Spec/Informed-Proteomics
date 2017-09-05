using System.Collections.Generic;

namespace InformedProteomics.Backend.MathAndStats
{
    /// <summary>
    /// Venn Diagram calculation
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class VennDiagram<T>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="set1"></param>
        /// <param name="set2"></param>
        public VennDiagram(ISet<T> set1, ISet<T> set2)
        {
            Set1 = set1;
            Set2 = set2;
            ComputeVennDiagram();
        }

        /// <summary>
        /// Set 1 data
        /// </summary>
        public ISet<T> Set1 { get; }

        /// <summary>
        /// Set 2 data
        /// </summary>
        public ISet<T> Set2 { get; }

        /// <summary>
        /// Intersection of Set 1 and set 2
        /// </summary>
        public ISet<T> Intersection { get; private set; }

        /// <summary>
        /// Set 1 - set 2
        /// </summary>
        public ISet<T> Set1Only { get; private set; }

        /// <summary>
        /// Set 2 - set 1
        /// </summary>
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

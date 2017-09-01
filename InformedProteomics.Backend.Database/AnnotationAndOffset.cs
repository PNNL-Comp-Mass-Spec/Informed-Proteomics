namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Annotation and offset data
    /// </summary>
    public class AnnotationAndOffset
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="offset"></param>
        /// <param name="sequence"></param>
        public AnnotationAndOffset(long offset, string sequence)
        {
            Offset = offset;
            Annotation = sequence;
        }

        /// <summary>
        /// Sequence offset
        /// </summary>
        public long Offset { get; }

        /// <summary>
        /// Sequence annotation
        /// </summary>
        public string Annotation { get; }
    }
}

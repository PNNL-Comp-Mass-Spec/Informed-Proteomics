using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// This class cataluges all possible modification instances
    /// </summary>
    public class ModificationParams
    {
        private ModificationInstance[] _modificationInstances;
        private int _numMaxDynMods;
        private Dictionary<int, Modification>[] _map;

        /// <summary>
        /// Storing all possible combinations of modifications up to numMaxDynMods
        /// </summary>
        /// <param name="modifications"></param>
        /// <param name="numMaxDynMods"></param>
        public ModificationParams(Modification[] modifications, int numMaxDynMods)
        {
            int numModifications = modifications.Count();
            int numModificationInstances = SimpleMath.NChooseK(numModifications + numMaxDynMods, numMaxDynMods);
            _modificationInstances = new ModificationInstance[numModificationInstances];
            GenerateModificationInstances(modifications.ToArray(), numMaxDynMods);
        }

        /// <summary>
        /// Gets the modification instance of the specified index
        /// </summary>
        /// <param name="index">modification instance index</param>
        /// <returns>the modification instance</returns>
        public ModificationInstance GetModificationInstance(int index)
        {
            return _modificationInstances[index];
        }

        /// <summary>
        /// Apply modification to prevModIns and returns the modification instance.
        /// Modification instances are referred by their indices.
        /// </summary>
        /// <param name="prevModIndex">the index of the modification instance to apply modification</param>
        /// <param name="modification">modification to apply</param>
        /// <returns>the index of the modification instance of prevModIns + modification</returns>
        public int GetModificationInstanceIndex(int prevModIndex, Modification modification)
        {
            throw new System.NotImplementedException();
        }

        /// <summary>
        /// Apply modification to prevModIns and returns the modification instance
        /// </summary>
        /// <param name="prevModIns">modification instance to apply modification</param>
        /// <param name="modification">modification to apply</param>
        /// <returns>the modification instance of prevModIns + modification</returns>
        public ModificationInstance GetModificationInstance(ModificationInstance prevModIns, Modification modification)
        {
            throw new System.NotImplementedException();
        }

        /// <summary>
        /// Gets the number of all possible modification instances
        /// </summary>
        /// <returns>the number of modification instances</returns>
        public int GetNumModificationInstances()
        {
            return _modificationInstances.Length;
        }

        private void GenerateModificationInstances(Modification[] modifications, int numMaxDynMods)
        {
            var modInsList = new List<ModificationInstance>[numMaxDynMods+1];
            modInsList[0] = new List<ModificationInstance>
                {
                    new ModificationInstance(new[] {Modification.NoModification})
                };

            _map = new Dictionary<int, Modification>[_modificationInstances.Length];

            // No modification

            for (int numMods = 1; numMods <= numMaxDynMods; numMods++)
            {
                //modInsList.AddRange(GetModificationInstances(modInsList));
            }
        }

        private IEnumerable<ModificationInstance> GetModificationInstances(IEnumerable<ModificationInstance> modificationInstances)
        {
            return null;
        }


    }
}

namespace InformedProteomics.Backend.Data.Spectrometry
{
	public class Fragment
	{
		/// <summary>
		/// The charge state of the fragment.
		/// </summary>
		public int ChargeState { get; set; }

		/// <summary>
		/// The m/z of the fragment.
		/// </summary>
		public double Mz { get; set; }

		/// <summary>
		/// The mass of the fragment.
		/// </summary>
		public double Mass { get; set; }

		/// <summary>
		/// If this fragment is y6, the FragmentIonClassBase is y
		/// </summary>
		public string IonType { get; set; }

		/// <summary>
		/// The index of the residue this fragment breaks the peptide at.
		/// </summary>
		public int ResidueNumber { get; set; }

		/// <summary>
		/// The full ion symbol of the fragment. e.g. y6 or y6++ or y6-H2O++
		/// </summary>
		public string IonSymbol { get; set; }

		public bool Equals(Fragment other)
		{
			if (ReferenceEquals(null, other)) return false;
			if (ReferenceEquals(this, other)) return true;
			return Equals(other.IonSymbol, IonSymbol);
		}

		public override bool Equals(object obj)
		{
			if (ReferenceEquals(null, obj)) return false;
			if (ReferenceEquals(this, obj)) return true;
			if (obj.GetType() != typeof(Fragment)) return false;
			return Equals((Fragment)obj);
		}

		public override int GetHashCode()
		{
			return IonSymbol.GetHashCode();
		}

		public static bool operator ==(Fragment left, Fragment right)
		{
			return Equals(left, right);
		}

		public static bool operator !=(Fragment left, Fragment right)
		{
			return !Equals(left, right);
		}

		public override string ToString()
		{
			return string.Format("ChargeState: {0}, Mz: {1}, AveragineMass: {2}, ResidueNumber: {3}, FragmentIonClassBase: {4}, IonSymbol: {5}", ChargeState, Mz, Mass, ResidueNumber, IonType, IonSymbol);
		}
	}
}

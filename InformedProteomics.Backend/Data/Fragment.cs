namespace InformedProteomics.Backend.Data
{
	public class Fragment
	{
		public int ChargeState { get; set; }
		public double Mz { get; set; }
		public double Mass { get; set; }
		public string IonType { get; set; }
		public int ResidueNumber { get; set; }
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
	}
}

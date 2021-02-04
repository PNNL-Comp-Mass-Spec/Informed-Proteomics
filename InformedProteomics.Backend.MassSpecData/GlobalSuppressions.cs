// This file is used by Code Analysis to maintain SuppressMessage
// attributes that are applied to this project.
// Project-level suppressions either have no target or are given
// a specific target and scoped to a namespace, type, member, etc.

using System.Diagnostics.CodeAnalysis;

[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Ignore errors here", Scope = "member", Target = "~M:InformedProteomics.Backend.MassSpecData.ProteoWizardReaderImplementation.FindPwizPath~System.String")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.Backend.MassSpecData.PbfLcMsRun.GetPrecursorChromatogramRange(System.Double,System.Double)~InformedProteomics.Backend.Data.Spectrometry.Xic")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.Backend.MassSpecData.PbfLcMsRun.GetPrecursorExtractedIonChromatogram(System.Double,System.Double)~InformedProteomics.Backend.Data.Spectrometry.Xic")]
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "Leave as-is for clarity", Scope = "member", Target = "~M:InformedProteomics.Backend.MassSpecData.ProteoWizardReaderImplementation.ReadSpectrum(System.Int32,System.Boolean)~InformedProteomics.Backend.Data.Spectrometry.Spectrum")]
[assembly: SuppressMessage("Usage", "RCS1236:Use exception filter.", Justification = "Leave as-is", Scope = "member", Target = "~M:InformedProteomics.Backend.MassSpecData.PbfLcMsRun.ConvertToPbf(System.String,InformedProteomics.Backend.MassSpecData.IMassSpecDataReader,System.Double,System.Double,System.String,System.IProgress{PRISM.ProgressData})~System.String")]

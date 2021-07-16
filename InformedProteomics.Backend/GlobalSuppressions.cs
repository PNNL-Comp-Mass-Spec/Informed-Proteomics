// This file is used by Code Analysis to maintain SuppressMessage
// attributes that are applied to this project.
// Project-level suppressions either have no target or are given
// a specific target and scoped to a namespace, type, member, etc.

using System.Diagnostics.CodeAnalysis;

[assembly: SuppressMessage("Readability", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Backend.Data.Spectrometry.MzComparerWithBinning.#ctor(System.Int32)")]
[assembly: SuppressMessage("Readability", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Backend.Data.Spectrometry.MzComparerWithBinning.GetTolerance~InformedProteomics.Backend.Data.Spectrometry.Tolerance")]
[assembly: SuppressMessage("Readability", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~P:InformedProteomics.Backend.Data.Spectrometry.MzComparerWithBinning.Ppm")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.Backend.Data.Spectrometry.AbstractFragmentScorer.FindMatchedPeaks(InformedProteomics.Backend.Data.Composition.Composition,System.Double,System.Double)~System.Collections.Generic.IEnumerable{InformedProteomics.Backend.Data.Spectrometry.DeconvolutedPeak}")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.Backend.Data.Spectrometry.AbstractFragmentScorer.FindMostIntensePeak(InformedProteomics.Backend.Data.Composition.Composition,System.Double,System.Double,System.Int32@,System.Double@,System.Double@)~InformedProteomics.Backend.Data.Spectrometry.Peak[]")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.Backend.Data.Spectrometry.AbstractFragmentScorer.GetMinMaxChargeRange(System.Double)~InformedProteomics.Backend.Utils.IntRange")]

/* Unmerged change from project 'InformedProteomics.Backend (netstandard2.0)'
Added:
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "<Pending>", Scope = "member", Target = "~M:InformedProteomics.Backend.FeatureFindingResults.Ms1FtEntry.WriteToFile(System.String,System.Collections.Generic.IEnumerable{InformedProteomics.Backend.FeatureFindingResults.Ms1FtEntry},System.Boolean)")]
*/
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "<Pending>", Scope = "member", Target = "~M:InformedProteomics.Backend.SearchResults.DatabaseSearchResultData.WriteResultsToFile(System.String,System.Collections.Generic.IEnumerable{InformedProteomics.Backend.SearchResults.DatabaseSearchResultData},System.Boolean)")]
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "<Pending>", Scope = "member", Target = "~M:InformedProteomics.Backend.FeatureFindingResults.Ms1FtEntry.WriteToFile(System.String,System.Collections.Generic.IEnumerable{InformedProteomics.Backend.FeatureFindingResults.Ms1FtEntry},System.Boolean)")]

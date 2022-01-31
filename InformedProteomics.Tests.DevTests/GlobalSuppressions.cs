// This file is used by Code Analysis to maintain SuppressMessage
// attributes that are applied to this project.
// Project-level suppressions either have no target or are given
// a specific target and scoped to a namespace, type, member, etc.

using System.Diagnostics.CodeAnalysis;

[assembly: SuppressMessage("Performance", "RCS1096:Convert 'HasFlag' call to bitwise operation (or vice versa).", Justification = "Leave as-is for readability", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TestMisc.TestPathUtils")]
[assembly: SuppressMessage("Roslynator", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.Obsolete.TestEdrn.ComputeSpikedInPeptideMzHist")]
[assembly: SuppressMessage("Roslynator", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TestLcMsCaching.TestMs1Filtering")]
[assembly: SuppressMessage("Roslynator", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TestLcMsCaching.TestNominalMassErrors")]
[assembly: SuppressMessage("Roslynator", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TestMisc.CreatePeptideAbundanceTableWithSkyline")]
[assembly: SuppressMessage("Roslynator", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TopDownAnalysis.TestLikelihoodScorer.GetMatchStatistics(InformedProteomics.Backend.Data.Spectrometry.ProductSpectrum,InformedProteomics.Backend.Data.Sequence.Sequence,System.Int32,System.IO.StreamWriter)")]
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "Leave as-is for clarity", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TestLcMsCaching.FilteringEfficiencyQcShew")]
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "Leave as-is for clarity", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TestUtex.CopyUTEX")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.Tests.DevTests.TopDownAnalysis.TestLikelihoodScorer.GetMatchStatistics(InformedProteomics.Backend.Data.Spectrometry.ProductSpectrum,InformedProteomics.Backend.Data.Sequence.Sequence,System.Int32,System.IO.StreamWriter)")]

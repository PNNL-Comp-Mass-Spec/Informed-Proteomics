﻿// This file is used by Code Analysis to maintain SuppressMessage
// attributes that are applied to this project.
// Project-level suppressions either have no target or are given
// a specific target and scoped to a namespace, type, member, etc.

using System.Diagnostics.CodeAnalysis;

[assembly: SuppressMessage("Performance", "RCS1096:Convert 'HasFlag' call to bitwise operation (or vice versa).", Justification = "Leave as-is for readability", Scope = "member", Target = "~M:InformedProteomics.TopDown.Execution.IcTopDownLauncher.RunSearch(System.Double,System.Nullable{System.Threading.CancellationToken},System.IProgress{PRISM.ProgressData})~System.Boolean")]
[assembly: SuppressMessage("Performance", "RCS1096:Convert 'HasFlag' call to bitwise operation (or vice versa).", Justification = "Leave as-is for readability", Scope = "member", Target = "~M:InformedProteomics.TopDown.Execution.MzidResultsWriter.WriteResultsToMzid(System.Collections.Generic.IEnumerable{InformedProteomics.Backend.SearchResults.DatabaseSearchResultData},System.String)")]
[assembly: SuppressMessage("Simplification", "RCS1180:Inline lazy initialization.", Justification = "Leave as-is for readability", Scope = "member", Target = "~M:InformedProteomics.TopDown.Execution.IcTopDownLauncher.RunGeneratingFunction(System.Collections.Generic.IReadOnlyList{System.Collections.Generic.SortedSet{InformedProteomics.TopDown.Execution.DatabaseSequenceSpectrumMatch}},System.Nullable{System.Threading.CancellationToken},System.IProgress{PRISM.ProgressData})~System.Collections.Generic.Dictionary{System.Int32,System.Collections.Generic.List{InformedProteomics.TopDown.Execution.DatabaseSequenceSpectrumMatch}}")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.TopDown.Execution.IcTopDownLauncher.SearchForMatches(InformedProteomics.Backend.Database.AnnotationAndOffset,InformedProteomics.Backend.Data.Spectrometry.ISequenceFilter,System.Collections.Generic.SortedSet{InformedProteomics.TopDown.Execution.DatabaseSequenceSpectrumMatch}[],System.Int32,System.Boolean,System.Nullable{System.Threading.CancellationToken})")]
[assembly: SuppressMessage("Usage", "RCS1246:Use element access.", Justification = "Prefer to use .First()", Scope = "member", Target = "~M:InformedProteomics.TopDown.Scoring.MatchedPeakPostScorer.GetHyperGeometricScore~System.Double")]
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "Leave as-is for clarity", Scope = "member", Target = "~M:InformedProteomics.TopDown.Execution.MsAlignRescorer.Rescore(System.String,System.String)")]
[assembly: SuppressMessage("Usage", "RCS1146:Use conditional access.", Justification = "Leave as-is for clarity", Scope = "member", Target = "~M:InformedProteomics.TopDown.SequenceTag.IdentifiedSequenceTag.GenerateSequence(System.String,System.String)~InformedProteomics.Backend.Data.Sequence.Sequence")]
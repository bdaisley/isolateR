## Update: 2024-08-11 (isolateR 1.0.1)
-Updated <code>sanger_assembly</code> function (aka "sanger_consensus" prior to v1.0.1) to allow assembly of paired Sanger sequences
  - This function loads in the CSV results table from isoQC and merges related sequences based on user input. 
  - Original file names before isoQC step need to have a common prefix and differentiating suffixes. (e.g. SAMPLE_01_F.ab1, SAMPLE_01_R.ab1).
  - After aligning paired sequences, the consensus sequence is extracted and priority is given to the read with higher quality. Phred quality scores are reassigned in the final output table in a basic way by taking the mean of both input sequences.
  - Note: This function is designed to be used after the isoQC step and before the isoTAX step.
-Fixed warning issue during isoTAX step caused by  usage of <code>girafe</code> function from ggiraph package.


## Update: 2024-07-21 (isolateR 1.0.0)
- Added feature to <code>isoTAX</code> function to allow input of FASTA format files (relevant to addressing Issue #9). 
  - FASTA file type is automatically detected if first character of file is ">". 
  - File extension (.fasta/.fa/.fna) does not impact processing.
  - Can handle single- or multi-line FASTA format. 
  - If input is a FASTA file, the sequence(s) will be converted and saved as an isoQC-formatted output file in the current working directory ("isolateR_output/01_isoQC_mock_table.csv"). Sequence date, name, length, and number of ambiguous bases (Ns) will be calculated from the input file and used to populate the relevant columns. Phred quality scores (phred_trim) will be set to the maximum value (60) and the remaining columns will be populated with mock data to allow compatibility with the isoTAX function. The main purpose of this output file is for flexibility and to allow users to edit/modify the sequence metadata before continuing with subsequent steps.
- Fixed error with <code>isoLIB</code> function not producing an HTML output file when only one group exists due to all sequences in a given run being identical or very similar in terms of pairwise nucleotide identity.
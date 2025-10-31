## Update: 2025-10-31 (isolateR 1.0.3)
- Fixed issue related to #18 (Error: object 'TreeLine' is not exported by 'namespace:DECIPHER') with internal dependency on the deprecated <code>TreeLine</code> function in later versions of the DECIPHER package.
  - Changes made: <code>TreeLine</code> -> <code>Treeline</code> (Line 164 and Line 234) replacement in "isoLIB_function.R"
- Fixed issue related to #16 causing working directory to be reset and isoQC/isoTAX/isoLIB requiring absolute paths.
  - Changes made: Added <code>normalizePath()</code> functions to isoQC/isoTAX/isoLIB functions, allowing both absolute and relative paths to be speicifed as input.

## Update: 2025-01-20 (isolateR 1.0.2)
- Updated <code>isoALL</code> function to allow input of custom databases (e.g. isoALL(.., db_path=/path/to/databases, ...)
  - This was already available using the isoTAX function but was not implemented through isoALL.
  - The path should be a FASTA-formatted database sequence file. Ignored if 'db' parameter is set to anything other than NULL or "custom".
  - Expects a semicolon (;) delimited FASTA header as follows: Accession_no;d__Domain;p__Phylum;c__Class; o__Order;f__Family;g__Genus;s__Species. See <code>get_db</code> function documentation for examples and details on automatically generating custom databaes for offline use or within an HPC cluster environment.
- Fixed an issue causing the Phylum-level taxonomy to be shifted (e.g. the new Kingdom "Bacillati" incorrectly labelled as phylum in place of "Bacillota"). This issue arised from the recent subdivision of the domain Bacteria into the kingdoms Bacillati, Fusobacteriati, Pseudomonadati and Thermotogati. See [Göker and Oren (2024, International Journal of Systematic and Evolutionary Microbiology)](https://doi.org/10.1099/ijsem.0.006242) for more details (DOI:10.1099/ijsem.0.006242).
  - The isoTAX and isoALL functions have both been updated to return the correct phylum names.

## Update: 2024-08-11 (isolateR 1.0.1)
- Updated <code>sanger_assembly</code> function (aka "sanger_consensus" prior to v1.0.1) to allow assembly of paired Sanger sequences
  - This function loads in the CSV results table from isoQC and merges related sequences based on user input. 
  - Original file names before isoQC step need to have a common prefix and differentiating suffixes. (e.g. SAMPLE_01_F.ab1, SAMPLE_01_R.ab1).
  - After aligning paired sequences, the consensus sequence is extracted and priority is given to the read with higher quality. Phred quality scores are reassigned in the final output table in a basic way by taking the mean of both input sequences.
  - Note: This function is designed to be used after the isoQC step and before the isoTAX step.
- Fixed warning issue during isoTAX step caused by  usage of <code>girafe</code> function from ggiraph package.


## Update: 2024-07-21 (isolateR 1.0.0)
- Added feature to <code>isoTAX</code> function to allow input of FASTA format files (relevant to addressing Issue #9). 
  - FASTA file type is automatically detected if first character of file is ">". 
  - File extension (.fasta/.fa/.fna) does not impact processing.
  - Can handle single- or multi-line FASTA format. 
  - If input is a FASTA file, the sequence(s) will be converted and saved as an isoQC-formatted output file in the current working directory ("isolateR_output/01_isoQC_mock_table.csv"). Sequence date, name, length, and number of ambiguous bases (Ns) will be calculated from the input file and used to populate the relevant columns. Phred quality scores (phred_trim) will be set to the maximum value (60) and the remaining columns will be populated with mock data to allow compatibility with the isoTAX function. The main purpose of this output file is for flexibility and to allow users to edit/modify the sequence metadata before continuing with subsequent steps.
- Fixed error with <code>isoLIB</code> function not producing an HTML output file when only one group exists due to all sequences in a given run being identical or very similar in terms of pairwise nucleotide identity.
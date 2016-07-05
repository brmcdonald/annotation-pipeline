# annotation-pipeline
JSON-based customizable genome annotation pipeline for bacterial genomes

Work in progress. Features that need to be added:
  parsing genbak files as input
  generating genbank files from JSON files
  adding blast and last parsers
  
Known bugs:
  Labeling of the contig a feature is found on is inconsistent between RNA and CDS features.
  Coordinates for RNA features on the reverse strand are not recorded properly
  CDS features with antifam hits may not be properly removed

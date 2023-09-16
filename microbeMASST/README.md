# Search for microbial metabolites with microbeMASST

For each compound you will find these files. The HTML files are for a visual inspection. The tsv files more for parsing.

| File 	| Format	| File present?	| Description	|
| ----	| ------	| -------------	| -----------	|
| matches.tsv	| tsv	| Always generated if fasstMASST search succeeded. Empty if no matches	| fastMASST GNPS dataset matches, one match per file with highest cosine similarity	|
| library.tsv	| tsv	| Only if matches to datasets > 0	| fastMASST GNPS library matches	|
| datasets.tsv	| tsv	| Only if matches to datasets > 0	| Frequency of matches in datasets	|
| microbe.html	| html,JS,CSS	| Only if matches to microbe metadata >0	| Interactive sharable html version of tables and json-tree 	|
| microbe.json	| json	| Only if matches to microbe metadata >0	| json-tree for NCBI taxonomy and enriched with dataset matches	|
| counts_microbe.tsv	| tsv	| Only if matches to microbe metadata >0	| For each NCBI: number of matching files	|
| food.html	| html,JS,CSS	| Only if matches to food metadata >0	| Interactive sharable html version of tables and json-tree 	|
| food.json	| json	| Only if matches to food metadata >0	| json-tree for global foodomics (GFOP) taxonomy and enriched with dataset matches	|
| counts_food.tsv	| tsv	| Only if matches to food metadata >0	| For each GFOP id: number of matching files	|

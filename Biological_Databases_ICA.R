#Setting the working directory 
setwd("/home/student_number/ICA/") 
 
################################################################################ 
#                 Part 1: Obtaining the data from the databases                # 
################################################################################ 
#Installing and loading the necessary packages 
install.packages("tidyverse") 
install.packages("UniprotR") 
library(tidyverse) #<- used to manipulate and organise tables 
library(httr) #<- used to work with URLs 
library(jsonlite) #<- allows us to convert between R and JSON objects 
library(rentrez) #<- allows for access to NCBI databases 
library(biomaRt) #<- allows for access to Biomart/Ensembl directly from the R terminal 
library(UniprotR) #<- allows for access to UniProt directly from the R terminal 
library(RMySQL) #<- allows us to interact with the MySQL database 
 
################################################################################ 
#Database 1: ENSEMBL 
#The following code was adapted from the code made available in lecture BD4_AccessingDatabasesUsingWebServices_2per.pdf, pages 6-7 
#creating a string variable containing the base url for ENSEMBL 
server <- "https://rest.ensembl.org"  
#creating a string variable representing a gene by its specific ID 
ext <- "/lookup/id" 
#Using POST() to conduct the REST API request  
ensembl_search <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = '{ "ids" : ["ENSMUSG00000036061", "ENSMUSG00000000555", "ENSMUSG00000023055", "ENSMUSG00000075394", "ENSMUSG00000001655", "ENSMUSG00000022485", "ENSMUSG00000001657", "ENSMUSG00000001661", "ENSMUSG00000076010", "ENSMUSG00000023048"] }') 
#Checking the connection between the http and the R terminal.  
stop_for_status(ensembl_search)  
#Inspecting the data frame 
data.frame(t(sapply(content(ensembl_search),c))) 
res <-content(ensembl_search) 
#Convert each element of the res list into a data frame, and then bind each data frame into one large data frame "ensembl" 
ensembl <- do.call(rbind, lapply(res, as.data.frame)) 
#Establishing the order of the column headers 
header <- c("id","assembly_name","biotype", "canonical_transcript", "db_type", "description", "display_name", "logic_name", "object_type", "seq_region_name", "source", "species",  "start", "end", "strand", "version") 
#Ensuring the correct columns are assigned to the appropriate column headers 
ensembl <- ensembl[, header] 
#Renaming the column containing gene names (relevant in Part 2) 
names(ensembl)[names(ensembl) == 'display_name'] <- 'gene_name' 
#Reordering the columns (relevant in Part 2) 
ensembl <- select(ensembl, "gene_name", "id", "assembly_name", "biotype", "canonical_transcript", 
       "db_type", "description", "logic_name", "object_type", "seq_region_name",  
       "source", "species", "start", "end", "strand", "version") 
#Organsising the rows alphabetically based on gene name 
ensembl_data <- ensembl[order(ensembl$gene_name),] 
#Exporting the data in "ensembl_data" to a tab-delimited text file "ensembl_output.txt" 
Ensemble_table <- write.table(ensembl_data, file = "ensembl_output.txt", sep = "\t", row.names = FALSE) 
 
#for loop in R 
#1:length(res) generates a sequence of numbers from 1 to the length of object res. x takes on each value in this sequence one at a time during each iteration of the loop. 
for (x in 1:length(res)) 
{ 
  #print the element of res at the position of x, and continue to do so until the loop has printed out all the elements of rest one by one. 
  print(res[x]) 
} 
 
#Sending our data to MySQL. 
db <- dbConnect(RMySQL::MySQL(), 
                user='<student_number>',  
                password='<password>', 
                dbname='<student_number>') #<- connecting to the MySQL database. 
 
#Create the table in MySQL, specifying the column names. 
#Note: when checking the table in MySQL, you need to use the command 'SELECT * from ensembl_ICA' to properly display the columns. 
ensembl_query <- dbSendQuery(db,"CREATE TABLE ENSEMBL_data( 
                      gene_name VARCHAR(50), 
                      id VARCHAR(50), 
                      assembly_name VARCHAR(50), 
                      biotype VARCHAR(50), 
                      canonical_transcript VARCHAR(50), 
                      db_type VARCHAR(50), 
                      description TEXT, 
                      logic_name TEXT, 
                      object_type VARCHAR(50), 
                      seq_region_name VARCHAR(50), 
                      source VARCHAR(50), 
                      species VARCHAR(50), 
                      start INT, 
                      end INT, 
                      strand INT, 
                      version INT);") 
 
#exporting the data from "ensemble_output.txt" into the table ensembl_ICA 
ensembl_query2 <- dbSendQuery(db, "LOAD DATA LOCAL INFILE './ensembl_output.txt' INTO TABLE ENSEMBL_data IGNORE 1 LINES;") 
 
################################################################################ 
#Database 2: NCBI 'Gene' 
#The following code was adapted from the code made available in BD5_AccessingDatabases_UsingCustomAPIs_2per.pdf, pages 9-11. 
#Access the rentrez() manual: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#rate-limiting-and-api-keys 
#List the available NCBI databases 
entrez_dbs() 
 
#We are going to be using the 'Gene' database, and want to explore it. 
entrez_db_summary("gene") 
entrez_db_searchable("gene") 
 
#Querying the gene database for the given Ensembl IDs. We're searching in the WORD field, which contains the "free text associated with record", and we are searching for the Ensembl IDs, which are not the actual gene names. 
#We asked for this to be returned as a json file format, because this allows us to isolate the Gene IDs 
result <- entrez_search(db   = "gene", 
                        term = "(ENSMUSG00000036061[WORD])  
                       OR (ENSMUSG00000000555[WORD])  
                       OR (ENSMUSG00000023055[WORD])  
                       OR (ENSMUSG00000075394[WORD])  
                       OR (ENSMUSG00000001655[WORD])  
                       OR (ENSMUSG00000022485[WORD])  
                       OR (ENSMUSG00000001657[WORD]) 
                       OR (ENSMUSG00000001661[WORD]) 
                       OR (ENSMUSG00000076010[WORD]) 
                       OR (ENSMUSG00000023048[WORD])", 
                        retmode="json") 
 
#Obtaining the gene IDs from our previous query and putting them into a variable 
gene_ids <- result$file$esearchresult$idlist 
print(gene_ids) 
 
#Creating a data frame from the information we can get from the obtained Gene IDs 
gene_data <- XML::xmlToDataFrame(entrez_fetch(db="gene", gene_ids, rettype="xml", parsed=TRUE)) 
#Renaming the column containing gene names (relevant in Part 2) 
names(gene_data)[names(gene_data) == 'Entrezgene_gene'] <- 'gene_name' 
#Reorganising column names (relevant in Part 2) 
gene_data <- select(gene_data, "gene_name", "Entrezgene_track-info", "Entrezgene_type",             
                    "Entrezgene_source", "Entrezgene_prot", "Entrezgene_summary", 
                    "Entrezgene_location", "Entrezgene_gene-source",      
                    "Entrezgene_locus", "Entrezgene_properties",       
                    "Entrezgene_comments", "Entrezgene_unique-keys", "Entrezgene_xtra-index-terms") 
#Organising the rows in alphabetical order based on gene name (relevant for Part 2) 
ncbi_gene <- gene_data[order(gene_data$gene_name),] 
 
#Creating a table of our genes  
ncbi_table <- write.table(ncbi_gene, file = "entrez_gene_output.txt", sep = "\t") 
 
#Connecting to MySQL 
db <- dbConnect(RMySQL::MySQL(), 
                user='<student_number>',  
                password='<password>', 
                dbname='<student_number>') #<- connecting to the MySQL database. 
 
#Create the table in MySQL, specifying the column names. 
ncbi_query <- dbSendQuery(db,"CREATE TABLE NCBI_gene( 
    Row_num VARCHAR(50), 
    gene_name TEXT, 
    Entrezgene_track_info VARCHAR(50), 
    Entrezgene_type VARCHAR(50), 
    Entrezgene_source TEXT, 
    Entrezgene_prot TEXT, 
    Entrezgene_summary TEXT,  
    Entrezgene_location VARCHAR(50), 
    Entrezgene_gene_source VARCHAR(50), 
    Entrezgene_locus TEXT, 
    Entrezgene_properties TEXT, 
    Entrezgene_comments TEXT,  
    Entrezgene_unique_keys VARCHAR(100), 
    Entrezgene_xtra_index_terms VARCHAR(50));") 
 
#Exporting the table to MySQL 
ncbi_query2 <- dbSendQuery(db, "LOAD DATA LOCAL INFILE './entrez_gene_output.txt'  
                  INTO TABLE NCBI_gene IGNORE 1 LINES;") 
 
################################################################################ 
#Database 3: UniProt 
#The most reliable way to search for proteins in UniProt is by using primary accession numbers. 
#To obtain the relevant primary accession numbers from the given ENSEMBL IDs, we must access the Biomart database, which can be done using the necessary package. 
#How to use the biomaRt package: https://www.ensembl.org/info/data/biomart/biomart_r_package.html 
#Lines 168-182 were adapted from the code available on https://www.biostars.org/p/429062/ 
#Loading and viewing the available data sets to determine which one you need to use for this task. 
biomart <- useMart('ENSEMBL_MART_ENSEMBL') 
datasets <- listDatasets(biomart) 
view(datasets) 
 
#From looking at the previous two databases, we know that we need to use the mmusculus_gene_ensembl dataset 
biomart <- useDataset("mmusculus_gene_ensembl", biomart) 
 
#Make a new data frame to obtain the relevant ENSEMBL IDs and UniProt Primary Access Codes 
ensemble_to_uniprot_table <- getBM(mart = biomart, attributes = c( 
    "ensembl_gene_id", "external_gene_name", "uniprot_gn_id"), 
  filter = 'ensembl_gene_id', 
  values = c("ENSMUSG00000036061", "ENSMUSG00000000555", "ENSMUSG00000023055",  
             "ENSMUSG00000075394", "ENSMUSG00000001655", "ENSMUSG00000022485",  
             "ENSMUSG00000001657", "ENSMUSG00000001661", "ENSMUSG00000076010", 
             "ENSMUSG00000023048"), uniqueRows = TRUE) 
view(ensemble_to_uniprot_table) 
 
#When you view the above data frame, you will notice that it contains too many objects 
#This is because in some cases, the system is also classifying the entry name as a Primary Access Code 
#Create a new data frame only containing the rows with the Primary Accession Codes 
new_table <- ensemble_to_uniprot_table[c(1, 3, 4, 6, 7, 9, 12, 15, 18, 19),] 
view(new_table) 
#Create a variable only containing the Primary Accession Codes 
accession_numbers <- new_table$uniprot_gn_id 
print(accession_numbers) 
 
#Access the manual for the UniprotR package: https://cran.r-project.org/web/packages/UniprotR/UniprotR.pdf 
#In the UniProtR package, we cannot submit a blind request to obtain information regarding the relevant proteins 
#We need to use specific functions listed in the UniProtR package to obtain specific information as data frames 
names <- GetNamesTaxa(accession_numbers) 
general_info <- GetGeneral_Information(accession_numbers) 
sequences <- GetSequences(accession_numbers) 
length <- GetSeqLength(accession_numbers) 
structure <- GetStructureInfo(accession_numbers) 
functions <- GetProteinFunction(accession_numbers) 
subcellular_info <- GetSubcellular_location(accession_numbers) 
misc <- GetMiscellaneous(accession_numbers) 
 
#Create one data frame containing all the information from the previous 8  data frames 
joined = left_join(rownames_to_column(names), rownames_to_column(general_info), by = join_by(rowname)) 
joined = joined %>%  
  left_join(., rownames_to_column(sequences),by = join_by(rowname)) %>%  
  left_join(., rownames_to_column(length),by = join_by(rowname)) %>%  
  left_join(., rownames_to_column(structure),by = join_by(rowname)) %>%  
  left_join(., rownames_to_column(functions),by = join_by(rowname)) %>% 
  left_join(., rownames_to_column(subcellular_info),by = join_by(rowname)) %>%  
  left_join(., rownames_to_column(misc), by = join_by(rowname)) 
view(joined) 
 
#Cleaning the data frame 
alltogether <- joined[ , colSums(is.na(joined))==0] 
#Renaming the column containing gene names (relevant in Part 2) 
names(alltogether)[names(alltogether) == 'Gene.Names..primary.'] <- 'gene_name' 
names(alltogether)[names(alltogether) == 'rowname'] <- 'primary_accession_number' 
#Reorganising column names (relevant in Part 2) 
alltogether <- select(alltogether, "gene_name", "primary_accession_number", "Entry",  
                      "Entry.Name", "Gene.Names", "Organism", "Organism..ID.", 
                      "Protein.names", "Proteomes", "Taxonomic.lineage",                  
                      "Date.of.creation", "Date.of.last.modification",          
                      "Date.of.last.sequence.modification", "Entry.version",                      
                      "Alternative.products..isoforms.","Length.x",                           
                      "Mass", "Sequence", "Sequence.version", "Length.y",                           
                      "Function..CC.", "Subcellular.location..CC.", "Annotation", 
                      "Keywords", "Features", "Protein.existence", "Reviewed") 
#Organising the rows alphabetically based on gene name (relevant for Part 2) 
uniprot <- alltogether[order(alltogether$gene_name),] 
 
#Connecting to MySQL 
db <- dbConnect(RMySQL::MySQL(), 
                user='<student_number>',  
                password='<password>', 
                dbname='<student_number>') #<- connecting to the MySQL database. 
 
#Create the table in MySQL, specifying the column names. 
uniprot_query <- dbSendQuery(db,"CREATE TABLE UniProt( 
                      gene_name VARCHAR(50), 
                      primary_accession_number VARCHAR(50), 
                      entry VARCHAR(50), 
                      entry_names VARCHAR(50), 
                      gene_names VARCHAR(50), 
                      organism VARCHAR(50), 
                      organims_id VARCHAR(50), 
                      protein_names TEXT, 
                      proteomes VARCHAR(100), 
                      taxonomic_lineage TEXT, 
                      date_of_creation DATE, 
                      date_of_last_modification DATE, 
                      date_of_last_sequence_modification DATE, 
                      entry_version VARCHAR(50), 
                      alternative_products_isoforms VARCHAR(250), 
                      length INT,  
                      mass INT,  
                      SEQUENCE TEXT, 
                      sequence_version INT, 
                      length_y INT, 
                      protein_function TEXT, 
                      subcellular_location TEXT, 
                      annotation INT, 
                      keywords TEXT,  
                      features TEXT,  
                      protein_existence VARCHAR(250), 
                      reviewed VARCHAR(50) 
                      );") 
 
uniprot_table <- write.table(uniprot, file = "uniprot_output.txt", sep = "\t", row.names = FALSE) 
 
#exporting the data from "ensemble_output.txt" into the table uniprot_ICA 
uniprot_query2 <- dbSendQuery(db, "LOAD DATA LOCAL INFILE './uniprot_output.txt' INTO TABLE UniProt IGNORE 1 LINES;") 
 
################################################################################ 
#                        Part 2: Collating the 3 tables                        # 
################################################################################ 
#If you look at the ncbi_gene data frame, you will notice that the gene names are listed in full, not the abbreviated format.  
#This can prove a bit of an issue when using left_join.  
#Since we organised all the rows in alphabetical order according to gene_name, we know that the column contents of ensembl_data and ncbi_gene should correspond. 
#Knowing this, we will overwrite the gene_names column in ncbi_gene with the contents in the same column from ensembl_data. 
ncbi_gene$gene_name = ensembl_data$gene_name 
#Join the ensembl_data and ncbi_gene data frames. 
ensembl_and_ncbi =  left_join(ensembl_data, ncbi_gene, by = "gene_name") 
 
#Collate the previously made dataframe with the uniprot data frame. 
big3 = left_join(ensembl_and_ncbi, uniprot, by = "gene_name") 
 
#Clean the data to get rid of duplicate information 
kept_columns <- big3[c("gene_name", "id", "assembly_name", "biotype", "canonical_transcript", 
                 "description", "object_type", "seq_region_name", "source", 
                 "species", "start", "end", "strand", "version", "Entrezgene_track-info", 
                 "Entrezgene_type", "Entrezgene_prot", "Entrezgene_summary", 
                 "Entrezgene_location", "Entrezgene_gene-source", "Entrezgene_locus",                
                 "Entrezgene_properties", "Entrezgene_comments", "Entrezgene_unique-keys", 
                 "primary_accession_number", "Protein.names", "Entry.version", 
                 "Alternative.products..isoforms.", "Length.x", "Mass",                               
                 "Sequence", "Sequence.version", "Length.y", "Function..CC.",                     
                 "Subcellular.location..CC.", "Annotation", "Features", "Protein.existence")] 
 
#Create a txt. file of the cleaned table 
final_table <- write.table(kept_columns, file = "databases_ICA.txt", sep = "\t", row.names = FALSE) 
 
#Creating a final table in MySQL 
penultimate_query <- dbSendQuery(db,"CREATE TABLE FINAL_TABLE( 
                      gene_name VARCHAR(50), 
                      id VARCHAR(50), 
                      assembly_name VARCHAR(50), 
                      biotype VARCHAR(50), 
                      canonical_transcript VARCHAR(100), 
                      description TEXT, 
                      object_type VARCHAR(50), 
                      seq_region_name VARCHAR(50), 
                      source VARCHAR(50), 
                      species VARCHAR(50), 
                      start INT, 
                      end INT, 
                      strand INT, 
                      version INT, 
                      Entrezgene_track_info VARCHAR(50), 
                      Entrezgene_type VARCHAR(50), 
                      Entrezgene_prot TEXT, 
                      Entrezgene_summary TEXT, 
                      Entrezgene_location VARCHAR(50), 
                      Entrezgene_gene_source VARCHAR(50), 
                      Entrezgene_locus TEXT, 
                      Entrezgene_properties TEXT, 
                      Entrezgene_comments TEXT, 
                      Entrezgene_unique_keys VARCHAR(250), 
                      primary_accession_number VARCHAR(50), 
                      Protein_names TEXT, 
                      Entry_version INT, 
                      Alternative_products_isoforms TEXT, 
                      Length_x INT, 
                      Mass INT, 
                      Sequence TEXT, 
                      Sequence_version INT, 
                      Length_y INT, 
                      Protein_function TEXT, 
                      Subcellular_location TEXT, 
                      Annotation INT, 
                      Features TEXT, 
                      Protein_existence VARCHAR(50));") 
 
#Exporting the final table to MySQL 
final_query <- dbSendQuery(db, "LOAD DATA LOCAL INFILE './databases_ICA.txt' INTO TABLE FINAL_TABLE IGNORE 1 LINES;") 

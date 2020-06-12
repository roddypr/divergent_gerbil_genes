
mrna <- read.table("tmp/psammomys_obesus.mrna")
mrna <- as.character(mrna$V1)

mrna <- strsplit(mrna, split = ";")


# ID=Pob_R000007;genename=Med7;source_id=ENSMUSP00000020665-D1;

mrna_id   <- sapply(mrna, function(x) x[grep("ID=", x)])
mrna_id   <- gsub("ID=", "", mrna_id)

mrna_gene <- sapply(mrna, function(x) x[grep("genename=", x)])
mrna_gene <- gsub("genename=", "", mrna_gene)
mrna_gene[mrna_gene == "NULL"] <- ""

index_df <- data.frame(species = "psammomys_obesus",
                       gene = mrna_id,
                       tx = mrna_id,
                       pep = mrna_id,
                       gene_name = mrna_gene)

write.table(index_df, file="tmp/psammomys_obesus.index",
           sep = " ", col.names=FALSE, row.names=FALSE, quote=FALSE)

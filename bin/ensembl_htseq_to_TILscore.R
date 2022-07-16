# read a single HTseq count file which contains ensembl ids and convert them to Symbols
suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
library(dplyr)
library(readr)
library(singscore)
set.seed(42)

opts = list(
  make_option("--htseqcountsfile", action="store", default="", type="character", help="HTSeq counts file that needs converted to Symbols (default: NA)"),
  make_option("--symbollengthfile", action="store", default="", type="character", help="The Symbols lengths file (default: NA)"),
  make_option("--tilscorefile", action="store", default="", type="character", help="The TILscore219.csv file (default: NA)")
)
opts = parse_args(OptionParser(option_list = opts))

if (opts$htseqcountsfile == "" | opts$symbollengthfile == "" | opts$tilscorefile == "") {
  cat("\n")
  message("Error: See script usage (--help)")
  quit()
}

sample = sub("\\.counts\\.txt", "", basename(opts$htseqcountsfile))



counts <- read_tsv(opts$htseqcountsfile, col_names = c("Name", "counts"))
head(counts,n=2)
tail(counts)

#remove htseq specific features (starts with "__")
counts = counts[grep("^__", counts$Name, invert = TRUE),]
tail(counts)

#remove version numbers from ends of ensembl gene IDs to match the Symbol genelengths file
counts$Name <- sub('\\..*$', '', counts$Name)
grep("PAR", counts$Name, value = TRUE)
head(counts,n=2)


# Read Symbol lengths file (contains mapping from EnsemblID to Symbol)
symbol_lengths = read_tsv(opts$symbollengthfile)
head(symbol_lengths,n=2)
dim(symbol_lengths)

# remove duplicate rows  - this can happen because removing the EnsemblID suffix generates the same EnsemblID value
example_dupid = symbol_lengths[duplicated(symbol_lengths),][1,1] %>% pull()
symbol_lengths %>% filter(Geneid == example_dupid)

symbol_lengths_dedup = symbol_lengths[!duplicated(symbol_lengths),]
dim(symbol_lengths_dedup)
length(unique(symbol_lengths_dedup$external_gene_name))  # 59427 unique symbols


# merge
symbolcounts = counts %>%
  inner_join(symbol_lengths_dedup, by = c("Name"="Geneid")) %>%
  select(external_gene_name, counts)
head(symbolcounts,n=2)
dim(symbolcounts)


outfilename = paste0(sample, ".symbolcounts.tsv")
write_tsv(symbolcounts, outfilename)


# sum counts across duplicate Symbols.  Result is a table with unique Symbols n=59427 and raw counts
bound = symbolcounts %>% group_by(external_gene_name) %>% summarize(sum_counts = sum(counts))
head(bound,n=2)

outfilenameBound = paste0(sample, ".symbolboundcounts.tsv")
write_tsv(bound, outfilenameBound)



#
# Calculate TPM
#
# In order to calculate TPM on each Symbol, need to remove duplicate Symbols from Symbols_lengths_dedup table
dup_external_gene_name = !duplicated(symbol_lengths_dedup$external_gene_name)
length <- symbol_lengths_dedup[dup_external_gene_name, ]
length <- arrange(length,external_gene_name) %>% select(-Geneid) #length <- length[,-1]
length <- na.omit(length)
head(length,n=2)


# TPM function
tpm<- function(counts, gene_length){
  # Transcripts per Million. See also
  # https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
  # Comparing TPM across samples makes sense IF the total number of transcripts
  # stays the same. This condition may be far from real.
  # 
  stopifnot(length(counts) == length(gene_length))
  stopifnot(is.numeric(counts))
  stopifnot(gene_length > 0)
  rate<- counts / gene_length  # Normalize gene count by length
  denom<- sum(rate)            # Total transcription
  xtpm<- rate / denom * 1e6    # Proportion each gene contributes to tot transcription;
  return(xtpm)                 # times 1M to have it has "count of mRNAs if the total
  # n of transcripts in the cell was 1M."
}

bound_wlength = bound %>% inner_join(length)
head(bound_wlength,n=2)

bound_tpm = bound_wlength %>%
  mutate(tpm = tpm(counts = bound_wlength$sum_counts, gene_length = bound_wlength$Length)) %>%
  select(external_gene_name, tpm)
sum(bound_tpm$tpm)   # check that it is equal to 1M
head(bound_tpm, n=2)

outfilenameTPM = paste0(sample, ".symbolboundtpm.tsv")
write_tsv(bound_tpm, outfilenameTPM)



#
# Calculate TIL scores
#
tilscore = read_csv(opts$tilscorefile)
head(tilscore,n=2)
tilscore <- na.omit(unique(tilscore$TILscore))
sigs <- list(tilscore)
names(sigs) <- c("TILscore")

## Run singscore
#Options for "UP"set and "DOWN"set (we are only using UPset option, but for multiple direction signatures this may be useful)
#rank per sample column of expressional value
genes = data.frame(bound_tpm)
rownames(genes) = genes$external_gene_name
genes = genes[,-1, drop = FALSE]
colnames(genes) = sample
head(genes,n=2)

rankData <- rankGenes(genes)
head(rankData)
#run singscore for a single signature
ss_tilscore <- simpleScore(rankData = rankData, upSet = tilscore)

## Clean and calculator probability of 10% reactive TIL cultures with UM92G38 equation
ss_tilscore <- ss_tilscore %>% rename("TILscore" = "TotalScore")
ss_tilscore$TS_probability <- (1/(1+(exp(-((ss_tilscore$TILscore*16.637663)-4.514455)))))
ss_tilscore$TS_zscore <- (ss_tilscore$TILscore - 0.2359621)/0.08796414
#ss_tilscore$probability <- (exp((ss_tilscore$TILscore*21.40348)-5.67883))/(1+(exp((ss_tilscore$TILscore*21.40348)-5.67883
ss_tilscore$sample <- rownames(ss_tilscore)
ss_tilscore <- ss_tilscore %>% relocate(sample)
ss_tilscore <- ss_tilscore[,-3]   # remove TotalDispersion
ss_tilscore

outfilenameTIL = paste0(sample, ".tilscore.tsv")
write_tsv(ss_tilscore, outfilenameTIL)



sessionInfo()

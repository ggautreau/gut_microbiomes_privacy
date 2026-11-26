#!/usr/bin/env Rscript
suppressMessages(library(ExperimentHub))
OUT <- "local/crocodeel_valles/matrices"
load("local/sm.rda"); META <- sampleMetadata               # sampleMetadata of curatedMetagenomicData, saved with save(sampleMetadata, file = "local/sm.rda")
eh <- ExperimentHub(); cmd <- query(eh, "curatedMetagenomicData")
tit <- cmd$title; ids <- names(cmd)
studies <- sort(unique(META$study_name))
done <- sub("\\.species_relab\\.tsv$","", list.files(OUT, pattern="species_relab.tsv$"))
todo <- setdiff(studies, done)
cat("à extraire:", length(todo), "sur", length(studies), "\n")
for (s in todo) {
  sel <- grep(paste0("^[0-9-]+\\.", s, "\\.relative_abundance$"), tit)
  if (!length(sel)) { cat("ABSENT ", s, "\n"); next }
  sel <- sel[order(tit[sel], decreasing=TRUE)][1]
  m <- tryCatch(suppressMessages(cmd[[ ids[sel] ]]), error=function(e) NULL)
  if (is.null(m)) { cat("ECHEC  ", s, "\n"); next }
  sp <- grep("\\|s__", rownames(m)); m <- m[sp, , drop=FALSE]
  rownames(m) <- sub(".*\\|s__", "", rownames(m))
  write.table(data.frame(species=rownames(m), m, check.names=FALSE),
              file.path(OUT, paste0(s, ".species_relab.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  md <- META[META$study_name==s, ]
  keep <- intersect(c("sample_id","subject_id","body_site","family","age","gender",
                      if (any(!is.na(md$family))) "infant_age",                          # the 15 studies with family records carry infant_age too
                      "days_from_first_collection","number_reads","NCBI_accession","country",
                      "sequencing_platform","DNA_extraction_kit"), colnames(md))
  write.table(md[, keep, drop=FALSE], file.path(OUT, paste0(s, ".metadata.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("OK %-26s %4d esp. x %5d ech.\n", s, nrow(m), ncol(m)))
}
cat("EXTRACTION_FINIE\n")

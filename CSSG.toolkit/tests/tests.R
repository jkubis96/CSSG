library(testthat)

test_that("loading GTF", {
  GTF <- GTF.tool::load_annotation("test_anno.gtf")
  expect_true(all(c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9") %in% colnames(GTF)))
  expect_gt(nrow(GTF), 0)
})

GTF <- GTF.tool::load_annotation("test_anno.gtf")

test_that("create_GTF_df", {
  GTF2 <- GTF.tool::create_GTF_df(GTF, optimize = TRUE, shift = 100000)
  expect_true(all(c(
    "chr", "source", "annotationType", "start", "end", "score",
    "strand", "phase", "gene_id", "gene_name", "transcript_name",
    "transcript_id", "gene_type", "metadata"
  ) %in% colnames(GTF2)))
  expect_gt(nrow(GTF2), 0)
})

GTF2 <- GTF.tool::create_GTF_df(GTF, optimize = TRUE, shift = 100000)

test_that("add_CDS", {
  GTF2 <- GTF.tool::create_GTF_df(GTF, optimize = TRUE, shift = 100000)
  GTF2 <- GTF2[!GTF2$annotationType %in% c("CDS", "transcript", "mRNA"), ]
  GTF3 <- GTF.tool::add_CDS(GTF2)
  expect_true(all(c(
    "chr", "source", "annotationType", "start", "end", "score",
    "strand", "phase", "gene_id", "gene_name", "transcript_name",
    "transcript_id", "gene_type", "metadata"
  ) %in% colnames(GTF3)))
  expect_gt(nrow(GTF3), 0)
  expect_true(any(GTF3$annotationType == "CDS"))
  expect_true(any(GTF3$source == "JBIO-predicted"))
})

test_that("add_introns", {
  GTF4 <- GTF.tool::add_introns(GTF2)
  expect_true(all(c(
    "chr", "source", "annotationType", "start", "end", "score",
    "strand", "phase", "gene_id", "gene_name", "transcript_name",
    "transcript_id", "gene_type", "metadata"
  ) %in% colnames(GTF4)))
  expect_gt(nrow(GTF4), 0)
  expect_true(any(GTF4$annotationType == "intron"))
  expect_true(any(GTF4$source == "JBIO-predicted"))
})


test_that("add_UTR", {
  GTF5 <- GTF.tool::add_UTR(GTF2, five_prime_utr = 300, three_prime_utr = 800)
  expect_true(all(c(
    "chr", "source", "annotationType", "start", "end", "score",
    "strand", "phase", "gene_id", "gene_name", "transcript_name",
    "transcript_id", "gene_type", "metadata"
  ) %in% colnames(GTF5)))
  expect_gt(nrow(GTF5), 0)
  expect_true(any(GTF5$annotationType %in% c("five_prime_UTR", "three_prime_UTR")))
  expect_true(any(GTF5$source == "JBIO-predicted"))
})

test_that("create_full_GTF", {
  GTF6 <- GTF.tool::create_full_GTF(GTF2)
  expect_gt(nrow(GTF6), 0)
  expect_true(all(c(
    "chr", "source", "annotationType", "start", "end", "score",
    "strand", "phase", "combine"
  ) %in% colnames(GTF6)))
})


test_that("create_reduced_GTF", {
  GTF7 <- GTF.tool::create_reduced_GTF(GTF2)
  expect_gt(nrow(GTF7), 0)
  expect_true(all(c(
    "chr", "start", "end", "strand", "transcript_name",
    "transcript_id", "gene_name", "gene_id", "gene_type",
    "annotationType", "transcriptType"
  ) %in% colnames(GTF7)))
})

test_that("refflat_create", {
  ref_flat <- GTF.tool::refflat_create(GTF2)
  expect_gt(nrow(ref_flat), 0)
  expect_true(all(c(
    "geneName", "name", "chrom", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"
  ) %in% colnames(ref_flat)))
  expect_true(any(!grepl("e+", ref_flat$exonEnds)))
})

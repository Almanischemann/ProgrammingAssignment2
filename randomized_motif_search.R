# Randomized Motif Search implementation with pseudocounts

count_with_pseudocounts <- function(motifs) {
  k <- nchar(motifs[1])
  counts <- matrix(1, nrow = 4, ncol = k,
                   dimnames = list(c("A", "C", "G", "T"), NULL))
  for (motif in motifs) {
    bases <- strsplit(motif, split = "")[[1]]
    for (i in seq_len(k)) {
      base <- bases[i]
      counts[base, i] <- counts[base, i] + 1
    }
  }
  counts
}

profile_with_pseudocounts <- function(motifs) {
  counts <- count_with_pseudocounts(motifs)
  counts / colSums(counts)
}

consensus <- function(motifs) {
  k <- nchar(motifs[1])
  counts <- count_with_pseudocounts(motifs)
  apply(counts, 2, function(col) {
    names(which.max(col))
  }) |> paste0(collapse = "")
}

score_motifs <- function(motifs) {
  cons <- consensus(motifs)
  motif_matrix <- do.call(rbind, strsplit(motifs, split = ""))
  mismatches <- sum(motif_matrix != rep(strsplit(cons, split = "")[[1]],
                                        each = nrow(motif_matrix)))
  mismatches
}

profile_most_probable_kmer <- function(text, k, profile) {
  best_prob <- -1
  best_kmer <- substr(text, 1, k)
  n <- nchar(text)
  for (start in 1:(n - k + 1)) {
    kmer <- substr(text, start, start + k - 1)
    bases <- strsplit(kmer, split = "")[[1]]
    probs <- mapply(function(base, idx) profile[base, idx], bases, seq_along(bases))
    prob <- prod(probs)
    if (prob > best_prob) {
      best_prob <- prob
      best_kmer <- kmer
    }
  }
  best_kmer
}

randomized_motif_search_once <- function(Dna, k, t) {
  motifs <- character(t)
  for (i in seq_len(t)) {
    text <- Dna[i]
    n <- nchar(text)
    start <- sample.int(n - k + 1, 1)
    motifs[i] <- substr(text, start, start + k - 1)
  }
  best_motifs <- motifs

  repeat {
    profile <- profile_with_pseudocounts(motifs)
    motifs <- vapply(Dna, profile_most_probable_kmer, character(1), k = k, profile = profile)
    if (score_motifs(motifs) < score_motifs(best_motifs)) {
      best_motifs <- motifs
    } else {
      return(best_motifs)
    }
  }
}

randomized_motif_search <- function(Dna, k, t, iterations = 1000) {
  best_motifs <- randomized_motif_search_once(Dna, k, t)
  for (i in seq_len(iterations - 1)) {
    motifs <- randomized_motif_search_once(Dna, k, t)
    if (score_motifs(motifs) < score_motifs(best_motifs)) {
      best_motifs <- motifs
    }
  }
  best_motifs
}

run <- function() {
  tokens <- scan(file = "stdin", what = character(), quiet = TRUE)
  if (length(tokens) < 3) {
    stop("Input must contain k, t, and at least one DNA string")
  }
  k <- as.integer(tokens[1])
  t <- as.integer(tokens[2])
  Dna <- tokens[-c(1, 2)]
  if (length(Dna) != t) {
    stop("Number of DNA strings does not match t")
  }
  result <- randomized_motif_search(Dna, k, t, iterations = 1000)
  cat(paste(result, collapse = " "))
}

if (identical(environment(), globalenv())) {
  run()
}

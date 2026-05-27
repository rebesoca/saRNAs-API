infoGenProm <- function(geneVector, upstream = 0){
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(httr)
    library(jsonlite)
  })
  
  # Load EPD promoters data
  data_promo <- read.table(
    file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"Hs_EPDnew.bed"),
    quote = "\"",
    comment.char = ""
  )
  
  colnames(data_promo) <- c("Chromosome", "Start", "End", "ID_promoter", "Score", "Strand", "StartP", "EndP")
  result <- data.frame()
  
  for (i in 1:length(geneVector)) {
    
    gene <- geneVector[i]
    
    # Skip NAs
    if (is.na(gene)) {
      next
    }
    
    x <- 1
    found <- FALSE
    
    repeat {
      gene_x <- sprintf("%s_%d", gene, x)
      r_x <- data_promo[data_promo$ID_promoter == gene_x, ]
      
      if (nrow(r_x) == 0) {
        break
      } else {
        r_x$Source <- "EPD"
        result <- rbind(result, r_x)
        found <- TRUE
        x <- x + 1
      }
    }
    
    # If not in EPD, search in Ensembl TSS
    if(!found){
      
      message(paste("Gene", gene, "not found in EPD. Searching in Ensembl..."))
      
      tryCatch({
        
        # Ensembl query
        server <- "https://rest.ensembl.org"
        ext <- paste0("/lookup/symbol/human/", gene, "?content-type=application/json")
        r <- GET(paste(server, ext, sep = ""))
        stop_for_status(r) 
        
        info <- fromJSON(content(r, as = "text", encoding = "UTF-8"))
        
        chr <- paste0("chr", info$seq_region_name)
        strand <- info$strand
        
        # Calculate promoter coordinates
        if(strand == 1) {
          startP <-  info$start - upstream
          endP <- info$start
        } else {
          startP <- info$end
          endP <- info$end + upstream
        }
        
        row <- data.frame(
          Chromosome = chr,
          Start = startP,
          End = endP,
          ID_promoter = gene,
          Score = NA,
          Strand = strand,
          StartP = startP,
          EndP = endP,
          Source = "Ensembl"
        )
        result <- rbind(result, row)
        
      }, error = function(e){
        message(paste("Gene", gene, "not found in EPD or Ensembl:", e$message))
      })
      
    }
    
  }
  
  return(result)
}

getPromoSeq <- function(x, upstring = 2000, downstring = 0){
  
  suppressPackageStartupMessages({
    library(httr)
    library(jsonlite)
    library(xml2)
  })
  
  result <- data.frame()
  
  for (i in 1:nrow(x)) {
    tryCatch({
      
      start <- as.numeric(x$Start[[i]])
      end <- as.numeric(x$End[[i]])
      chr <- x$Chromosome[[i]]
      chr2 <- gsub("^chr(.*)$", "\\1", chr)
      
      # Extend the region if the user ask for it
      if (upstring != 0) start <- start - upstring
      if (downstring != 0) end <- end + downstring
      
      # Ensembl query
      server <- "https://rest.ensembl.org"
      ext <- paste0("/sequence/region/human/", chr2, ":", start, "..", end, "?")
      r <- GET(paste(server, ext, sep = ""), content_type("text/x-fasta"))
      stop_for_status(r)
      
      # Sequence extraction
      r_split <- strsplit(content(r), split = "\n")
      symbol <- x$ID_promoter[[i]]
      chain <- paste0(r_split[[1]][-1], collapse = "")
      
      fila <- data.frame(ID_Promoter = symbol, Chain = chain)
      result <- rbind(result, fila)
    }, error = function(e) {
      message(paste("Error retrieving sequence for:", x$ID_promoter[[i]], "-", e$message))
    })
  }
  
  return(result)
}

getSeqT <- function(transcriptID){
  
  suppressPackageStartupMessages({
    library(httr)
    library(stringr)
    library(jsonlite)
  })
  
  # Ensembl query
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", transcriptID, "?type=cdna")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  
  chain <- fromJSON(content(r, as = "text", encoding = "UTF-8"))
  seq <- chain$seq
  
  
  return(seq)
}

alignSequences <- function(geneSeq, promoterSeq, identityThreshold = 0.6){
  
  message(paste("Started at:", format(Sys.time(), "%H:%M:%S")))
  start_time <- Sys.time()
  
  results <- data.frame()
  
  # Create tmp for fasta files
  tmp_dir <- tempdir()
  
  # Write gene sequence in fasta file
  gene_fasta <- file.path(tmp_dir, "gene.fasta")
  writeLines(c(">Gene", geneSeq), gene_fasta)
  
  for (i in 1:nrow(promoterSeq)) {
    
    promoter_id <- promoterSeq$ID_Promoter[i]
    promoter_seq <- promoterSeq$Chain[i]
    
    message(paste("Processing promoter:", promoter_id))
    
    tryCatch({
      
      # Write promoter fasta file
      promoter_fasta <- file.path(tmp_dir, paste0(promoter_id, ".fasta"))
      writeLines(c(paste0(">", promoter_id), promoter_seq), promoter_fasta)
      
      # Output file
      output_file <- file.path(tmp_dir, paste0(promoter_id, "_result.txt"))
      
      # Matcher
      system2("matcher", 
              args = c(
                "-asequence", gene_fasta,
                "-bsequence", promoter_fasta,
                "-outfile", output_file,
                "-aformat3", "markx0"
              ),
              stdout = F, stderr = F)
      
      # Parse output
      if (file.exists(output_file)){
        
        output <- readLines(output_file)
        #print(output)
        
        #Extract identity %
        identity_line <- grep("^# Identity:", output, value = T)[1]
        identity_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", identity_line))
        
        simi_line <- grep("^# Similarity:", output, value = T)[1]
        simi_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", simi_line))
        
        gaps_line <- grep("^# Gaps:", output, value = T)[1]
        gaps_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", gaps_line))
        
        # Extract length
        length_line <- grep("^# Length:", output, value = T)[1]
        alig_length <- as.numeric(gsub("# Length:\\s+", "", length_line))
        
        # Extract score
        score_line <- grep("^# Score:", output, value = T)[1]
        alig_score <- as.numeric(gsub("# Score:\\s+", "", score_line))
        
        # Extract aligned sequences
        gene_lines <- grep("^\\s*Gene\\s", output, value = T)
        promoter_lines <- grep(paste0("^\\s*", promoter_id, "\\s"), output, value = T)
        
        if (length(promoter_lines) == 0){
          short_id <- substr(promoter_id, 1, 6)
          promoter_lines <- grep(paste0("^\\s*", short_id, "\\s"), output, value = T)
        }
        
        # Extract sequences part for each line
        extract_seq <- function(lines){
          seqs <- sapply(lines, function(line){
            # Keep sequences
            parts <- strsplit(trimws(line), "\\s+")[[1]]
            parts[length(parts)]
          })
          paste0(seqs, collapse = "")
        }
        
        gene_alig <- extract_seq(gene_lines)
        promoter_alig <- extract_seq(promoter_lines)
        
        # Clean gaps
        gene_alig_clean <- gsub("-", "", gene_alig)
        promoter_alig_clean <- gsub("-", "", promoter_alig)
        
        if (!is.na(identity_pct) && identity_pct >= identityThreshold * 100){
          results <- rbind(results, data.frame(
            Gene = sub("_\\d+$", "", promoter_id),
            ID_Promoter = promoter_id,
            Identity_pct = identity_pct,
            Similarity_pct = simi_pct,
            Gaps_pct = gaps_pct,
            Alignment_length = alig_length,
            Score = alig_score,
            Gene_aligned = gene_alig_clean,
            Promoter_aligned = gsub("T", "U", promoter_alig_clean),
            Gene_aligned_gaps = gene_alig,
            Promoter_aligned_gaps = promoter_alig
          ))
        }
      }
    }, error = function(e){
      message(paste("Error procesing promoter:", promoter_id, "-", e$message))
    })
  }
  
  end_time <- Sys.time()
  total_time <- round(difftime(end_time, start_time, units = "mins"), 2)
  message(paste("Finished at:", format(Sys.time(), "%H:%M:%S")))
  message(paste("Total time:", total_time, "minutes"))
  
  if (nrow(results) == 0){
    message(paste("No aligments found with identity >=", identityThreshold *100,
                  "%. Consider lowering the threshold."))
    return(data.frame())
  }
  
  # Sort by identity
  results <- results[order(-results$Identity_pct),]
  rownames(results) <- NULL
  
  return(results)
}

extractCandidates <- function(alignmentsResults){
  
  suppressPackageStartupMessages(library(Biostrings))
  
  candidates <- data.frame()
  
  for (i in 1:nrow(alignmentsResults)){
    
    gene_alig_gaps <- gsub("U", "T", alignmentsResults$Gene_aligned_gaps[i])
    promoter_alig_gaps <- alignmentsResults$Promoter_aligned_gaps[i]
    promoter_id <- alignmentsResults$ID_Promoter[i]
    gene_name <- alignmentsResults$Gene[i]
    
    # Conserve global alignment metrics
    global_identity <- alignmentsResults$Identity_pct[i]
    global_similarity <- alignmentsResults$Similarity_pct[i]
    global_gaps <- alignmentsResults$Gaps_pct[i]
    global_score <- alignmentsResults$Score[i]
    global_length <- alignmentsResults$Alignment_length[i]
    
    if (nchar(gene_alig_gaps) < 19) next
    
    for (win_size in c(19, 20, 21)) {
      
      # Slide window over gapped alignment
      pos <- 1
      non_gap_count <- 0
      window_starts <- c()
      
      for (start in 1:(nchar(gene_alig_gaps) - win_size +1)) {
        
        window_gapped <- substr(gene_alig_gaps, start, start + win_size - 1)
        window_clean <- gsub("-", "", window_gapped)
        
        if (nchar(window_clean) == win_size){
          
          promoter_window_gapped <- substr(promoter_alig_gaps, start, start + win_size - 1)
          promoter_window_clean <- gsub("-", "", promoter_window_gapped)
          
          if (nchar(promoter_window_clean) < 1) next
          
          window_rc <- as.character(reverseComplement(DNAString(window_clean)))
          
          # Calculate identity against promoter window
          min_len <- min(nchar(window_rc), nchar(promoter_window_clean))
          rc_chars <- strsplit(substr(window_rc, 1, min_len), "")[[1]]
          prom_chars <- strsplit(substr(promoter_window_clean, 1, min_len), "")[[1]]
          matches <- sum(rc_chars == prom_chars)
          identity <- matches / win_size
          
          candidates <- rbind(candidates, data.frame(
            Gene = gene_name,
            ID_Promoter = promoter_id,
            Global_identity_pct = global_identity,
            Global_similarity_pct = global_similarity,
            Global_gaps_pct = global_gaps,
            Global_score = global_score,
            Global_alignment_length = global_length,
            Window_size = win_size,
            Start_position = start,
            saRNA_sequence = gsub("T", "U", window_clean),
            saRNA_RC = gsub("T", "U", window_rc),
            Promoter_region = gsub("T", "U", promoter_window_clean),
            Window_identity_pct = round(identity * 100, 2)
          ))
        }
      }
    }
  }
  
  if (nrow(candidates) == 0){
    message("Candidates not found.")
    return(data.frame())
  }
  
  candidates <- candidates[order(-candidates$Global_score), ]
  rownames(candidates) <- NULL
  
  return(candidates)
}

interpretAlignment <- function(candidates){
  
  suppressPackageStartupMessages({
    library(Biostrings)
  })
  
  if(nrow(candidates) == 0) return(data.frame())
  
  message("Calculating seed region identity and folding energies...")
  message("Processing ", nrow(candidates), " candidates...")
  
  for (i in 1:nrow(candidates)){
    
    saRNA_seq <- gsub("U", "T", candidates$saRNA_sequence[i])
    saRNA_rc <- gsub("U", "T", candidates$saRNA_RC[i])
    promoter_region <- gsub("U", "T", candidates$Promoter_region[i])
    
    # B: Seed region analysis (Check if seed region has perfect complementarity)
    if (nchar(saRNA_rc) >= 8 & nchar(promoter_region) >= 8){
      seed_saRNA <- substr(saRNA_rc, 2, 8)
      seed_promoter <- substr(promoter_region, 2, 8)
      seed_matches <- sum(strsplit(seed_saRNA, "")[[1]] == strsplit(seed_promoter, "")[[1]])
      score_seed <- round((seed_matches / 7) * 100, 2)
    } else {
      score_seed <- NA
    }
    
    # Delta1: Secondary structure RNA alone
    tryCatch({
      rnafold_input <- gsub("T", "U", saRNA_seq)
      rnafold_output <- system2("RNAfold",
                                input = rnafold_input,
                                stdout = TRUE,
                                stderr = FALSE)
      # Get output
      dg1_line <- rnafold_output[2]
      dg1 <- as.numeric(gsub(".*\\(\\s*(-?[0-9.]+)\\s*\\).*", "\\1", dg1_line))
      structure_saRNA <- trimws(gsub("\\s*\\(.*", "", dg1_line))
    }, error = function(e){
      dg1 <- NA
      structure_saRNA <- NA
    })
    
    # Delta2: Hybridization energy saRNA-promoter (RNAcofold)
    tryCatch({
      rnacofold_input <- paste0(
        gsub("T", "U", saRNA_rc),
        "&",
        gsub("T", "U", promoter_region)
      )
      rnacofold_output <- system2("RNAcofold",
                                  input = rnacofold_input,
                                  stdout = TRUE,
                                  stderr = FALSE)
      # Get output
      dg2_line <- rnacofold_output[2]
      dg2 <- as.numeric(gsub(".*\\(\\s*(-?[0-9.]+)\\s*\\).*", "\\1", dg2_line))
      structure_duplex <- trimws(gsub("\\s*\\(.*", "", dg2_line))
    }, error = function(e){
      dg2 <- NA
      structure_duplex <- NA
    })
    
    #candidates$Gen <- candidates$Gen
    #candidates$ID_Promoter <- candidates$ID_Promoter
    #candidates$Window_size <- candidates$Window_size
    #candidates$Start_position <- candidates$Start_position
    candidates$Structure_saRNA[i] <- structure_saRNA
    candidates$Structure_duplex[i] <- structure_duplex
    candidates$Score_seed[i] <- score_seed
    candidates$DeltaG1_saRNA[i] <- dg1
    candidates$DeltaG2_duplex[i] <- dg2
    #candidates$Global_score <- candidates$Global_score
    #candidates$Global_similarity <- candidates$Global_similarity
  }
  
  # Sort by DeltaG2
  candidates <- candidates[order(-candidates$Global_score, candidates$DeltaG2_duplex), ]
  rownames(candidates) <- NULL
  
  return(candidates)
}

filterCandidates <- function(candidates, max_deltaG1 = -1,
                             max_gaps_pct = 10, min_score_seed = 30.0,
                             min_window_identity = 30.0){
  
  if (nrow(candidates) == 0){
    message("No candidates to filter.")
    return(data.frame())
  }
  
  n1 <- nrow(candidates)
  message("Candidates before filtering:", n1)
  
  # Filter 1
  f1 <- candidates[!is.na(candidates$DeltaG1_saRNA) & 
                     candidates$DeltaG1_saRNA > max_deltaG1, ]
  
  # Filter 2
  f2 <- f1[!is.na(f1$Global_gaps_pct) & 
             f1$Global_gaps_pct < max_gaps_pct, ]
  
  # Filter 3
  f3 <- f2[!is.na(f2$Score_seed) & 
             f2$Score_seed > min_score_seed, ]
  # Filter 4
  f4 <- f3[!is.na(f3$Window_identity_pct) & 
             f3$Window_identity_pct > min_window_identity, ]
  
  if (nrow(f4) == 0){
    message("No candidates pased all filters. Consider relaxing thresholds.")
  }
  
  message(paste("Final candidates:", nrow(f4), 
                "(", round(nrow(f4)/n1*100, 1), "% of original)"))
  
  rownames(f4) <- NULL
  return(f4)
  
}

#
# This is a Plumber API. In RStudio 1.2 or newer you can run the API by
# clicking the 'Run API' button above.
#
# In RStudio 1.1 or older, see the Plumber documentation for details
# on running the API.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)
# Source
source("functions.R")

#* @apiTitle saRNA Prediction API
#* @apiDescription API for predicting saRNA sequences from genes
#* @apiContact list(name = "API support", email = "r.sotocampillo@um.es")
#* @apiVersion 1.0.0

#* Hello world
#* @get /hello-world
function(){
  "Hello world!
  saRNA Prediction API is working"
}

#* Predict candidate saRNA sequences
#* @param gene Gen symbol (e.g. ATOH1)
#* @param transcript_id Optional: Ensembl transcript ID (e.g. ENST00000306011). If empty the API will search automatically.
#* @param similarity Functional similarity threshold between genes (default: 0.8)
#* @param min_similarity Minimum functional similarity accepted (default: 0.5)
#* @param identity_threshold Sequence identity threshold for alignment (default: 0.6)
#* @param upstream Length in base pairs of the promoter region upstream of TSS (default: 2000)
#* @param max_delta_saRNAs Maximum DeltaG1 allowed for saRNA secondary structure. Values close to 0 indicate accessible sequences for AGO2 loading (default: -0.1)
#* @param max_gapspct Maximum percentage of gaps allowed in the global alignment. High gap percentage indicates poor alignment quality (default: 10)
#* @param min_score_seed Minimum identity percentage required in the seed region (positions 2-8). Critical for saRNA activity according to Meng et al. 2016 (default: 30)
#* @param min_window Minimum identity percentage required in the candidate window. Low values indicate insufficient complementarity for stable duplex formation (default: 30)
#* @get /predict
#* 


function(gene, transcript_id = NULL, similarity = 0.8, 
         min_similarity = 0.5, identity_threshold = 0.6, upstream = 2000 ,max_delta_saRNAs = -1,
         max_gapspct =  10, min_score_seed = 30, min_window = 30){
  
  start_time <- Sys.time()
  
  suppressPackageStartupMessages({
    library(httr)
    library(jsonlite)
    library(reticulate)
    reticulate::py_require("pymongo")
  })
  
  # Convert parameters to correct types
  similarity <- as.numeric(similarity)
  min_similarity <- as.numeric(min_similarity)
  identity_threshold <- as.numeric(identity_threshold)
  upstream <- as.numeric(upstream)
  max_delta_saRNAs  <- as.numeric(max_delta_saRNAs)
  max_gapspct <- as.numeric(max_gapspct)
  min_score_seed <- as.numeric(min_score_seed)
  min_window <- as.numeric(min_window)
  
  # Parameter validation
  if (similarity < 0 | similarity > 1){
    stop("Similarity value must be between 0 and 1")
  }
  
  if (min_similarity < 0 | min_similarity > similarity){
    stop("Minimun functional similarity accepted must be between 0 and similarity value.")
  }
  
  message(paste("All the parameters have been received correctly."))
  
  if(is.null(transcript_id) | length(transcript_id) == 0 | identical(transcript_id, "")){
    transcript_id <- NULL
  }
  
  # 1: Get the transcript ID if not provided
  if(is.null(transcript_id)){
    message(paste("Getting the transcript ID for", gene))
    tryCatch({
      
      # Ensembl query
      server <- "https://rest.ensembl.org"
      ext <- paste0("/lookup/symbol/human/", gene, "?content-type=application/json")
      r <- GET(paste(server, ext, sep = ""))
      stop_for_status(r)
      info <- fromJSON(content(r, as = "text", encoding = "UTF-8"))
      
      # Removing version suffix
      transcript_id <- sub("\\..*$", "", info$canonical_transcript)
      message(paste("Transcript ID found:", transcript_id))
    
    }, error = function(e){
      stop(paste("Could not find transcript ID for gene", gene, ":", e$message))
    })
  } else {
    
    # Clean transcript ID if provided manually
    transcript_id <- sub("\\..*$", "", transcript_id)
    message(paste("Using provided transcript ID:", transcript_id))
    
  }
  
  # 2: Funtional similarity search in MongoDB
  message(paste("Searching functionally similar genes for", gene))
  tryCatch({
    py_script <- reticulate::py_run_file("similarity.py")
    similar_genes <- py_script$calculate_similarity(
      gene = gene, similarity = similarity,
      min_similarity = min_similarity
    )
    
    if (length(similar_genes) == 0){
      stop(paste("No functionally similar genes found for", gene))
    }
    
    gene_names <- sapply(similar_genes, function(x) x$gene_name)
    message(paste("Found", length(gene_names), "similar genes:", paste(gene_names, collapse = ", ")))
    
  }, error = function(e){
    stop(paste("Error in similarity search:", e$message))
  })
  
  # 3: Get promoter coordinates
  message("Getting promoter coordinates...")
  
  tryCatch({
    positions <- infoGenProm(gene_names)
    
    if (nrow(positions) == 0){
      stop("No promoter sequences found for similar genes")
    }
    
    message(paste("Found promoters for", nrow(positions), "sequences"))
    
  }, error = function(e){
    stop(paste("Error getting promoter coordinates:", e$message))
  })
  
  # 4: Get promoter sequences
  message("Getting promoter sequences...")
  
  tryCatch({
    promoter_seqs <- getPromoSeq(positions, upstring = upstream)
    
    if (nrow(promoter_seqs) == 0){
      stop("Could not retrieve promoter sequences")
    }
    message(paste("Retrieved", nrow(promoter_seqs), "promoter sequences"))
    
  }, error = function(e){
    stop(paste("Error getting promoter sequences:", e$message))
  })
  
  # 5: Get transcript sequence
  message(paste("Getting transcript sequence for", transcript_id))
  
  tryCatch({
    truncated_seq <- getSeqT(transcript_id)
    message(paste("Sequence length:", nchar(truncated_seq), "bases"))
    
  }, error = function(e){
    stop(paste("Error getting transcript sequence:", e$message))
  })
  
  # 6: Alignment and saRNA prediction
  message("Running alignment and predicting saRNAs...")
  
  tryCatch({
    sarnas <- alignSequences(truncated_seq, promoter_seqs, 
                             identityThreshold = identity_threshold)
    
    if (nrow(sarnas) == 0){
      return(list(
        status = "success",
        gene = gene,
        transcript_id = transcript_id,
        similar_genes = gene_names,
        message = "No candidate saRNAs found. Consider lowering identity_threshold.",
        results = list()
      ))
    }
    
  }, error = function(e){
    stop(paste("Error in alignment:", e$message))
  })
  
  # 7: Extracting candidates...
  message("Extracting saRNAs...")
  
  tryCatch({
    sarnas2 <- extractCandidates(sarnas)
    
  }, error = function(e){
    stop(paste("Error in extraction:", e$message))
  })
  
  # 8: Interpret alignment results...
  message("Interpreting alignment...")
  
  tryCatch({
    sarnas3 <- interpretAlignment(sarnas2)
  }, error = function(e){
    stop(paste("Error in alignment interpretation:", e$message))
  })
  
  # 9: Filtering...
  message("Filtering...")
  
  tryCatch({
    sarnas4 <- filterCandidates(sarnas3, max_deltaG1 = max_delta_saRNAs,
                                max_gaps_pct = max_gapspct, min_score_seed = min_score_seed,
                                min_window_identity = min_window)
  }, error = function(e){
    stop(paste("Error in filtering:", e$message))
  })
  
  end_time <- Sys.time()
  time <- round(difftime(end_time, start_time, units = "secs"), 2)
  message(paste("Total pipeline time", time, "seconds"))
  
  # Return results
  return(list(
    status = "success",
    gene = gene,
    transcript_id = transcript_id,
    similar_genes = gene_names,
    n_candidates = nrow(sarnas4),
    results = sarnas4
  ))
  
}








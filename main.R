# Francesco Cabiddu, CabidduF@cardiff.ac.uk
# Study 1 - PhD, Cardiff University

# load libraries ----------------------------------------------------------
lib <- c(
  "magrittr", "tidyverse",
  "stringr",
  "beepr", "fastmatch",
  "mailR", "childesr", "viridis"
  #"stringdist", "data.table"
)
lapply(lib, require, character.only = TRUE)
rm(lib)


# homemade funs -----------------------------------------------------------
convert_phonemes <- function(df) {
  # create a column of converted phonemes to single characters for each word
  df %>%
    mutate(phon_converted = phon %>%
             str_split("_") %>%
             sapply(function(i) {
               paste0(phonemes_converted[i], collapse="")
             }))
}

adj_poss <- function(x) {
  # all possible adjacent combinations of phonemes in a word, including the word itself
  n <- length(x)
  if(n == 1L) return(NA) # if a word is monophonemic return NA
  idx <- expand.grid(start = 1L:n, len = 2L:(n))
  idx$end <- idx$start + idx$len - 1L
  idx <- idx[idx$end <= n, ]
  Map(function(start, end) x[start:end], idx$start, idx$end)
}

subset_utterances <- function(df, corpora_age, extension_random) {
  # Note. this function matches N utterances and MLU but does not consider MLU SD.
  df %>%
    filter(speaker_role != "Target_Child") %>%
    filter(target_child_age_year %in% corpora_age) %>%
    group_by(target_child_age_year) %>%
    sample_n(corpora$all %>%
               filter(speaker_role != "Target_Child") %>%
               group_by(target_child_age_year) %>%
               summarise(N_utt = n()) %>%
               ungroup %>%
               filter(target_child_age_year %in% 1:4) %>%
               {.$N_utt} %>%
               min %>%
               # first subset each year number of utterances similar to minimum
               {. + extension_random}) %>%
    ungroup %>%
    mutate(LU = phon %>%
             str_split(" ") %>%
             sapply(length)) %>%
    arrange(target_child_age_year, desc(LU)) %>%
    group_by(target_child_age_year) %>%
    group_split %>%
    lapply(function(DF) {
      lu_reference <- corpora$all %>%
        filter(speaker_role != "Target_Child") %>%
        filter(target_child_age_year %in% 1) %>%
        mutate(LU = phon %>%
                 str_split(" ") %>%
                 sapply(length)) %>%
        arrange(target_child_age_year, desc(LU)) %>% 
        summarise(MLU = mean(LU)) %>%
        {.$MLU}
      
      while (mean(DF$LU) > lu_reference) {
        # then exclude utterances from longest to shorter to match MLU
        DF <- DF[-1, ]
      }
      
      DF 
    }) %>%
    bind_rows
}

save_input <- function(utt, file_name, segmented) {
  path <- "/Users/francesco/Dropbox/phd/Incrementality/Segmentation/input/"
  
  if (segmented == TRUE) {
    write_lines(utt, paste(path, file_name, sep = ""))
  } else if (segmented == FALSE) {
    write_lines(utt %>%
                  str_remove_all("[|]"), paste(path, file_name, sep = ""))
  }
}

# import corpora ----------------------------------------------------------
corpora <- list() %>%
  (function(x) {
    x$all <- get_utterances(
      # database version used: '2018.1'
      collection = "Eng-UK", 
      corpus = c("Belfast", "Fletcher", "Manchester", 
                 "Thomas", "Tommerdahl", "Wells",
                 "Forrester", "Lara")
    ) %>%
      # only select columns of interest
      select(corpus_id = corpus_name,
             transcript_id,
             target_child_id,
             target_child_name,
             target_child_sex,
             target_child_age,
             speaker_role,
             speaker_name, 
             speaker_id,
             speaker_code,
             utt_id = utterance_order,
             utt_type = type, 
             gloss,
             stem,
             part_of_speech) %>%
      arrange(transcript_id) %>%
      # symbols + and _ are used interchangeably for compounds, 
      # so convert + to _ for consistency
      mutate(gloss = gloss %>%
               str_replace_all("\\+", "_")) %>%
      # every word in gloss and stem columns to lowercase
      mutate_at(.funs = list(~tolower(.)), .vars = vars(gloss:stem)) 
    
    return(x)
  })

# import CMU and words transcription ----------------------------
load("~/Dropbox/phd/Incrementality/analysis_R/CMU_DICT.RData")

CMU_DICT %<>%
  (function(list_dfs) {
    additional_compounds <- list_dfs$words %>%
      # extract some compounds that have the _ sign, 
      # and are still included in the single words dataframe!
      filter(str_detect(word, "_"))
    
    list_dfs$words %<>%
      # modify words for consistency with corpora  
      filter(!str_detect(word, "_")) %>%
      mutate(word = word %>%
               tolower %>%
               str_replace_all("@", "'"))
    
    list_dfs$compounds %<>%
      # add additional compounds to compound dataframe
      # and modify words for consistency with corpora
      rbind(
        additional_compounds
      ) %>%
      mutate(word = word %>%
               tolower %>%
               str_replace_all("[+]{1}", "_") %>%
               str_replace_all("@", "'"))
    
    list_dfs
  })

corpora$converted <- tibble(word_type = corpora$all$gloss %>%
                        str_split(" ") %>%
                        unlist %>%
                        unique %>%
                        sort) %>%
  (function(x) {
    # create a similar CMU list, with dataframes for single words and compounds
    list(single_word = filter(x, !str_detect(word_type, "[_+]+")) %>%
           # add a column with phonetic transcription
           # if no trascription is available, NA value is assigned
           left_join(., CMU_DICT$words, by = c("word_type" = "word")),
         multi_word = filter(x, str_detect(word_type, "[_+]+")) %>%
           left_join(., CMU_DICT$compounds, by = c("word_type" = "word")) %>%
           # convert each single word within a compound
           # if one sigle words is not translated assign NA to compound
           (function(df) {
             df_ready <- df %>%
               # leave out compounds already translated with previous step
               filter(!is.na(phon))
             
             df_not_ready <- df %>%
               filter(is.na(phon)) %>% 
               mutate(phon = word_type %>%
                        str_split("_") %>%
                        sapply(function(y) {
                          CMU_DICT$words$phon[fmatch(y,
                                                     CMU_DICT$words$word)] %>%
                            paste(collapse = "_")
                        })) %>%
               mutate(phon = ifelse(str_detect(phon, "NA"), NA, phon))
             
             df_ready %>%
               rbind(df_not_ready)
           }))
  }) %>%
  # reunite single words and compounds into a single dataframe, unpacking the list
  bind_rows %>%
  arrange(word_type)

corpora$all %<>%
  # assign phon variable with transcribed utterances (takes ~5 min)
  mutate(phon = gloss %>%
           str_split(" ") %>%
           lapply(function(utt) {
             corpora$converted$phon[fmatch(utt, 
                                           corpora$converted$word_type)] %>%
               paste(collapse = " ")
           })) %>%
  # if a word in an utterance is not transcribed, assign NA to whole utterance
  mutate(phon = ifelse(str_detect(phon, "NA"), NA, phon) %>%
           unlist) %>%
  # filter out utterances that are not transcribed correctly
  filter(!is.na(phon)) %>%
  # create a year age variable 
  mutate(target_child_age_year = ifelse(target_child_age < 24, 1,
                                        ifelse(target_child_age >= 24 & target_child_age < 36, 2,
                                               ifelse(target_child_age >= 36 & target_child_age < 48, 3,
                                                      ifelse(target_child_age >= 48 & target_child_age < 60, 4,
                                                             ifelse(target_child_age >= 60 & target_child_age < 72, 5,
                                                                    ifelse(target_child_age >= 72 & target_child_age < 84, 6,
                                                                           ifelse(target_child_age >= 84 & target_child_age < 96, 7,
                                                                                  ifelse(target_child_age >= 96 & target_child_age < 108, 8, 9))))))))) %>%
  
  rename(target_child_age_month = target_child_age) %>%
  select(corpus_id:target_child_age_month, target_child_age_year, speaker_role:phon) 

# corpora preparation  -----------------------------
# convert phonemes into single characters 
phonemes_converted <- c(letters, LETTERS)[1:39] %>%
  # convert each possible phoneme in a single letter (for word segmentation algorithms)
  (function(x) {
    names(x) <- corpora$all$phon %>%
      str_split(" ") %>% 
      unlist %>%
      str_split("_") %>%
      unlist %>%
      unique %>%
      sort
    
    x
  })

corpora$converted %<>%
  # create a variable of converted phonetic words
  convert_phonemes %>%
  mutate(phon_converted = ifelse(str_detect(phon_converted, "NA"), NA, phon_converted)) 

segmentation <- list() %>%
  # match corpora by age (sampling equal number of utterances)
  (function(LIST) {
    LIST$input <- list()
    
    set.seed(18854)
    LIST$input$basic <- corpora$all %>%
      # select only utterances spoken by speakers who are not the child
      # and sample number of utterances as less populated year
      filter(speaker_role != "Target_Child") %>%
      filter(target_child_age_year %in% 1:4) %>%
      group_by(target_child_age_year) %>%
      sample_n(corpora$all %>%
                 filter(speaker_role != "Target_Child") %>%
                 group_by(target_child_age_year) %>%
                 summarise(N_utt = n()) %>%
                 ungroup %>%
                 filter(target_child_age_year %in% 1:4) %>%
                 {.$N_utt} %>%
                 min) %>%
      ungroup %>%
      arrange(target_child_age_year) %>%
      # created a variable of converted phonetic utterances (word separator = "|")
      mutate(phon_converted = phon %>%
               str_split(" ") %>%
               lapply(function(word) {
                 corpora$converted$phon_converted[fmatch(word, corpora$converted$phon)] %>%
                   paste(collapse = "|")
               })) %>%
      # add separator at the end of utterance
      mutate(phon_converted = phon_converted %>%
               str_replace("$", "|"))
    
    return(LIST)
  })

segmentation$input$mlu <- corpora$all %>%
  # match corpora by age and MLU (takes ~20 min)
  (function(DF) {
    set.seed(219900)
    
    DF %>%
      filter(speaker_role != "Target_Child") %>%
      filter(target_child_age_year %in% 1) %>%
      mutate(LU = phon %>%
               str_split(" ") %>%
               sapply(length)) %>%
      arrange(target_child_age_year, desc(LU)) %>%
      rbind(
        subset_utterances(corpora$all, 
                          corpora_age = 2, 
                          extension_random = 6000)
      ) %>%
      rbind(
        subset_utterances(corpora$all, 
                          corpora_age = 3, 
                          extension_random = 11000)
      ) %>%
      rbind(
        subset_utterances(corpora$all, 
                          corpora_age = 4, 
                          extension_random = 13000)
      ) %>%
      # created a variable of converted phonetic utterances (word separator = "|")
      mutate(phon_converted = phon %>%
               str_split(" ") %>%
               lapply(function(word) {
                 corpora$converted$phon_converted[fmatch(word, corpora$converted$phon)] %>%
                   paste(collapse = "|")
               })) %>%
      # add separator at the end of utterance
      mutate(phon_converted = phon_converted %>%
               str_replace("$", "|")) %>%
      # reshuffle the utterances at each year 
      # to avoid descending order of length of utterance
      group_by(target_child_age_year) %>%
      sample_n(n()) %>%
      ungroup
  })

# save input to workspace

#save_input(segmentation$input$basic$phon_converted,
#           "input_basic_segmented.txt",
#           segmented = TRUE)

#save_input(segmentation$input$basic$phon_converted,
#           "input_basic_unsegmented.txt",
#           segmented = FALSE)

#save_input(segmentation$input$mlu$phon_converted,
#           "input_mlu_segmented.txt",
#           segmented = TRUE)

#save_input(segmentation$input$mlu$phon_converted,
#           "input_mlu_unsegmented.txt",
#           segmented = FALSE)

# spok_bnc, plurals, lexical measures ----------------------------
# import Spoken BNC types with associated lexical measures
load("~/Dropbox/phd/Incrementality/analysis_R/spok_bnc_lexical.RData")
# import Spoken BNC tokens
load("~/Dropbox/phd/Incrementality/analysis_R/spok_bnc_tokens.RData")

spok_bnc <- "empty" %>%
  (function(empty) {
    lexical_measures <- measures %>%
      # only select spok_bnc measures
      select(phon, contains("spok_bnc")) %>%
      mutate(phonemic_len = phon %>%
               str_split("_") %>%
               sapply(length)) %>%
      (function(x) {
        colnames(x) <- str_remove(colnames(x), "spok_bnc_")
        x
      }) %>%
      filter(!is.na(freq))
    
    list(
      all = SPOK_BNC,
      lexical_measures = lexical_measures,
      plurals = SPOK_BNC %>%
        # create table of spoken bnc plurals
        filter(str_detect(c5, "NN2")) %>% 
        select(word) %>% 
        mutate(word = tolower(word) %>%
                 str_replace_all("\\+", "_") %>%
                 str_replace_all("@", "'")) %>%
        distinct(word),
      types = SPOK_BNC %>%
        # create spoken bnc types table and exclude plurals
        distinct(phon) %>%
        inner_join(., lexical_measures, by = c("phon")) # it will exclude plurals
    )
  }) ; rm(measures, SPOK_BNC)

corpora$types <- corpora$all %>%
  # create a corpora type table with lexical measures associated
  filter(speaker_role != "Target_Child") %>%
  mutate(gloss = gloss %>%
           str_split(" "),
         phon = phon %>%
           str_split(" ")) %>%
  group_by(corpus_id) %>%
  group_split() %>%
  lapply(function(sub_df) {
    tibble(corpus_id = sub_df$corpus_id[1],
           gloss = sub_df$gloss %>% unlist,
           phon = sub_df$phon %>% unlist)
  }) %>%
  bind_rows %>%
  # exclude spoken bnc plurals
  filter(!gloss %in% spok_bnc$plurals$word) %>%
  group_by(corpus_id) %>%
  distinct(phon) %>%
  select(corpus_id, phon) %>%
  # assign measures
  left_join(., spok_bnc$lexical_measures, by = "phon")

randomisation <- "empty" %>%
  (function(empty) {
    number_to_sample <- corpora$all %>%
      # minimum number of tokens in a corpus
      filter(speaker_role != "Target_Child") %>%
      group_by(corpus_id) %>%
      summarise(N_tok = gloss %>%
                  str_split(" ") %>%
                  unlist %>%
                  length) %>%
                  {.$N_tok} %>%
      min
    
    list(
      types =
        corpora$all %>%
          # create type tables of corpora after sampling by minimum number of tokens in a corpus
          filter(speaker_role != "Target_Child") %>%
          mutate(gloss = gloss %>%
                   str_split(" "),
                 phon = phon %>%
                   str_split(" ")) %>%
          (function(DF) {
            DF %>%
              group_by(corpus_id) %>%
              group_split %>%
              lapply(function(corpus_df) {
                
                
                random_samples <- list()
                
                set.seed(48593)
                for (i in 1:10) {
                  random_samples[[i]] <- tibble(corpus_id = corpus_df$corpus_id[1],
                                                random_sample = i,
                                                word = corpus_df$gloss %>% unlist,
                                                phon = corpus_df$phon %>% unlist) %>%
                    sample_n(number_to_sample) %>%
                    filter(!word %in% spok_bnc$plurals$word) %>%
                    distinct(phon, .keep_all = TRUE) %>%
                    left_join(., spok_bnc$lexical_measures, by = "phon") %>%
                    filter(!is.na(freq))
                }
                
                random_samples %>%
                  bind_rows
              })
          }) %>%
          bind_rows %>%
          select(-word) %>%
        rbind(
          spok_bnc$all %>%
            # create type tables of spoken bnc after sampling by minimum number of tokens in a child corpus
            (function(DF) {
              random_samples <- list()
              
              set.seed(48200)
              for (i in 1:10) {
                random_samples[[i]] <- DF %>%
                  sample_n(number_to_sample) %>%
                  distinct(phon) %>%
                  inner_join(., spok_bnc$lexical_measures, by = c("phon")) %>%
                  mutate(random_sample = i)
              }
              
              random_samples %>%
                bind_rows
            }) %>%
            mutate(corpus_id = "Spoken BNC") %>%
            select(corpus_id, random_sample, phon, everything())
      )
    )
  }) 

# performance - tokens ------------------------------------------------------
segmentation$output <- "empty" %>%
  # import output algorithms
  (function(empty) {
    import_output <- function(name_file, name_input) {
      read_delim(name_file, 
                 col_names = c("utt_id", "segmented_utt", "model_utt", 
                               "word_tp", "word_fp", "word_fn", "word_acc", "word_rec"),
                 delim = " ") %>%
        (function(DF) {
          name_input %>%
            cbind(
              DF %>%
                select(-utt_id, -segmented_utt)
            )
        })
    }
    
    path <- "/Users/francesco/Dropbox/phd/Incrementality/Segmentation/output/performance/"
    
    list(puddle_basic = import_output(paste(path, "PUDDLE_basic_alignment.txt", sep = ""),
                                      segmentation$input$basic),
         puddle_mlu = import_output(paste(path, "PUDDLE_mlu_alignment.txt", sep = ""),
                                    segmentation$input$mlu),
         ftp_basic = import_output(paste(path, "FTP_basic_alignment.txt", sep = ""),
                                      segmentation$input$basic),
         ftp_mlu = import_output(paste(path, "FTP_mlu_alignment.txt", sep = ""),
                                    segmentation$input$mlu)
    )
  })

segmentation$output %<>%
  # add a colum with logical values 
  # to indicate if a segment is a legal multiword sequence
  # takes ~30min
  lapply(function(DF) {
    DF %>%
      mutate(is_multiword_tp = sapply(1:nrow(.), function(i) {
        model_utt[i] %>%
          str_remove("\\|$") %>%
          str_split("\\|") %>%
          unlist %>%
          fmatch(
            adj_poss(phon_converted[i] %>%
                       str_remove("\\|$") %>%
                       str_split("\\|") %>%
                       unlist) %>%
              sapply(paste, collapse = "")) %>%
              {!is.na(.)}
      }))
  })

segmentation$output %<>%
  # logical columns of multiword false alarms: 
  # all undersegmented words or multi-word sequences
  lapply(function(DF) {
    DF %>%
      mutate(is_multiword_fp = sapply(1:nrow(.), function(i) {
        model_utt[i] %>%
          str_remove("\\|$") %>%
          str_split("\\|") %>%
          unlist %>%
          {.[!unlist(is_multiword_tp[i])]}
      })) %>%
      mutate(is_multiword_fp = sapply(1:nrow(.), function(i) {
        !str_detect(phon_converted[i], is_multiword_fp[[i]])
      }))
  })

segmentation$output %<>%
  # count number of legal multiword sequences in each utterance
  lapply(function(DF) {
    DF %>%
      mutate(multiword_tp = is_multiword_tp %>%
               sapply(sum))
  })

segmentation$output %<>%
  # 1) count multiword false alarms
  # 2) accuracy legal multiword learnt for each utterance
  lapply(function(DF) {
    DF %>%
      mutate(multiword_fp = is_multiword_fp %>%
               sapply(sum),
             multiword_acc = multiword_tp / (multiword_tp + multiword_fp))
  })

segmentation$input$unpacked <- "empty" %>%
  # unpack utterances: create a column with a word for each observation
  # takes ~10 min
  (function(empty) {
    unpack_utterances <- function(df) {
      df %>%
        mutate(phon = phon %>%
                 str_split(" "),
               part_of_speech = part_of_speech %>%
                 str_split(" "),
               row_id = 1:n(),
               MLU = sapply(phon, length)) %>% 
        group_by(target_child_age_year, row_id) %>%
        group_split() %>%
        lapply(function(sub_df) {
          tibble(row_id = sub_df$row_id[1],
                 target_child_age_year = sub_df$target_child_age_year[1],
                 MLU = sub_df$MLU[1],
                 phon = sub_df$phon %>% unlist)
        }) %>%
        bind_rows
    }
    
    list(puddle_basic = unpack_utterances(
      segmentation$output$puddle_basic
    ),
    puddle_mlu = unpack_utterances(
      segmentation$output$puddle_mlu
    ),
    ftp_basic = unpack_utterances(
      segmentation$output$ftp_basic
    ),
    ftp_mlu = unpack_utterances(
      segmentation$output$ftp_mlu
    ))
  })

segmentation$output %<>%
  # create consecutive stages every 1k utterances
  (function(LIST) {
    LIST$puddle_basic %<>%
      mutate(stage = c(rep(1:217, each = 1000), rep(217, 96))) %>%
      {.[1:217000, ]}
    
    LIST$puddle_mlu %<>%
      mutate(stage = c(rep(1:216, each = 1000), rep(216, 356))) %>%
      {.[1:216000, ]}
    
    LIST$ftp_basic %<>%
      mutate(stage = c(rep(1:217, each = 1000), rep(217, 96))) %>%
      {.[1:217000, ]}
    
    LIST$ftp_mlu %<>%
      mutate(stage = c(rep(1:216, each = 1000), rep(216, 356))) %>%
      {.[1:216000, ]}
    
    return(LIST)
  })

# performance - types  ------------------------------------------------
spok_bnc$gram <- spok_bnc$all %>%
  # most frequent grammatical category for each phonemic word
  group_by(phon) %>%
  count(pos) %>%
  group_split %>%
  lapply(function(sub_df) {
    sub_df %>%
      arrange(desc(n)) %>%
      {.[1, ]}
  }) %>%
  bind_rows %>%
  select(-n)

segmentation$output$unpacked_chunks <- "empty" %>%
  (function(empty) {
    chunk_extraction <- function(DF) {
       DF %>%
        # extract puddle_basic chunks segmented
        mutate(row_id = 1:n()) %>%
        group_by(row_id) %>%
        group_split %>%
        lapply(function(sub_df) {
          tibble(target_child_age_year = sub_df$target_child_age_year[1],
                 utt_num = sub_df$row_id[1],
                 utt_type = sub_df$utt_type[1],
                 stage = sub_df$stage[1],
                 phon_converted = sub_df$phon_converted[1],
                 model_chunk = sub_df$model_utt[1] %>%
                   str_remove("\\|$") %>%
                   str_split("\\|") %>%
                   unlist,
                 is_multiword = sub_df$is_multiword_tp %>% unlist)
        }) %>%
        bind_rows %>%
        # put the segmentation within the multiword chunks
        mutate(phon_converted = phon_converted %>%
                 str_remove("\\|$")) %>%
        mutate(model_chunk_converted = sapply(1:n(), function(i) {
          phon_converted[i] %>% 
            str_split("\\|") %>%
            unlist %>%
            adj_poss %>%
            (function(x) {
              c(
                sapply(x, paste, collapse = "|"),
                sapply(x, paste, collapse = "")
              )
            }) %>%
            matrix(ncol = 2) %>%
            {.[, 1][fmatch(.[, 2], model_chunk[i], nomatch = FALSE) %>%
                      as.logical]} %>%
                      {.[1]}
        })) %>%
        # create a variable with phonetic transcription of chunk if it is a word
        rowwise %>%
        mutate(is_word = str_detect(phon_converted, 
                                    paste("\\|", model_chunk, "\\|", sep = "") %>%
                                      paste(paste("^", model_chunk, "$", sep = ""), sep = "|") %>%
                                      paste(paste("\\|", model_chunk, "$", sep = ""), sep = "|") %>%
                                      paste(paste("^", model_chunk, "\\|", sep = ""), sep = "|"))) %>% 
        ungroup %>%
        rowwise %>%
        mutate(is_oversegmented = str_detect(phon_converted, 
                                            paste("\\|", model_chunk, "[a-zA-Z]+\\|", sep = "") %>%
                                              paste(paste("\\|[a-zA-Z]+", model_chunk, "\\|", sep = ""), sep = "|") %>%
                                              paste(paste("\\|[a-zA-Z]+", model_chunk, "[a-zA-Z]+\\|", sep = ""), sep = "|") %>%
                                              paste(paste("^", model_chunk, "[a-zA-Z]+\\|", sep = ""), sep = "|") %>%
                                              paste(paste("^[a-zA-Z]+", model_chunk, "\\|", sep = ""), sep = "|") %>%
                                              paste(paste("^[a-zA-Z]+", model_chunk, "[a-zA-Z]+\\|", sep = ""), sep = "|") %>%
                                              paste(paste("\\|[a-zA-Z]+", model_chunk, "$", sep = ""), sep = "|") %>%
                                              paste(paste("\\|", model_chunk, "[a-zA-Z]+$", sep = ""), sep = "|") %>%
                                              paste(paste("\\|[a-zA-Z]+", model_chunk, "[a-zA-Z]+$", sep = ""), sep = "|") %>%
                                              paste(paste("^", model_chunk, "[a-zA-Z]+$", sep = ""), sep = "|") %>%
                                              paste(paste("^[a-zA-Z]+", model_chunk, "$", sep = ""), sep = "|") %>%
                                              paste(paste("^[a-zA-Z]+", model_chunk, "[a-zA-Z]+$", sep = ""), sep = "|"))) %>% 
        ungroup %>%
        # logical: undersegmented chunk
        mutate(is_undersegmented = ifelse(is_multiword == FALSE & is_word == FALSE & is_oversegmented == FALSE, 
                                          TRUE,
                                          FALSE))
    } 
    
    segmentation$output %>%
      lapply(chunk_extraction)
  })

segmentation$output$unpacked_chunks$types <- "empty" %>%
  (function(empty) {
    chunk_types <- function(DF) {
      DF %>%
        # calculate proportion of chunks by stage over number of chunks learned
        mutate(rowid = 1:n()) %>%
        group_by(model_chunk, is_multiword, is_word, is_undersegmented, is_oversegmented) %>%
        filter(stage == min(stage)) %>%
        filter(rowid == min(rowid)) %>%
        ungroup %>% 
        filter(!(is_word == TRUE & is_oversegmented == TRUE)) %>%
        group_by(target_child_age_year, stage) %>%
        summarise(multiword_prop = sum(is_multiword),
                  word_prop = sum(is_word),
                  undersegmented_prop = sum(is_undersegmented),
                  oversegmented_prop = sum(is_oversegmented),
                  n = n()) %>%
        ungroup %>%
        mutate(multiword_prop = cumsum(multiword_prop),
               word_prop = cumsum(word_prop),
               undersegmented_prop = cumsum(undersegmented_prop),
               oversegmented_prop = cumsum(oversegmented_prop),
               n = cumsum(n)) %>%
        rowwise %>%
        mutate(multiword_prop = multiword_prop / n,
               word_prop = word_prop / n,
               undersegmented_prop = undersegmented_prop / n,
               oversegmented_prop = oversegmented_prop / n) %>%
        ungroup
    }
    
    segmentation$output$unpacked_chunks %>%
      lapply(chunk_types)
  })

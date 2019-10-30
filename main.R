# Francesco Cabiddu, CabidduF@cardiff.ac.uk

# load libraries ----------------------------------------------------------
lib <- c(
  "magrittr", "tidyverse",
  "stringr",
  "beepr", "fastmatch"#,
  #"mailR", "stringdist", "data.table"
)
lapply(lib, require, character.only = TRUE)
rm(lib)

# homemade funs -----------------------------------------------------------
# import and fist cleaning of text
first_wash <- function(DF, text_set) {
  DF %>%
    lapply(function(file) {
      print(file)
      
      output <- read_lines(paste(path_to_corpora, 
                                 file, 
                                 sep = "/")) %>%
        str_replace_all("\t", " ") %>%
        enframe(name = NULL) %>%
        (function(x) {
          # fix problem with cut utterances
          for (i in seq_along(x$value)) {
            if (str_detect(x$value[i], "^ ")) {
              if (str_detect(x$value[i-1], "^[^ ]{1}")) {
                x$value[i-1] <- paste(x$value[i-1], x$value[i], sep = "")
              } else if (str_detect(x$value[i-2], "^[^ ]{1}")) {
                x$value[i-2] <- paste(x$value[i-2], x$value[i], sep = "")
              } else if (str_detect(x$value[i-3], "^[^ ]{1}")) {
                x$value[i-3] <- paste(x$value[i-3], x$value[i], sep = "")
              } else if (str_detect(x$value[i-4], "^[^ ]{1}")) {
                x$value[i-4] <- paste(x$value[i-4], x$value[i], sep = "")
              } else if (str_detect(x$value[i-5], "^[^ ]{1}")) {
                x$value[i-5] <- paste(x$value[i-5], x$value[i], sep = "")
              } else if (str_detect(x$value[i-6], "^[^ ]{1}")) {
                x$value[i-6] <- paste(x$value[i-6], x$value[i], sep = "")
              } else if (str_detect(x$value[i-7], "^[^ ]{1}")) {
                x$value[i-7] <- paste(x$value[i-7], x$value[i], sep = "")
              } else if (str_detect(x$value[i-8], "^[^ ]{1}")) {
                x$value[i-8] <- paste(x$value[i-8], x$value[i], sep = "")
              }
            }
          }
          
          x
        }) %>%
        filter(!str_detect(value, "^ |^@")) 
      
      if (text_set %in% c("wells")) {
        output %<>%
          mutate(value = value %>% 
                   str_remove_all("[(]{1}[0-9]+[.]{1}[)]{1} "))
      }
      
      if (text_set %in% c("thomas", "tommerdahl")) {
        output %<>%
          mutate(value = value %>% 
                   str_remove_all(" \025[0-9]+_[0-9]+\025$"))
      }
      
      if (text_set %in% c("general", "thomas", "tommerdahl", "wells")) {
        output %>%
        {.$value %>% 
            paste(collapse = "") %>%
            str_extract_all("[*]{1}[A-Z]{3}:[^%*]+%mor[^%]+%gra[^*]+")} %>%
            {.[[1]]} %>%
          enframe(name = NULL) %>%
          separate(value, c("id_utt", "mor_gra"), sep = "%mor: ") %>% 
          separate(mor_gra, c("mor", "gra"), sep = "%gra: ") %>%
          mutate(id = str_match(id_utt, "^[*]{1}([A-Z]{3}):")[,2],
                 id_utt = str_remove(id_utt, "^[*]{1}[A-Z]{3}: ")) %>%
          rename(utt = id_utt) %>%
          select(id, everything(), -gra) %>%
          mutate_at(.funs = list(~str_split(., " ")), .vars = vars(utt:mor)) %>%
          mutate_at(.funs = list(length = ~sapply(., length)), .vars = vars(utt:mor)) %>%
          filter(utt_length == mor_length) %>%
          select(-c(utt_length, mor_length)) %>%
          mutate(utt_id = 1:n(),
                 corpus_id = file) %>% 
          group_by(utt_id) %>%
          group_split() %>%
          lapply(function(sub_df) {
            tibble(corpus_id = sub_df$corpus_id,
                   utt_id = sub_df$utt_id,
                   id = sub_df$id,
                   utt = sub_df$utt %>% unlist,
                   mor = sub_df$mor %>% unlist)
          }) %>%
          bind_rows %>%
          # delete special charachters and punctuation
          mutate(utt = tolower(utt) %>%
                   str_remove("@[a-z]+$") %>%
                   str_remove_all("[^a-z'+_]+") %>%
                   str_replace_all("[+]+", "_")) %>% 
          filter(str_detect(utt, "[a-z]+")) %>%
          group_by(utt_id) %>%
          group_split %>%
          lapply(function(sub_df) {
            tibble(corpus_id = sub_df$corpus_id[1],
                   utt_id = sub_df$utt_id[1],
                   id = sub_df$id[1],
                   utt = sub_df$utt %>% paste(collapse = " "),
                   mor = sub_df$mor %>% paste(collapse = " "))
          }) %>%
          bind_rows %>%
          mutate(utt = utt %>%
                   str_replace_all("' | '", " "))
      } else if (text_set == "forrester") {
        output %>%
          filter(!str_detect(value, "^%gpx|^%com")) %>%
          filter(!str_detect(value, "xx")) %>%
          mutate(value = str_replace(value, "\\: ", "_9_9_9_")) %>% 
          separate(value, c("id", "utt"), sep = "_9_9_9_") %>% 
          mutate(id = id %>% str_remove("^\\*")) %>%
          mutate(utt = tolower(utt)) %>%
          mutate(utt = utt %>%
                   str_remove_all("\\[% [a-zA-Z]+\\] | \\[% [^\\]]+\\]") %>%
                   str_remove_all("&=[a-zA-Z]+ | &=[a-zA-Z]+")) %>%
          mutate(utt = utt %>%
                   str_remove_all("[^ a-zA-Z'+_]+")) %>%
          (function(z) {
            z %>%
              mutate(corpus_id = file,
                     utt_id = 1:nrow(z),
                     mor = NA)
          }) %>% 
          select(corpus_id, utt_id, id, utt, mor) %>%
          mutate(utt = utt %>%
                   str_replace_all("' | '", " "))
      }
    })
}

extract_age <- function(file) {
  print(file)
  
  # the output is a dataframe with years, months, and days
  output <- read_lines(paste(path_to_corpora, 
                             file, 
                             sep = "/")) %>%
    str_replace_all("\t", " ") %>%
    enframe(name = NULL) %>%
    filter(str_detect(value, "^@ID")) %>%
    unlist %>%
    str_extract("[|]{1}([0-9]+;[0-9+]+[.]{1}[0-9]+)[|]{1}|[|]{1}([0-9]+;[0-9+]+[.]{1})[|]{1}|[|]{1}([0-9]+;)[|]{1}") %>%
    na.omit %>%
    str_split(";|[.]{1}")
  
  if (is_empty(output)) {
    output <- rep(NA, 3)
  } else {
    output %<>%
    {.[[1]]}
  }
  
  output %>%
    str_remove_all("[|]{1}") %>%
    as.numeric %>%
    (function(x) {
      length(x) <- 3
      
      x %<>%
        as.matrix %>%
        t
      
      colnames(x) <- c("year", "month", "day")
      x
    }) %>%
    as_tibble %>%
    mutate(corpus_id = file) %>%
    select(corpus_id, everything())
}

convert_phonemes <- function(df) {
  # create a column of converted phonemes to single characters for each word
  df %>%
    mutate(phon_converted = phon %>%
             str_split("_") %>%
             sapply(function(i) {
               paste0(phonemes_converted[i], collapse="")
             }))
}

# import corpora ----------------------------------------------------------
path_to_corpora <- 
  "/Users/francesco/Dropbox/phd/Incrementality/corpora/CDS/TXT_format"

filelist <- list.files(
  path = path_to_corpora,
  pattern = "txt$",
  all.files = TRUE,
  recursive = TRUE)

# leave out forrester because has a different pattern to match
filelist_not_forr_tho_to_we <- filelist[filelist %>% str_detect("Forrester|Thomas|Tommerdahl|Wells") %>% {!.}]
filelist_forr <- filelist[filelist %>% str_detect("Forrester")]
filelist_tho <- filelist[filelist %>% str_detect("Thomas")]
filelist_to <- filelist[filelist %>% str_detect("Tommerdahl")]
filelist_we <- filelist[filelist %>% str_detect("Wells")]

time1 <- Sys.time()

corpora_not_forr_tho_to_we <- filelist_not_forr_tho_to_we %>%
  first_wash(text_set = "general")
corpora_forr <- filelist_forr %>%
  first_wash(text_set = "forrester")
corpora_tho <- filelist_tho %>%
  first_wash(text_set = "thomas")
corpora_to <- filelist_to %>%
  first_wash(text_set = "tommerdahl")
corpora_we <- filelist_we %>%
  first_wash(text_set = "wells")

corpora <- corpora_not_forr_tho_to_we %>%
  bind_rows %>%
  rbind(
    corpora_forr %>%
      bind_rows
  ) %>%
  rbind(
    corpora_tho %>%
      bind_rows
  ) %>%
  rbind(
    corpora_to %>%
      bind_rows
  ) %>%
  rbind(
    corpora_we %>%
      bind_rows
  )

corpora %<>%
  mutate(utt = utt %>%
           str_replace_all("' | '|'$|^'", " ") %>%
           str_remove_all("[ ]{2,}") %>%
           str_remove_all("[_]{2,}") %>%
           str_remove_all("_$|_ $|^_|^_ | _ ") %>%
           str_remove_all("[+]{1} | [+]{1}|^\\+|\\+$") %>%
           str_remove_all("^ | $") %>%
           str_remove_all("[ ]{2,}"))
           
# compare corpora lexical measures ----------------------------------------
corpora %<>%
  mutate(utt = str_split(utt, " "),
         mor = str_split(mor, " "))

corpora %<>%
  mutate(row_id = 1:length(id)) %>%
  group_by(row_id) %>%
  group_split %>%
  lapply(function(sub_df) {
    tibble(corpus_id = sub_df$corpus_id,
           utt_id = sub_df$utt_id,
           id = sub_df$id,
           utt = sub_df$utt %>% unlist,
           mor = sub_df$mor %>% unlist)
  }) %>% bind_rows

corpora_typ <- corpora %>%
  distinct(utt) %>%
  arrange(utt) %>%
  (function(x) {
    list(single_word = filter(x, !str_detect(utt, "[_+]+")),
         multi_word = filter(x, str_detect(utt, "[_+]+")))
  })

time2 <- Sys.time()
time2 - time1

# import dictonary
load("~/Dropbox/phd/Incrementality/analysis_R/CMU_DICT.RData")

CMU_DICT %<>%
  (function(list_dfs) {
    additional_compounds <- list_dfs$words %>%
      filter(str_detect(word, "_"))
    
    list_dfs$words %<>%
      filter(!str_detect(word, "_")) %>%
      mutate(word = word %>%
               tolower %>%
               str_replace_all("@", "'"))
    
    list_dfs$compounds %<>%
      rbind(
        additional_compounds
      ) %>%
      mutate(word = word %>%
               tolower %>%
               str_replace_all("[+]{1}", "_") %>%
               str_replace_all("@", "'"))
    
    list_dfs
  })

corpora_typ$single_word %<>%
  left_join(., CMU_DICT$words, by = c("utt" = "word"))

corpora_typ$multi_word %<>%
  left_join(., CMU_DICT$compounds, by = c("utt" = "word"))

corpora_typ$multi_word %<>%
  filter(is.na(phon)) %>%
  mutate(phon = utt %>%
           str_split("_") %>%
           sapply(function(x) {
             CMU_DICT$words$phon[fmatch(x,
                                        CMU_DICT$words$word)] %>%
               paste(collapse = "_")
           })) %>%
  mutate(phon = ifelse(str_detect(phon, "NA"), NA, phon))

corpora_typ %<>%
  bind_rows %>%
  arrange(utt)
  
  # convert single and multi-word sequences
corpora %<>%
  left_join(., corpora_typ, by = "utt")

  # get rid of utterances where conversion is missing for a word
corpora %<>%
  group_by(corpus_id, utt_id) %>%
  filter(all(!is.na(phon))) %>%
  ungroup

# assign measures
load("~/Dropbox/phd/Incrementality/analysis_R/spok_bnc_lexical.RData")
load("~/Dropbox/phd/Incrementality/analysis_R/spok_bnc_tokens.RData")

spok_bnc_typ_measures <- measures %>%
  select(phon, contains("spok_bnc")) %>%
  mutate(phonemic_len = phon %>%
           str_split("_") %>%
           sapply(length)) %>%
  (function(x) {
    colnames(x) <- str_remove(colnames(x), "spok_bnc_")
    x
  }) %>%
  filter(!is.na(freq)) ; rm(measures)
  
corpora %<>%
    left_join(., spok_bnc_typ_measures, by = "phon")

spok_bnc_typ <- SPOK_BNC %>%
  filter(!str_detect(c5, "NN2")) %>%
  left_join(., spok_bnc_typ_measures, by = "phon") %>%
  distinct(phon, .keep_all = TRUE)

# prepare dfs for segmentation algorithms  -------------------------------------------------------------------
age_table <- filelist_not_forr_tho_to_we %>%
  c(filelist_forr, filelist_tho, filelist_to, filelist_we) %>%
  lapply(extract_age) %>%
  bind_rows %>%
  arrange(year, month, day)

# exclude files that don't have any age assigned
age_table %<>%
  na.omit

corpora %<>%
  left_join(., age_table, by = "corpus_id")

number_to_sample <- corpora %>%
  filter(!id %in% c("CHI")) %>%
  arrange(year, month) %>%
  filter(!is.na(year) & !is.na(month) & !is.na(day)) %>%
  group_by(year, corpus_id) %>%
  summarise(N_utterance = length(unique(utt_id))) %>%
  group_by(year) %>%
  summarise(N_utterance = sum(N_utterance)) %>%
  ungroup %>%
  filter(year %in% 1:4) %>%
  filter(N_utterance == min(N_utterance)) %>%
  {.$N_utterance[1]}

# convert phonemes into single characters 
phonemes_converted <- c(letters, LETTERS)[1:39] %>%
  (function(x) {
    names(x) <- corpora$phon %>%
      str_split("_") %>%
      unlist %>%
      unique %>%
      sort
    
    x
  })

corpora_input <- corpora %>%
  filter(!id %in% c("CHI")) %>%
  filter(year %in% 1:4) %>%
  convert_phonemes %>%
  group_by(corpus_id, utt_id) %>%
  group_split %>%
  lapply(function(sub_df) {
    tibble(corpus_id = sub_df$corpus_id[1],
           utt_id = sub_df$utt_id[1],
           id = sub_df$id[1],
           year = sub_df$year[1],
           month = sub_df$month[1],
           day = sub_df$day[1],
           utt_orth = paste(sub_df$utt, collapse = " "),
           mor = paste(sub_df$mor, collapse = " "),
           utt_phon = paste(sub_df$phon, collapse = " "),
           utt_phon_converted = paste(sub_df$phon_converted, collapse = "|"))
  }) %>%
    bind_rows

set.seed(12776)
corpora_input %<>%
  group_by(year) %>%
  sample_n(number_to_sample) %>%
  ungroup

corpora_input %<>%
  arrange(year, month)
  
#vowels <- c("AA", "AE", "AH", "AO", "AW", "AY", "EH", "ER", 
#            "EY", "IH", "IY", "OW", "OY", "UH", "UW") 


#write_lines(corpora_input$utt_phon_converted %>%
#              str_replace("$", "|"), "input_PUDDLE_segmented.txt")
#write_lines(corpora_input$utt_phon_converted %>%
#              str_remove_all("[|]"), "input_PUDDLE_unsegmented.txt")


# performance PUDDLE ------------------------------------------------------
puddle <- read_tsv("results_temp.txt", col_names = FALSE) %>%
  separate(X1, c("utt_id", "segmented_utt", "model_utt", 
                 "tp", "fp", "fn", "acc", "rec"), 
           sep =" ") %>%
  mutate(stage = rep(1:27, each = 1000),
         year = corpora_input$year[27000])



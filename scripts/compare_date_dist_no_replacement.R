

library(reshape2)
library(tidyverse)
library(dplyr)

# data
genetic_file <- "" # symmetrical distance matrix with matching file names
date_file <- "" # two columns (1- sample name 2 - date)
output_file <- ""

genetic <- read.delim(sep="\t", header = T, file = genetic_file)
date <- read.delim(sep="\t", header = T, file = date_file)
output <- output_file

# adjust dates entries to remove unwanted characters
date2 <- date
date2$isolate  <- gsub(pattern = "^X", replacement = "", x = date$isolate, perl=T)
date2$isolate <- gsub(pattern = "\\|", replacement = "\\.", x = date2$isolate, perl=T)
date2$isolate <- gsub(pattern = "\\-", replacement = "\\.", x = date2$isolate, perl=T)
date2$isolate <- gsub(pattern = "\\/", replacement = "\\.", x = date2$isolate, perl=T)

# make dates into ditsance matrix
date_dm <- as.data.frame(abs(outer(date$year, date$year, '-')))
colnames(date_dm) <- date$isolate
date_dm$X <- date$isolate

# melt
g_pw <- melt(genetic, id = c("X"))
d_pw <- melt(date_dm, id = c("X"))

# remove NAs
#g_pw <- g_pw[!(is.na(g_pw[,3])),] 
#d_pw <- g_pw[!(is.na(d_pw[,3])),] 

# set column names
colnames(g_pw) <- c("A", "B", "genetic")
colnames(d_pw) <- c("A", "B", "date")

# remove unwanted charactersfrom second column 
g_pw[,2] <- gsub(pattern = "^X", replacement = "", x = g_pw[,2], perl=T)
d_pw[,2] <- gsub(pattern = "^X", replacement = "", x = d_pw[,2], perl=T)
#g_pw[,2] <- gsub(pattern = "\\.", replacement = "|", x = g_pw[,1], perl=T)
#d_pw[,2] <- gsub(pattern = "\\.", replacement = "|", x = d_pw[,1], perl=T)
g_pw[,1] <- gsub(pattern = "\\|", replacement = "\\.", x = g_pw[,1], perl=T)
d_pw[,1] <- gsub(pattern = "\\|", replacement = "\\.", x = d_pw[,1], perl=T)
d_pw[,2] <- gsub(pattern = "\\|", replacement = "\\.", x = d_pw[,2], perl=T)
g_pw[,1] <- gsub(pattern = "\\-", replacement = "\\.", x = g_pw[,1], perl=T)
d_pw[,1] <- gsub(pattern = "\\-", replacement = "\\.", x = d_pw[,1], perl=T)
d_pw[,2] <- gsub(pattern = "\\-", replacement = "\\.", x = d_pw[,2], perl=T)
g_pw[,1] <- gsub(pattern = "\\/", replacement = "\\.", x = g_pw[,1], perl=T)
d_pw[,1] <- gsub(pattern = "\\/", replacement = "\\.", x = d_pw[,1], perl=T)
d_pw[,2] <- gsub(pattern = "\\/", replacement = "\\.", x = d_pw[,2], perl=T)

# check unique
setdiff(d_pw$A, d_pw$B)
u1 <- unique(c(d_pw$A, d_pw$B))
u2 <- unique(c(g_pw$A, g_pw$B))
ui <- unique(u1, u2)

unique_values  <- sort(unique(c(d_pw$date)))
for (uv in unique_values ){ 
  
  # uv <- 1 # test value
  pw <- left_join(g_pw, d_pw) %>% filter(date>=uv)
  
  # unique samples
  u <- unique(c(pw$A, pw$B))
  
  if ( length(u) > 1 ){
    
    # make output
    result_list <- list()
    result_count <- 0
    
    # check for replacement in B
    replace_anc <- data.frame( id = u, inc = rep(0, length(u))) 
    replace_mod <- data.frame( id = u, inc = rep(0, length(u))) 
    
    tc <- 0
    
    for (i in 1:length(u)){
      
      tc=tc+1
      
      #i<-79
      test <- u[i]
      test_df <- as.data.frame(pw %>% filter(A == test))
      test_df$yearA <- NULL
      test_df$yearB <- NULL
      
      # polarise comparisons
      for (j in 1:length(test_df[,1])){
        
        #j<- 1
        val1 <- date2$year[date2$isolate == as.character(test_df$A[j])]
        val2 <- date2$year[date2$isolate == as.character(test_df$B[j])]
        
        test_df$yearA[j] <- val1
        test_df$yearB[j] <- val2
        
        # make column A ancestral
        if ( val1 > val2 ){
          temp <- test_df$A[j]
          test_df$A[j] <- test_df$B[j]
          test_df$B[j] <- temp
          
          test_df$yearA[j] <- val2
          test_df$yearB[j] <- val1
        }
          
      }
      
      # filter table for samples to remove
      include_anc <- replace_anc$id[replace_anc$inc==0]
      include_mod <- replace_mod$id[replace_mod$inc==0]
      test_df <- test_df %>% filter(A %in% include_anc) %>% filter(B %in% include_mod)
      
      # remove same same 
      test_df <- test_df[!(test_df$A==test_df$B),]
      
      # select lowest genetic distance 
      test_df <- test_df %>% arrange(genetic)
      
      if (length(test_df[,1]) > 0){
        
        result_count <- result_count + 1
        
        # select fist value
        test_include <- test_df[1,]
        
        tc_out = sprintf("%s\t%s", result_count, tc)
        # print(tc_out) # test
         
          # store result 
          names(test_include) <- colnames(test_include)
          result_list[[result_count]] <- test_include
        
          # exclude isolates from the relevant comparisons
          replace_anc$inc[as.character(replace_anc$id) == as.character(test_include$A[1]) ] <- 1
          replace_mod$inc[as.character(replace_mod$id) == as.character(test_include$B[1]) ] <- 1
          
      }else{
        print("no pair")
      }
    
    }
    
    # make df
    result <- result_list %>% 
      map_df(as_tibble)
    #result <- data.frame(A = rep( 0, length(u) ), B = rep( 0, length(u) ), genetic = rep( 0, length(u) ) , date = rep( 0, length(u) ) )
    
    # write to file
    write.table(result, sprintf("%s.%i.tsv", output, uv), sep = "\t", quote = F, row.names=F, col.names = T)
    
  }

}


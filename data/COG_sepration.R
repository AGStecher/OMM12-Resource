# Load libraries
library(dplyr)
library(tidyr)
library(data.table)

# 1. Import your CSV file

df <- fread("C:/Users/ge24biy/Documents/OMM12_website/data/acutalibacter_muris_kb18_cog_annotations.txt", sep="\t", header= TRUE, stringsAsFactors = FALSE)

# 2. Split multi-letter COG categories into separate rows
df2 <- df %>%
  separate_rows(COG_category, sep = "") %>%   # split into individual characters
  filter(COG_category != "")                  # remove empty values created by splitting

# 3. View the result
print(df2)

# 4. (Optional) Save the output to a new CSV
write.table(df2, "C:/Users/ge24biy/Documents/OMM12_website/data/acutalibacter_muris_kb18_cog.tsv",sep="\t",quote = FALSE, row.names = FALSE)

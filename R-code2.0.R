# Read the expression matrix data from a CSV file without a header
dataor <- read.csv("D:\\Users\\Lenovo\\Desktop\\expression_matrix.csv", header = FALSE)

# Load the 'biomaRt' library
library(biomaRt)

# List available Ensembl archives
listEnsemblArchives()

# List available BioMart datasets on the specified host
listMarts(host = 'jan2020.archive.ensembl.org')

# Use the 'biomaRt' package to connect to the Ensembl database
ensembl99 <- useMart(host = 'jan2020.archive.ensembl.org',
                     biomart = 'ENSEMBL_MART_ENSEMBL',
                     dataset = 'hsapiens_gene_ensembl')

# Read gene data from a CSV file
gene <- read.csv("D:\\Users\\Lenovo\\Desktop\\mart_export.txt")

# Read data from another CSV file with headers
data <- read.csv("D:\\Users\\Lenovo\\Desktop\\change_expression_matrix.csv", header = TRUE)

# Read genes data from a CSV file
genes <- read.csv("D:\\Users\\Lenovo\\Desktop\\d2.csv")

# Perform a left join of 'genes' and 'data' based on the 'ensembl_gene_id' column
dnew2 <- left_join(genes, data, by = 'ensembl_gene_id')

# Duplicate 'dnew' variable, but the actual assignment is missing
d_m <- dnew

# Write 'd_m' data to a CSV file
write.csv(d_m, "D:\\Users\\Lenovo\\Desktop\\merge.csv")

# Merge 'data' and 'genes' data frames based on 'ensembl_gene_id' column
dnew3 <- merge(data, genes, by = 'ensembl_gene_id')

# Remove duplicate rows from 'dnew3' and store in 'd_mer'
d_mer <- distinct(dnew3)

# Read data from another CSV file
datan <- read.csv("D:\\Users\\Desktop\\MData.csv")

# Remove duplicate rows from 'datan' and store in 'data1'
data1 <- datan[!duplicated(datan), ]

# Define a function 'filter_rows' that filters rows based on certain conditions
filter_rows <- function(df) {
  num_cols <- ncol(df)
  min_nonzero <- num_cols * 0.8
  at_least_one_gt_1 <- rowSums(df > 1) >= 1
  at_least_80_percent_gt_0 <- rowSums(df > 0) >= min_nonzero
  df[at_least_one_gt_1 & at_least_80_percent_gt_0, ]
}

# Use the 'filter_rows' function to filter the 'data1' dataset
f <- filter_rows(data1)

# Read data from another CSV file
data4 <- read.csv("D:\\Users\\Desktop\\rows_metadata.csv")

# Merge 'data4' and 'f' data frames based on 'ensembl_gene_id' column
pinjie <- merge(data4, f, by = 'ensembl_gene_id')

# Select columns from 'pinjie' data frame
new <- pinjie[, c(4, 6:529)]

# Write 'new' data to a CSV file
write.csv(new, "D:\\USERS\\Desktop\\S.Data.csv")

# Transpose the 'new' data frame
dz <- t(new)

# Write the transposed data to a CSV file
write.csv(dz, "D:\\Users\\Desktop\\T.Data.csv")  # Open the file here and delete the first row

# Read data from another CSV file
dv <- read.csv("D:\\Users\\Desktop\\T.Data.csv")

# Read data from another CSV file
dy <- read.csv("D:\\USERS\\Desktop\\columns_metadata.csv")

# Combine 'dy' and 'dv' data frames column-wise
dx <- cbind(dy, dv)

# Write 'dx' data to a CSV file
write.csv(dx, "D:\\USERS\\Desktop\\S.Final.csv")

# Select specific rows from 'dx'
dw <- dx[31:524, ]

# Read data from another CSV file
data8 <- read.csv("D:\\USERS\\Desktop\\DD.Final.csv")

# Split 'data8' into different datasets based on the 'structure_acronym' column
split_datasets <- split(data8, data8$structure_acronym)

# Specify the export path
export_path <- "D:\\USERS\\Desktop\\"

# Create a named list to store datasets
named_datasets <- list()

# Process and export each dataset
for (name in names(split_datasets)) {
  processed_data <- split_datasets[[name]][, -c(2:9)]
  transposed_data <- t(processed_data)
  colnames(transposed_data) <- transposed_data[1, ]
  transposed_data <- transposed_data[-1, ]
  transposed_data <- cbind(buding2, transposed_data)  # Unclear variable 'buding2'
  file_path <- paste0(export_path, name, ".csv")
  write.csv(transposed_data, file = file_path, row.names = FALSE)
  named_datasets[[name]] <- transposed_data
}

# Read data from another CSV file
ENSG <- read.csv("D:\\Users\\Desktop\\Biostatistics\\rows_metadata.csv")

# Read data from another CSV file
dB_old <- read.csv("D:\\Users\\Desktop\\Biostatistics\\OLD_B_CLASS.csv")

# Merge 'dB_old' and 'ENSG' data frames based on 'gene' column
B_ENSG_old <- merge(dB_old, ENSG, by = 'gene')

# Extract a specific column from 'B_ENSG_old'
OLDB <- B_ENSG_old[, 5]

# Write 'OLDB' data to a CSV file
write.csv(OLDB, "D:\\USERS\\Desktop\\OLDB.csv")

# Read data from another CSV file
dD_old <- read.csv("D:\\Users\\Desktop\\Biostatistics\\OLD_D_CLASS.csv")

# Merge 'dD_old' and 'ENSG' data frames based on 'gene' column
D_ENSG_old <- merge(dD_old, ENSG, by = 'gene')

# Extract a specific column from 'D_ENSG_old'
OLDD <- D_ENSG_old[, 5]

# Write 'OLDD' data to a CSV file
write.csv(OLDD, "D:\\USERS\\Desktop\\OLDD.csv")

# Read data from another CSV file
dB_new <- read.csv("D:\\Users\\Desktop\\Biostatistics\\NEW_E_CLASS.csv")

# Merge 'dB_new' and 'ENSG' data frames based on 'gene' column
B_ENSG_new <- merge(dB_new, ENSG, by = 'gene')

# Extract a specific column from 'B_ENSG_new'
NEWB <- B_ENSG_new[, 6]

# Write 'NEWB' data to a CSV file
write.csv(NEWB, "D:\\USERS\\Desktop\\NEWEnum.csv")

# Read data from another CSV file
dB_new <- read.csv("D:\\Users\\Desktop\\Biostatistics\\ANEW.csv")

# Merge 'dB_new' and 'ENSG' data frames based on 'gene' column
B_ENSG_new <- merge(dB_new, ENSG, by = 'gene')

# Extract a specific column from 'B_ENSG_new'
NEWB <- B_ENSG_new[, 5]

# Write 'NEWB' data to a CSV file
write.csv(NEWB, "D:\\USERS\\Desktop\\NEWAnum.csv")

# Read data from another CSV file
dD_new1 <- read.csv("D:\\Users\\Desktop\\Biostatistics\\DNEW.csv")

# Merge 'dD_new1' and 'ENSG' data frames based on 'gene' column
D_ENSG_new <- merge(dD_new1, ENSG, by = 'gene')

# Extract a specific column from 'D_ENSG_new'
NEWD1 <- D_ENSG_new[, 5]

# Write 'NEWD1' data to a CSV file
write.csv(NEWD1, "D:\\USERS\\Desktop\\NEWDnum.csv")

# Read data from another CSV file
dCnew <- read.csv("D:\\Users\\Desktop\\Biostatistics\\CNEW.csv")

# Merge 'dCnew' and 'ENSG' data frames based on 'gene' column
C_ENSG_new <- merge(dCnew, ENSG, by = 'gene')

# Extract a specific column from 'C_ENSG_new'
NEWC <- C_ENSG_new[, 5]

# Write 'NEWC' data to a CSV file
write.csv(NEWC, "D:\\USERS\\Desktop\\NEWCnum.csv")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read data from a CSV file
data <- read.csv("D:/USERS/Desktop/OLDD.csv")

# Define colors
COLS <- c("BP" = "blue", "MF" = "green", "CC" = "red")

# Simplify category labels
data$class <- recode(data$class,
                     "Biological Process" = "BP",
                     "Molecular Function" = "MF",
                     "Cellular Component" = "CC")

# Calculate the number of genes based on Description and class
gene_counts <- data %>%
  group_by(Description, class) %>%
  summarise(GeneNumber = sum(num)) %>%
  arrange(class, desc(GeneNumber)) %>%
  group_by(class) %>%
  slice_head(n = 10)

# Create a horizontal bar plot, grouped by class, and set different font sizes and colors for various text elements. Use facet_wrap to group categories together.
bar_plot <- ggplot(data = gene_counts, aes(x = GeneNumber, y = Description, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = COLS) +
  theme_minimal() +
  ylab("Description") +
  xlab("Num of Genes") +
  labs(title = "") +
  theme(axis.text.y = element_text(face = "bold", color = "black", hjust = 1, family = "Times New Roman", size = 32),
        axis.text.x = element_text(size = 24, color = "black"),  # Font size and color for gene number scale
        axis.title = element_text(face = "bold", color = "black", family = "Times New Roman", size = 40),
        legend.title = element_blank(),
        plot.title = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 40)) +  # Legend text size
  facet_wrap(~class, ncol = 1, scales = "free_y")  # Group categories together

# Display the horizontal bar plot
print(bar_plot)

# Load necessary libraries
library(ggplot2)
library(extrafont)

# Read data from a CSV file with specific file encoding
result <- read.csv("D:\\Users\\Desktop\\PLOTT.csv", fileEncoding = "GBK")

# Install and import the Times New Roman font
font_import(pattern = "Times New Roman")
loadfonts()

# Define custom colors
custom_colors <- c('#ED0000E5', '#42B540E5', '#00468BE5')

# Create a basic bar plot
p <- ggplot(data = result, aes(x = Region, y = count)) +
  geom_bar(aes(fill = Subset), stat = "identity") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  scale_x_discrete(labels = 1:19) +
  labs(x = "Region", y = "Count", fill = "Subset") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 48, color = "black", family = "Times New Roman"),
    axis.text.y = element_text(size = 48, color = "black", family = "Times New Roman"),
    axis.title = element_text(size = 50, color = "black", family = "Times New Roman"),
    legend.text = element_text(size = 42, color = "black", family = "Times New Roman"),
    legend.title = element_text(size = 48, color = "black", family = "Times New Roman"),
    plot.title = element_text(size = 40, color = "black", family = "Times New Roman")
  ) +
  scale_fill_manual(values = custom_colors)

# Add y-axis labels
p + annotate("text", x = 8, y = 500, label = "Down-regulated Genes", size = 24, family = "Times New Roman") +
  annotate("text", x = 8, y = -100, label = "Up-regulated Genes", size = 24, family = "Times New Roman", vjust = 2)

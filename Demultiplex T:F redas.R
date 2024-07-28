# Load necessary libraries
library(ggplot2)
library(reshape2)

# Sample data based on your table (just a small subset for demonstration)
data <- data.frame(
  Barcode = factor(1:10),
  'FALSE' = c(231650, 44654, 51912, 45768, 34088, 36252, 59338, 40705, 28700, 45664),
  'TRUE' = c(412439, 103222, 59462, 68540, 75633, 80385, 103512, 69630, 55380, 115684)
)

# Melt the data for ggplot2
data_melted <- melt(data, id.vars = "Barcode", variable.name = "Assignment", value.name = "Count")

# Stacked Bar Plot
ggplot(data_melted, aes(x = Barcode, y = Count, fill = Assignment)) +
  geom_bar(stat = "identity") +
  labs(title = "Read Assignment per Barcode", x = "Barcode", y = "Read Count") +
  theme_minimal()

# Grouped Bar Plot
ggplot(data_melted, aes(x = Barcode, y = Count, fill = Assignment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Read Assignment per Barcode", x = "Barcode", y = "Read Count") +
  theme_minimal()

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Full data based on your table
data <- data.frame(
  Barcode = factor(1:96),
  `FALSE` = c(231650, 44654, 51912, 45768, 34088, 36252, 59338, 40705, 28700, 45664,
              49396, 75738, 35611, 24085, 33485, 34629, 32997, 42579, 49926, 29178,
              34050, 40476, 42717, 62662, 53108, 70700, 39882, 43001, 82326, 52512,
              22181, 64427, 67327, 39943, 226561, 47137, 52660, 62128, 82059, 45273,
              155490, 40545, 26193, 141153, 95783, 42218, 34560, 25834, 54950, 27613,
              38190, 31783, 44209, 41402, 33946, 22113, 23537, 34537, 27056, 14279,
              147431, 53196, 42266, 24076, 79087, 62521, 30939, 53008, 20896, 51667,
              23311, 35702, 38830, 58804, 26006, 28194, 26695, 27123, 54454, 20870,
              24016, 36165, 20249, 48867, 29719, 60333, 49565, 64673, 54106, 72502,
              31406, 33948, 17581, 40773, 21678, 20040),
  `TRUE` = c(412439, 103222, 59462, 68540, 75633, 80385, 103512, 69630, 55380, 115684,
             101780, 105626, 81781, 26235, 66922, 72206, 73579, 89742, 67773, 41814,
             60957, 101775, 82977, 51296, 167507, 231555, 105795, 116257, 161255, 99058,
             54646, 165915, 150581, 57576, 380975, 92473, 173814, 134687, 196771, 155012,
             287777, 94199, 62336, 211869, 264562, 116763, 65705, 82555, 78969, 62035,
             66967, 58427, 39064, 86816, 56889, 43592, 45127, 76911, 57021, 26613,
             902, 101941, 78404, 55338, 116825, 143998, 62199, 101058, 44211, 126127,
             43954, 63003, 78423, 124652, 64069, 63689, 48736, 68703, 54281, 43433,
             59098, 61411, 32378, 93229, 73612, 136886, 105046, 96805, 100881, 116273,
             54141, 49108, 33490, 90001, 42007, 31182)
)

# Melt the data for ggplot2
data_melted <- melt(data, id.vars = "Barcode", variable.name = "Assignment", value.name = "Count")

# Stacked Bar Plot
ggplot(data_melted, aes(x = Barcode, y = Count, fill = Assignment)) +
  geom_bar(stat = "identity") +
  labs(title = "Read Assignment per Barcode", x = "Barcode", y = "Read Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Grouped Bar Plot
ggplot(data_melted, aes(x = Barcode, y = Count, fill = Assignment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Read Assignment per Barcode", x = "Barcode", y = "Read Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

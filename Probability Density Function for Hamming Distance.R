

# Load necessary libraries
library(ggplot2)
library(viridis)

# Define a function to calculate the PDF of the Hamming distance
calculate_pdf <- function(lambda, k) {
  return((lambda^k * exp(-lambda)) / factorial(k))
}

# Define the Hamming distances (as a vector)
hamming_distances <- 0:8  # Limited to 0 to 8

# Hardcoded example for lambda values calculated previously
# Replace this with the actual lambda values from your earlier script if needed.
lambda_values <- c(0.001, 0.01, 0.1, 1, 5, 10)  # Example lambda values; adjust accordingly

# Create a data frame to store results
pdf_data <- data.frame()

# Calculate the PDF for each combination of lambda and Hamming distance
for (lambda in lambda_values) {
  # Calculate the PDF for the current lambda
  pdf_values <- sapply(hamming_distances, calculate_pdf, lambda = lambda)
  
  # Create a temporary data frame for the current lambda
  temp_data <- data.frame(Hamming_Distance = hamming_distances, 
                          PDF = pdf_values, 
                          Lambda = as.factor(lambda))
  
  # Combine the temporary data frame with the main data frame
  pdf_data <- rbind(pdf_data, temp_data)
}

# Plot the results with the Viridis "D" color palette
ggplot(pdf_data, aes(x = Hamming_Distance, y = PDF, color = Lambda)) +
  geom_line(size = 0.5) +
  labs(title = "Probability Density Function of Hamming Distance",
       x = "Hamming Distance (k)",
       y = "Probability Density Function P(k ≥ 0)",
       color = "Lambda (λ)") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +  # Set y-axis limit to 0-1
  scale_x_continuous(breaks = 0:8, limits = c(0, 8), expand = c(0, 0)) +  # Show x-axis from 0 to 8
  scale_color_viridis_d(option = "D") +  # Use the viridis D color palette
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"))  # Set axes to black

### END ###

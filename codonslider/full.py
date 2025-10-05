import matplotlib.pyplot as plt
from codonslider import maxslider
from codonslider import geneslider

# Inputs
steps = 3
codon_of_interest = str(input("Which codon should be processed? ").upper())
initial_window_size = int(input("Initial size of the sliding window? "))
max_window_size = int(input("Maximum size of the sliding window? "))
path_to_file = 'CCDS_nucleotide.current.fasta'

# Run the codon slider: for longest codon stretch per sliding window size, use 'maxslider'. For identifying the genes, use 'geneslider'
final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot = maxslider(path_to_file=path_to_file,
                                                                                steps=steps, codon_of_interest=codon_of_interest,
                                                                                initial_window_size=initial_window_size,
                                                                                max_window_size=max_window_size)

# Save the final DataFrame to a CSV file
final_result_df.to_csv("codon_stretches.csv", index=False)

# Display and save the graph
plt.plot(window_sizes_for_plot, max_codon_frequencies_for_plot, marker='o')
plt.title(codon_of_interest)
plt.xlabel('Window Size (nt)')
plt.ylabel('Highest Codon Density')
plt.grid(False)
plt.savefig("Highest codon frequency across all records.tif", format="tif")
plt.show()
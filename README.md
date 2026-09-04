# Codonslider -> codon stretch analyzer

codonslider is a new package implemented in python to determine the percentage of codon stretches within certain region of the gene.
It works by applying the sliding window principle which scans the CCDS sequence and records the percentage of the codon of interest.
The sliding window can grow overtime and record the codon accumulation per window size. 
It is a useful tool for people who study protein translation and tRNA biology.

- For more information about the gtAI: http://willbeadded!

## Python Support

Python >=3.13 is required.

## Dependencies

1. Biopython
2. pandas
3. datetime

## Installation
**Using pip**

```python
pip install codonslider
```
## Contribution Guidelines

We welcome contributions to the software!

For reporting bugs or making suggestions, the most effective method is to raise an issue on the GitHub issue tracker. 
GitHub allows you to classify your issues, letting us know whether it's a bug report, feature request, or feedback for the authors.

If you would like to contribute changes to the code, please submit a [pull request](https://github.com/AliYoussef96/gtAI/pulls). 
For guidance on creating a pull request, please refer to the [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests).

## Usage

```python
# Import "maxslider" for identifying the longest stretch in the transcriptome
# Import "geneslider" if you would like to identify the genes and the stretch associated with them
# The insturctions are identical for both functions

from codonslider import maxslider
import matplotlib.pyplot as plt

path_to_file = "CCDS.fasta"
steps = 3 # it should always be 3 if analysing codons!
codon_of_interest = "AGA"
initial_window_size = 3
max_window_size = 12

final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot = maxslider(path_to_file = path_to_file,
                                                                                steps = steps, codon_of_interest = codon_of_interest,
                                                                                initial_window_size = initial_window_size,
                                                                                max_window_size = max_window_size)
# Save the final DataFrame to a CSV file
final_result_df.to_csv("codon_stretches.csv", index=False)

# Display and save the graph
plt.plot(window_sizes_for_plot, max_codon_frequencies_for_plot, marker='o')
plt.title(codon_of_interest)
plt.xlabel('Window Size (nt)')
plt.ylabel('Highest Codon Density')
plt.grid(False)
plt.savefig("codon_stretches.tif", format="tif")
plt.show()
```

Where:

```
path_to_file: a CCDS fasta file which contains all sequences that needs to be analyzed
steps = How many nucleotides each slide moves. It should always be 3 since each codon consists of 3 nucleotides
codon_of_interest = a codon you would like to analyse
initial_window_size = the initial size of the sliding window. Minimum is 3
max_window_size = the final size of the growing sliding window
```

Returns:

```
final_result_df = a dataframe with a list of all genes and their codon frequency
window_sizes_for_plot = the progression of the window expansion
max_codon_frequencies_for_plot = maximum codon accumulation recorded for the current window size
```


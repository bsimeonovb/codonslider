from Bio import SeqIO
from datetime import datetime
import pandas as pd

def slidwin(elements, current_window_size, step, codon, max_window_size, identifier):
    if len(elements) <= current_window_size:
        return {'Window_Size': current_window_size, 'Codon_Frequency': 0, 'Sequence': elements}

    max_window_size = min(max_window_size, len(elements))  # Ensure the maximum window size is within the sequence length
    max_codon_freq = 0  # Track the maximum codon frequency across all records for the current window size
    sequence_with_max_codon_freq = ''

    for a in range(0, len(elements) - current_window_size + 1, step):
        win = elements[a:a + current_window_size]

        # Count the codon of interest in the current window
        if len(win) % 3 == 0:
            codons = [win[b:b + 3] for b in range(0, len(win), 3)]
            codon_freq = codons.count(codon) / (len(win) / step)

            # Update the maximum codon frequency information if needed
            if codon_freq >= max_codon_freq:
                max_codon_freq = codon_freq
                sequence_with_max_codon_freq = win

    return {'Window_Size': current_window_size, 'Codon_Frequency': max_codon_freq, 'Sequence': sequence_with_max_codon_freq, 'Gene': identifier}


def maxslider(path_to_file, steps, codon_of_interest, initial_window_size, max_window_size):
    """
    Reads a FASTA file, iterates through sequence records, and applies the
    sliding window analysis (slidwin) for a range of window sizes.

    Args:
        path_to_file (str): Path to the input CCDS FASTA file.
        steps (int): The step size for the slidwin function.
        codon_of_interest (str): The codon to count (e.g., 'ATG').
        initial_window_size (int): The starting size for the sliding window.
        max_window_size (int): The maximum size for the sliding window.

    Returns:
        tuple: A tuple containing (DataFrame of results, list of window sizes, list of max frequencies).
    """

    # Initialize variables to keep track of the highest density information
    all_results = []
    current_window_size = initial_window_size
    window_sizes_for_plot = []
    max_codon_frequencies_for_plot = []

    # storing the current time in the variable
    date = datetime.now()
    # Displays Time
    current_time = date.strftime('%H:%M:%S')
    print("Starting time:", date)
    # Loop over different window sizes
    while current_window_size <= max_window_size:
        # NOTE: Opening the file inside the loop is necessary since SeqIO.parse
        # consumes the generator. You must open/close on each iteration.
        with open(path_to_file, mode='r') as handle:
            max_codon_freq_across_records = 0
            sequence_with_max_codon_freq_across_records = ''

            for record in SeqIO.parse(handle, 'fasta'):
                identifier = record.id
                sequence = record.seq

                # Skip records with sequences shorter than the current window size
                if len(sequence) < current_window_size:
                    # print(f"Skipping record {identifier} with sequence length {len(sequence)}")
                    continue

                # Use the imported slidwin function
                active_info = slidwin(str(sequence), current_window_size, steps, codon_of_interest,
                                      max_window_size, identifier)

                # Update the maximum codon frequency information across all records
                if active_info['Codon_Frequency'] >= max_codon_freq_across_records:
                    max_codon_freq_across_records = active_info['Codon_Frequency']
                    sequence_with_max_codon_freq_across_records = active_info['Sequence']

            # Append the maximum codon frequency across all records for the current window size
            all_results.append({'Window_Size': current_window_size,
                                'Codon_Frequency': max_codon_freq_across_records,
                                'Sequence': sequence_with_max_codon_freq_across_records})

            # Append data for plotting
            window_sizes_for_plot.append(current_window_size)
            max_codon_frequencies_for_plot.append(max_codon_freq_across_records)

            print(
                f"Highest codon frequency across all records for window size {current_window_size} / Frequency: {max_codon_freq_across_records} / Gene: {identifier} / finished at: {datetime.now().strftime('%H:%M:%S')}")

        # Increase the window size by 3 units
        current_window_size += 3

    final_result_df = pd.DataFrame(all_results)

    # Return the results needed by the remaining script
    return final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot


def geneslider(path_to_file, steps, codon_of_interest, initial_window_size, max_window_size):
    """
    Reads a FASTA file, iterates through sequence records, and applies the
    sliding_growing_window analysis for a range of window sizes.

    Args:
        path_to_file (str): Path to the input FASTA file.
        steps (int): The step size for the sliding_growing_window function.
        codon_of_interest (str): The codon to count (e.g., 'ATG').
        initial_window_size (int): The starting size for the sliding window.
        max_window_size (int): The maximum size for the sliding window.

    Returns:
        tuple: A tuple containing (DataFrame of results, list of window sizes, list of max frequencies).
    """

    # Initialize variables to keep track of the highest density information
    all_results = []
    current_window_size = initial_window_size
    window_sizes_for_plot = []
    max_codon_frequencies_for_plot = []

    # storing the current time in the variable
    date = datetime.now()
    # Displays Time
    current_time = date.strftime('%H:%M:%S')
    print("Starting time:", date)

    # Loop over different window sizes
    while current_window_size <= max_window_size:
        # Variables for tracking highest frequency in the current window size iteration
        frequencies_in_current_window = []

        # NOTE: The file must be opened inside the loop because SeqIO.parse consumes the generator
        with open(path_to_file, mode='r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                identifier = record.id
                sequence = str(record.seq)  # Ensure sequence is a string for slicing

                # Skip sequences shorter than the current window size
                if len(sequence) < current_window_size:
                    # print(f"Skipping protein {identifier} with sequence length {len(sequence)}")
                    continue

                # Apply the sliding window function
                active_info = slidwin(sequence, current_window_size, steps, codon_of_interest,
                                                     max_window_size, identifier)

                # Append result for the current sequence/window size combination
                all_results.append(active_info)

                # Keep track of the frequencies just for the current window size iteration
                frequencies_in_current_window.append(active_info['Codon_Frequency'])

        # --- End of iteration for the current_window_size ---

        # Update plotting lists with the highest frequency found in this window size
        if frequencies_in_current_window:
            max_freq_for_current_window = max(frequencies_in_current_window)
            window_sizes_for_plot.append(current_window_size)
            max_codon_frequencies_for_plot.append(max_freq_for_current_window)

            # Print information for the current window size
            print(
                f"Highest codon frequency across all proteins for window size {current_window_size} / Frequency: {max_freq_for_current_window} / Gene: {identifier} / finished at: {datetime.now().strftime('%H:%M:%S')}")
        else:
            print(f"No records processed for window size {current_window_size}.")

        # Increase the window size by 3 units
        current_window_size += 3

    # Create the final DataFrame
    final_result_df = pd.DataFrame(all_results)

    # Return the results needed by the remaining script (for saving/plotting)
    return final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot

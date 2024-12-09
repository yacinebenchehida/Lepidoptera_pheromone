import os
import sys

# Define function for filtering lines from file
def filter_lines(file_path, threshold_col4=300, threshold_col5=0.20, threshold_col7=100):
    """Process the file and filter lines based on conditions."""
    cdef dict pairs = {}  # Store the pairs as seen in your awk script
    cdef list columns
    cdef float col4, col5, col7
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.split()
            # Ensure we have at least 7 columns for proper filtering
            if len(columns) >= 7:
                try:
                    # Correct column indexing based on your input
                    col4 = float(columns[3])  # 4th column (index 3)
                    col5 = float(columns[4])  # 5th column (index 4)
                    col7 = float(columns[6])  # 7th column (index 6)
                except ValueError as e:
                    # Print and skip any lines with non-numeric values in the relevant columns
                    print(f"Skipping line (invalid number): {line.strip()}")
                    continue

                # Apply the filters for all three conditions
                if col5 > threshold_col5 and col7 < threshold_col7 and col4 > threshold_col4:
                    key1 = columns[0]
                    key2 = columns[1]
                    # Store the pairs, appending the values if they exist
                    if key1 not in pairs:
                        pairs[key1] = []
                    if key2 not in pairs:
                        pairs[key2] = []
                    pairs[key1].append(key2)
                    pairs[key2].append(key1)
    
    return pairs

# Define function to print groups of linked keys (pairs)
def print_groups(pairs):
    """Print the groups of linked keys (pairs)."""
    cdef set visited = set()  # Use set to track visited items
    cdef list group, to_visit, links
    cdef str current, neighbor, key

    for key in pairs:
        if key not in visited:
            group = [key]
            visited.add(key)
            to_visit = [key]
            while to_visit:
                current = to_visit.pop()
                for neighbor in pairs[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        group.append(neighbor)
                        to_visit.append(neighbor)
            # Print the group
            print(' '.join(group))

# Main function to run the script
def main():
    """Main function to run the script."""
    if len(sys.argv) != 5:
        print("Usage: python script.py <directory_path> <threshold_col4> <threshold_col5> <threshold_col7>")
        sys.exit(1)

    # Get the directory containing txt files and the thresholds from arguments
    directory_path = sys.argv[1]
    threshold_col4 = float(sys.argv[2])
    threshold_col5 = float(sys.argv[3])
    threshold_col7 = float(sys.argv[4])

    cdef dict all_pairs = {}
    cdef bint found_files = False  # Using 'bint' for boolean value
    
    # Traverse through all files in the directory
    for root, _, files in os.walk(directory_path):  # Use os.walk to traverse subdirectories
        for file_name in files:
            if file_name.endswith(".txt"):
                found_files = True
                file_path = os.path.join(root, file_name)
                pairs = filter_lines(file_path, threshold_col4, threshold_col5, threshold_col7)
                # Merge the pairs into a single dictionary
                for key, value in pairs.items():
                    if key not in all_pairs:
                        all_pairs[key] = []
                    all_pairs[key].extend(value)

    if not found_files:
        print("No .txt files found in the specified directory.")
        sys.exit(1)

    print_groups(all_pairs)

# This is needed to call the main function when running the script
if __name__ == "__main__":
    main()

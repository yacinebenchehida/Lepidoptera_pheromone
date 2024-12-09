import os
import sys

def filter_lines(file_path, threshold_col4=300, threshold_col5=0.20, threshold_col7=100):
    """Process the file and filter lines based on conditions."""
    pairs = {}  # Store the pairs as seen in your awk script
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

def print_groups(pairs):
    """Print the groups of linked keys (pairs)."""
    visited = set()
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

    all_pairs = {}
    found_files = False  # To check if we find any valid files

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

if __name__ == "__main__":
    main()

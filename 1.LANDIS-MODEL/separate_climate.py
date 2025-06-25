import pandas as pd
from pathlib import Path

# Load the dataset as a plain text file
input_file = Path(__file__).parent / "LANDIS_run" / "ca_klamath_miroc585.csv"
output_file = Path(__file__).parent / "LANDIS_run" / "ca_klamath_miroc585_NJTedits.csv"

# Initialize variables
current_table = []
table_headers = None
current_table_name = None
all_tables = []  # List to hold all processed tables

# Open the file and process line by line
with open(input_file, "r") as f:
    for line in f:
        line = line.strip()

        # Check if the line is a separator
        if line.startswith("#"):
            # Process and save the previous table if it exists
            if current_table:
                # Create a DataFrame from the current table
                table_df = pd.DataFrame(current_table, columns=table_headers)

                # Retain only the specified rows and columns
                selected_rows = pd.concat(
                    [table_df.iloc[[0]], table_df.iloc[1828:5480]]
                )
                selected_columns = selected_rows.iloc[:, [0, 86, 256, 426]]
                # Insert a row at the top with the table name
                header_row = pd.DataFrame(
                    [
                        [f"#{current_table_name}"]
                        + [""] * (selected_columns.shape[1] - 1)
                    ],
                    columns=selected_columns.columns,
                )
                column_headers = pd.DataFrame(
                    [selected_columns.columns.tolist()],
                    columns=selected_columns.columns,
                )
                selected_columns = pd.concat(
                    [header_row, column_headers, selected_columns], ignore_index=True
                )
                all_tables.append(selected_columns)

            # Reset for the new table
            current_table = []
            table_headers = None
            current_table_name = line[1:].strip().replace(",", "")  # Extract table name
            print(f"Reading {current_table_name} table...")
        else:
            # Parse the line into columns
            row = line.split(",")  # Adjust delimiter if necessary

            # Determine if it's a header or a data row
            if table_headers is None:
                table_headers = row  # First row after separator is the header
            else:
                current_table.append(row)

# Process the last table if it exists
if current_table:
    table_df = pd.DataFrame(current_table, columns=table_headers)
    selected_rows = pd.concat([table_df.iloc[[0]], table_df.iloc[1828:5480]])
    selected_columns = selected_rows.iloc[:, [0, 86, 256, 426]]
    header_row = pd.DataFrame(
        [[f"#{current_table_name}"] + [""] * (selected_columns.shape[1] - 1)],
        columns=selected_columns.columns,
    )
    column_headers = pd.DataFrame(
        [selected_columns.columns.tolist()], columns=selected_columns.columns
    )
    selected_columns = pd.concat(
        [header_row, column_headers, selected_columns], ignore_index=True
    )
    all_tables.append(selected_columns)

# Combine all processed tables into one CSV file
with open(output_file, "w", newline="") as f:  # Ensure no extra newline is added
    for table in all_tables:
        table.to_csv(
            f, index=False, header=False, line_terminator="\n"
        )  # Avoid extra blank rows

print(f"Recombined dataset saved to {output_file}")

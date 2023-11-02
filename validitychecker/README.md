# Search Scheme Validator

This Python script is designed to check the validity of a search scheme specified in a text file. A valid search scheme consists of a set of searches, where each error distribution for at most k errors and p parts is covered.

## Requirements

- Python 3.x

## Usage

1. Clone this repository or download the script to your local environment.

2. Ensure you have Python 3.x installed on your system.

3. Run the script from the command line as follows:

   ```bash
   python validitychecker.py <file> <k>
   ```

   - `validitychecker.py` is the name of the Python script.
   - vile is the path to the file containing the search scheme.
   - `k` is an integer representing the value for 'k'.

## File Format for Search Scheme

The search scheme file  must adhere to the following format:

- Each line represents a single search definition, with three lists enclosed in curly braces. The lists are separated by whitespace.
- Each list should contain a comma-separated list of integers.
- The three sets must represent the following:
  1. The permutation (`pi`) of the parts.
  2. The lower bounds (`L`) for each part.
  3. The upper bounds (`U`) for each part.
- The sets should not have any spaces within the curly braces.

Example of a valid search definition line:

```plaintext
{0,1,2,3} {0,0,0,1} {0,1,3,3}
```

## Example

Suppose you have a search scheme file `search_scheme.txt` for `k = 2` with the following content:

```plaintext
{0,1,2,3} {0,0,0,2} {0,1,2,2}
{1,2,3,0} {0,0,1,1} {0,1,1,2}
{2,3,1,0} {0,0,0,0} {0,0,2,2}
```

You can check its validity with the following command:

```bash
python validitychecker.py search_scheme.txt 2
```

## Output

- If the search scheme is valid, the script will display "No errors were raised. Search scheme is valid!" on the console.

- If the search scheme is invalid, the script will raise an `InvalidSearchSchemeError` or a `ValueError` exception and provide information about the error.

## Author

- Luca Renders

## License

This project is licensed under the AGPL-3.0 license - see the [LICENSE.md](../LICENSE.md) file for details.

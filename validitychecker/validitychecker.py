from typing import List, Tuple
from itertools import combinations
import os

# Function to calculate n choose k
def n_choose_k(n: int, k: int) -> int:
    """
    Calculate the number of combinations "n choose k."

    Args:
        n (int): Total number of items.
        k (int): Number of items to choose.

    Returns:
        int: The number of combinations.
    """
    res = 1
    for i in range(1, k + 1):
        res *= (n - k + i) / i
    return int(res + 0.01)

# Function to generate the lexicographic m-subset successor of a set
def successor(A: List[int], m: int, n: int) -> bool:
    """
    Generate the lexicographic m-subset successor of a given set A.

    Args:
        A (List[int]): The input set.
        m (int): Size of the subset.
        n (int): Total number of elements.

    Returns:
        bool: True if a successor exists, False if not.
    """
    i = m - 1
    while A[i] == n - m + i:
        i -= 1
        if i < 0:
            return False

    for j in range(m - 1, i - 1, -1):
        A[j] = A[i] + 1 + j - i
    return True

# Function to generate all error distributions with at most k errors divided over p parts
def generate_error_distributions(k: int, p: int) -> List[int]:
    """
    Generate all possible error distributions with at most k errors divided over p parts.

    Args:
        k (int): Maximum number of errors.
        p (int): Number of parts to divide the errors.

    Returns:
        List[int]: List of error distributions.
    """
    A = list(range(p))
    A.append(k + p)

    distributions = []

    while True:
        distribution = [A[i + 1] - A[i] - 1 for i in range(p)]
        distributions.append(distribution)

        if not successor(A, p, p + k):
            return distributions

# Custom exception class for invalid search schemes
class InvalidSearchSchemeError(Exception):
    """
    Custom exception class for handling invalid search schemes.
    """
    def __init__(self, message):
        super().__init__(message)

# Class to represent a search
class Search:
    def __init__(self, line):
        """
        Initialize a Search object from a given line.

        Args:
            line (str): Input line containing search information.
        """
        self.pi = []
        self.L = []
        self.U = []
        self.p = 0

        # Split the input line into three parts using curly braces
        parts = line.strip().split()
        if len(parts) != 3:
            raise ValueError("Invalid input format. Must have three sets in curly braces.")

        for i in range(3):
            elements = parts[i].strip('{}').split(',')
            elements = [int(e) for e in elements]

            if i == 0:
                self.pi = elements
                self.p = len(self.pi)

                # Check if pi is a permutation of 0..p-1
                if sorted(self.pi) != list(range(self.p)):
                    raise ValueError("pi must be a permutation of 0..p-1")

                # Check if elements in pi are connected
                for j in range(1, self.p):
                    if self.pi[j] != min(self.pi[0:j]) - 1 and self.pi[j] != max(self.pi[0:j]) + 1:
                        raise ValueError("Elements in pi must be connected.")
            elif i == 1:
                self.L = elements
            elif i == 2:
                self.U = elements

        # Check if all three lists have the same size 'p'
        if len(self.L) != self.p or len(self.U) != self.p:
            raise ValueError("All three lists should have the same size 'p'.")

        # Check if L[i] is no greater than U[i]
        for i in range(self.p):
            if self.L[i] > self.U[i]:
                raise ValueError("L[i] must not be greater than U[i].")

        # Check if elements in U and L are non-decreasing
        if any(self.U[i] < self.U[i - 1] for i in range(1, self.p)) or any(self.L[i] < self.L[i - 1] for i in range(1, self.p)):
            raise ValueError("Elements in U and L should be non-decreasing.")

    # Function to check if a distribution covers a search
    def covers(self, distribution: List) -> bool:
        """
        Check if a given error distribution covers the search.

        Args:
            distribution (List[int]): Error distribution.

        Returns:
            bool: True if the distribution covers the search, False otherwise.
        """
        cumulError = 0
        for i in range(self.p):
            cumulError += distribution[self.pi[i]]
            if self.L[i] > cumulError or self.U[i] < cumulError:
                return False
        return True

    def __str__(self):
        pi_str = f"pi: {self.pi}"
        L_str = f"L: {self.L}"
        U_str = f"U: {self.U}"
        return f"[{pi_str},{L_str},{U_str}]"

    def __repr__(self):
        return self.__str()

# Class to represent a search scheme
class SearchScheme:
    def __init__(self, filename, k):
        """
        Initialize a SearchScheme from a file and a value for 'k'.

        Args:
            filename (str): Path to the search scheme file.
            k (int): Value for 'k'.
        """
        self.searches = []
        self.k = k

        if not os.path.isfile(filename):
            raise FileNotFoundError(f"The file '{filename}' does not exist.")

        with open(filename, 'r') as file:
            for line_number, line in enumerate(file, start=1):
                line = line.strip()
                if not line:
                    continue

                try:
                    search = Search(line)
                    self.searches.append(search)
                except Exception as error:
                    raise ValueError(f"Problem with line {line_number}: {line}\n{error}")

        if not self.same_p_check():
            raise InvalidSearchSchemeError(f"Invalid scheme, not all searches have the same number of parts!")

        # check if U or L don't go over k
        for search in self.searches:
            if search.L[self.p - 1] > k or search.U[self.p - 1] > k:
                raise InvalidSearchSchemeError(f"Bounds must not be greater than k (={self.k})\nViolated in search: {search}")
        try:
            self.check_coverage()
        except Exception as error:
            raise error

    # Function to check if all searches have the same 'p'
    def same_p_check(self):
        if not self.searches:
            return True

        first_p = self.searches[0].p
        self.p = first_p
        return all(search.p == first_p for search in self.searches)

    # Function to check if a distribution is covered by any of the searches
    def covers(self, distribution: List[int]) -> bool:
        """
        Check if a given error distribution is covered by any of the searches.

        Args:
            distribution (List[int]): Error distribution.

        Returns:
            bool: True if the distribution is covered, False otherwise.
        """
        return any(search.covers(distribution) for search in self.searches)

    # Function to check coverage for generated error distributions
    def check_coverage(self):
        print(f"Generating {n_choose_k(self.p + self.k, self.k)} error distributions")
        distributions = generate_error_distributions(self.k, self.p)
        print(f"Checking coverage of {len(distributions)} error distributions...")

        for distribution in distributions:
            if not (self.covers(distribution)):
                raise InvalidSearchSchemeError(f"The provided scheme is invalid.\nDistribution {distribution} is not covered.")

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check the validity of a search scheme file.")
    parser.add_argument("filename", help="Path to the search scheme file.")
    parser.add_argument("k", type=int, help="Value for 'k'.")

    args = parser.parse_args()
    filename = args.filename
    k = args.k

    search_scheme = SearchScheme(filename, k)
    print("No errors were raised. Search scheme is valid!")

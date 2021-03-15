/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2021 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/
#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include <cstdio>
#include <iostream>
#include <vector>

typedef uint32_t length_t;

// ============================================================================
// CLASS EDIT DISTANCE
// ============================================================================

// The edit Distance class has two fields, an unsigned integer that stores the
// actual edit distance and a boolean flagged, signifiying if this ED originated
// from a flagged path

// the data is a 32 unsigned bit and the most significant bit is reserved for
// the flag while the other 31 bits are for the actual edit distance
class EditDistance {
  private:
    length_t data;

  public:
    EditDistance() {
        data = 0u;
    }

    EditDistance(int edit) {
        data = edit;
    }

    EditDistance(int edit, bool flag) {
        data = (flag) ? (1u << 31) | edit : edit;
    }

    void setED(int edit) {
        data = (data & 31) | edit;
    }

    unsigned int getED() const {
        unsigned int ret = (data << 1) >> 1;
        return ret;
    }

    bool isFlagged() const {
        return ((data >> 31) > 0);
    }

    /**
     * Operator overloading creates a new edit distance with the added distance
     * @param toAdd the value to add to this edit distance
     * @returns a new EditDistance instance with the increased edit distance and
     * same flag as the current one
     */
    EditDistance operator+(unsigned int toAdd) const {

        return EditDistance(getED() + toAdd, isFlagged());
    }

    /**
     * Operator overloading. Comparing two edit distances goes as follows: if
     * their ED's are equal then a flagged instance is smaller than a non
     * flagged. Else the one with the smallest ED is the smallest
     * @param argument the EditDistance to compare to this
     * @returns a bool indicating wheter this is smaller
     */
    bool operator<=(const EditDistance& argument) const {
        // first, look at equal edit distance

        if (getED() == argument.getED()) {
            // both flagged or not flagged => equal => true
            // if this flagged and other not => true
            // if this is not flagged and other is => false
            return (isFlagged() || !argument.isFlagged());
        }
        return getED() < argument.getED();
    }

    /**
     * Operator overloading.
     * This returns a bool indicating wheter this is smaller then some integer
     * @param argument the integer to compare to
     * @returns true if this is not flagged and smaller or equal, else it will
     * return false
     */
    bool operator<=(unsigned int argument) const {
        // a flagged reported instance is never smaller then an argument
        return !isFlagged() && getED() <= argument;
    }

    /**
     * Makes a string representation of this.
     * @returns the string representation of this EditDistance
     */
    std::string to_string() {
        std::string flagString = (isFlagged() ? "*" : "");
        return std::to_string(getED()) + flagString;
    }
};

// ============================================================================
// CLASS BANDED MATRIX
// ============================================================================

// The band matrix class can best be understood as m horizontal bands each
// with a width of 2W+1 elements. It is allocated as a single m*(2W+1) array.
// For example, for W = 2 and m = 6.
// XX|XXX...|
//  X|XXXX..|
//   |XXXXX.|
//   |.XXXXX|
//   |..XXXX|X
//   |...XXX|XX
// The actual matrix is between |.| The Xs left and right of the |.| are
// allocated but should in principle not be addressed.
// The storage order is row-major.

class BandMatrix {
  private:
    std::vector<length_t> matrix;
    length_t W; // The off diagonal width of this bandmatrix
    length_t m; // number of rows
    length_t n; // number of columns
    const static int Wprod = 16;

    int finalCellFirstCol; // the final row of the zeroth column
    int rowsPerColumn; // the number of rows per column (starting from the Wth
                       // column)

    void initializeMatrix(const std::vector<int>& eds, const int& increase) {

        // initialize the first column
        length_t row = 0;
        for (auto rIt = eds.begin(); rIt != eds.end(); ++rIt) {
            operator()(row++, 0) = *rIt + increase;
        }
        for (; operator()(row - 1, 0) <= eds[0] + increase + W; row++) {
            operator()(row, 0) = operator()(row - 1, 0) + 1;
        }

        // intialize the first row
        for (length_t i = 1; i <= W; i++) {
            operator()(0, i) = i + eds[0] + increase;
        }
        // initialize the cells to the right of band
        for (length_t c = W + 1; c < m; c++) {
            operator()(c - (W + 1), c) = W + eds[0] + increase + 1;
        }

        // initialize the cells to the left of band
        for (length_t r = 1; r + row <= m; r++) {
            operator()(r + row - 1, r) = W + eds[0] + increase + 1;
        }
        finalCellFirstCol = row - 2;
    }

    void initializeMatrix(length_t startValue) {
        finalCellFirstCol = W;
        for (length_t i = 0; i <= W + 1; i++) {
            operator()(0, i) = i + startValue;
            operator()(i, 0) = i + startValue;
        }
        // set max elements at sides
        // first the elements on rows [1, W]
        for (length_t i = 1; i <= W; i++) {
            // right of band
            operator()(i, i + W + 1) = W + 1 + startValue;
        }

        // then the elements on rows [W + 1, x]

        for (length_t i = W + 1; i + W + 1 < n; i++) {
            // right of band
            operator()(i, i + W + 1) = W + 1 + startValue;
            // left of band
            operator()(i, i - (W + 1)) = W + 1 + startValue;
        }

        for (length_t i = std::max<int>((int)n - (W + 1), W + 1); i < m; i++) {
            // left of band
            operator()(i, i - (W + 1)) = W + 1 + startValue;
        }
    }

  public:
    /**
     * Constructor
     * @param pieceSize, the size of the piece to match, this will initialize
     * the top row
     * @param W, the width needed for this matrix
     * @param startValue, the value found at the startmatch, this should be the
     * minimum value of the vector eds if this vector is not empty
     * @param eds, a vector to initialize the first column, if an empty vector
     * is provided the first column will be initialized starting from startvalue
     * and up to W
     */
    BandMatrix(length_t pieceSize, int W, int startValue,
               const std::vector<int>& eds)
        : W(W), n(pieceSize + 1) {

        if (eds.empty()) {
            m = pieceSize + W + 1;
            matrix.resize(m * Wprod);
            rowsPerColumn = 2 * W + 1;
            initializeMatrix(startValue);
            return;
        }
        // minimum value should be put at startmatch...

        m = pieceSize + eds.size() + (W + eds[0] - eds.back());
        matrix.resize(m * Wprod);
        rowsPerColumn = m - pieceSize + W;

        initializeMatrix(eds, startValue);
    }

    /**
     * Constructor
     * @param pieceSize, the size of the piece to match, this will initialize
     * the top row
     * @param W, the width needed for this matrix
     * @param startValue, the value found at the startmatch, this will be put at
     * the origin
     */
    BandMatrix(length_t pieceSize, int W, int startValue)
        : W(W), m(pieceSize + W + 1), n(pieceSize + 1) {
        matrix.resize(m * Wprod);
        rowsPerColumn = 2 * W + 1;
        initializeMatrix(startValue);
        return;
    }

    /**
     * Constructor
     * @param m Number of rows
     * @param W Number of off-diagonal elements (one sided)
     */
    BandMatrix(length_t m, int W) : W(W), m(m) {
        matrix.resize(m * Wprod);
        n = m - W;
        initializeMatrix(0);
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Element at position (i, j)
     */
    length_t operator()(length_t i, int j) const {
        return matrix[i * Wprod + j - i + W];
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    length_t& operator()(length_t i, int j) {
        return matrix[i * Wprod + j - i + W];
    }

    /**
     * Get the band width
     * @return The band width
     */
    const length_t getWidth() const {
        return W;
    }

    void printMatrix(length_t maxRow = 500) const {
        length_t mRow = std::min<length_t>(maxRow + 1, m);
        for (length_t i = 0; i < mRow; i++) {
            length_t firstCol = 0;
            length_t lastCol = n - 1;
            std::string rowNumber = ((i < 10) ? "0" : "") + std::to_string(i);
            std::string row = "row " + rowNumber + ": ";

            for (length_t j = firstCol; j <= lastCol; j++) {
                int number = operator()(i, j);
                row += std::to_string(number) + " ";
            }
            for (length_t j = lastCol + 1; j < 2 * W + 1; j++) {
                row += ".\t";
            }
            std::cout << row << std::endl;
        }
        std::cout << "----------------------------------------------\n";
    }

    /**
     * Update the matrix by calculating the element at position row, column.
     * @param match whether the character was a match
     * @param row the row of the elment to update
     * @param collumn the collumn of the element to update
     * @returns the new value at row, column
     */
    length_t updateMatrix(bool notMatch, unsigned int row,
                          unsigned int column) {
        length_t diag = operator()(row - 1, column - 1) + notMatch;
        length_t gapX = operator()(row, column - 1) + 1;
        length_t gapY = operator()(row - 1, column) + 1;

        length_t returnValue =
            std::min<length_t>(diag, std::min<length_t>(gapX, gapY));

        operator()(row, column) = returnValue;
        return returnValue;
    }

    /**
     * Retrieves the first column that needs to be filled in for the row
     * @param row the row to fill in
     * @returns the first column to fill in
     */
    const int getFirstColumn(int row) const {
        // leftmost cell of band
        return std::max(1, row - finalCellFirstCol);
    }
    /**
     * Retrieves the last column that needs to be filled in for the row
     * @param row the row to fill in
     * @returns the last column to fill in
     */
    const int getLastColumn(int row) const {
        // rightmost cell of band
        return std::min(n - 1, W + row);
    }

    const int getNumberOfRows() const {
        return m - 1;
    }

    const int getSizeOfFinalColumn() const {
        if (n > W) {
            return rowsPerColumn;
        }
        return finalCellFirstCol + n;
    }

    const int getLastColumn() const {
        return n - 1;
    }
};

// ============================================================================
// CLASS EDIT MATRIX
// ============================================================================

// This is a band matrix but it elements are of class
// EditDistance. This checks whether an EditDistance originates from a path that
// has a leading gap and will flag as such

class EditMatrix {
  private:
    std::vector<EditDistance> matrix; // the matrix
    length_t W;                       // the width of the matrix
    const static int Wprod = 16;

    /** Helper function for update matrix, it checks which path to take, the
     * diagonal, horizontal or verical path.
     * @param diag the value for the element if the diagonal path is chosen
     * @param gapX the value for the element if the horizontal path is chosen
     * @param gapY the value for the elemnt if the vertical path is chosen
     *
     * @returns the best element out of the three provided.
     */
    EditDistance chooseBestElement(EditDistance& diag, EditDistance& gapX,
                                   EditDistance& gapY) {
        // need to check if for the optimal value one of the paths is flagged
        // so we need to find the smallest

        if ((diag <= gapX) && (diag <= gapY)) {
            return diag;
        }
        if (gapX <= gapY) {
            return gapX;
        }
        return gapY;
    }

    /**
     * Helper function for initializing the matrix. The top row and left most
     * collumn are filled in with as value their row/column number plus a
     * startvalue. If the noLeadingGaps argument is true than the leftmost
     * column will be flagged (except the origin).
     * @param noLeadingGaps a bool to indicate if leading gaps are allowed
     * @param startValue the value of the origin, default is 0
     */
    void initializeMatrix(bool noLeadingGapsText, bool noLeadingGapsPattern,
                          int startValue = 0) {
        for (length_t i = 0; i <= W; i++) {
            operator()(i, 0) = EditDistance(i + startValue, noLeadingGapsText);
            operator()(0, i) =
                EditDistance(i + startValue, noLeadingGapsPattern);
        }
        operator()(0, 0) = EditDistance(startValue, false);
    }

  public:
    /**
     * Constructor
     * @param m Number of rows
     * @param W Number of off-diagonal elements (one sided)
     */
    EditMatrix(length_t m, int W, int startValue, bool leadingGapsAllowedText,
               bool leadingGapsAllowedPattern)
        : W(W) {
        matrix.resize(m * Wprod);
        initializeMatrix(!leadingGapsAllowedText, !leadingGapsAllowedPattern,
                         startValue);
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Element at position (i, j)
     */
    EditDistance operator()(length_t i, int j) const {
        return matrix[i * Wprod + j - i + W];
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    EditDistance& operator()(length_t i, int j) {
        return matrix[i * Wprod + j - i + W];
    }

    /**
     * Get the band width
     * @return The band width
     */
    int getWidth() const {
        return W;
    }

    /**
     * Prints the matrix, for debugging purposes
     */
    void printMatrix() {
        length_t m = matrix.size() / Wprod;

        for (length_t i = 0; i < m; i++) {
            length_t firstCol = std::max<int>(1, i - W);
            length_t lastCol = std::min<int>(m, i + W);
            std::string row = "";
            for (length_t j = 0; j < firstCol; j++) {
                row += ".\t";
            }
            for (length_t j = firstCol; j <= lastCol; j++) {
                row += (operator()(i, j)).to_string() + "\t";
            }
            for (length_t j = lastCol + 1; j < 2 * W + 1; j++) {
                row += ".\t";
            }
            std::cout << row << std::endl;
        }
    }

    void printRow(int i) {
        length_t m = matrix.size() / Wprod;
        length_t firstCol = std::max<int>(1, i - W);
        length_t lastCol = std::min<int>(m, i + W);
        std::string row = "";
        for (length_t j = 0; j < firstCol; j++) {
            row += ".\t";
        }
        for (length_t j = firstCol; j <= lastCol; j++) {
            row += (operator()(i, j)).to_string() + "\t";
        }
        for (length_t j = lastCol + 1; j < 2 * W + 1; j++) {
            row += ".\t";
        }
        std::cout << row << std::endl;
    }

    /**
     * Update the matrix by calculating the element at position row, column.
     * @param match whether the character was a match
     * @param row the row of the elment to update
     * @param collumn the collumn of the element to update
     */
    void updateMatrix(bool notMatch, unsigned int row, unsigned int collumn) {
        EditDistance diag = operator()(row - 1, collumn - 1) + notMatch;

        // watch out for integer overflow
        unsigned int diffRowAndW = (row >= W) ? row - W : 0;
        EditDistance gapX =
            (collumn > diffRowAndW) ? operator()(row, collumn - 1) + 1
                                    : diag + 1;
        EditDistance gapY =
            (collumn < row + W) ? operator()(row - 1, collumn) + 1 : diag + 1;

        operator()(row, collumn) = chooseBestElement(diag, gapX, gapY);
    }
};

#endif

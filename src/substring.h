/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Luca Renders <luca.renders@ugent.be> and        *
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

#ifndef SUBSTRING_H
#define SUBSTRING_H

#include "definitions.h"

#include <string>

// ============================================================================
// CLASS SUBSTRING
// ============================================================================

/**
 * A class that represents a substring of a string.
 */
class Substring {
  private:
    const std::string* text; // pointer to the string this is a substring of
    length_t startIndex;     // the startIndex of this substring in the text
    length_t
        endIndex; // the endIndex of this substring in the text (non-inclusive)
    Direction d;  // The direction of this substring

    // function pointer for the [] operator depending on d
    using AccessorFunc = char (Substring::*)(length_t) const;

    AccessorFunc accessor;

    /**
     * Accessor if the direction is forward
     * @param i the index of the character to get
     */
    char forwardAccessor(length_t i) const {
        return (*text)[startIndex + i];
    }

    /**
     * Accessor if the direction is backward
     * @param i the index of the character to get
     */
    char reverseAccessor(length_t i) const {
        return (*text)[endIndex - i - 1];
    }

  public:
    /**
     * Constructor, the start and end index default of 0 and the size of the
     * text
     * @param t the text to point to
     * @param dir the direction (defaults to FORWARD)
     */
    Substring(const std::string& t, Direction dir = FORWARD)
        : text(&t), startIndex(0), endIndex(t.size()), d(dir),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
    }

    /**
     * Constructor, the start and end index default of 0 and the size of the
     * text
     * @param t the text to point to
     * @param dir the direction (defaults to FORWARD)
     */
    Substring(const std::string* t, Direction dir = FORWARD)
        : text(t), startIndex(0), endIndex(t->size()), d(dir),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
    }

    /**
     * Constructor, the direction defaults to FORWARD
     * @param t the text to point to
     * @param start the start index of this substring in t
     * @param end the end index of this substring in t (non-inclusive)
     * @param dir the direction (defaults to FORWARD)
     */
    Substring(const std::string* t, length_t start, length_t end,
              Direction dir = FORWARD)
        : text(t), startIndex(start), endIndex(end), d(dir),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
        check();
    }

    /**
     * Constructs a substring of the text another substring points to
     * @param s pointer to the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     */
    Substring(const Substring* s, length_t start, length_t end)
        : text(s->text), startIndex(start), endIndex(end), d(s->d),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
        check();
    }

    /**
     * Constructs a substring of the text another substring points and sets the
     * direction
     * @param s pointer to the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     * @param dir the direction of this new substring
     */
    Substring(const Substring* s, length_t start, length_t end, Direction dir)
        : text(s->text), startIndex(start), endIndex(end), d(dir),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
        check();
    }
    /**
     * Constructs a substring of the text another substring points to
     * @param s the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     */
    Substring(const Substring& s, length_t start, length_t end)
        : text(s.text), startIndex(start), endIndex(end), d(s.d),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
        check();
    }

    /**
     * Constructs a substring of the text another substring points and sets the
     * @param s the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     * @param dir the direction of this new substring
     */
    Substring(const Substring& s, length_t start, length_t end, Direction dir)
        : text(s.text), startIndex(start), endIndex(end), d(dir),
          accessor(d == FORWARD ? &Substring::forwardAccessor
                                : &Substring::reverseAccessor) {
        check();
    }

    /**
     * Set the direction of this substring
     * @param nd the new direction
     */
    void setDirection(Direction nd) {
        d = nd;
        accessor = (d == FORWARD ? &Substring::forwardAccessor
                                 : &Substring::reverseAccessor);
    }

    /**
     * Get the direction of this substring
     */
    Direction getDirection() const {
        return d;
    }

    /**
     * Get the character at index i of this substring
     * @param i the index to get the character from
     * @returns the character at index i
     */
    char operator[](length_t i) const {
        return (this->*accessor)(i);
    }

    /**
     * Get the size of this substring
     * @returns the size of this substring
     */
    size_t size() const {
        if (empty()) {
            return 0;
        }
        return endIndex - startIndex;
    }

    /**
     * Get the length of this substring (equals the size)
     * @returns the length of this substring
     */
    size_t length() const {
        return size();
    }

    /**
     * Check if this substring is empty
     * @returns a bool that indicates whether the substring was empty
     */
    bool empty() const {
        return endIndex <= startIndex;
    }

    void check() {

        if (endIndex > text->size()) {
            endIndex = text->size();
        }
    }

    /**
     * Get the end of this substring
     * @returns the end Index of this substring (non-inclusive)
     */
    length_t end() const {
        return endIndex;
    }

    /**
     * Get the begin of this substring
     * @returns the begin index of the substring
     */
    length_t begin() const {
        return startIndex;
    }

    Substring& operator=(const Substring& other) {
        this->text = other.text;
        this->startIndex = other.begin();
        this->endIndex = other.end();
        this->d = other.d;

        return *this;
    }

    /**
     * Converts the Substring to a c++ string
     * @returns, the Substring as a c++ string
     */
    std::string tostring() const {
        if (empty()) {
            return "";
        }
        return text->substr(startIndex, endIndex - startIndex);
    }

    void setEnd(length_t newEnd) {
        endIndex = newEnd;
    }

    void setBegin(length_t n) {
        startIndex = n;
    }

    void incrementEnd() {
        endIndex++;
    }

    void decrementBegin() {
        startIndex--;
    }

    /**
     * Check if this substring contains an N
     */
    bool containsN() const {
        for (length_t i = startIndex; i < endIndex; i++) {
            if (text->at(i) == 'N') {
                return true;
            }
        }
        return false;
    }


    /**
     * @brief Checks if this substring equals a given string
     * @param s the string to compare to
     * @returns true if the substring equals the string
     */
    bool equals(const std::string& s) const {
        if (s.size() != size()) {
            return false;
        }
        for (length_t i = 0; i < size(); i++) {
            if ((*this)[i] != s[i]) {
                return false;
            }
        }
        return true;}
};

#endif

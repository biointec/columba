#ifndef SUBSTRING_H
#define SUBSTRING_H

// ============================================================================
// ENUMS
// ============================================================================

/**
 * An enum for the direction of the search
 */
enum Direction { FORWARD, BACKWARD };

// ============================================================================
// CLASS SUBSTRING
// ============================================================================

class Substring {
  private:
    const std::string* text; // pointer to the string this is a substring of
    unsigned int startIndex; // the startIndex of this substring in the text
    unsigned int
        endIndex; // the endIndex of this substring in the text (non-inclusive)
    Direction d;  // The direction of this substring

  public:
    /**
     * Constructor, the start and end index default of 0 and the size of the
     * text
     * @param t, the text to point to
     * @param dir, the direction (defaults to FORWARD)
     */
    Substring(const std::string& t, Direction dir = FORWARD)
        : startIndex(0), endIndex(t.size()), d(dir) {
        text = &t;
    }

    /**
     * Constructor, the direction defaults to FORWARD
     * @param t, the text to point to
     * @param start, the start index of this substring in t
     * @param end, the end index of this substring in t (non-inclusive)
     * @param dir, the direction (defaults to FORWARD)
     */
    Substring(const std::string* t, unsigned int start, unsigned int end,
              Direction dir = FORWARD)
        : text(t), startIndex(start), endIndex(end), d(dir) {
    }

    /**
     * Constructs a substring of the text another substring points to
     * @param s, pointer to the other substring
     * @param start, the start index of this new substring in the original text
     * @param end, the end index of this new stubstring
     */
    Substring(const Substring* s, unsigned int start, unsigned int end)
        : text(s->text), startIndex(start), endIndex(end), d(s->d) {
    }

    /**
     * Constructs a substring of the text another substring points and sets the
     * direction
     * @param s, pointer to the other substring
     * @param start, the start index of this new substring in the original text
     * @param end, the end index of this new stubstring
     */
    Substring(const Substring* s, unsigned int start, unsigned int end,
              Direction dir)
        : text(s->text), startIndex(start), endIndex(end), d(dir) {
    }
    /**
     * Constructs a substring of the text another substring points to
     * @param s, the other substring
     * @param start, the start index of this new substring in the original text
     * @param end, the end index of this new stubstring
     */
    Substring(const Substring& s, unsigned int start, unsigned int end)
        : text(s.text), startIndex(start), endIndex(end), d(s.d) {
    }

    /**
     * Constructs a substring of the text another substring points and sets the
     * @param s, the other substring
     * @param start, the start index of this new substring in the original text
     * @param end, the end index of this new stubstring
     */
    Substring(const Substring& s, unsigned int start, unsigned int end,
              Direction dir)
        : text(s.text), startIndex(start), endIndex(end), d(dir) {
    }

    /**
     * Set the direction of this substring
     * @param nd, the new direction
     */
    void setDirection(Direction nd) {
        d = nd;
    }

    /**
     * Creates a substring of a substring, skipping the first skip characters
     * (relative to the direction)
     * @param skip, the number of characters to skip
     */
    const Substring getSubPiece(unsigned int skip) const {
        if (d == FORWARD) {
            return Substring(this, startIndex + skip, endIndex);
        } else {
            return Substring(this, startIndex, endIndex - skip);
        }
    }

    /**
     * Get the character at index i of this substring
     * @param i the index to get the character from
     * @returns the character at index i
     */
    char operator[](unsigned int i) const {
        return (d == FORWARD) ? text->at(startIndex + i)
                              : text->at(endIndex - i - 1);
    }

    /**
     * Get the size of this substring
     * @returns the size of this substring
     */
    unsigned int size() const {
        if (empty()) {
            return 0;
        }
        return endIndex - startIndex;
    }

    /**
     * Get the length of this substring (equals the size)
     * @returns the length of this substring
     */
    unsigned int length() const {
        return size();
    }

    /**
     * Check if this substring is empty
     * @returns a bool that indicates whether the substring was empty
     */
    bool empty() const {
        return endIndex <= startIndex;
    }

    /**
     * Get the end of this substring
     * @returns the end Index of this substring (non-inclusive)
     */
    unsigned int end() const {
        return endIndex;
    }

    /**
     * Get the begin of this substring
     * @returns the begin index of the substring
     */
    unsigned int begin() const {
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

    void setEnd(unsigned int newEnd) {
        endIndex = newEnd;
    }

    void setBegin(unsigned int n) {
        startIndex = n;
    }

    void incrementEnd() {
        endIndex++;
    }

    void decrementBegin() {
        startIndex--;
    }
};

#endif

!    * A_TO_I4 returns the index of an alphabetic character.
!    * B4_IEEE_TO_R4 converts a 4 byte IEEE word into an R4.
!    * B4_IEEE_TO_SEF converts an IEEE real word to S * 2**E * F format.
!    * B4_ULTRIX_TO_R4 converts a 4 byte ULTRIX word into an R4.
!    * B4_VMS_TO_R4 converts a 4 byte VMS word into an R4.
!    * BASE_TO_I4 returns the value of an integer represented in some base.
!    * BINARY_TO_I4 converts a binary representation into an integer value.
!    * BINARY_TO_R4 converts a binary representation into an R4.
!    * BINARY_TO_R8 converts a binary representation into an R4.
!    * BITS_TO_I4 converts a bit string into a 32 bit integer.
!    * BITS_TO_R4 converts a bit string into an R4.
!    * CH_COUNT_CHVEC_ADD adds a character vector to a character count.
!    * CH_COUNT_FILE_ADD adds characters in a file to a character count.
!    * CH_COUNT_HISTOGRAM_PRINT prints a histogram of a set of character counts.
!    * CH_COUNT_INIT initializes a character count.
!    * CH_COUNT_PRINT prints a set of character counts.
!    * CH_COUNT_S_ADD adds a character string to a character histogram.
!    * CH_EXTRACT extracts the next nonblank character from a string.
!    * CH_INDEX is the first occurrence of a character in a string.
!    * CH_INDEX_LAST is the last occurrence of a character in a string.
!    * CH_INDEXI is the (case insensitive) first occurrence of a character in a string.
!    * CH_IS_ALPHA is TRUE if CH is an alphabetic character.
!    * CH_IS_ALPHANUMERIC is TRUE if CH is alphanumeric.
!    * CH_IS_CONTROL is TRUE if a character is a control character.
!    * CH_IS_DIGIT is TRUE if a character is a decimal digit.
!    * CH_IS_FORMAT_CODE is TRUE if a character is a FORTRAN format code.
!    * CH_IS_LOWER is TRUE if a character is a lower case letter.
!    * CH_IS_PRINTABLE is TRUE if C is printable.
!    * CH_IS_SPACE is TRUE if a character is a whitespace character.
!    * CH_IS_UPPER is TRUE if CH is an upper case letter.
!    * CH_LOW lowercases a single character.
!    * CH_NEXT reads the next character from a string, ignoring blanks and commas.
!    * CH_NOT_CONTROL = CH is NOT a control character.
!    * CH_ROMAN_TO_I4 returns the integer value of a single Roman digit.
!    * CH_SCRABBLE returns the character on a given Scrabble tile.
!    * CH_SCRABBLE_FREQUENCY returns the Scrabble frequency of a character.
!    * CH_SCRABBLE_POINTS returns the Scrabble point value of a character.
!    * CH_SCRABBLE_SELECT selects a character with the Scrabble probability.
!    * CH_SWAP swaps two characters.
!    * CH_TO_AMINO_NAME converts a character to an amino acid name.
!    * CH_TO_BRAILLE converts an ASCII character to a Braille character string.
!    * CH_TO_CH3_AMINO converts a 1 character code to a 3 character code for an amino acid.
!    * CH_TO_DIGIT returns the integer value of a base 10 digit.
!    * CH_TO_DIGIT_BIN returns the integer value of a binary digit.
!    * CH_TO_DIGIT_HEX returns the integer value of a hexadecimal digit.
!    * CH_TO_DIGIT_OCT returns the integer value of an octal digit.
!    * CH_TO_EBCDIC converts a character to EBCDIC.
!    * CH_TO_MILITARY converts an ASCII character to a Military code word.
!    * CH_TO_MORSE converts an ASCII character to a Morse character string.
!    * CH_TO_ROT13 converts a character to its ROT13 equivalent.
!    * CH_TO_SCRABBLE returns the Scrabble index of a character.
!    * CH_TO_SOUNDEX converts an ASCII character to a Soundex character.
!    * CH_TO_SYM returns a printable symbol for any ASCII character.
!    * CH_UNIFORM returns a random character in a given range.
!    * CH3_TO_CH_AMINO converts a 3 character code to a 1 character code for an amino acid.
!    * CH4_TO_I4 converts a four character string to an integer.
!    * CH4_TO_R4 converts a 4 character string to an R4.
!    * CH4VEC_TO_I4VEC converts an string of characters into an array of integers.
!    * CHR4_TO_8 replaces pairs of hexadecimal digits by a character.
!    * CHR8_TO_4 replaces characters by a pair of hexadecimal digits.
!    * CHRA_TO_S replaces control characters by printable symbols.
!    * CHRASC converts a vector of ASCII codes into character strings.
!    * CHRASS "understands" an assignment statement of the form LHS = RHS.
!    * CHRCTF reads an integer or rational fraction from a string.
!    * CHRCTG reads an integer, decimal fraction or a ratio from a string.
!    * CHRCTI2 finds and reads an integer from a string.
!    * CHRCTP reads a parenthesized complex number from a string.
!    * CHRS_TO_A replaces all control symbols by control characters.
!    * CHVEC_PERMUTE permutes a character vector in place.
!    * CHVEC_PRINT prints a character vector.
!    * CHVEC_REVERSE reverses the elements of a character vector.
!    * CHVEC_TO_S converts a character vector to a string.
!    * CHVEC2_PRINT prints two vectors of characters.
!    * COMMA moves commas left through blanks in a string.
!    * DEC_TO_S_LEFT returns a left-justified representation of IVAL * 10**JVAL.
!    * DEC_TO_S_RIGHT returns a right justified representation of IVAL * 10**JVAL.
!    * DEGREES_TO_RADIANS converts an angle from degrees to radians.
!    * DIGIT_BIN_TO_CH returns the character representation of a binary digit.
!    * DIGIT_HEX_TO_CH returns the character representation of a hexadecimal digit.
!    * DIGIT_INC increments a decimal digit.
!    * DIGIT_OCT_TO_CH returns the character representation of an octal digit.
!    * DIGIT_TO_CH returns the character representation of a decimal digit.
!    * EBCDIC_TO_CH converts an EBCDIC character to ASCII.
!    * EBCDIC_TO_S converts a string of EBCDIC characters to ASCII.
!    * FILE_NAME_INC increments a partially numeric filename.
!    * FILLCH writes a string into a subfield of a string.
!    * FILLIN writes an integer into a subfield of a string.
!    * FILLRL writes a real into a subfield of a string.
!    * FLT_TO_S returns a representation of MANT * 10**IEXP.
!    * FORCOM splits a FORTRAN line into "fortran" and "comment".
!    * HEX_TO_I4 converts a hexadecimal string to its integer value.
!    * HEX_TO_S converts a hexadecimal string into characters.
!    * I2_BYTE_SWAP swaps bytes in an 8-byte word.
!    * I4_BYTE_SWAP swaps bytes in a 4-byte word.
!    * I4_GCD finds the greatest common divisor of I and J.
!    * I4_EXTRACT "extracts" an I4 from the beginning of a string.
!    * I4_INPUT prints a prompt string and reads an I4 from the user.
!    * I4_LENGTH computes the number of characters needed to print an I4.
!    * I4_NEXT "reads" I4's from a string, one at a time.
!    * I4_NEXT_READ finds and reads the next I4 in a string.
!    * I4_RANGE_INPUT reads a pair of I4's from the user, representing a range.
!    * I4_SQZ compresses the integer information in I4VEC into JVEC.
!    * I4_SWAP swaps two I4's.
!    * I4_TO_A returns the I-th alphabetic character.
!    * I4_TO_AMINO_CODE converts an integer to an amino code.
!    * I4_TO_BASE represents an integer in any base up to 16.
!    * I4_TO_BINARY produces the binary representation of an I4.
!    * I4_TO_BINHEX returns the I-th character in the BINHEX encoding.
!    * I4_TO_BITS converts an I4 to a string of 32 bits.
!    * I4_TO_CH4 converts an I4 to a 4 character string.
!    * I4_TO_HEX produces the hexadecimal representation of an I4.
!    * I4_TO_MONTH_NAME returns the name of a given month.
!    * I4_TO_NUNARY produces the "base -1" representation of an integer.
!    * I4_TO_OCT produces the octal representation of an integer.
!    * I4_TO_S_LEFT converts an I4 to a left-justified string.
!    * I4_TO_S_RIGHT converts an I4 to a right justified string.
!    * I4_TO_S_RIGHT_COMMA converts an I4 to a right justified string with commas.
!    * I4_TO_S_ROMAN converts an I4 to a string of Roman numerals.
!    * I4_TO_S_ZERO converts an I4 to a string, with zero padding.
!    * I4_TO_S32 converts an I4 to a 32 character string.
!    * I4_TO_UNARY produces the "base 1" representation of an I4.
!    * I4_TO_UUDECODE returns the I-th character in the UUDECODE encoding.
!    * I4_TO_XXDECODE returns the I-th character in the XXDECODE encoding.
!    * I4VEC_TO_CH4VEC converts an I4VEC into a string.
!    * IC_TO_IBRAILLE converts an ASCII integer code to a Braille code.
!    * IC_TO_IEBCDIC converts an ASCII character code to an EBCDIC code.
!    * IC_TO_IMORSE converts an ASCII integer code to a Morse integer code.
!    * IC_TO_ISOUNDEX converts an ASCII integer code to a Soundex integer code.
!    * IEBCDIC_TO_IC converts an EBCDIC character code to ASCII.
!    * ISTRCMP compares two strings, returning +1, 0, or -1.
!    * ISTRNCMP compares the start of two strings, returning +1, 0, or -1.
!    * LEN_NONNULL returns the length of a string up to the last non-null character.
!    * LOWER returns a lowercase version of a string.
!    * MALPHNUM2 is TRUE if a string contains only alphanumerics and underscores.
!    * MILITARY_TO_CH converts a Military code word to an ASCII character.
!    * MONTH_NAME_TO_I4 returns the month number of a given month
!    * NAMEFL replaces "lastname, firstname" by "firstname lastname".
!    * NAMELF replaces "firstname lastname" by "lastname, firstname".
!    * NAMELS reads a NAMELIST line, returning the variable name and value.
!    * NEXCHR returns the next nonblank character from a string.
!    * NEXSTR returns the next nonblank characters from a string.
!    * NUMBER_INC increments the integer represented by a string.
!    * OCT_TO_I4 converts an octal string to its integer value.
!    * PERM_CHECK checks that a vector represents a permutation.
!    * PERM_INVERSE3 produces the inverse of a given permutation.
!    * PERM_UNIFORM selects a random permutation of N objects.
!    * R4_TO_B4_IEEE converts an R4 to a 4 byte IEEE word.
!    * R4_TO_B4_ULTRIX converts an R4 to a 4 byte ULTRIX word.
!    * R4_TO_BINARY represents an R4 as a string of binary digits.
!    * R4_TO_BITS converts an R4 to a string of 32 bits.
!    * R4_TO_CH4 converts an R4 to a 4 character string.
!    * R4_TO_FLT computes the scientific representation of an R4.
!    * R4_TO_S_LEFT writes an R4 into a left justified character string.
!    * R4_TO_S_RIGHT writes an R4 into a right justified character string.
!    * R4_TO_S32 encodes an R4 as 32 characters.
!    * R4_TO_SEF represents an R4 as R = S * 2**E * F.
!    * R8_EXTRACT "extracts" an R8 from the beginning of a string.
!    * R8_INPUT prints a prompt string and reads an R8 from the user.
!    * R8_NEXT "reads" R8's from a string, one at a time.
!    * R8_PRINT prints a scalar, vector or array or type R8.
!    * R8_TO_BINARY represents an R8 as a string of binary digits.
!    * R8_TO_S_LEFT writes an R8 into a left justified string.
!    * R8_TO_S_LEFT writes an R8 into a right justified string.
!    * R8VEC_TO_S "writes" an R8VEC into a string.
!    * RANGER "understands" a range defined by a string like '4:8'.
!    * RAT_TO_S_LEFT returns a left-justified representation of IVAL/JVAL.
!    * RAT_TO_S_RIGHT returns a right-justified representation of IVAL/JVAL.
!    * RUBOUT deletes the pair "character" + Backspace from a string.
!    * S_ADJUSTL flushes a string left.
!    * S_ADJUSTR flushes a string right.
!    * S_AFTER_SS_COPY copies a string after a given substring.
!    * S_ALPHA_LAST returns the location of the last alphabetic character.
!    * S_ANY_ALPHA is TRUE if a string contains any alphabetic character.
!    * S_ANY_CONTROL is TRUE if a string contains any control characters.
!    * S_B2U replaces interword blanks by underscores.
!    * S_BEFORE_SS_COPY copies a string up to a given substring.
!    * S_BEGIN is TRUE if one string matches the beginning of the other.
!    * S_BEHEAD_SUBSTRING "beheads" a string, removing a given substring.
!    * S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!    * S_BLANKS_DELETE replaces consecutive blanks by one blank.
!    * S_BLANKS_INSERT inserts blanks into a string, sliding old characters over.
!    * S_CAP replaces any lowercase letters by uppercase ones in a string.
!    * S_CAT concatenates two strings to make a third string.
!    * S_CAT1 concatenates two strings, with a single blank separator.
!    * S_CENTER centers the non-blank portion of a string.
!    * S_CENTER_INSERT inserts one string into the center of another.
!    * S_CH_BLANK replaces each occurrence of a particular character by a blank.
!    * S_CH_DELETE removes all occurrences of a character from a string.
!    * S_CH_LAST returns the last nonblank character in a string.
!    * S_CHOP "chops out" a portion of a string, and closes up the hole.
!    * S_CONTROL_BLANK replaces control characters with blanks.
!    * S_CONTROL_COUNT returns the number of control characters in a string.
!    * S_CONTROL_DELETE removes all control characters from a string.
!    * S_COPY copies one string into another.
!    * S_DETAG removes from a string all substrings marked by angle brackets.
!    * S_EQI is a case insensitive comparison of two strings for equality.
!    * S_EQIDB compares two strings, ignoring case and blanks.
!    * S_ESCAPE_TEX de-escapes TeX escape sequences.
!    * S_FILL overwrites every character of a string by a given character.
!    * S_FIRST_NONBLANK returns the location of the first nonblank.
!    * S_GEI = ( S1 is lexically greater than or equal to S2 ).
!    * S_GTI = S1 is lexically greater than S2.
!    * S_INDEX seeks the first occurrence of a substring.
!    * S_INDEX_SET searches a string for any of a set of characters.
!    * S_INDEXI is a case-insensitive INDEX function.
!    * S_INDEX_LAST finds the LAST occurrence of a given substring.
!    * S_INDEX_LAST_C finds the LAST occurrence of a given character.
!    * S_I_APPEND appends an integer to a string.
!    * S_INC "increments" a string.
!    * S_INPUT prints a prompt string and reads a string from the user.
!    * S_IS_ALPHA returns TRUE if the string contains only alphabetic characters.
!    * S_IS_ALPHANUMERIC = string contains only alphanumeric characters.
!    * S_IS_DIGIT returns TRUE if a string contains only decimal digits.
!    * S_IS_F77_NAME = input string represent a legal FORTRAN-77 identifier.
!    * S_IS_F90_NAME = input string represent a legal FORTRAN 90 identifier.
!    * S_IS_I is TRUE if a string represents an integer.
!    * S_IS_R is TRUE if a string represents a real number.
!    * S_LEFT_INSERT inserts one string flush left into another.
!    * S_LEI = ( S1 is lexically less than or equal to S2 ).
!    * S_LOW replaces all uppercase letters by lowercase ones.
!    * S_LTI = ( S1 is lexically less than S2 ).
!    * S_NEQI compares two strings for non-equality, ignoring case.
!    * S_NO_CONTROL = string contains no control characters.
!    * S_OF_I converts an integer to a left-justified string.
!    * S_ONLY_ALPHAB checks if a string is only alphabetic and blanks.
!    * S_ONLY_DIGITB returns TRUE if the string contains only digits or blanks.
!    * S_OVERLAP determines the overlap between two strings.
!    * S_PAREN_CHECK checks the parentheses in a string.
!    * S_PLOT plots a character string onto a graphics aimage.
!    * S_R_APPEND appends a real number to a string.
!    * S_REPLACE_CH replaces all occurrences of one character by another.
!    * S_REPLACE_ONE replaces the first occurrence of SUB1 with SUB2.
!    * S_REPLACE_REC is a recursive replacement of one string by another.
!    * S_REPLACE replaces all occurrences of SUB1 by SUB2 in a string.
!    * S_REPLACE_I replaces all occurrences of SUB1 by SUB2 in a string.
!    * S_REVERSE reverses the characters in a string.
!    * S_RIGHT_INSERT inserts a string flush right into another.
!    * S_ROMAN_TO_I4 converts a Roman numeral to an integer.
!    * S_S_DELETE removes all occurrences of a substring from a string.
!    * S_S_DELETE2 recursively removes a substring from a string.
!    * S_S_INSERT inserts a substring into a string.
!    * S_SET_DELETE removes any characters in one string from another string.
!    * S_SHIFT_CIRCULAR circular shifts the characters in a string to the right.
!    * S_SHIFT_LEFT shifts the characters in a string to the left and blank pads.
!    * S_SHIFT_RIGHT shifts the characters in a string to the right and blank pads.
!    * S_SKIP_SET finds the first entry of a string that is NOT in a set.
!    * S_SPLIT divides a string into three parts, given the middle.
!    * S_SWAP swaps two strings.
!    * S_TAB_BLANK replaces each TAB character by one space.
!    * S_TAB_BLANKS replaces TAB characters by 6 spaces.
!    * S_TO_C4 reads a complex number from a string.
!    * S_TO_CHVEC converts a string to a character vector.
!    * S_TO_DATE converts the F90 date string to a more usual format.
!    * S_TO_DEC reads a number from a string, returning a decimal result.
!    * S_TO_EBCDIC converts a character string from ASCII to EBCDIC.
!    * S_TO_FORMAT reads a FORTRAN format from a string.
!    * S_TO_HEX replaces a character string by a hexadecimal representation.
!    * S_TO_I4 reads an integer value from a string.
!    * S_TO_I4VEC reads an integer vector from a string.
!    * S_TO_L reads a logical value from a string.
!    * S_TO_R4 reads an R4 value from a string.
!    * S_TO_R4VEC reads an R4VEC from a string.
!    * S_TO_R8 reads an R8 value from a string.
!    * S_TO_R8VEC reads an R8VEC from a string.
!    * S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.
!    * S_TO_SOUNDEX computes the Soundex code of a string.
!    * S_TO_W reads the next blank-delimited word from a string.
!    * S_TOKEN_EQUAL checks whether a string is equal to any of a set of strings.
!    * S_TOKEN_MATCH matches the beginning of a string and a set of tokens.
!    * S_TRIM_ZEROS removes trailing zeros from a string.
!    * S_U2B replaces underscores by blanks.
!    * S_WORD_APPEND appends a word to a string.
!    * S_WORD_CAP capitalizes the first character of each word in a string.
!    * S_WORD_COUNT counts the number of "words" in a string.
!    * S_WORD_EXTRACT_FIRST extracts the first word from a string.
!    * S_WORD_FIND finds the word of a given index in a string.
!    * S_WORD_INDEX finds the word of a given index in a string.
!    * S_WORD_NEXT "reads" words from a string, one at a time.
!    * S_WORD_PERMUTE permutes the words in a string.
!    * S32_TO_I4 returns an I4 equivalent to a 32 character string.
!    * S32_TO_R4 converts a 32-character variable into an R4.
!    * SEF_TO_B4_IEEE converts SEF information to a 4 byte IEEE real word.
!    * SEF_TO_B4_IEEE converts SEF information to a 4 byte IEEE real word.
!    * SEF_TO_R4 converts SEF information to an R4 = S * 2.0**E * F.
!    * SQZ_I4 uncompresses JVEC to extract integers stored there.
!    * STATE_ID returns the 2 letter Postal Code for one of the 50 states.
!    * STATE_NAME returns the name of one of the 50 states.
!    * SVEC_LAB makes an index array for an array of (repeated) strings.
!    * SVEC_MERGE_A merges two ascending sorted string arrays.
!    * SVEC_PERMUTE permutes a string vector in place.
!    * SVEC_REVERSE reverses the elements of a string vector.
!    * SVEC_SEARCH_BINARY_A searches an ascending sorted string vector.
!    * SVEC_SORT_HEAP_A ascending sorts an SVEC using heap sort.
!    * SVEC_SORT_HEAP_A_INDEX: case-sensitive indexed heap sort of an SVEC.
!    * SVEC_SORTED_UNIQUE: number of unique entries in a sorted SVEC.
!    * SVECI_SEARCH_BINARY_A searches an ascending sorted SVEC of implicitly capitalized strings.
!    * SVECI_SORT_HEAP_A heap sorts an SVEC of implicitly capitalized strings.
!    * SVECI_SORT_HEAP_A_INDEX index heap sorts an SVEC of implicitly capitalized strings.
!    * SYM_TO_CH returns the character represented by a symbol.
!    * TOKEN_EXPAND makes sure certain tokens have spaces surrounding them.
!    * TOKEN_EXTRACT "extracts" a token from the beginning of a string.
!    * TOKEN_INDEX finds the N-th FORTRAN variable name in a string.
!    * TOKEN_NEXT finds the next FORTRAN variable name in a string.
!    * UPPER returns an uppercase version of a string.
!    * WORD_BOUNDS returns the start and end of each word in a string.
!    * WORD_LAST_READ returns the last word from a string.
!    * WORD_NEXT finds the next (blank separated) word in a string.
!    * WORD_NEXT_READ "reads" words from a string, one at a time.
!    * WORD_NEXT2 returns the first word in a string.
!    * WORD_SWAP swaps two words in a given string.

function a_to_i4 ( ch )

  !*****************************************************************************80
  !
  !! A_TO_I4 returns the index of an alphabetic character.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Example:
  !
  !    CH  A_TO_I4
  !
  !    'A'   1
  !    'B'   2
  !    ...
  !    'Z'  26
  !    'a'  27
  !    'b'  28
  !    ...
  !    'z'  52
  !    '$'   0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 February 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, a character.
  !
  !    Output, integer ( kind = 4 ) A_TO_I4, is the alphabetic index of the
  !    character, between 1 and 26 if the character is a capital letter,
  !    between 27 and 52 if it is lower case, and 0 otherwise.
  !
  implicit none
  integer ( kind = 4 ) a_to_i4
  integer ( kind = 4 ), parameter :: cap_shift = 64
  character            ch
  integer ( kind = 4 ), parameter :: low_shift = 96
  if ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) then
     a_to_i4 = iachar ( ch ) - cap_shift
  else if ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) then
     a_to_i4 = iachar ( ch ) - low_shift + 26
  else
     a_to_i4 = 0
  end if
  return
end function a_to_i4


subroutine b4_ieee_to_r4 ( word, r )

  !*****************************************************************************80
  !
  !! B4_IEEE_TO_R4 converts a 4 byte IEEE word into an R4.
  !
  !  Discussion:
  !
  !    An "R4" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 4 )".
  !
  !    This routine does not seem to work reliably for unnormalized data.
  !
  !    The word containing the real value may be interpreted as:
  !
  !      /SEEEEEEE/EFFFFFFF/FFFFFFFF/FFFFFFFF/
  !
  !      /33222222/22222222/22222100/00000000/
  !      /10987654/32109876/54321098/76543210/  <-- Bit numbering
  !
  !    where
  !
  !      S is the sign bit,
  !      E are the exponent bits,
  !      F are the mantissa bits.
  !
  !    The mantissa is usually "normalized"; that is, there is an implicit
  !    leading 1 which is not stored.  However, if the exponent is set to
  !    its minimum value, this is no longer true.
  !
  !    The exponent is "biased".  That is, you must subtract a bias value
  !    from the exponent to get the true value.
  !
  !    If we read the three fields as integers S, E and F, then the
  !    value of the resulting real number R can be determined by:
  !
  !    * if E = 255
  !        if F is nonzero, then R = NaN;
  !        if F is zero and S is 1, R = -Inf;
  !        if F is zero and S is 0, R = +Inf;
  !    * else if 0 < E then R = (-1)**(S) * 2**(E-127) * (1 + (F/2**24))
  !    * else if E = 0
  !        if F is nonzero, R = (-1)**(S) * 2**(E-126) * (F/2**24)
  !        if F is zero and S is 1, R = -0;
  !        if F is zero and S is 0, R = +0;
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 November 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    IEEE Standards Committee 754,
  !    IEEE Standard for Binary Floating Point Arithmetic,
  !    ANSI/IEEE Standard 754-1985,
  !    SIGPLAN Notices,
  !    Volume 22, Number 2, 1987, pages 9-25.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) WORD, the word to be decoded.
  !
  !    Output, real ( kind = 4 ) R, the value of the real number.
  !
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) i
  real    ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) word
  !
  !  Read the fields.
  !
  s = 0
  call mvbits ( word, 31, 1, s, 0 )

  e = 0
  call mvbits ( word, 23, 8, e, 0 )

  f = 0
  call mvbits ( word, 0, 23, f, 0 )
  !
  !  Don't bother trying to return NaN or Inf just yet.
  !
  if ( e == 255 ) then
     r = 0.0E+00
  else if ( 0 < e ) then
     r = ( -1.0E+00 )**s * 2.0E+00**(e-127-23) * real ( 8388608 + f, kind = 4 )
  else if ( e == 0 ) then
     r = ( -1.0E+00 )**s * 2.0E+00**(-126) * real ( f, kind = 4 )
     do i = 1, 23
        r = r / 2.0E+00
     end do
  end if

  return
end subroutine b4_ieee_to_r4
subroutine b4_ieee_to_sef ( word, s, e, f )

  !*****************************************************************************80
  !
  !! B4_IEEE_TO_SEF converts an IEEE real word to S * 2**E * F format.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 November 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) WORD, a word containing an IEEE real number.
  !
  !    Output, integer ( kind = 4 ) S, the sign bit:
  !    0, if R is nonnegative;
  !    1, if R is negative.
  !
  !    Output, integer ( kind = 4 ) E, the exponent base 2.
  !
  !    Output, integer ( kind = 4 ) F, the mantissa.
  !
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) e2
  integer ( kind = 4 ) f
  integer ( kind = 4 ) s
  integer ( kind = 4 ) word

  s = 0
  call mvbits ( word, 31, 1, s, 0 )

  e2 = 0
  call mvbits ( word, 23, 8, e2, 0 )

  if ( e2 == 255 ) then

     e = 128

     call mvbits ( word, 0, 23, f, 0 )

     if ( f == 0 ) then
        f = 0
     else
        f = 2**23 - 1
     end if

  else if ( 0 < e2 ) then

     e = e2 - 127 - 23
     f = 2**23
     call mvbits ( word, 0, 23, f, 0 )

     do while ( mod ( f, 2 ) == 0 )
        f = f / 2
        e = e + 1
     end do

  else if ( e2 == 0 ) then

     e = e2 - 127 - 23
     f = 0
     call mvbits ( word, 0, 23, f, 0 )

     if ( f == 0 ) then

        e = 0

     else

        do while ( 0 < f .and. mod ( f, 2 ) == 0 )
           f = f / 2
           e = e + 1
        end do

     end if

  end if

  return
end subroutine b4_ieee_to_sef
subroutine base_to_i4 ( s, base, i )

  !*****************************************************************************80
  !
  !! BASE_TO_I4 returns the value of an integer represented in some base.
  !
  !  Discussion:
  !
  !    BASE = 1 is allowed, in which case we allow the digits '1' and '0',
  !    and we simply count the '1' digits for the result.
  !
  !    Negative bases between -16 and -2 are allowed.
  !
  !    The base -1 is allowed, and essentially does a parity check on
  !    a string of 1's.
  !
  !  Example:
  !
  !        Input      Output
  !    -------------  ------
  !         S   BASE       I
  !    ------  -----  ------
  !      '101'     2       5
  !    '-1000'     3     -27
  !      '100'     4      16
  !   '111111'     2      63
  !   '111111'    -2      21
  !   '111111'     1       6
  !   '111111'    -1       0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string.  The elements of S are
  !    blanks, a plus or minus sign, and digits.  Normally, the digits
  !    are representations of integers between 0 and |BASE-1|.  In the
  !    special case of base 1 or base -1, we allow both 0 and 1 as digits.
  !
  !    Input, integer ( kind = 4 ) BASE, the base in which the representation
  !    is given.  Normally, 2 <= BASE <= 16.  However, there are two exceptions.
  !
  !    Output, integer ( kind = 4 ) I, the integer.
  !
  implicit none

  integer   ( kind = 4 ) base
  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) isgn
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  integer   ( kind = 4 ) state

  i = 0
  s_length = len_trim ( s )

  if ( base == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'BASE_TO_I4 - Serious error!'
     write ( *, '(a)' ) '  The input base is zero.'
     i = -1
     return
  end if

  if ( 16 < abs ( base ) ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'BASE_TO_I4 - Serious error!'
     write ( *, '(a)' ) '  The input base is greater than 16!'
     i = -1
     return
  end if

  state = 0
  isgn = 1
  ichr = 1

  do while ( ichr <= s_length )

     c = s(ichr:ichr)
     !
     !  Blank.
     !
     if ( c == ' ' ) then

        if ( state == 2 ) then
           exit
        end if
        !
        !  Sign, + or -.
        !
     else if ( c == '-' ) then

        if ( state /= 0 ) then
           exit
        end if

        state = 1
        isgn = -1

     else if ( c == '+' ) then

        if ( state /= 0 ) then
           exit
        end if

        state = 1

     else
        !
        !  Digit?
        !
        call hex_digit_to_i4 ( c, idig )

        if ( abs ( base ) == 1 .and. ( idig == 0 .or. idig == 1 ) ) then

           i = base * i + idig
           state = 2

        else if ( 0 <= idig .and. idig < abs ( base ) ) then

           i = base * i + idig
           state = 2

        else

           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'BASE_TO_I4 - Serious error!'
           write ( *, '(a)' ) '  Illegal digit = "' // c // '"'
           write ( *, '(a)' ) '  Conversion halted prematurely!'
           return

        end if

     end if

     ichr = ichr + 1

  end do
  !
  !  Once we're done reading information, we expect to be in state 2.
  !
  if ( state /= 2 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'BASE_TO_I4 - Serious error!'
     write ( *, '(a)' ) '  Unable to decipher input!'
     return
  end if
  !
  !  Account for the sign.
  !
  i = isgn * i

  return
end subroutine base_to_i4
subroutine binary_to_i4 ( s, i )

  !*****************************************************************************80
  !
  !! BINARY_TO_I4 converts a binary representation into an integer value.
  !
  !  Example:
  !
  !        S        I
  !
  !      '101'      5
  !    '-1000'     -8
  !        '1'      1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the binary representation.
  !
  !    Output, integer ( kind = 4 ) I, the integer whose representation was input.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) isgn
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  integer   ( kind = 4 ) state

  s_length = len_trim ( s )

  i = 0
  ichr = 1
  state = 0
  isgn = 1

  do while ( ichr <= s_length )

     c = s(ichr:ichr)
     !
     !  Blank.
     !
     if ( c == ' ' ) then

        if ( state == 2 ) then
           state = 3
        end if
        !
        !  Sign, + or -.
        !
     else if ( c == '-' ) then

        if ( state == 0 ) then
           state = 1
           isgn = -1
        else
           state = -1
        end if

     else if ( c == '+' ) then

        if ( state == 0 ) then
           state = 1
        else
           state = -1
        end if
        !
        !  Digit, 0 or 1.
        !
     else if ( c == '1' ) then

        i = 2 * i
        i = i + 1
        state = 2

     else if ( c == '0' ) then

        i = 2 * i
        state = 2
        !
        !  Illegal or unknown sign.
        !
     else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BINARY_TO_I4 - Serious error!'
        write ( *, '(a)' ) '  Illegal digit = "' // c // '"'
        write ( *, '(a)' ) '  Conversion halted prematurely!'
        return

     end if

     if ( state == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BINARY_TO_I4 - Serious error!'
        write ( *, '(a)' ) '  Unable to decipher input!'
        return
     end if

     if ( 3 <= state ) then
        exit
     end if

     ichr = ichr + 1

  end do
  !
  !  Apply the sign.
  !
  i = isgn * i

  return
end subroutine binary_to_i4
subroutine binary_to_r4 ( s, r )

  !*****************************************************************************80
  !
  !! BINARY_TO_R4 converts a binary representation into an R4.
  !
  !  Discussion:
  !
  !    An "R4" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 4 )".
  !
  !  Example:
  !
  !        S                         R
  !
  !    -1010.11                   -10.75
  !    0.011011                     0.4218750
  !    0.01010101010101010101010    0.3333333
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the binary representation.
  !
  !    Output, real ( kind = 4 ) R, the real number.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) power
  real      ( kind = 4 ) r
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  integer   ( kind = 4 ) state

  s_length = len_trim ( s )

  intval = 0
  ichr = 1
  state = 0
  isgn = 1
  r = 0.0E+00
  power = 0

  do while ( ichr <= s_length )

     c = s(ichr:ichr)
     !
     !  Blank.
     !
     if ( c == ' ' ) then

        if ( state == 4 ) then
           state = 5
        end if
        !
        !  Sign, + or -.
        !
     else if ( c == '-' ) then

        if ( state == 0 ) then
           state = 1
           isgn = -1
        else
           state = -1
        end if

     else if ( c == '+' ) then

        if ( state == 0 ) then
           state = 1
        else
           state = -1
        end if
        !
        !  Digit, 0 or 1.
        !
     else if ( c == '1' ) then

        intval = 2 * intval + 1

        if ( state == 0 .or. state == 1 ) then
           state = 2
        else if ( state == 3 ) then
           state = 4
        end if

        if ( state == 4 ) then
           power = power + 1
        end if

     else if ( c == '0' ) then

        intval = 2 * intval

        if ( state == 0 .or. state == 1 ) then
           state = 2
        else if ( state == 3 ) then
           state = 4
        end if

        if ( state == 4 ) then
           power = power + 1
        end if
        !
        !  Decimal point.
        !
     else if ( c == '.' ) then

        if ( state <= 2 ) then
           state = 3
        else
           state = -1
        end if
        !
        !  Illegal or unknown sign.
        !
     else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BINARY_TO_R4 - Serious error!'
        write ( *, '(a)' ) '  Illegal character = "' // c // '"'
        write ( *, '(a)' ) '  Conversion halted prematurely!'
        stop

     end if

     if ( state == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BINARY_TO_R4 - Serious error!'
        write ( *, '(a)' ) '  Unable to decipher input!'
        stop
     end if

     if ( 5 <= state ) then
        exit
     end if

     ichr = ichr + 1

  end do
  !
  !  Apply the sign and the scale factor.
  !
  r = real ( isgn * intval, kind = 4 ) / 2.0E+00**power

  return
end subroutine binary_to_r4
subroutine binary_to_r8 ( s, r )

  !*****************************************************************************80
  !
  !! BINARY_TO_R8 converts a binary representation into an R8.
  !
  !  Discussion:
  !
  !    An "R8" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 8 )".
  !
  !  Example:
  !
  !        S                         R
  !
  !    -1010.11                   -10.75
  !    0.011011                     0.4218750
  !    0.01010101010101010101010    0.3333333
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 June 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the binary representation.
  !
  !    Output, real ( kind = 8 ) R, the real number.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) power
  real      ( kind = 8 ) r
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  integer   ( kind = 4 ) state

  s_length = len_trim ( s )

  intval = 0
  ichr = 1
  state = 0
  isgn = 1
  r = 0.0D+00
  power = 0

  do while ( ichr <= s_length )

     c = s(ichr:ichr)
     !
     !  Blank.
     !
     if ( c == ' ' ) then

        if ( state == 4 ) then
           state = 5
        end if
        !
        !  Sign, + or -.
        !
     else if ( c == '-' ) then

        if ( state == 0 ) then
           state = 1
           isgn = -1
        else
           state = -1
        end if

     else if ( c == '+' ) then

        if ( state == 0 ) then
           state = 1
        else
           state = -1
        end if
        !
        !  Digit, 0 or 1.
        !
     else if ( c == '1' ) then

        intval = 2 * intval + 1

        if ( state == 0 .or. state == 1 ) then
           state = 2
        else if ( state == 3 ) then
           state = 4
        end if

        if ( state == 4 ) then
           power = power + 1
        end if

     else if ( c == '0' ) then

        intval = 2 * intval

        if ( state == 0 .or. state == 1 ) then
           state = 2
        else if ( state == 3 ) then
           state = 4
        end if

        if ( state == 4 ) then
           power = power + 1
        end if
        !
        !  Decimal point.
        !
     else if ( c == '.' ) then

        if ( state <= 2 ) then
           state = 3
        else
           state = -1
        end if
        !
        !  Illegal or unknown sign.
        !
     else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BINARY_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Illegal character = "' // c // '"'
        write ( *, '(a)' ) '  Conversion halted prematurely!'
        stop

     end if

     if ( state == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BINARY_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Unable to decipher input!'
        stop
     end if

     if ( 5 <= state ) then
        exit
     end if

     ichr = ichr + 1

  end do
  !
  !  Apply the sign and the scale factor.
  !
  r = real ( isgn * intval, kind = 8 ) / 2.0D+00**power

  return
end subroutine binary_to_r8
subroutine ch_cap ( ch )

  !*****************************************************************************80
  !
  !! CH_CAP capitalizes a single character.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
  !    which guarantee the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character CH, the character to capitalize.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
     ch = achar ( itemp - 32 )
  end if

  return
end subroutine ch_cap
subroutine ch_count_chvec_add ( n, chvec, count )

  !*****************************************************************************80
  !
  !! CH_COUNT_CHVEC_ADD adds a character vector to a character count.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
  !
  !    Input, character CHVEC(N), a vector of characters.
  !
  !    Input/output, integer ( kind = 4 ) COUNT(0:255), the character counts.
  !
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) count(0:255)
  character              chvec(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j

  do i = 1, n
     j = iachar ( chvec(i) )
     count(j) = count(j) + 1
  end do

  return
end subroutine ch_count_chvec_add
subroutine ch_count_file_add ( file_name, count )

  !*****************************************************************************80
  !
  !! CH_COUNT_FILE_ADD adds characters in a file to a character count.
  !
  !  Discussion:
  !
  !    Each line is counted up to the last nonblank.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) FILE_NAME, the name of the file to examine.
  !
  !    Output, integer ( kind = 4 ) COUNT(0:255), the character counts.
  !
  implicit none

  integer   ( kind = 4 )  count(0:255)
  character ( len = * )   file_name
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  iunit
  character ( len = 255 ) line
  !
  !  Open the file.
  !
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CH_COUNT_FILE_ADD - Fatal error!'
     write ( *, '(a)' ) '  Could not open the file:'
     write ( *, '(a)' ) '    ' // trim ( file_name )
     return
  end if

  do

     read ( iunit, '(a)', iostat = ios ) line

     if ( ios /= 0 ) then
        exit
     end if

     call ch_count_s_add ( trim ( line ), count )

  end do

  close ( unit = iunit )

  return
end subroutine ch_count_file_add
subroutine ch_count_histogram_print ( count, title )

  !*****************************************************************************80
  !
  !! CH_COUNT_HISTOGRAM_PRINT prints a histogram of a set of character counts.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) COUNT(0:255), the character counts.
  !
  !    Input, character ( len = * ) TITLE, a title to be printed.
  !
  implicit none

  character              c
  character ( len = 4 )  ch4(0:255)
  integer   ( kind = 4 ) count(0:255)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) percent
  integer   ( kind = 4 ) row
  character ( len = 4 )  s(0:255)
  character ( len = * )  title
  integer   ( kind = 4 ) total

  total = sum ( count )

  do i = 0, 255
     c = achar ( i )
     call ch_to_sym ( c, ch4(i) )
  end do

  do i = 0, 255

     if ( total == 0 ) then
        percent = 0
     else
        percent = nint ( real ( 100 * count(i), kind = 4 ) &
             / real ( total, kind = 4 ) )
     end if

     if ( percent == 0 ) then
        s(i) = '   .'
     else
        write ( s(i), '(i4)' ) percent
     end if

  end do

  if ( 0 < len_trim ( title ) ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Character Histogram (Percentages).'
  write ( *, '(a)' ) ' '

  do row = 1, 16
     ilo = ( row - 1 ) * 16
     ihi =   row       * 16 - 1
     write ( *, '(2x,i3,a4,i3,3x,16a4)' ) ilo, ' to ', ihi, ch4(ilo:ihi)
     write ( *, '(12x,16a4)' ) s(ilo:ihi)
  end do

  return
end subroutine ch_count_histogram_print
subroutine ch_count_init ( count )

  !*****************************************************************************80
  !
  !! CH_COUNT_INIT initializes a character count.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer ( kind = 4 ) COUNT(0:255), the character counts.
  !
  implicit none

  integer ( kind = 4 ) count(0:255)

  count(0:255) = 0

  return
end subroutine ch_count_init
subroutine ch_count_print ( count, title )

  !*****************************************************************************80
  !
  !! CH_COUNT_PRINT prints a set of character counts.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) COUNT(0:255), the character counts.
  !
  !    Input, character ( len = * ) TITLE, a title to be printed.
  !
  implicit none

  character              c
  character ( len = 4 )  ch4(0:255)
  integer   ( kind = 4 ) count(0:255)
  integer   ( kind = 4 ) i
  real      ( kind = 4 ) percent
  character ( len = * )  title
  integer   ( kind = 4 ) total

  total = sum ( count )

  do i = 0, 255
     c = achar ( i )
     call ch_to_sym ( c, ch4(i) )
  end do

  if ( 0 < len_trim ( title ) ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Char  Count  Percentages.'
  write ( *, '(a)' ) ' '
  do i = 0, 255
     if ( 0 < count(i) ) then
        if ( total == 0 ) then
           percent = 0.0E+00
        else
           percent = real ( 100 * count(i), kind = 4 ) / real ( total, kind = 4 )
        end if
        write ( *, '(2x,a4,2x,i8,2x,f6.3)' ) ch4(i), count(i), percent
     end if
  end do

  return
end subroutine ch_count_print
subroutine ch_count_s_add ( s, count )

  !*****************************************************************************80
  !
  !! CH_COUNT_S_ADD adds a character string to a character histogram.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string to be examined.
  !
  !    Input/output, integer ( kind = 4 ) COUNT(0:255), the character counts.
  !
  implicit none

  integer   ( kind = 4 ) count(0:255)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  character ( len = * )  s

  do i = 1, len ( s )
     j = iachar ( s(i:i) )
     count(j) = count(j) + 1
  end do

  return
end subroutine ch_count_s_add
function ch_eqi ( c1, c2 )

  !*****************************************************************************80
  !
  !! CH_EQI is a case insensitive comparison of two characters for equality.
  !
  !  Discussion:
  !
  !    CH_EQI ( 'A', 'a' ) is TRUE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character C1, C2, the characters to compare.
  !
  !    Output, logical CH_EQI, the result of the comparison.
  !
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical   ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
     ch_eqi = .true.
  else
     ch_eqi = .false.
  end if

  return
end function ch_eqi
subroutine ch_extract ( s, ch )

  !*****************************************************************************80
  !
  !! CH_EXTRACT extracts the next nonblank character from a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string.  On output, the
  !    first nonblank character has been removed, and the string
  !    has been shifted left.
  !
  !    Output, character CH, the leading character of the string.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) get
  character ( len = * )  s
  integer   ( kind = 4 ) s_len

  s_len = len_trim ( s )
  ch = ' '

  get = 1

  do while ( get <= s_len )

     if ( s(get:get) /= ' ' ) then
        ch = s(get:get)
        call s_shift_left ( s, get )
        exit
     end if

     get = get + 1

  end do

  return
end subroutine ch_extract
function ch_index_first ( s, ch )

  !*****************************************************************************80
  !
  !! CH_INDEX_FIRST is the first occurrence of a character in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character CH, the character to be searched for.
  !
  !    Output, integer ( kind = 4 ) CH_INDEX_FIRST, the location of the first
  !    occurrence of the character in the string, or -1 if it does not occur.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) ch_index_first
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  ch_index_first = - 1
  s_length = len_trim ( s )

  do i = 1, s_length

     if ( s(i:i) == ch ) then
        ch_index_first = i
        return
     end if

  end do

  return
end function ch_index_first
function ch_index_last ( s, ch )

  !*****************************************************************************80
  !
  !! CH_INDEX_LAST is the last occurrence of a character in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character CH, the character to be searched for.
  !
  !    Output, integer ( kind = 4 ) CH_INDEX_LAST, the location of the last
  !    occurrence of the character in the string, or -1 if it does not occur.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) ch_index_last
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  ch_index_last = -1
  s_length = len_trim ( s )

  do i = s_length, 1, -1

     if ( s(i:i) == ch ) then
        ch_index_last = i
        return
     end if

  end do

  return
end function ch_index_last
function ch_indexi ( s, ch )

  !*****************************************************************************80
  !
  !! CH_INDEXI: (case insensitive) first occurrence of a character in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character CH, the character to be searched for.
  !
  !    Output, integer ( kind = 4 ) CH_INDEXI, the location of the first
  !    occurrence of the character (upper or lowercase), or -1 if it does
  !    not occur.
  !
  implicit none

  character              ch
  logical                ch_eqi
  integer   ( kind = 4 ) ch_indexi
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  ch_indexi = -1
  s_length = len_trim ( s )

  do i = 1, s_length

     if ( ch_eqi ( s(i:i), ch ) ) then
        ch_indexi = i
        return
     end if

  end do

  return
end function ch_indexi
function ch_is_alpha ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_ALPHA is TRUE if CH is an alphabetic character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, a character to check.
  !
  !    Output, logical CH_IS_ALPHA is TRUE if CH is an alphabetic character.
  !
  implicit none

  character ch
  logical ch_is_alpha

  if ( ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) .or. &
       ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) ) then
     ch_is_alpha = .true.
  else
     ch_is_alpha = .false.
  end if

  return
end function ch_is_alpha
function ch_is_alphanumeric ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_ALPHANUMERIC is TRUE if CH is alphanumeric.
  !
  !  Discussion:
  !
  !    Alphanumeric characters are 'A' through 'Z', 'a' through 'z' and
  !    '0' through '9'.
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be checked.
  !
  !    Output, logical CH_IS_ALPHANUMERIC, is TRUE if the character is
  !    alphabetic or numeric.
  !
  implicit none

  character              ch
  logical                ch_is_alphanumeric
  integer   ( kind = 4 ) i

  i = iachar ( ch )

  if ( ( 65 <= i .and. i <= 90 ) .or. &
       ( 97 <= i .and. i <= 122 ) .or. &
       ( 48 <= i .and. i <= 57 ) ) then

     ch_is_alphanumeric = .true.

  else

     ch_is_alphanumeric = .false.

  end if

  return
end function ch_is_alphanumeric
function ch_is_control ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_CONTROL is TRUE if a character is a control character.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be tested.
  !
  !    Output, logical CH_IS_CONTROL, TRUE if the character is a control
  !    character, and FALSE otherwise.
  !
  implicit none

  character ch
  logical   ch_is_control

  if ( iachar ( ch ) <= 31 .or. 127 <= iachar ( ch ) ) then
     ch_is_control = .true.
  else
     ch_is_control = .false.
  end if

  return
end function ch_is_control
function ch_is_digit ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_DIGIT is TRUE if a character is a decimal digit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be analyzed.
  !
  !    Output, logical CH_IS_DIGIT, is TRUE if the character is a digit.
  !
  implicit none

  character ch
  logical   ch_is_digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
     ch_is_digit = .true.
  else
     ch_is_digit = .false.
  end if

  return
end function ch_is_digit
function ch_is_format_code ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_FORMAT_CODE is TRUE if a character is a FORTRAN format code.
  !
  !  Discussion:
  !
  !    The format codes accepted here are not the only legal format
  !    codes in FORTRAN90.  However, they are more than sufficient
  !    for my needs!
  !
  !  Table:
  !
  !    A  Character
  !    B  Binary digits
  !    D  Real number, exponential representation
  !    E  Real number, exponential representation
  !    F  Real number, fixed point
  !    G  General format
  !    I  Integer
  !    L  Logical variable
  !    O  Octal digits
  !    Z  Hexadecimal digits
  !    *  Free format
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 November 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be analyzed.
  !
  !    Output, logical CH_IS_FORMAT_CODE, is TRUE if the character is a
  !    FORTRAN format code.
  !
  implicit none

  character ch
  logical   ch_eqi
  logical   ch_is_format_code

  ch_is_format_code = .true.

  if ( ch_eqi ( ch, 'A' ) ) then
     return
  else if ( ch_eqi ( ch, 'B' ) ) then
     return
  else if ( ch_eqi ( ch, 'D' ) ) then
     return
  else if ( ch_eqi ( ch, 'E' ) ) then
     return
  else if ( ch_eqi ( ch, 'F' ) ) then
     return
  else if ( ch_eqi ( ch, 'G' ) ) then
     return
  else if ( ch_eqi ( ch, 'I' ) ) then
     return
  else if ( ch_eqi ( ch, 'L' ) ) then
     return
  else if ( ch_eqi ( ch, 'O' ) ) then
     return
  else if ( ch_eqi ( ch, 'Z' ) ) then
     return
  else if ( ch == '*' ) then
     return
  end if

  ch_is_format_code = .false.

  return
end function ch_is_format_code
function ch_is_lower ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_LOWER is TRUE if a character is a lower case letter.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 May 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be analyzed.
  !
  !    Output, logical CH_IS_LOWER, is TRUE if the character is a lower
  !    case letter.
  !
  implicit none

  character ch
  logical   ch_is_lower

  if ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) then
     ch_is_lower = .true.
  else
     ch_is_lower = .false.
  end if

  return
end function ch_is_lower
function ch_is_printable ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_PRINTABLE is TRUE if C is printable.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 July 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, a character to check.
  !
  !    Output, logical CH_IS_PRINTABLE is TRUE if C is a printable character.
  !
  implicit none

  character              ch
  logical                ch_is_printable
  integer   ( kind = 4 ) i

  i = iachar ( ch )

  if ( 32 <= i .and. i <= 126 ) then
     ch_is_printable = .true.
  else
     ch_is_printable = .false.
  end if

  return
end function ch_is_printable
function ch_is_space ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_SPACE is TRUE if a character is a whitespace character.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    A whitespace character is a space, a form feed, a newline,
  !    a carriage return, a tab, or a vertical tab.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 October 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, a character to check.
  !
  !    Output, logical CH_IS_SPACE is TRUE if the character is a whitespace
  !    character.
  !
  implicit none

  character ch
  logical ch_is_space

  if ( ch == ' ' ) then
     ch_is_space = .true.
  else if ( ch == achar ( 12 ) ) then
     ch_is_space = .true.
  else if ( ch == achar ( 10 ) ) then
     ch_is_space = .true.
  else if ( ch == achar ( 13 ) ) then
     ch_is_space = .true.
  else if ( ch == achar ( 9 ) ) then
     ch_is_space = .true.
  else if ( ch == achar ( 11 ) ) then
     ch_is_space = .true.
  else
     ch_is_space = .false.
  end if

  return
end function ch_is_space
function ch_is_upper ( ch )

  !*****************************************************************************80
  !
  !! CH_IS_UPPER is TRUE if CH is an upper case letter.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 May 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be analyzed.
  !
  !    Output, logical CH_IS_UPPER, is TRUE if CH is an upper case letter.
  !
  implicit none

  character ch
  logical ch_is_upper

  if ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) then
     ch_is_upper = .true.
  else
     ch_is_upper = .false.
  end if

  return
end function ch_is_upper
subroutine ch_low ( ch )

  !*****************************************************************************80
  !
  !! CH_LOW lowercases a single character.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
  !    which guarantee the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character CH, the character to be lowercased.
  !
  implicit none

  character ch
  integer ( kind = 4 ) i

  i = iachar ( ch )

  if ( 65 <= i .and. i <= 90 ) then
     ch = achar ( i + 32 )
  end if

  return
end subroutine ch_low
subroutine ch_next ( s, ch, done )

  !*****************************************************************************80
  !
  !! CH_NEXT reads the next character from a string, ignoring blanks and commas.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = ' A  B, C    DE  F'
  !
  !    Output:
  !
  !      'A', 'B', 'C', 'D', 'E', 'F', and then blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string of characters.  Blanks and
  !    commas are considered insignificant.
  !
  !    Output, character CH.  If DONE is FALSE, then the
  !    "next" character.  If DONE is TRUE, then a blank.
  !
  !    Input/output, logical DONE.
  !    On input with a fresh value of S, the user should set
  !    DONE to TRUE.
  !    On output, the routine sets DONE to FALSE if another character
  !    was read, or TRUE if no more characters could be read.
  !
  implicit none

  character ch
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  if ( done ) then
     next = 1
     done = .false.
  end if

  s_length = len_trim ( s )

  do i = next, s_length

     if ( s(i:i) /= ' ' .and. s(i:i) /= ',' ) then
        ch = s(i:i)
        next = i + 1
        return
     end if

  end do

  done = .true.
  next = 1
  ch = ' '

  return
end subroutine ch_next
function ch_not_control ( ch )

  !*****************************************************************************80
  !
  !! CH_NOT_CONTROL = CH is NOT a control character.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 January 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH is the character to be tested.
  !
  !    Output, logical CH_NOT_CONTROL, TRUE if CH is not a control character,
  !    and FALSE otherwise.
  !
  implicit none

  character ch
  logical ch_not_control

  if ( iachar ( ch ) <= 31 .or. 128 <= iachar ( ch ) ) then
     ch_not_control = .true.
  else
     ch_not_control = .false.
  end if

  return
end function ch_not_control
function ch_roman_to_i4 ( ch )

  !*****************************************************************************80
  !
  !! CH_ROMAN_TO_I4 returns the integer value of a single Roman digit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, a Roman digit.
  !
  !    Output, integer ( kind = 4 ) CH_ROMAN_TO_I4, the value of the Roman
  !    numeral.  If the Roman numeral was not recognized, 0 is returned.
  !
  implicit none

  character ch
  integer ( kind = 4 ) ch_roman_to_i4
  integer ( kind = 4 ) i

  if ( ch == 'M' .or. ch == 'm' ) then
     i = 1000
  else if ( ch == 'D' .or. ch == 'd' ) then
     i = 500
  else if ( ch == 'C' .or. ch == 'c' ) then
     i = 100
  else if ( ch == 'L' .or. ch == 'l' ) then
     i = 50
  else if ( ch == 'X' .or. ch == 'x' ) then
     i = 10
  else if ( ch == 'V' .or. ch == 'v' ) then
     i = 5
  else if ( ch == 'I' .or. ch == 'i' .or. &
       ch == 'J' .or. ch == 'j' ) then
     i = 1
  else
     i = 0
  end if

  ch_roman_to_i4 = i

  return
end function ch_roman_to_i4
function ch_scrabble ( tile )

  !*****************************************************************************80
  !
  !! CH_SCRABBLE returns the character on a given Scrabble tile.
  !
  !  Discussion:
  !
  !    The tiles are numbered 1 to 100, and are labeled 'A' through 'Z',
  !    plus two blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) TILE, the index of the desired Scrabble tile.
  !
  !    Output, character CH_SCRABBLE, the character on the given tile.
  !
  implicit none

  character ch_scrabble
  character, dimension ( 1 : 100 ) :: scrabble = (/ &
       'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B', &
       'B', 'C', 'C', 'D', 'D', 'D', 'D', 'E', 'E', 'E', &
       'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'F', &
       'F', 'G', 'G', 'G', 'H', 'H', 'I', 'I', 'I', 'I', &
       'I', 'I', 'I', 'I', 'I', 'J', 'K', 'L', 'L', 'L', &
       'L', 'M', 'M', 'N', 'N', 'N', 'N', 'N', 'N', 'O', &
       'O', 'O', 'O', 'O', 'O', 'O', 'O', 'P', 'P', 'Q', &
       'R', 'R', 'R', 'R', 'R', 'R', 'S', 'S', 'S', 'S', &
       'T', 'T', 'T', 'T', 'T', 'T', 'U', 'U', 'U', 'U', &
       'V', 'V', 'W', 'W', 'X', 'X', 'Y', 'Z', ' ', ' ' /)
  integer ( kind = 4 ) tile

  if ( 1 <= tile .and. tile <= 100 ) then
     ch_scrabble = scrabble(tile)
  else
     ch_scrabble = '?'
  end if

  return
end function ch_scrabble
function ch_scrabble_frequency ( ch )

  !*****************************************************************************80
  !
  !! CH_SCRABBLE_FREQUENCY returns the Scrabble frequency of a character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character.
  !
  !    Output, integer ( kind = 4 ) CH_SCRABBLE_FREQUENCY, the frequency of
  !    the character.
  !
  implicit none

  character ch
  integer ( kind = 4 ) ch_scrabble_frequency
  integer ( kind = 4 ) ch_to_scrabble
  integer ( kind = 4 ), dimension ( 27 ) :: frequency = (/ &
       9,  2,  2,  4, 12, &
       2,  3,  2,  9,  1, &
       1,  4,  2,  6,  8, &
       2,  1,  6,  4,  6, &
       4,  2,  2,  1,  2, &
       1,  2 /)
  integer ( kind = 4 ) ic
  !
  !  Convert character to a Scrabble character index.
  !
  ic = ch_to_scrabble ( ch )

  if ( 1 <= ic .and. ic <= 27 ) then
     ch_scrabble_frequency = frequency(ic)
  else
     ch_scrabble_frequency = 0
  end if

  return
end function ch_scrabble_frequency
function ch_scrabble_points ( ch )

  !*****************************************************************************80
  !
  !! CH_SCRABBLE_POINTS returns the Scrabble point value of a character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character.
  !
  !    Output, integer ( kind = 4 ) CH_SCRABBLE_POINTS, the point value of
  !    the character.
  !
  implicit none

  character ch
  integer ( kind = 4 ) ch_scrabble_points
  integer ( kind = 4 ) ch_to_scrabble
  integer ( kind = 4 ) ic
  integer ( kind = 4 ), dimension ( 27 ) :: points = (/ &
       1,  3,  3,  2,  1, &
       4,  2,  4,  1,  8, &
       5,  1,  3,  1,  1, &
       3, 10,  1,  1,  1, &
       1,  4,  4,  8,  4, &
       10,  0 /)
  !
  !  Convert character to a Scrabble character index.
  !
  ic = ch_to_scrabble ( ch )

  if ( 1 <= ic .and. ic <= 27 ) then
     ch_scrabble_points = points(ic)
  else
     ch_scrabble_points = 0
  end if

  return
end function ch_scrabble_points
function ch_scrabble_select ( seed )

  !*****************************************************************************80
  !
  !! CH_SCRABBLE_SELECT selects a character with the Scrabble probability.
  !
  !  Discussion:
  !
  !    There are 100 Scrabble tiles, including two blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) SEED, a seed for the random
  !    number generator.
  !
  !    Output, character CH_SCRABBLE_SELECT, the character on a randomly
  !    chosen Scrabble tile.
  !
  implicit none

  character              ch_scrabble
  character              ch_scrabble_select
  integer   ( kind = 4 ) i4_uniform
  integer   ( kind = 4 ) seed
  integer   ( kind = 4 ) tile
  !
  !  Choose a tile between 1 and 100.
  !
  tile = i4_uniform ( 1, 100, seed )
  !
  !  Retrieve the character on that tile.
  !
  ch_scrabble_select = ch_scrabble ( tile )

  return
end function ch_scrabble_select
subroutine ch_swap ( ch1, ch2 )

  !*****************************************************************************80
  !
  !! CH_SWAP swaps two characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character CH1, CH2.  On output, the values
  !    have been interchanged.
  !
  implicit none

  character ch1
  character ch2
  character ch3

  ch3 = ch1
  ch1 = ch2
  ch2 = ch3

  return
end subroutine ch_swap
subroutine ch_to_amino_name ( ch, amino_name )

  !*****************************************************************************80
  !
  !! CH_TO_AMINO_NAME converts a character to an amino acid name.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carl Branden, John Tooze,
  !    Introduction to Protein Structure,
  !    Garland Publishing, 1991.
  !
  !  Parameters:
  !
  !    Input, character CH, the one letter code for an amino acid.
  !    Lower and upper case letters are treated the same.
  !
  !    Output, character ( len = * ) AMINO_NAME, the full name of the
  !    corresponding amino acid.  The longest name is 27 characters.
  !    If the input code is not recognized, then AMINO_NAME will be set to '???'.
  !
  implicit none

  integer ( kind = 4 ), parameter :: n = 23

  character ( len = * ) amino_name
  character ( len = 27 ), dimension ( n ) :: amino_table = (/ &
       'Alanine                    ', &
       'Aspartic acid or Asparagine', &
       'Cysteine                   ', &
       'Aspartic acid              ', &
       'Glutamic acid              ', &
       'Phenylalanine              ', &
       'Glycine                    ', &
       'Histidine                  ', &
       'Isoleucine                 ', &
       'Lysine                     ', &
       'Leucine                    ', &
       'Methionine                 ', &
       'Asparagine                 ', &
       'Proline                    ', &
       'Glutamine                  ', &
       'Arginine                   ', &
       'Serine                     ', &
       'Threonine                  ', &
       'Valine                     ', &
       'Tryptophan                 ', &
       'Undetermined amino acid    ', &
       'Tyrosine                   ', &
       'Glutamic acid or Glutamine ' /)
  character ch
  logical ch_eqi
  character, dimension ( n ) :: ch_table = (/ &
       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
       'X', 'Y', 'Z' /)
  integer ( kind = 4 ) i

  do i = 1, n
     if ( ch_eqi ( ch, ch_table(i) ) ) then
        amino_name = amino_table(i)
        return
     end if
  end do

  amino_name = '???'

  return
end subroutine ch_to_amino_name
subroutine ch_to_braille ( ch, ncol, braille )

  !*****************************************************************************80
  !
  !! CH_TO_BRAILLE converts an ASCII character to a Braille character string.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 August 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the ASCII character.
  !
  !    Output, integer ( kind = 4 ) NCOL, the number of columns used to represent
  !    the character.
  !
  !    Output, character ( len = 6 ) BRAILLE(3), contains, in rows 1
  !    through 3 and character columns 1 through NCOL, either a '*' or a ' '.
  !
  implicit none

  integer ( kind = 4 ), parameter :: num_symbol = 37

  character ( len = 6 )  braille(3)
  character              ch
  logical                ch_is_digit
  logical                ch_is_upper
  integer   ( kind = 4 ) iascii
  integer   ( kind = 4 ) ic_to_ibraille
  integer   ( kind = 4 ) ibraille
  integer   ( kind = 4 ) ncol
  !
  !       space       Aa1        Bb2        Cc3        Dd4
  !       Ee5         Ff6        Gg7        Hh8        Ii9
  !       Jj0         Kk         Ll         Mm         Nn
  !       Oo          Pp         Qq         Rr         Ss
  !       Tt          Uu         Vv         Ww         Xx
  !       Yy          Zz         &          ,          ;
  !       :           .          !          ()         "?
  !       '           -
  !
  character ( len = 6 ), parameter, dimension ( num_symbol ) :: symbol = (/ &
       '      ',  '*     ',  '* *   ',  '**    ',  '** *  ', &
       '*  *  ',  '***   ',  '****  ',  '* **  ',  ' **   ', &
       ' ***  ',  '*   * ',  '* * * ',  '**  * ',  '** ** ', &
       '*  ** ',  '*** * ',  '***** ',  '* *** ',  ' ** * ', &
       ' **** ',  '*   **',  '* * **',  ' *** *',  '**  **', &
       '** ***',  '*  ***',  '*** **',  '  *   ',  '  * * ', &
       '  **  ',  '  ** *',  '  *** ',  '  ****',  '  * **', &
       '    * ',  '    **' /)

  ncol = 0

  braille(1)(1:6) = ' '
  braille(2)(1:6) = ' '
  braille(3)(1:6) = ' '
  !
  !  A space is treated specially.
  !
  if ( ch == ' ' ) then

     braille(1)(1:2) = '  '
     braille(2)(1:2) = '  '
     braille(3)(1:2) = '  '
     ncol = 2
     return

  end if
  !
  !  Get the ASCII numeric code of the character.
  !
  iascii = iachar ( ch )
  !
  !  Get the index of the Braille equivalent.
  !
  ibraille = ic_to_ibraille ( iascii )

  if ( 0 <= ibraille ) then
     !
     !  Upper case characters are preceded by a special mark.
     !
     if ( ch_is_upper ( ch ) ) then

        braille(1)(1:3) = '   '
        braille(2)(1:3) = '   '
        braille(3)(1:3) = ' * '

        ncol = 3
        !
        !  Digits are preceded by a special mark.
        !
     else if ( ch_is_digit ( ch ) ) then

        braille(1)(1:3) = ' * '
        braille(2)(1:3) = ' * '
        braille(3)(1:3) = '** '

        ncol = 3

     end if

     braille(1)(ncol+1:ncol+2) = symbol(ibraille)(1:2)
     braille(2)(ncol+1:ncol+2) = symbol(ibraille)(3:4)
     braille(3)(ncol+1:ncol+2) = symbol(ibraille)(5:6)

     ncol = ncol + 2
     !
     !  Add a trailing "half space".
     !
     braille(1)(ncol+1:ncol+1) = ' '
     braille(2)(ncol+1:ncol+1) = ' '
     braille(3)(ncol+1:ncol+1) = ' '

     ncol = ncol + 1

  end if

  return
end subroutine ch_to_braille
subroutine ch_to_ch3_amino ( ch, ch3 )

  !*****************************************************************************80
  !
  !! CH_TO_CH3_AMINO converts a 1 character to a 3 character code for amino acids.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carl Branden, John Tooze,
  !    Introduction to Protein Structure,
  !    Garland Publishing, 1991.
  !
  !  Parameters:
  !
  !    Input, character CH, the one letter code for an amino acid.
  !    Lower and upper case letters are treated the same.
  !
  !    Output, character ( len = 3 ) CH3, the three letter code for the
  !    amino acid.  If the input code is not recognized, then CH3 will be '???'.
  !
  implicit none

  integer ( kind = 4 ), parameter :: n = 23

  character ch
  logical ch_eqi
  character, parameter, dimension ( n ) :: ch_table = (/ &
       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
       'X', 'Y', 'Z' /)
  character ( len = 3 ) ch3
  character ( len = 3 ), parameter, dimension ( n ) :: ch3_table = (/ &
       'Ala', 'Asx', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ise', 'Lys', &
       'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Val', 'Trp', &
       'X  ', 'Tyr', 'Glx' /)
  integer ( kind = 4 ) i

  do i = 1, n
     if ( ch_eqi ( ch, ch_table(i) ) ) then
        ch3 = ch3_table(i)
        return
     end if
  end do

  ch3 = '???'

  return
end subroutine ch_to_ch3_amino
subroutine ch_to_digit ( ch, digit )

  !*****************************************************************************80
  !
  !! CH_TO_DIGIT returns the integer value of a base 10 digit.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Example:
  !
  !     CH  DIGIT
  !    ---  -----
  !    '0'    0
  !    '1'    1
  !    ...  ...
  !    '9'    9
  !    ' '    0
  !    'X'   -1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the decimal digit, '0' through '9' or blank
  !    are legal.
  !
  !    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
  !    If CH was 'illegal', then DIGIT is -1.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

     digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

     digit = 0

  else

     digit = - 1

  end if

  return
end subroutine ch_to_digit
subroutine ch_to_digit_bin ( ch, digit )

  !*****************************************************************************80
  !
  !! CH_TO_DIGIT_BIN returns the integer value of a binary digit.
  !
  !  Discussion:
  !
  !    This routine handles other traditional binary pairs of "digits"
  !    besides '0' and '1'.
  !
  !  Example:
  !
  !     CH  DIGIT
  !    ---  -----
  !    '0'    0
  !    '1'    1
  !    'T'    1
  !    'F'    0
  !    'Y'    1
  !    'N'    0
  !    '+'    1
  !    '-'    0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the binary digit.
  !
  !    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
  !    If CH was 'illegal', then DIGIT is -1.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  if ( ch == '0' .or. &
       ch == 'F' .or. &
       ch == 'f' .or. &
       ch == '-' .or. &
       ch == 'N' .or. &
       ch == 'n' ) then

     digit = 0

  else if ( ch == '1' .or. &
       ch == 'T' .or. &
       ch == 't' .or. &
       ch == '+' .or. &
       ch == 'Y' .or. &
       ch == 'y' ) then

     digit = 1

  else

     digit = -1

  end if

  return
end subroutine ch_to_digit_bin
subroutine ch_to_digit_oct ( ch, i )

  !*****************************************************************************80
  !
  !! CH_TO_DIGIT_OCT returns the integer value of an octal digit.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the octal digit, '0' through '7'.
  !
  !    Output, integer ( kind = 4 ) I, the corresponding integer value, or
  !    -1 if CH was illegal.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i

  i = iachar ( ch )

  if ( lle ( '0', ch ) .and. lle ( ch, '7' ) ) then

     i = i - 48

  else if ( ch == ' ' ) then

     i = 0

  else

     i = -1

  end if

  return
end subroutine ch_to_digit_oct
function ch_to_ebcdic ( ch )

  !*****************************************************************************80
  !
  !! CH_TO_EBCDIC converts a character to EBCDIC.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, which
  !    guarantee the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the ASCII character.
  !
  !    Output, character CH_TO_EBCDIC, the corresponding EBCDIC character, or a
  !    blank character if no correspondence holds.
  !
  implicit none

  character              ch
  character              ch_to_ebcdic
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ic_to_iebcdic

  i = ic_to_iebcdic ( iachar ( ch ) )

  if ( i /= -1 ) then
     ch_to_ebcdic = achar ( i )
  else
     ch_to_ebcdic = ' '
  end if

  return
end function ch_to_ebcdic
subroutine ch_to_military ( ch, military )

  !*****************************************************************************80
  !
  !! CH_TO_MILITARY converts an ASCII character to a Military code word.
  !
  !  Example:
  !
  !    'A'  'Alpha'
  !    'B'  'Bravo'
  !    'Z'  'Zulu'
  !    'a'  'alpha'
  !    '7'  '7'
  !    '%'  '%'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the ASCII character.
  !
  !    Output, character ( len = 8 ) MILITARY, the military code word.
  !    If CH is not an alphabetic letter, then MILITARY is simply set equal to CH.
  !
  implicit none

  integer ( kind = 4 ) a_to_i4
  character ch
  character ( len = 8 ), parameter, dimension ( 26 ) :: code = (/ &
       'alpha   ', 'bravo   ', 'charlie ', 'delta   ', 'echo    ', &
       'foxtrot ', 'golf    ', 'hotel   ', 'india   ', 'juliet  ', &
       'kilo    ', 'lima    ', 'mike    ', 'november', 'oscar   ', &
       'papa    ', 'quebec  ', 'romeo   ', 'sierra  ', 'tango   ', &
       'uniform ', 'victor  ', 'whiskey ', 'x-ray   ', 'yankee  ', &
       'zulu    ' /)
  integer ( kind = 4 ) i
  character ( len = * ) military

  if ( 'A' <= ch .and. ch <= 'Z' ) then
     i = a_to_i4 ( ch )
     military = code(i)
     call ch_cap ( military(1:1) )
  else if ( 'a' <= ch .and. ch <= 'z' ) then
     i = a_to_i4 ( ch ) - 26
     military = code(i)
  else
     military = ch
  end if

  return
end subroutine ch_to_military
subroutine ch_to_morse ( ch, morse )

  !*****************************************************************************80
  !
  !! CH_TO_MORSE converts an ASCII character to a Morse character string.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 September 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the ASCII character.
  !
  !    Output, character ( len = 6 ) MORSE, the Morse character string.
  !
  implicit none

  integer   ( kind = 4 ), parameter :: num_symbol = 45

  character              ch
  integer   ( kind = 4 ) iascii
  integer   ( kind = 4 ) ic_to_imorse
  integer   ( kind = 4 ) imorse
  character ( len = 6 )  morse
  character ( len = 6 ), parameter, dimension ( num_symbol ) :: msymbol = (/ &
       '      ', '.-    ', '-...  ', '-.-.  ', '-..   ', &
       '.     ', '..-.  ', '--.   ', '....  ', '..    ', &
       '.---  ', '-.-   ', '.-..  ', '--    ', '-.    ', &
       '---   ', '.--.  ', '--.-  ', '.-.   ', '...   ', &
       '-     ', '..-   ', '...-  ', '.--   ', '-..-  ', &
       '-.--  ', '--..  ', '.---- ', '..--- ', '...-- ', &
       '....- ', '..... ', '-.... ', '--... ', '---.. ', &
       '----. ', '----- ', '.-.-.-', '--..--', '---...', &
       '..--..', '.----.', '-....-', '-..-. ', '.-..-.' /)

  iascii = iachar ( ch )
  imorse = ic_to_imorse ( iascii )

  if ( imorse == -1 ) then
     morse = ' '
  else
     morse = msymbol ( imorse )
  end if

  return
end subroutine ch_to_morse
function ch_to_rot13 ( ch )

  !*****************************************************************************80
  !
  !! CH_TO_ROT13 converts a character to its ROT13 equivalent.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, which
  !    guarantees the ASCII collating sequence.
  !
  !    Two applications of CH_TO_ROT13 to a character will return the original.
  !
  !    As a further scrambling, digits are similarly rotated using
  !    a "ROT5" scheme.
  !
  !  Example:
  !
  !    Input:  Output:
  !
  !    a       n
  !    C       P
  !    J       W
  !    1       6
  !    5       0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be converted.
  !
  !    Output, character CH_TO_ROT13, the ROT13 equivalent of the character.
  !
  implicit none

  character              ch
  character              ch_to_rot13
  integer   ( kind = 4 ) itemp

  itemp = iachar ( ch )
  !
  !  [0:4] -> [5:9]
  !
  if ( 48 <= itemp .and. itemp <= 52 ) then
     itemp = itemp + 5
     !
     !  [5:9] -> [0:4]
     !
  else if ( 53 <= itemp .and. itemp <= 57 ) then
     itemp = itemp - 5
     !
     !  [A:M] -> [N:Z]
     !
  else if ( 65 <= itemp .and. itemp <= 77 ) then
     itemp = itemp + 13
     !
     !  [N:Z] -> [A:M]
     !
  else if ( 78 <= itemp .and. itemp <= 90 ) then
     itemp = itemp - 13
     !
     !  [a:m] -> [n:z]
     !
  else if ( 97 <= itemp .and. itemp <= 109 ) then
     itemp = itemp + 13
     !
     !  [n:z] -> [a:m]
     !
  else if ( 110 <= itemp .and. itemp <= 122 ) then
     itemp = itemp - 13
  end if

  ch_to_rot13 = achar ( itemp )

  return
end function ch_to_rot13
function ch_to_scrabble ( ch )

  !*****************************************************************************80
  !
  !! CH_TO_SCRABBLE returns the Scrabble index of a character.
  !
  !  Discussion:
  !
  !    'A' through 'Z' have indices 1 through 26, and blank is index 27.
  !    Case is ignored.  All other characters return index -1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character.
  !
  !    Output, integer ( kind = 4 ) CH_TO_SCRABBLE, the Scrabble index of
  !    the character.
  !
  implicit none

  integer   ( kind = 4 ) a_to_i4
  character              ch
  character              ch_copy
  integer   ( kind = 4 ) ch_to_scrabble
  integer   ( kind = 4 ) ic

  if ( ch == ' ' ) then
     ch_to_scrabble = 27
     return
  end if

  ch_copy = ch
  call ch_cap ( ch_copy )
  ic = a_to_i4 ( ch_copy )

  if ( 1 <= ic .and. ic <= 26 ) then
     ch_to_scrabble = ic
  else
     ch_to_scrabble = -1
  end if

  return
end function ch_to_scrabble
subroutine ch_to_soundex ( ch, soundex )

  !*****************************************************************************80
  !
  !! CH_TO_SOUNDEX converts an ASCII character to a Soundex character.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, which
  !    guarantees the ASCII collating sequence.
  !
  !    The soundex code is used to replace words by a code of up to four
  !    digits.  Similar sounding words will often have identical soundex
  !    codes.
  !
  !    Soundex  Letters
  !    -------  ---------------
  !       0     A E I O U Y H W
  !       1     B B P V
  !       2     C G J K Q S X Z
  !       3     D T
  !       4     L
  !       5     M N
  !       6     R
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 January 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the ASCII character.
  !
  !    Output, character SOUNDEX, the Soundex character, which is
  !    '0', '1', '2', '3', '4', '5', '6', or ' '.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) iascii
  integer   ( kind = 4 ) ic_to_isoundex
  integer   ( kind = 4 ) isoundex
  character              soundex

  iascii = iachar ( ch )
  isoundex = ic_to_isoundex ( iascii )

  if ( isoundex == -1 ) then
     soundex = ' '
  else
     soundex = achar ( isoundex )
  end if

  return
end subroutine ch_to_soundex
subroutine ch_to_sym ( ch, sym )

  !*****************************************************************************80
  !
  !! CH_TO_SYM returns a printable symbol for any ASCII character.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the character to be represented.
  !
  !    Output, character ( len = 4 ) SYM, is the printable symbol for CHR.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ), parameter :: i4_128 = 128
  integer   ( kind = 4 ) put
  character ( len = 4 )  sym

  i = iachar ( ch )

  sym = ' '

  put = 0
  !
  !  Characters 128-255 are symbolized with a ! prefix.
  !  Then shift them down by 128.
  !  Now all values of I are between 0 and 127.
  !
  if ( 128 <= i ) then
     i = mod ( i, i4_128 )
     put = put + 1
     sym(put:put) = '!'
  end if
  !
  !  Characters 0-31 are symbolized with a ^ prefix.
  !  Shift them up by 64.  Now all values of I are between 32 and 127.
  !
  if ( i <= 31 ) then
     i = i + 64
     put = put + 1
     sym(put:put) = '^'
  end if
  !
  !  Character 32 becomes SP.
  !  Characters 32 through 126 are themselves.
  !  Character 127 is DEL.
  !
  if ( i == 32 ) then
     put = put + 1
     sym(put:put+1) = 'SP'
  else if ( i <= 126 ) then
     put = put + 1
     sym(put:put) = achar ( i )
  else if ( i == 127 ) then
     put = put + 1
     sym(put:put+2) = 'DEL'
  end if

  return
end subroutine ch_to_sym
function ch_uniform ( clo, chi, seed )

  !*****************************************************************************80
  !
  !! CH_UNIFORM returns a random character in a given range.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
  !    which guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CLO, CHI, the minimum and maximum acceptable characters.
  !
  !    Input/output, integer ( kind = 4 ) SEED, a seed for the random
  !    number generator.
  !
  !    Output, character CH_UNIFORM, the randomly chosen character.
  !
  implicit none

  character              ch_uniform
  character              chi
  character              clo
  real      ( kind = 4 ) r4_uniform_01
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) seed

  ilo = iachar ( clo )
  ihi = iachar ( chi )

  i = ilo + int ( r4_uniform_01 ( seed ) * real ( ihi + 1 - ilo, kind = 4 ) )

  i = max ( i, ilo )
  i = min ( i, ihi )

  ch_uniform = achar ( i )

  return
end function ch_uniform
subroutine ch3_to_ch_amino ( ch3, ch )

  !*****************************************************************************80
  !
  !! CH3_TO_CH_AMINO converts a 3 character to a 1 character code for amino acids.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carl Branden, John Tooze,
  !    Introduction to Protein Structure,
  !    Garland Publishing, 1991.
  !
  !  Parameters:
  !
  !    Input, character ( len = 3 ) CH3, presumably the 3 letter code for an
  !    amino acid.  Lower and upper case letters are treated the same.
  !
  !    Output, character CH, the one letter code for the amino acid.
  !    If the input code is not recognized, then CH will be '?'.
  !
  implicit none

  integer ( kind = 4 ), parameter :: n = 23

  character ch
  character, parameter, dimension ( n ) :: ch_table = (/ &
       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
       'X', 'Y', 'Z' /)
  character ( len = 3 ) ch3
  character ( len = 3 ), parameter, dimension ( n ) :: ch3_table = (/ &
       'Ala', 'Asx', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ise', 'Lys', &
       'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Val', 'Trp', &
       'X  ', 'Tyr', 'Glx' /)
  integer ( kind = 4 ) i
  logical s_eqi

  do i = 1, n
     if ( s_eqi ( ch3, ch3_table(i) ) ) then
        ch = ch_table(i)
        return
     end if
  end do

  ch = '?'

  return
end subroutine ch3_to_ch_amino
subroutine ch4_to_i4 ( ch4, i4 )

  !*****************************************************************************80
  !
  !! CH4_TO_I4 converts a four character string to an integer.
  !
  !  Example:
  !
  !    Adam    1097097581
  !    Bill    1114205292
  !    Crow    1131573111
  !    Dave    1147237989
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = 4 ) CH4, the character value.
  !
  !    Output, integer ( kind = 4 ) I4, a corresponding integer value.
  !
  implicit none

  character              c1
  character              c2
  character              c3
  character              c4
  character ( len = 4 )  ch4
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) j1
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j3
  integer   ( kind = 4 ) j4

  read ( ch4, '(4a1)' ) c1, c2, c3, c4

  j1 = iachar ( c1 )
  j2 = iachar ( c2 )
  j3 = iachar ( c3 )
  j4 = iachar ( c4 )

  call mvbits ( j1, 0, 8, i4,  0 )
  call mvbits ( j2, 0, 8, i4,  8 )
  call mvbits ( j3, 0, 8, i4, 16 )
  call mvbits ( j4, 0, 8, i4, 24 )

  return
end subroutine ch4_to_i4
subroutine ch4_to_r4 ( ch4, r4 )

  !*****************************************************************************80
  !
  !! CH4_TO_R4 converts a 4 character string to an R4.
  !
  !  Discussion:
  !
  !    The MVBITS routine requires the two word arguments to be of the
  !    same arithmetic type, so we first need to use the TRANSFER
  !    function so that the data inside an integer word can be copied
  !    verbatin into a real.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = 4 ) CH4, the character value.
  !
  !    Output, real ( kind = 4 ) R4, a corresponding real value.
  !
  implicit none

  character              c1
  character              c2
  character              c3
  character              c4
  character ( len = 4 )  ch4
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) j
  real      ( kind = 4 ) r4

  read ( ch4, '(4a1)' ) c1, c2, c3, c4

  j = iachar ( c1 )
  call mvbits ( j, 0, 8, i4,  0 )

  j = iachar ( c2 )
  call mvbits ( j, 0, 8, i4,  8 )

  j = iachar ( c3 )
  call mvbits ( j, 0, 8, i4, 16 )

  j = iachar ( c4 )
  call mvbits ( j, 0, 8, i4, 24 )

  r4 = transfer ( i4, r4 )

  return
end subroutine ch4_to_r4
subroutine ch4vec_to_i4vec ( n, s, i4vec )

  !*****************************************************************************80
  !
  !! CH4VEC_TO_I4VEC converts an string of characters into an array of integers.
  !
  !  Discussion:
  !
  !    This routine can be useful when trying to write character data to an
  !    unformatted direct access file.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of sets of 4 characters
  !    in the string.
  !
  !    Input, character ( len = 4*N ) S, the string of characters.
  !    Each set of 4 characters is assumed to represent an integer.
  !
  !    Output, integer ( kind = 4 ) I4VEC(N), the integers encoded in the string.
  !
  implicit none

  integer   ( kind = 4 )  n

  integer   ( kind = 4 )  i
  integer   ( kind = 4 )  i4vec(n)
  integer   ( kind = 4 )  j
  character ( len = 4*n ) s

  do i = 1, n
     j = 4 * ( i - 1 ) + 1
     call ch4_to_i4 ( s(j:j+3), i4vec(i) )
  end do

  return
end subroutine ch4vec_to_i4vec
subroutine chr4_to_8 ( s1, s2 )

  !*****************************************************************************80
  !
  !! CHR4_TO_8 replaces pairs of hexadecimal digits by a character.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be decoded.
  !
  !    Output, character ( len = * ) S2, the output string.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) j1
  integer   ( kind = 4 ) k1
  integer   ( kind = 4 ) nchar2
  integer   ( kind = 4 ) nroom
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  !
  !  Set S1_LENGTH to the number of characters to be copied.
  !
  nchar2 = 0
  s1_length = len ( s1 )

  if ( mod ( s1_length, 2 ) == 1 ) then
     s1_length = s1_length - 1
  end if

  if ( s1_length <= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHR4_TO_8 - Serious error!'
     write ( *, '(a)' ) '  The input string has nonpositive length!'
     return
  end if
  !
  !  Make sure we have enough room.
  !
  nroom = len ( s2 )

  if ( 2 * nroom < s1_length ) then

     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHR4_TO_8 - Warning!'
     write ( *, '(a)' ) '  Not enough room in the output string.'
     write ( *, '(a,i8)' ) '  Positions available = ', nroom
     write ( *, '(a,i8)' ) '  Positions needed =    ', s1_length / 2
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) '  The program will drop excess characters.'

     s1_length = 2 * nroom

  end if

  do i = 1, s1_length, 2

     call hex_digit_to_i4 ( s1(i:i), j1 )
     call hex_digit_to_i4 ( s1(i+1:i+1), k1 )
     !
     !  Make sure that the values of J1 and K1 are legal.  If not,
     !  set I1 so that it returns a blank character.
     !
     if ( ( 0 <= j1 .and. j1 <= 15) .and. ( 0 <= k1 .and. k1 <= 15) ) then
        i1 = 16 * j1 + k1
     else
        i1 = 0
     end if

     nchar2 = nchar2 + 1
     s2(nchar2:nchar2) = achar ( i1 )

  end do

  return
end subroutine chr4_to_8
subroutine chr8_to_4 ( s1, s2 )

  !*****************************************************************************80
  !
  !! CHR8_TO_4 replaces characters by a pair of hexadecimal digits.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    Unprintable characters (0 through 31, or 127 through 255)
  !    can be displayed.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be replaced.
  !
  !    Output, character ( len = * ) S2, the output string.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j1
  integer   ( kind = 4 ) k1
  integer   ( kind = 4 ) nroom
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  !
  !  Set S1_LENGTH to the number of characters to be copied.
  !
  s1_length = len ( s1 )

  if ( s1_length <= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHR8_TO_4 - Serious error!'
     write ( *, '(a)' ) '  The input string has nonpositive length!'
     return
  end if
  !
  !  Make sure we have enough room.
  !
  nroom = len ( s2 )

  if ( nroom < 2 * s1_length ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHR8_TO_4 - Warning!'
     write ( *, '(a)' ) '  The output string isn''t long enough to hold'
     write ( *, '(a)' ) '  all the information!'
     write ( *, '(a,i8)' ) '  Positions available: ', nroom
     write ( *, '(a,i8)' ) '  Positions needed:    ', 2 * s1_length
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) '  We will do a partial conversion.'

     s1_length = nroom / 2

  end if

  j = 0

  do i = 1, s1_length

     c = s1(i:i)

     i1 = iachar ( c )
     !
     !  Compute J1 and K1 so that I1 = J1*16+K1.
     !
     j1 = i1 / 16
     k1 = i1 - 16 * j1

     j = j + 1
     call i4_to_hex_digit ( j1, s2(j:j) )

     j = j + 1
     call i4_to_hex_digit ( k1, s2(j:j) )

  end do

  return
end subroutine chr8_to_4
subroutine chra_to_s ( s1, s2 )

  !*****************************************************************************80
  !
  !! CHRA_TO_S replaces control characters by printable symbols.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Table:
  !
  !    IACHAR(c)    Symbol
  !    --------    ------
  !      0          ^@
  !      1          ^A
  !    ...         ...
  !     31          ^_
  !     32         (space)
  !    ...         ...
  !    126         ~
  !    127         DEL
  !    128         !^@
  !    ...         ...
  !    159         !^_
  !    160         !(space)
  !    ...         ...
  !    254         !~
  !    255         !DEL
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be operated on.
  !
  !    Output, character ( len = * ) S2, a copy of S1, except that each
  !    control character has been replaced by a symbol.
  !
  implicit none

  logical                ch_is_control
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  integer   ( kind = 4 ) lsym
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  character ( len = 4 )  sym

  s1_length = len_trim ( s1 )
  s2 = ' '

  put = 1

  do get = 1, s1_length

     if ( ch_is_control ( s1(get:get) ) ) then

        call ch_to_sym ( s1(get:get), sym )
        lsym = len_trim ( sym )

        s2(put:put+lsym-1) = sym(1:lsym)
        put = put + lsym

     else

        s2(put:put) = s1(get:get)
        put = put + 1

     end if

  end do

  return
end subroutine chra_to_s
subroutine chrasc ( iascii, nascii, string )

  !*****************************************************************************80
  !
  !! CHRASC converts a vector of ASCII codes into character strings.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IASCII(NASCII), a vector presumed to
  !    contain entries between 0 and 255, the ASCII codes of
  !    individual characters.
  !
  !    Input, integer ( kind = 4 ) NASCII, the number of ASCII codes input.
  !
  !    Output, character ( len = * ) STRING(*).  STRING is assumed to be
  !    a vector of sufficient size to contain the information
  !    input in IASCII.
  !
  !    The length of the strings is determined via the
  !    LEN function.  The entries in IASCII are converted and
  !    stored into the characters of STRING(1), and when that is
  !    full, into STRING(2) and so on until all the entries have
  !    been converted.
  !
  !    If any entry of IASCII is less than 0, or greater than
  !    255, it is handled as though it were 0.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iascii(*)
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) itemp
  integer   ( kind = 4 ) ix
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) nascii
  integer   ( kind = 4 ) nchar
  character ( len = * )  string(*)

  nchar = len ( string(1) )

  ix = 0
  ihi = ( (nascii-1) / nchar ) + 1

  do i = 1, ihi
     do j = 1, nchar

        ix = ix + 1

        if ( nascii <= ix ) then
           return
        end if

        itemp = iascii ( ix )

        if ( itemp < 0 .or. 255 < itemp ) then
           itemp = 0
        end if

        string(i)(j:j) = achar ( itemp )

     end do
  end do

  return
end subroutine chrasc
subroutine chrass ( s, lhs, rhs )

  !*****************************************************************************80
  !
  !! CHRASS "understands" an assignment statement of the form LHS = RHS.
  !
  !  Discussion:
  !
  !    CHRASS returns a string containing the left hand side, and another
  !    string containing the right hand side.
  !
  !    Leading and trailing spaces are removed from the right hand side
  !    and the left hand side.
  !
  !  Example:
  !
  !    S                            Rhs               Lhs
  !
  !    'a = 1.0'                    'a'               '1.0'
  !    'n = -17'                      'n'               '-17'
  !    'scale = +5.3E-2'            'scale'           '+5.3E-2'
  !    'filename = myprog.f'        'filename'        'myprog.f'
  !    '= A pot of gold'            ' '               'A pot of gold'
  !    'Fred'                       'Fred'            ' '
  !    '= Bob'                      ' '               'Bob'
  !    '1=2, 2=3, 3=4'              '1'               '2, 2=3, 3=4'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the assignment statement to be broken up.
  !
  !    Output, character ( len = * ) LHS.
  !
  !    LHS contains the left hand side of the assignment statement.
  !
  !    Normally, this will be the name of a variable, which is
  !    assumed to be whatever appears before the first equals
  !    sign in the string.
  !
  !    If the input line was blank, then LHS will equal ' '.
  !
  !    If the input line contains an equal sign, but nothing
  !    before the equals sign except blanks, then LHS will be ' '.
  !
  !    If the input line does not contain an "=" sign, then
  !    NAME will contain the text of the whole line.
  !
  !    If an error occurred while trying to process the
  !    input line, NAME will contain the text of the line..
  !
  !    If the line began with "#", then NAME will contain the
  !    text of the line.
  !
  !    If the line equals "end-of-input", then NAME will contain
  !    the text of the line.
  !
  !    Output, character ( len = * ) RHS.
  !
  !    RHS contains the right hand side of the assignment statement.
  !
  !    RHS is whatever appears on the right hand side of the
  !    first equals sign in the string.
  !
  !    If S is blank, then RHS is ' '.
  !
  !    If the string contains no equals sign, then RHS is ' '.
  !
  !    If the string contains nothing to the right of the first equals
  !    sign, but blanks, then RHS is ' '.
  !
  !    The user may read the data in RHS by
  !
  !      calling S_TO_R8 to read real ( kind = 8 ) data,
  !      calling CHRCTR to read real data,
  !      calling CHRCTI to read integer data,
  !      calling CHRCTL to read logical data,
  !      calling CHRCTC to read complex data.
  !
  implicit none

  integer   ( kind = 4 ) first
  integer   ( kind = 4 ) iequal
  character ( len = * )  lhs
  character ( len = * )  rhs
  character ( len = * )  s
  integer   ( kind = 4 ) s_first_nonblank
  integer   ( kind = 4 ) s_length
  !
  !  Set default values
  !
  lhs = ' '
  rhs = ' '
  !
  !  Find the last nonblank.
  !
  s_length = len_trim ( s )

  if ( s_length <= 0 ) then
     return
  end if
  !
  !  Look for the first equals sign.
  !
  iequal = index ( s, '=' )
  !
  !  If no equals sign, then LHS = S and return.
  !
  if ( iequal == 0 ) then
     first = s_first_nonblank ( s )
     lhs = s(first:s_length)
     return
  end if
  !
  !  Otherwise, copy LHS = S(1:IEQUAL-1), RHS = S(IEQUAL+1:).
  !
  lhs = s(1:iequal-1)

  if ( iequal + 1 <= s_length ) then
     rhs = s(iequal+1:)
  end if
  !
  !  Now shift the strings to the left.
  !
  lhs = adjustl ( lhs )
  rhs = adjustl ( rhs )

  return
end subroutine chrass
subroutine chrctf ( s, itop, ibot, ierror, length )

  !*****************************************************************************80
  !
  !! CHRCTF reads an integer or rational fraction from a string.
  !
  !  Discussion:
  !
  !    The integer may be in real format, for example '2.25'.  The routine
  !    returns ITOP and IBOT.  If the input number is an integer, ITOP
  !    equals that integer, and IBOT is 1.  But in the case of 2.25,
  !    the program would return ITOP = 225, IBOT = 100.
  !
  !    Legal input is:
  !
  !      blanks,
  !      initial sign,
  !      blanks,
  !      integer ( kind = 4 ) part,
  !      decimal point,
  !      fraction part,
  !      'E' or 'e' or 'D' or 'd', exponent marker,
  !      exponent sign,
  !      exponent integer part,
  !      blanks,
  !      final comma or semicolon.
  !
  !    with most quantities optional.
  !
  !  Example:
  !
  !    S               ITOP      IBOT
  !
  !    '1'               1         1
  !    '     1   '       1         1
  !    '1A'              1         1
  !    '12,34,56'        12        1
  !    '  34 7'          34        1
  !    '-1E2ABCD'        -100      1
  !    '-1X2ABCD'        -1        1
  !    ' 2E-1'           2         10
  !    '23.45'           2345      100
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate when no more characters
  !    can be read to form a legal integer.  Blanks, commas,
  !    or other nonnumeric data will, in particular, cause
  !    the conversion to halt.
  !
  !    Output, integer ( kind = 4 ) ITOP, the integer read from the string,
  !    assuming that no negative exponents or fractional parts
  !    were used.  Otherwise, the 'integer' is ITOP/IBOT.
  !
  !    Output, integer ( kind = 4 ) IBOT, the integer divisor required to
  !    represent numbers which are in real format or have a
  !    negative exponent.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0 if no errors,
  !    Value of IHAVE when error occurred otherwise.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
  !    the string to form the number.
  !
  implicit none

  character              c
  logical                ch_eqi
  integer   ( kind = 4 ) ibot
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) itop
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) ndig
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  ierror = 0
  length = -1
  isgn = 1
  itop = 0
  ibot = 1
  jsgn = 1
  jtop = 0
  ihave = 1
  iterm = 0

  do while ( length < s_length )

     length = length + 1
     c = s(length+1:length+1)
     !
     !  Blank.
     !
     if ( c == ' ' ) then

        if ( ihave == 2 ) then

        else if ( ihave == 6 .or. ihave == 7 ) then
           iterm = 1
        else if ( 1 < ihave ) then
           ihave = 11
        end if
        !
        !  Comma.
        !
     else if ( c == ',' .or. c == ';' ) then

        if ( ihave /= 1 ) then
           iterm = 1
           ihave = 12
           length = length + 1
        end if
        !
        !  Minus sign.
        !
     else if ( c == '-' ) then

        if ( ihave == 1 ) then
           ihave = 2
           isgn = -1
        else if ( ihave == 6 ) then
           ihave = 7
           jsgn = -1
        else
           iterm = 1
        end if
        !
        !  Plus sign.
        !
     else if ( c == '+' ) then

        if ( ihave == 1 ) then
           ihave = 2
        else if ( ihave == 6 ) then
           ihave = 7
        else
           iterm = 1
        end if
        !
        !  Decimal point.
        !
     else if ( c == '.' ) then

        if ( ihave < 4 ) then
           ihave = 4
        else
           iterm = 1
        end if
        !
        !  Exponent marker.
        !
     else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

        if ( ihave < 6 ) then
           ihave = 6
        else
           iterm = 1
        end if
        !
        !  Digit.
        !
     else if ( lle ( '0', c ) .and. lle ( c, '9' ) .and. ihave < 11 ) then

        if ( ihave <= 2 ) then
           ihave = 3
        else if ( ihave == 4 ) then
           ihave = 5
        else if ( ihave == 6 .or. ihave == 7 ) then
           ihave = 8
        end if

        call ch_to_digit ( c, ndig )

        if ( ihave == 3 ) then
           itop = 10 * itop + ndig
        else if ( ihave == 5 ) then
           itop = 10 * itop + ndig
           ibot = 10 * ibot
        else if ( ihave == 8 ) then
           jtop = 10 * jtop + ndig
        end if
        !
        !  Anything else is regarded as a terminator.
        !
     else
        iterm = 1
     end if

     if ( iterm == 1 ) then
        exit
     end if

  end do

  if ( iterm /= 1 .and. length + 1 == s_length ) then
     length = s_length
  end if
  !
  !  Number seems to have terminated.  Have we got a legal number?
  !
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
     ierror = ihave
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHRCTF - Serious error!'
     write ( *, '(a)' ) '  Illegal input:' // trim ( s )
     return
  end if
  !
  !  Number seems OK.  Form it.
  !
  if ( jsgn == 1 ) then
     itop = itop * 10**jtop
  else
     ibot = ibot * 10**jtop
  end if

  itop = isgn * itop

  return
end subroutine chrctf
subroutine chrctg ( s, itop, ibot, ierror, length )

  !*****************************************************************************80
  !
  !! CHRCTG reads an integer, decimal fraction or a ratio from a string.
  !
  !  Discussion:
  !
  !    CHRCTG returns an equivalent ratio (ITOP/IBOT).
  !
  !    If the input number is an integer, ITOP equals that integer, and
  !    IBOT is 1.   But in the case of 2.25, the program would return
  !    ITOP = 225, IBOT = 100.
  !
  !    A ratio is either
  !      a number
  !    or
  !      a number, "/", a number.
  !
  !    A "number" is defined as:
  !
  !      blanks,
  !      initial sign,
  !      integer ( kind = 4 ) part,
  !      decimal point,
  !      fraction part,
  !      E,
  !      exponent sign,
  !      exponent integer part,
  !      blanks,
  !      final comma or semicolon,
  !
  !    Examples of a number:
  !
  !      15, 15.0, -14E-7, E2, -12.73E-98, etc.
  !
  !    Examples of a ratio:
  !
  !      15, 1/7, -3/4.9, E2/-12.73
  !
  !  Example:
  !
  !    S               ITOP      IBOT
  !
  !    '1'               1         1
  !    '     1   '       1         1
  !    '1A'              1         1
  !    '12,34,56'        12        1
  !    '  34 7'          34        1
  !    '-1E2ABCD'        -100      1
  !    '-1X2ABCD'        -1        1
  !    ' 2E-1'           2         10
  !    '23.45'           2345      100
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate when no more characters
  !    can be read to form a legal integer.  Blanks, commas,
  !    or other nonnumeric data will, in particular, cause
  !    the conversion to halt.
  !
  !    Output, integer ( kind = 4 ) ITOP, the integer read from the string,
  !    assuming that no negative exponents or fractional parts
  !    were used.  Otherwise, the 'integer' is ITOP/IBOT.
  !
  !    Output, integer ( kind = 4 ) IBOT, the integer divisor required to
  !    represent numbers which are in decimal format or have a
  !    negative exponent.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0 if no errors,
  !    Value of IHAVE in CHRCTF when error occurred otherwise.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4_gcd
  integer   ( kind = 4 ) ibot
  integer   ( kind = 4 ) ibotb
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) itemp
  integer   ( kind = 4 ) itop
  integer   ( kind = 4 ) itopb
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) length2
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  itop = 0
  ibot = 1
  length = 0

  call chrctf ( s, itop, ibot, ierror, length )

  if ( ierror /= 0) then
     return
  end if
  !
  !  The number is represented as a fraction.
  !  If the next nonblank character is "/", then read another number.
  !
  s_length = len_trim ( s )

  do i = length + 1, s_length - 1

     if ( s(i:i) == '/' ) then

        call chrctf ( s(i+1:), itopb, ibotb, ierror, length2 )

        if ( ierror /= 0 ) then
           return
        end if

        itop = itop * ibotb
        ibot = ibot * itopb

        itemp = i4_gcd ( itop, ibot )

        itop = itop / itemp
        ibot = ibot / itemp

        length = i + length2

        return

     else if ( s(i:i) /= ' ' ) then

        return

     end if

  end do

  return
end subroutine chrctg
subroutine chrcti2 ( s, intval, ierror, length )

  !*****************************************************************************80
  !
  !! CHRCTI2 finds and reads an integer from a string.
  !
  !  Discussion:
  !
  !    The routine is given a string which may contain one or more integers.
  !    Starting at the first character position, it looks for the first
  !    substring that could represent an integer.  If it finds such a string,
  !    it returns the integer's value, and the position of the last character
  !    read.
  !
  !  Example:
  !
  !    S               INTVAL      LENGTH
  !
  !    'Apollo 13'       13          9
  !    '     1   '       1           6
  !    '1A'              1           1
  !    '12,34,56'        12          2
  !    'A1A2A3'          1           2
  !    '-1E2ABCD'        -1          2
  !    '-X20ABCD'        20          4
  !    '23.45'           23          2
  !    ' N = 34, $'      34          7
  !    'Oops!'           0           0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 September 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be read.
  !    Reading will begin at position 1 and terminate at the end of the
  !    string, or when no more characters can be read to form a legal integer.
  !    Blanks, commas, or other nonnumeric data will, in particular,
  !    cause the conversion to halt.
  !
  !    Output, integer ( kind = 4 ) INTVAL, the integer read from the string,
  !    or 0 if there was an error.
  !
  !    Output, integer ( kind = 4 ) IERROR, 0 an integer was found,
  !    1 if no integer found.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )
  ierror = 0
  i = 0
  isgn = 1
  intval = 0
  ihave = 0
  iterm = 0
  !
  !  Examine the next character.
  !
  do while ( iterm /= 1 )

     i = i + 1

     if ( s_length < i ) then

        iterm = 1

     else

        c = s(i:i)
        !
        !  Minus sign.
        !
        if ( c == '-' ) then

           if ( ihave == 0 ) then
              ihave = 1
              isgn = -1
           else
              iterm = 1
           end if
           !
           !  Plus sign.
           !
        else if ( c == '+' ) then

           if ( ihave == 0 ) then
              ihave = 1
           else
              iterm = 1
           end if
           !
           !  Digit.
           !
        else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

           ihave = 2

           call ch_to_digit ( c, idig )

           intval = 10 * intval + idig
           !
           !  Blank or TAB.
           !
        else

           if ( ihave == 2 ) then
              iterm = 1
           else
              ihave = 0
           end if

        end if

     end if

  end do

  if ( ihave == 2 ) then
     length = i - 1
     intval = isgn * intval
  else
     ierror = 0
     length = 0
     intval = 0
  end if

  return
end subroutine chrcti2
subroutine chrctp ( s, cval, ierror, length )

  !*****************************************************************************80
  !
  !! CHRCTP reads a parenthesized complex number from a string.
  !
  !  Discussion:
  !
  !    The routine will read as many characters as possible until it reaches
  !    the end of the string, or encounters a character which cannot be
  !    part of the number.
  !
  !    Legal input is:
  !
  !       1 blanks,
  !
  !       2 left parenthesis, REQUIRED
  !
  !       3 blanks
  !       4 '+' or '-' sign,
  !       5 blanks
  !       6 integer part,
  !       7 decimal point,
  !       8 fraction part,
  !       9 'E' or 'e' or 'D' or 'd', exponent marker,
  !      10 exponent sign,
  !      11 exponent integer part,
  !      12 exponent decimal point,
  !      13 exponent fraction part,
  !      14 blanks,
  !
  !      15 comma, REQUIRED
  !
  !      16 blanks
  !      17 '+' or '-' sign,
  !      18 blanks
  !      19 integer part,
  !      20 decimal point,
  !      21 fraction part,
  !      22 'E' or 'e' or 'D' or 'd', exponent marker,
  !      23 exponent sign,
  !      24 exponent integer part,
  !      25 exponent decimal point,
  !      26 exponent fraction part,
  !      27 blanks,
  !
  !      28 right parenthesis, REQUIRED
  !
  !  Example:
  !
  !    S                   CVAL      IERROR     LENGTH
  !
  !    '(1, 1)'              1 + 1 i   0           5
  !    '( 20 , 99 )'        20+99i     0          11
  !    '(-1.2E+2, +30E-2)'  -120+0.3i  0          17
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate at the end of the string, or when no more
  !    characters can be read to form a legal real.  Blanks,
  !    commas, or other nonnumeric data will, in particular,
  !    cause the conversion to halt.
  !
  !    Output, complex ( kind = 4 ) CVAL, the value read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    1, the string was empty.
  !    2, Did not find left parenthesis.
  !    3, Could not read A correctly.
  !    4, Did not find the comma.
  !    5, Could not read B correctly.
  !    6, Did not find right parenthesis.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read.
  !
  implicit none

  real      ( kind = 4 ) aval
  real      ( kind = 4 ) bval
  character              c
  complex   ( kind = 4 ) cval
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) length
  character ( len = * )  s
  !
  !  Initialize the return arguments.
  !
  ierror = 0
  aval = 0
  bval = 0
  cval = cmplx ( aval, bval, kind = 4 )
  length = 0
  !
  !  Get the length of the line, and if it's zero, return.
  !
  if ( len_trim ( s ) <= 0 ) then
     ierror = 1
     return
  end if
  !
  !  Is the next character a left parenthesis, like it must be?
  !
  call nexchr ( s, ichr, c )

  if ( c /= '(' ) then
     ierror = 2
     return
  end if

  length = ichr
  !
  !  Is the next character a comma?  Then a = 0.
  !
  call nexchr ( s(length+1:), ichr, c )

  if ( c == ',' ) then
     aval = 0
     length = length + ichr
     !
     !  Read the A value.
     !
  else

     call s_to_r4 ( s(length+1:), aval, ierror, ichr )

     if ( ierror /= 0 ) then
        ierror = 3
        length = 0
        return
     end if

     length = length + ichr
     !
     !  Expect to read the comma
     !
     if ( s(length:length) /= ',' ) then
        ierror = 4
        length = 0
        return
     end if

  end if
  !
  !  Is the next character a left parenthesis?  Then b = 0.
  !
  call nexchr ( s(length+1:), ichr, c )

  if ( c == ')' ) then
     bval = 0
     length = length + ichr
     !
     !  Read the B value.
     !
  else

     call s_to_r4 ( s(length+1:), bval, ierror, ichr )

     if ( ierror /= 0 ) then
        ierror = 5
        length = 0
        return
     end if

     length = length + ichr
     !
     !  Expect to read the right parenthesis.
     !
     call nexchr ( s(length+1:), ichr, c )

     if ( c /= ')' ) then
        ierror = 6
        length = 0
        return
     end if

  end if

  length = length + ichr

  cval = cmplx ( aval, bval, kind = 4 )

  return
end subroutine chrctp
subroutine chrs_to_a ( s1, s2 )

  !*****************************************************************************80
  !
  !! CHRS_TO_A replaces all control symbols by control characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1 is the string to be operated on.
  !
  !    Output, character ( len = * ) S2 is a copy of S1, except that each
  !    control symbol has been replaced by a control character.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) put
  integer   ( kind = 4 ) nchar2
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2

  s1_length = len_trim ( s1 )
  nchar2 = len ( s2 )

  ihi = 0
  put = 0

  do

     if ( s1_length <= ihi ) then
        return
     end if

     ilo = ihi + 1

     call sym_to_ch ( s1(ilo:), c, ihi )

     put = put + 1

     if ( nchar2 < put ) then
        exit
     end if

     s2(put:put) = c

  end do

  return
end subroutine chrs_to_a
subroutine chvec_permute ( n, a, p )

  !*****************************************************************************80
  !
  !! CHVEC_PERMUTE permutes a character vector in place.
  !
  !  Discussion:
  !
  !    This routine permutes an array of character "objects", but the same
  !    logic can be used to permute an array of objects of any arithmetic
  !    type, or an array of objects of any complexity.  The only temporary
  !    storage required is enough to store a single object.  The number
  !    of data movements made is N + the number of cycles of order 2 or more,
  !    which is never more than N + N/2.
  !
  !  Example:
  !
  !    Input:
  !
  !      N = 5
  !      P = (   2,   4,   5,   1,   3 )
  !      A = (  'B', 'D', 'E', 'A', 'C' )
  !
  !    Output:
  !
  !      A    = ( 'A', 'B', 'C', 'D', 'E' ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of objects.
  !
  !    Input/output, character A(N), the array to be permuted.
  !
  !    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
  !    that the I-th element of the output array should be the J-th
  !    element of the input array.  P must be a legal permutation
  !    of the integers from 1 to N, otherwise the algorithm will
  !    fail catastrophically.
  !
  implicit none

  integer   ( kind = 4 ) n

  character              a(n)
  character              a_temp
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  integer   ( kind = 4 ) istart
  integer   ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHVEC_PERMUTE - Fatal error!'
     write ( *, '(a)' ) '  The input array does not represent'
     write ( *, '(a)' ) '  a proper permutation.  In particular, the'
     write ( *, '(a,i8)' ) '  array is missing the value ', ierror
     stop
  end if
  !
  !  Search for the next element of the permutation that has not been used.
  !
  do istart = 1, n

     if ( p(istart) < 0 ) then

        cycle

     else if ( p(istart) == istart ) then

        p(istart) = -p(istart)
        cycle

     else

        a_temp = a(istart)
        get = istart
        !
        !  Copy the new value into the vacated entry.
        !
        do

           put = get
           get = p(get)

           p(put) = -p(put)

           if ( get < 1 .or. n < get ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'CHVEC_PERMUTE - Fatal error!'
              write ( *, '(a)' ) '  "get" character is out of bounds.'
              stop
           end if

           if ( get == istart ) then
              a(put) = a_temp
              exit
           end if

           a(put) = a(get)

        end do

     end if

  end do
  !
  !  Restore the signs of the entries.
  !
  p(1:n) = - p(1:n)

  return
end subroutine chvec_permute
subroutine chvec_print ( n, a, title )

  !*****************************************************************************80
  !
  !! CHVEC_PRINT prints a character vector.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of components of the vector.
  !
  !    Input, character A(N), the vector to be printed.
  !
  !    Input, character ( len = * ) TITLE, a title to be printed first.
  !    TITLE may be blank.
  !
  implicit none

  integer   ( kind = 4 ) n

  character              a(n)
  logical                ch_is_printable
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) j
  character ( len = 80 ) string
  character ( len = * )  title

  if ( title /= ' ' ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do ilo = 1, n, 80
     ihi = min ( ilo + 79, n )
     string = ' '
     do i = ilo, ihi
        j = i + 1 - ilo
        if ( ch_is_printable ( a(i) ) ) then
           string(j:j) = a(i)
        end if
     end do

     write ( *, '(a)' ) trim ( string )

  end do

  return
end subroutine chvec_print
subroutine chvec_reverse ( n, x )

  !*****************************************************************************80
  !
  !! CHVEC_REVERSE reverses the elements of a character vector.
  !
  !  Example:
  !
  !    Input:
  !
  !      N = 4, X = ( 'L', 'I', 'V', 'E' ).
  !
  !    Output:
  !
  !      X = ( 'E', 'V', 'I', 'L' ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in the array.
  !
  !    Input/output, character X(N), the array to be reversed.
  !
  implicit none

  integer   ( kind = 4 ) n

  character              cval
  integer   ( kind = 4 ) i
  character              x(n)

  do i = 1, n/2
     cval = x(i)
     x(i) = x(n+1-i)
     x(n+1-i) = cval
  end do

  return
end subroutine chvec_reverse
subroutine chvec_to_s ( n, chvec, s )

  !*****************************************************************************80
  !
  !! CHVEC_TO_S converts a character vector to a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of characters to convert.
  !
  !    Input, character CHVEC(N), a vector of characters.
  !
  !    Output, character ( len = * ) S, a string of characters.
  !
  implicit none

  integer   ( kind = 4 ) n

  character              chvec(n)
  integer   ( kind = 4 ) i
  character ( len = * )  s

  do i = 1, min ( n, len ( s ) )
     s(i:i) = chvec(i)
  end do

  return
end subroutine chvec_to_s
subroutine chvec2_print ( m, a, n, b, title )

  !*****************************************************************************80
  !
  !! CHVEC2_PRINT prints two vectors of characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the length of the first sequence.
  !
  !    Input, character A(M), the first sequence.
  !
  !    Input, integer ( kind = 4 ) N, the length of the second sequence.
  !
  !    Input, character B(N), the second sequence.
  !
  !    Input, character ( len = * ), a title.
  !
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  character              a(m)
  character              ai
  character              b(n)
  character              bi
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, max ( m, n )

     if ( i <= m ) then
        ai = a(i)
     else
        ai = ' '
     end if

     if ( i <= n ) then
        bi = b(i)
     else
        bi = ' '
     end if

     write ( *, '(i3,2x,a1,2x,a1)' ) i, ai, bi

  end do

  return
end subroutine chvec2_print
subroutine comma ( s )

  !*****************************************************************************80
  !
  !! COMMA moves commas left through blanks in a string.
  !
  !  Example:
  !
  !    Input:                    Output:
  !    -----                     ------
  !    "To Henry , our dog"      "To Henry,  our dog"
  !    " , , ,"                  ",,,  "
  !    "  14.0   ,"              "  14.0,  "
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a string in which the
  !    commas are to be shifted left through blanks.
  !
  implicit none

  integer   ( kind = 4 ) iblank
  integer   ( kind = 4 ) icomma
  character ( len = * )  s

  icomma = len_trim ( s )

  do while ( 1 < icomma )

     if ( s(icomma:icomma) == ',' ) then

        iblank = icomma

        do while ( 1 < iblank )
           if ( s(iblank-1:iblank-1) /= ' ' ) then
              exit
           end if
           iblank = iblank - 1
        end do

        if ( icomma /= iblank ) then
           s(icomma:icomma) = ' '
           s(iblank:iblank) = ','
        end if

     end if

     icomma = icomma - 1

  end do

  return
end subroutine comma
subroutine dec_to_s_left ( ival, jval, s )

  !*****************************************************************************80
  !
  !! DEC_TO_S_LEFT returns a left-justified representation of IVAL * 10**JVAL.
  !
  !  Example:
  !
  !    IVAL     JVAL       S
  !    ----     ----       ------
  !       0        0       0
  !      21        3       21000
  !      -3        0       -3
  !     147       -2       14.7
  !      16       -5       0.00016
  !      34       30       Inf
  !     123      -21       0.0000000000000000012
  !      34      -30       0.0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    13 September 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IVAL, JVAL, integers which represent
  !    the decimal.
  !
  !    Output, character ( len = * ) S, the representation of the value.
  !    The string is 'Inf' or '0.0' if the value was too large
  !    or small to represent with a fixed point format.
  !
  implicit none

  character ( len = 22 ) chrrep
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) get1
  integer   ( kind = 4 ) get2
  integer   ( kind = 4 ) put1
  integer   ( kind = 4 ) put2
  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) jval
  integer   ( kind = 4 ) ndigit
  integer   ( kind = 4 ) nleft
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s = ' '

  if ( ival == 0 ) then
     s = '0'
     return
  end if

  s_length = len ( s )
  !
  !  Store a representation of IVAL in CHRREP.
  !
  write ( chrrep, '(i22)' ) ival
  call s_blank_delete ( chrrep )
  ndigit = len_trim ( chrrep )
  !
  !  Inf if JVAL is positive, and S_LENGTH < NDIGIT + JVAL.
  !
  if ( 0 < jval ) then
     if ( s_length < ndigit + jval ) then
        s = 'Inf'
        return
     end if
  end if
  !
  !  Underflow if JVAL is negative, and S_LENGTH < 3 + NDIGIT - JVAL.
  !
  if ( jval < 0 ) then
     if ( 0 < ival ) then
        if ( s_length < 3 - ndigit - jval ) then
           s = '0.0'
           return
        end if
     else
        if ( s_length < 5 - ndigit - jval ) then
           s = '0.0'
           return
        end if
     end if
  end if
  !
  !  If JVAL is nonnegative, insert trailing zeros.
  !
  if ( 0 <= jval ) then

     s(1:ndigit) = chrrep(1:ndigit)

     do i = ndigit + 1, ndigit + jval
        s(i:i) = '0'
     end do

  else if ( jval < 0 ) then

     put2 = 0
     get2 = 0
     !
     !  Sign.
     !
     if ( ival < 0 ) then
        put1 = 1
        put2 = 1
        get2 = 1
        s(put1:put2) = '-'
        ndigit = ndigit - 1
     end if
     !
     !  Digits of the integral part.
     !
     if ( 0 < ndigit + jval ) then
        put1 = put2 + 1
        put2 = put1 + ndigit + jval -1
        get1 = get2 + 1
        get2 = get1 + ndigit+jval - 1
        s(put1:put2) = chrrep(get1:get2)
     else
        put1 = put2 + 1
        put2 = put1
        s(put1:put2) = '0'
     end if
     !
     !  Decimal point.
     !
     put1 = put2 + 1
     put2 = put1
     s(put1:put2) = '.'
     !
     !  Leading zeroes.
     !
     do i = 1, - jval - ndigit
        put1 = put2 + 1
        put2 = put1
        s(put1:put2) = '0'
     end do

     nleft = min ( -jval, ndigit )
     nleft = min ( nleft, s_length - put2 )
     put1 = put2 + 1
     put2 = put1 + nleft - 1
     get1 = get2 + 1
     get2 = get1 + nleft - 1
     s(put1:put2) = chrrep(get1:get2)

  end if

  return
end subroutine dec_to_s_left
subroutine dec_to_s_right ( ival, jval, s )

  !*****************************************************************************80
  !
  !! DEC_TO_S_RIGHT returns a right justified representation of IVAL * 10**JVAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IVAL, JVAL, the integers which represent the
  !    decimal fraction.
  !
  !    Output, character ( len = * ) S, a right justified string
  !    containing the representation of the decimal fraction.
  !
  implicit none

  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) jval
  character ( len = * )  s

  call dec_to_s_left ( ival, jval, s )
  call s_adjustr ( s )

  return
end subroutine dec_to_s_right
subroutine digit_bin_to_ch ( i, ch )

  !*****************************************************************************80
  !
  !! DIGIT_BIN_TO_CH returns the character representation of a binary digit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the integer, between 0 and 1.
  !
  !    Output, character CH, the character representation of the integer.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i

  if ( i == 0 ) then
     ch = '0'
  else if ( i == 1 ) then
     ch = '1'
  else
     ch = '*'
  end if

  return
end subroutine digit_bin_to_ch
subroutine digit_inc ( ch )

  !*****************************************************************************80
  !
  !! DIGIT_INC increments a decimal digit.
  !
  !  Example:
  !
  !    Input  Output
  !    -----  ------
  !    '0'    '1'
  !    '1'    '2'
  !    ...
  !    '8'    '9'
  !    '9'    '0'
  !    'A'    'A'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character CH, a digit to be incremented.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  call ch_to_digit ( ch, digit )

  if ( digit == -1 ) then
     return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
     digit = 0
  end if

  call digit_to_ch ( digit, ch )

  return
end subroutine digit_inc
subroutine digit_oct_to_ch ( i, ch )

  !*****************************************************************************80
  !
  !! DIGIT_OCT_TO_CH returns the character representation of an octal digit.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the integer, between 0 and 7.
  !
  !    Output, character CH, the character representation of the integer.
  !
  character              ch
  integer   ( kind = 4 ) i

  if ( 0 <= i .and. i <= 7 ) then
     ch = achar ( i + 48 )
  else
     ch = '*'
  end if

  return
end subroutine digit_oct_to_ch
subroutine digit_to_ch ( digit, ch )

  !*****************************************************************************80
  !
  !! DIGIT_TO_CH returns the character representation of a decimal digit.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Example:
  !
  !    DIGIT   CH
  !    -----  ---
  !      0    '0'
  !      1    '1'
  !    ...    ...
  !      9    '9'
  !     17    '*'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
  !
  !    Output, character CH, the corresponding character.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

     ch = achar ( digit + 48 )

  else

     ch = '*'

  end if

  return
end subroutine digit_to_ch
function ebcdic_to_ch ( e )

  !*****************************************************************************80
  !
  !! EBCDIC_TO_CH converts an EBCDIC character to ASCII.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character E, the EBCDIC character.
  !
  !    Output, character EBCDIC_TO_CH, the corresponding ASCII
  !    character, or a blank character if no correspondence holds.
  !
  implicit none

  character              e
  character              ebcdic_to_ch
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iebcdic_to_ic

  i = iebcdic_to_ic ( iachar ( e ) )

  if ( i /= -1 ) then
     ebcdic_to_ch = achar ( i )
  else
     ebcdic_to_ch = ' '
  end if

  return
end function ebcdic_to_ch
subroutine ebcdic_to_s ( s )

  !*****************************************************************************80
  !
  !! EBCDIC_TO_S converts a string of EBCDIC characters to ASCII.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.
  !    On input, the EBCDIC string.
  !    On output, the ASCII string.
  !
  implicit none

  character              ebcdic_to_ch
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  do i = 1, s_length
     s(i:i) = ebcdic_to_ch ( s(i:i) )
  end do

  return
end subroutine ebcdic_to_s
subroutine fillch ( s1, field, s2 )

  !*****************************************************************************80
  !
  !! FILLCH writes a string into a subfield of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S1, a string which is presumed
  !    to contain, somewhere, a substring that is to be filled in.
  !    The substring might be '?', for instance.
  !
  !    On output, the substring has been overwritten.
  !
  !    Input, character ( len = * ) FIELD, a substring to be searched for in
  !    S, which denotes the spot where the value should be placed.
  !    Trailing blanks are ignored.
  !
  !    Input, character ( len = * ) S2, the character string to be written
  !    into the subfield.  Trailing blanks are ignored.
  !
  implicit none

  character ( len = * )  field
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  integer   ( kind = 4 ) s_indexi
  character ( len = * )  s1
  character ( len = * )  s2

  i = s_indexi ( s1, field )

  if ( i /= 0 ) then

     lenc = len_trim ( field )
     call s_chop ( s1, i, i+lenc-1 )

     lenc = len_trim ( s2 )
     call s_s_insert ( s1, i, s2(1:lenc) )

  end if

  return
end subroutine fillch
subroutine fillin ( s, field, ival )

  !*****************************************************************************80
  !
  !! FILLIN writes an integer into a subfield of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a string which is presumed
  !    to contain, somewhere, a substring that is to be filled in.
  !    The substring might be '?', for instance.
  !
  !    On output, the substring has been overwritten by the value of IVAL.
  !
  !    Input, character ( len = * ) FIELD, a substring to be searched for in
  !    S, which denotes the spot where the value should be placed.
  !    Trailing blanks are ignored.
  !
  !    Input, integer ( kind = 4 ) IVAL, the value to be written
  !    into the subfield.
  !
  implicit none

  character ( len = * )  field
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) lenc
  integer   ( kind = 4 ) s_indexi
  character ( len = * )  s
  character ( len = 14 ) sval

  i = s_indexi ( s, field )

  if ( i /= 0 ) then

     lenc = len_trim ( field )
     call s_chop ( s, i, i+lenc-1 )

     call i4_to_s_left ( ival, sval )

     lenc = len_trim ( sval )
     call s_s_insert ( s, i, sval(1:lenc) )

  end if

  return
end subroutine fillin
subroutine fillrl ( s, field, r )

  !*****************************************************************************80
  !
  !! FILLRL writes a real into a subfield of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a string which is presumed
  !    to contain, somewhere, a substring that is to be filled in.
  !    The substring might be '?', for instance.
  !    On output, the substring has been overwritten by the value.
  !
  !    Input, character ( len = * ) FIELD, a substring to be searched for in
  !    S, which denotes the spot where the value should be placed.
  !    Trailing blanks are ignored.
  !
  !    Input, real  ( kind = 4 ) R, the value to be written into the subfield.
  !
  implicit none

  character ( len = * )  field
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  real      ( kind = 4 ) r
  character ( len = * )  s
  integer   ( kind = 4 ) s_indexi
  character ( len = 10 ) sval

  i = s_indexi ( s, field )

  if ( i /= 0 ) then

     lenc = len_trim ( field )

     call s_chop ( s, i, i+lenc-1 )

     call r4_to_s_right ( r, sval )
     call s_blank_delete ( sval )
     lenc = len_trim ( sval )

     call s_s_insert ( s, i, sval(1:lenc) )

  end if

  return
end subroutine fillrl
subroutine flt_to_s ( mant, iexp, ndig, s )

  !*****************************************************************************80
  !
  !! FLT_TO_S returns a representation of MANT * 10**IEXP.
  !
  !  Example:
  !
  !    MANT   IEXP   S
  !
  !       1      2   100
  !     101     -1   10.1
  !      23     -3   0.023
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) MANT, the mantissa of the representation.
  !    This is an integer whose magnitude is between 0 and
  !    10**NDIG, that is, 0 <= MANT < 10**NDIG.
  !
  !    Input, integer ( kind = 4 ) IEXP, the exponent of 10 that multiplies MULT.
  !
  !    Input, integer ( kind = 4 ) NDIG, the number of digits of accuracy
  !    in the representation.
  !
  !    Output, character ( len = * ) S, the representation of the quantity.
  !
  implicit none

  integer   ( kind = 4 ) iexp
  integer   ( kind = 4 ) jexp
  integer   ( kind = 4 ) mant
  integer   ( kind = 4 ) ndig
  character ( len = * )  s
  !
  !  Get the length of the string, and set it all to blanks.
  !
  s = ' '
  !
  !  If the mantissa is zero, the number is zero, and we have
  !  a special case: S = '0'.
  !
  if ( mant == 0 ) then
     s = '0'
     return
  else if ( 0 < mant ) then
     s(1:2) = '  '
  else if ( mant < 0 ) then
     s(1:2) = '- '
  end if
  !
  !  Now write the mantissa into S in positions 3 to NDIG+2.
  !
  call i4_to_s_left ( abs ( mant ), s(3:ndig+2) )
  !
  !  Insert a decimal place after the first digit.
  !
  s(2:2) = s(3:3)
  s(3:3) = '.'
  !
  !  Place the "e" representing the exponent.
  !
  s(ndig+3:ndig+3) = 'e'
  !
  !  Write the exponent.
  !
  jexp = 0

  do while ( 10**jexp <= abs ( mant ) )
     jexp = jexp + 1
  end do

  jexp = jexp + iexp - 1

  call i4_to_s_zero ( jexp, s(ndig+4:ndig+6) )
  !
  !  Remove all blanks, effectively shifting the string left too.
  !
  call s_blank_delete ( s )

  return
end subroutine flt_to_s
subroutine forcom ( s, fortran, comment )

  !*****************************************************************************80
  !
  !! FORCOM splits a FORTRAN line into "fortran" and "comment".
  !
  !  Discussion:
  !
  !    The "comment" portion is everything following the first occurrence
  !    of an exclamation mark (and includes the exclamation mark).
  !
  !    The "fortran" portion is everything before the first exclamation
  !    mark.
  !
  !    Either or both the data and comment portions may be blank.
  !
  !  Example:
  !
  !    S                             FORTRAN           COMMENT
  !
  !    '      x = 1952   ! Wow'      '      x = 1952'  '! Wow'
  !    '      continue'              '      continue'  ' '
  !    '! Hey, Abbott!'              ' '               '! Hey, Abbott!'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be analyzed.
  !
  !    Output, character ( len = * ) FORTRAN, the initial portion of the string,
  !    containing a FORTRAN statement.
  !
  !    Output, character COMMENT, the final portion of the string,
  !    containing a comment.
  !
  implicit none

  character ( len = * )  comment
  character ( len = * )  fortran
  integer   ( kind = 4 ) i
  character ( len = * )  s

  i = index ( s, '!' )

  if ( i == 0 ) then
     fortran = s
     comment = ' '
  else if ( i == 1 ) then
     fortran = ' '
     comment = s
  else
     fortran = s ( 1:i-1)
     comment = s ( i: )
  end if

  return
end subroutine forcom
subroutine get_unit ( iunit )

  !*****************************************************************************80
  !
  !! GET_UNIT returns a free FORTRAN unit number.
  !
  !  Discussion:
  !
  !    A "free" FORTRAN unit number is an integer between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer ( kind = 4 ) IUNIT.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5 and 6).
  !
  !    Otherwise, IUNIT is an integer between 1 and 99, representing a
  !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  !    are special, and will never return those values.
  !
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical              lopen

  iunit = 0

  do i = 1, 99

     if ( i /= 5 .and. i /= 6 ) then

        inquire ( unit = i, opened = lopen, iostat = ios )

        if ( ios == 0 ) then
           if ( .not. lopen ) then
              iunit = i
              return
           end if
        end if

     end if

  end do

  return
end subroutine get_unit
subroutine hex_digit_to_i4 ( ch, i )

  !*****************************************************************************80
  !
  !! HEX_DIGIT_TO_I4 converts a hexadecimal digit to an I4.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 August 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character CH, the hexadecimal digit, '0'
  !    through '9', or 'A' through 'F', or also 'a' through 'f'
  !    are allowed.
  !
  !    Output, integer ( kind = 4 ) I, the corresponding integer, or -1 if
  !    CH was illegal.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i

  i = iachar ( ch )

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

     i = i - 48

  else if ( 65 <= i .and. i <= 70 ) then

     i = i - 55

  else if ( 97 <= i .and. i <= 102 ) then

     i = i - 87

  else if ( ch == ' ' ) then

     i = 0

  else

     i = -1

  end if

  return
end subroutine hex_digit_to_i4
subroutine hex_to_binary_digits ( hex_digit, binary_digits )

  !*****************************************************************************80
  !
  !! HEX_TO_BINARY_DIGITS converts a hexadecimal digit to 4 binary digits.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 August 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character HEX_DIGIT, the hexadecimal digit.
  !
  !    Output, character ( len = 4 ) BINARY_DIGITS, the binary digits.
  !
  implicit none

  character ( len = 4 ) binary_digits
  character             hex_digit

  if ( hex_digit == '0' ) then
     binary_digits = '0000'
  else if ( hex_digit == '1' ) then
     binary_digits = '0001'
  else if ( hex_digit == '2' ) then
     binary_digits = '0010'
  else if ( hex_digit == '3' ) then
     binary_digits = '0011'
  else if ( hex_digit == '4' ) then
     binary_digits = '0100'
  else if ( hex_digit == '5' ) then
     binary_digits = '0101'
  else if ( hex_digit == '6' ) then
     binary_digits = '0110'
  else if ( hex_digit == '7' ) then
     binary_digits = '0111'
  else if ( hex_digit == '8' ) then
     binary_digits = '1000'
  else if ( hex_digit == '9' ) then
     binary_digits = '1001'
  else if ( hex_digit == 'A' .or. hex_digit == 'a' ) then
     binary_digits = '1010'
  else if ( hex_digit == 'B' .or. hex_digit == 'b' ) then
     binary_digits = '1011'
  else if ( hex_digit == 'C' .or. hex_digit == 'c' ) then
     binary_digits = '1100'
  else if ( hex_digit == 'D' .or. hex_digit == 'd' ) then
     binary_digits = '1101'
  else if ( hex_digit == 'E' .or. hex_digit == 'e' ) then
     binary_digits = '1110'
  else if ( hex_digit == 'F' .or. hex_digit == 'f' ) then
     binary_digits = '1111'
  else
     binary_digits = '    '
  end if

  return
end subroutine hex_to_binary_digits
subroutine hex_to_i4 ( s, i4 )

  !*****************************************************************************80
  !
  !! HEX_TO_I4 converts a hexadecimal string to its integer value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string of hexadecimal digits.
  !
  !    Output, integer ( kind = 4 ) I4, the corresponding integer value.
  !
  implicit none

  integer   ( kind = 4 ) first
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) j
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )
  !
  !  Determine if there is a plus or minus sign.
  !
  isgn = 1

  first = s_length + 1

  do j = 1, s_length

     if ( s(j:j) == '-' ) then
        isgn = -1
     else if ( s(j:j) == '+' ) then
        isgn = + 1
     else if ( s(j:j) /= ' ' ) then
        first = j
        exit
     end if

  end do
  !
  !  Read the numeric portion of the string.
  !
  i4 = 0

  do j = first, s_length
     call hex_digit_to_i4 ( s(j:j), idig )
     i4 = i4 * 16 + idig
  end do

  i4 = isgn * i4

  return
end subroutine hex_to_i4
subroutine hex_to_s ( hex, s )

  !*****************************************************************************80
  !
  !! HEX_TO_S converts a hexadecimal string into characters.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Example:
  !
  !    Input:
  !
  !      '414243'
  !
  !    Output:
  !
  !      'ABC'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) HEX, a string of pairs of hexadecimal values.
  !
  !    Output, character ( len = * ) S, the corresponding character string.
  !
  implicit none

  character ( len = * )  hex
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) ndo
  integer   ( kind = 4 ) nhex
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )
  nhex = len_trim ( hex )
  ndo = min ( nhex / 2, s_length )

  s = ' '

  do i = 1, ndo
     j = 2 * i - 1
     call hex_to_i4 ( hex(j:j+1), intval )
     s(i:i) = achar ( intval )
  end do

  return
end subroutine hex_to_s
subroutine i2_byte_swap ( iword, bytes )

  !*****************************************************************************80
  !
  !! I2_BYTE_SWAP swaps bytes in an 8-byte word.
  !
  !  Discussion:
  !
  !    This routine uses the MVBITS routines to carry out the swaps.  The
  !    relationship between the bits in the word (0 through 63) and the
  !    bytes (1 through 8) is machine dependent.  That is, byte 1 may
  !    comprise bits 0 through 7, or bits 56 through 63.  So some
  !    experimentation may be necessary the first time this routine
  !    is used.
  !
  !    This routine was originally written simply to take the drudgery
  !    out of swapping bytes in a VAX word that was to be read by
  !    another machine.
  !
  !    The statement
  !
  !      call i2_byte_swap ( IWORD, (/ 1, 2, 3, 4, 5, 6, 7, 8 /) )
  !
  !    will do nothing to IWORD, and
  !
  !      call i2_byte_swap ( IWORD, (/ 8, 7, 6, 5, 4, 3, 2, 1 /) )
  !
  !    will reverse the bytes in IWORD, and
  !
  !      call i2_byte_swap ( IWORD, (/ 2, 2, 2, 2, 2, 2, 2, 2 /) )
  !
  !    will replace IWORD with a word containing byte(2) repeated 8 times.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) IWORD, the word whose bits are to
  !    be swapped.
  !
  !    Input, integer ( kind = 4 ) BYTES(8), indicates which byte in the input
  !    word should overwrite each byte of the output word.
  !
  implicit none

  integer ( kind = 4 ), parameter :: bytes_num = 8

  integer ( kind = 4 ) bytes(bytes_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) jword

  jword = iword

  do i = 1, bytes_num

     if ( bytes(i) < 1 .or. bytes_num < bytes(i) ) then
        cycle
     end if

     if ( bytes(i) == i ) then
        cycle
     end if

     call mvbits ( jword, (bytes(i)-1)*8, 8, iword, 0 )

  end do

  return
end subroutine i2_byte_swap
subroutine i4_byte_swap ( iword, bytes )

  !*****************************************************************************80
  !
  !! I4_BYTE_SWAP swaps bytes in a 4-byte word.
  !
  !  Discussion:
  !
  !    This routine uses the MVBITS routines to carry out the swaps.  The
  !    relationship between the bits in the word (0 through 31) and the
  !    bytes (1 through 4) is machine dependent.  That is, byte 1 may
  !    comprise bits 0 through 7, or bits 24 through 31.  So some
  !    experimentation may be necessary the first time this routine
  !    is used.
  !
  !    This routine was originally written simply to take the drudgery
  !    out of swapping bytes in a VAX word that was to be read by
  !    another machine.
  !
  !    The statement
  !
  !      call i4_byte_swap ( IWORD, (/ 1, 2, 3, 4 /) )
  !
  !    will do nothing to IWORD, and
  !
  !      call i4_byte_swap ( IWORD, (/ 4, 3, 2, 1 /) )
  !
  !    will reverse the bytes in IWORD, and
  !
  !      call i4_byte_swap ( IWORD, (/ 2, 2, 2, 2 /) )
  !
  !    will replace IWORD with a word containing byte(2) repeated 4 times.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) IWORD, the word whose bits are to
  !    be swapped.
  !
  !    Input, integer ( kind = 4 ) BYTES(4), indicates which byte in the
  !    input word should overwrite each byte of the output word.
  !
  implicit none

  integer ( kind = 4 ), parameter :: NUM_BYTES = 4

  integer ( kind = 4 ), parameter :: bit_length = 8
  integer ( kind = 4 ) bytes(NUM_BYTES)
  integer ( kind = 4 ) from_pos
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) jword
  integer ( kind = 4 ) to_pos

  jword = iword

  do i = 1, NUM_BYTES

     if ( bytes(i) < 1 .or. NUM_BYTES < bytes(i) ) then
        cycle
     end if

     if ( bytes(i) == i ) then
        cycle
     end if

     from_pos = 8 * ( bytes(i) - 1 )
     to_pos = 8 * ( i - 1 )
     call mvbits ( jword, from_pos, bit_length, iword, to_pos )

  end do

  return
end subroutine i4_byte_swap
subroutine i4_extract ( s, i, ierror )

  !*****************************************************************************80
  !
  !! I4_EXTRACT "extracts" an I4 from the beginning of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S; on input, a string from
  !    whose beginning an integer is to be extracted.  On output,
  !    the integer, if found, has been removed.
  !
  !    Output, integer ( kind = 4 ) I.  If IERROR is 0, then I contains the
  !    next integer read from S; otherwise I is 0.
  !
  !    Output, integer ( kind = 4 ) IERROR.
  !    0, no error.
  !    nonzero, an integer could not be extracted from the beginning of the
  !    string.  I is 0 and S is unchanged.
  !
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0

  call s_to_i4 ( s, i, ierror, length )

  if ( ierror /= 0 .or. length == 0 ) then
     ierror = 1
     i = 0
  else
     call s_shift_left ( s, length )
  end if

  return
end subroutine i4_extract
function i4_gcd ( i, j )

  !*****************************************************************************80
  !
  !! I4_GCD finds the greatest common divisor of I and J.
  !
  !  Discussion:
  !
  !    Note that only the absolute values of I and J are
  !    considered, so that the result is always nonnegative.
  !
  !    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
  !
  !    If I and J have no common factor, I4_GCD is returned as 1.
  !
  !    Otherwise, using the Euclidean algorithm, I_GCD is the
  !    largest common factor of I and J.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, J, two numbers whose greatest common
  !    divisor is desired.
  !
  !    Output, integer ( kind = 4 ) I4_GCD, the greatest common divisor of
  !    I and J.
  !
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j

  i4_gcd = 1
  !
  !  Return immediately if either I or J is zero.
  !
  if ( i == 0 ) then
     i4_gcd = max ( 1, abs ( j ) )
     return
  else if ( j == 0 ) then
     i4_gcd = max ( 1, abs ( i ) )
     return
  end if
  !
  !  Set IP to the larger of I and J, IQ to the smaller.
  !  This way, we can alter IP and IQ as we go.
  !
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
  !
  !  Carry out the Euclidean algorithm.
  !
  do

     ir = mod ( ip, iq )

     if ( ir == 0 ) then
        exit
     end if

     ip = iq
     iq = ir

  end do

  i4_gcd = iq

  return
end function i4_gcd
function i4_huge ( )

  !*****************************************************************************80
  !
  !! I4_HUGE returns a "huge" I4.
  !
  !  Discussion:
  !
  !    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
  !    bit pattern should be
  !
  !     01111111111111111111111111111111
  !
  !    In this case, its numerical value is 2147483647.
  !
  !    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
  !    use I4_HUGE() and HUGE interchangeably.
  !
  !    However, when using the G95, the values returned by HUGE were
  !    not equal to 2147483647, apparently, and were causing severe
  !    and obscure errors in my random number generator, which needs to
  !    add I4_HUGE to the seed whenever the seed is negative.  So I
  !    am backing away from invoking HUGE, whereas I4_HUGE is under
  !    my control.
  !
  !    Explanation: because under G95 the default integer type is 64 bits!
  !    So HUGE ( 1 ) = a very very huge integer indeed, whereas
  !    I4_HUGE ( ) = the same old 32 bit big value.
  !
  !    An I4 is an integer ( kind = 4 ) value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 January 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
  !
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end function i4_huge
subroutine i4_input ( string, value, ierror )

  !*****************************************************************************80
  !
  !! I4_INPUT prints a prompt string and reads an I4 from the user.
  !
  !  Discussion:
  !
  !    If the input line starts with a comment character ('#') or is
  !    blank, the routine ignores that line, and tries to read the next one.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 March 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) STRING, the prompt string.
  !
  !    Output, integer ( kind = 4 ) VALUE, the value input by the user.
  !
  !    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero
  !    if no error occurred.
  !
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  character ( len = 80 ) line
  character ( len = * ) string
  integer ( kind = 4 ) value

  ierror = 0
  value = huge ( value )
  !
  !  Write the prompt.
  !
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

     read ( *, '(a)', iostat = ierror ) line

     if ( ierror /= 0 ) then
        return
     end if
     !
     !  If the line begins with a comment character, go back and read the next line.
     !
     if ( line(1:1) == '#' ) then
        cycle
     end if

     if ( len_trim ( line ) == 0 ) then
        cycle
     end if
     !
     !  Extract integer information from the string.
     !
     call s_to_i4 ( line, value, ierror, last )

     if ( ierror /= 0 ) then
        value = huge ( value )
        return
     end if

     exit

  end do

  return
end subroutine i4_input
function i4_length ( i4 )

  !*****************************************************************************80
  !
  !! I4_LENGTH computes the number of characters needed to print an I4.
  !
  !  Example:
  !
  !        I4    I4_LENGTH
  !
  !         0       1
  !         1       1
  !        -1       2
  !      1952       4
  !    123456       6
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, the integer whose length is desired.
  !
  !    Output, integer ( kind = 4 ) I4_LENGTH, the number of characters required
  !    to print the integer.
  !
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_copy
  integer ( kind = 4 ) i4_length

  if ( i4 < 0 ) then
     i4_length = 1
  else if ( i4 == 0 ) then
     i4_length = 1
     return
  else if ( 0 < i4 ) then
     i4_length = 0
  end if

  i4_copy = abs ( i4 )

  do while ( i4_copy /= 0 )

     i4_length = i4_length + 1

     i4_copy = i4_copy / 10

  end do

  return
end function i4_length
subroutine i4_next ( s, ival, done )

  !*****************************************************************************80
  !
  !! I4_NEXT "reads" I4's from a string, one at a time.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string, presumably containing
  !    integer ( kind = 4 )s.  These may be separated by spaces or commas.
  !
  !    Output, integer ( kind = 4 ) IVAL.  If DONE is FALSE, then IVAL contains
  !    the "next" integer read.  If DONE is TRUE, then IVAL is zero.
  !
  !    Input/output, logical DONE.
  !    On input with a fresh string, the user should set DONE to TRUE.
  !    On output, the routine sets DONE to FALSE if another integer
  !    was read, or TRUE if no more integers could be read.
  !
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s

  ival = 0

  if ( done ) then
     next = 1
     done = .false.
  end if

  if ( len ( s ) < next ) then
     done = .true.
     return
  end if

  call s_to_i4 ( s(next:), ival, ierror, length )

  if ( ierror /= 0 .or. length == 0 ) then
     done = .true.
     next = 1
  else
     done = .false.
     next = next + length
  end if

  return
end subroutine i4_next
subroutine i4_next_read ( s, intval, ierror )

  !*****************************************************************************80
  !
  !! I4_NEXT_READ finds and reads the next I4 in a string.
  !
  !  Discussion:
  !
  !    This routine can be used to extract, one at a time, the integers in
  !    a string.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = 'Data set #12 extends from (5,-43) and is worth $4.56'
  !      IERROR = -1
  !
  !    Output:
  !
  !      (on successive calls)
  !
  !      INTVAL  IERROR
  !      ------  ------
  !           1       0
  !           2       0
  !           5       0
  !         -43       0
  !           4       0
  !          56       0
  !           0       1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 August 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, integer ( kind = 4 ) INTVAL, the next integer in the string, or 0
  !    if no integer could be found.
  !
  !    Input/output, integer ( kind = 4 ) IERROR.
  !    On the first call for a given string, set IERROR = -1.
  !    Thereafter, the routine will return IERROR = 0 if another
  !    integer ( kind = 4 ) was found, or 1 if no more integers were found.
  !
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) intval
  integer ( kind = 4 ), save :: istart = 1
  integer ( kind = 4 ) length
  character ( len = * ) s

  if ( ierror == -1 ) then
     istart = 1
  end if

  ierror = 0
  intval = 0

  if ( len_trim ( s ) < istart ) then
     ierror = 1
     return
  end if

  call chrcti2 ( s(istart:), intval, ierror, length )

  if ( ierror == 0 .and. 0 < length ) then
     istart = istart + length
  else
     ierror = 1
  end if

  return
end subroutine i4_next_read
subroutine i4_range_input ( string, value1, value2, ierror )

  !*****************************************************************************80
  !
  !! I4_RANGE_INPUT reads a pair of I4's from the user, representing a range.
  !
  !  Discussion:
  !
  !    If the input line starts with a comment character ('#') or is blank,
  !    the routine ignores that line, and tries to read the next one.
  !
  !    The pair of integers may be separated by spaces or a comma or both.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 March 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) STRING, the prompt string.
  !
  !    Output, integer ( kind = 4 ) VALUE1, VALUE2, the values entered by
  !    the user.
  !
  !    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero
  !    if no error occurred.
  !
  implicit none

  character, parameter :: comma = ','
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) last2
  character ( len = 80 ) line
  character, parameter :: space = ' '
  character ( len = * ) string
  integer ( kind = 4 ) value1
  integer ( kind = 4 ) value2

  ierror = 0
  value1 = huge ( value1 )
  value2 = huge ( value2 )
  !
  !  Write the prompt.
  !
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

     read ( *, '(a)', iostat = ierror ) line

     if ( ierror /= 0 ) then
        return
     end if
     !
     !  If the line begins with a comment character, go back and read the next line.
     !
     if ( line(1:1) == '#' ) then
        cycle
     end if

     if ( len_trim ( line ) == 0 ) then
        cycle
     end if
     !
     !  Replace commas by spaces.
     !
     call s_replace_ch ( line, comma, space )
     !
     !  Extract integer information from the string.
     !
     call s_to_i4 ( line, value1, ierror, last )

     if ( ierror /= 0 ) then
        value1 = huge ( value1 )
        return
     end if

     call s_to_i4 ( line(last+1:), value2, ierror, last2 )

     if ( ierror /= 0 ) then
        value2 = huge ( value2 )
        return
     end if

     exit

  end do

  return
end subroutine i4_range_input
subroutine i4_swap ( i, j )

  !*****************************************************************************80
  !
  !! I4_SWAP swaps two I4's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 November 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
  !    J have been interchanged.
  !
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end subroutine i4_swap
function i4_to_a ( i )

  !*****************************************************************************80
  !
  !! I4_TO_A returns the I-th alphabetic character.
  !
  !  Discussion:
  !
  !    Instead of CHAR, we now use the ACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Example:
  !
  !    I  I4_TO_A
  !
  !   -8  ' '
  !    0  ' '
  !    1  'A'
  !    2  'B'
  !   ..
  !   26  'Z'
  !   27  'a'
  !   52  'z'
  !   53  ' '
  !   99  ' '
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 February 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the index of the letter to be returned.
  !    0 is a space;
  !    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
  !    27 through 52 requests 'a' through 'z', (ASCII 97:122);
  !
  !    Output, character I4_TO_A, the requested alphabetic letter.
  !
  implicit none

  integer ( kind = 4 ), parameter :: cap_shift = 64
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ), parameter :: low_shift = 96

  if ( i <= 0 ) then
     i4_to_a = ' '
  else if ( 1 <= i .and. i <= 26 ) then
     i4_to_a = achar ( cap_shift + i )
  else if ( 27 <= i .and. i <= 52 ) then
     i4_to_a = achar ( low_shift + i - 26 )
  else if ( 53 <= i ) then
     i4_to_a = ' '
  end if

  return
end function i4_to_a
subroutine i4_to_amino_code ( i, ch )

  !*****************************************************************************80
  !
  !! I4_TO_AMINO_CODE converts an integer to an amino code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carl Branden, John Tooze,
  !    Introduction to Protein Structure,
  !    Garland Publishing, 1991.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the index of an amino acid, between 1
  !    and 23.
  !
  !    Output, character CH, the one letter code for an amino acid.
  !
  implicit none

  integer ( kind = 4 ), parameter :: n = 23

  character ch
  character, dimension ( n ) :: ch_table = (/ &
       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
       'X', 'Y', 'Z' /)
  integer ( kind = 4 ) i

  if ( 1 <= i .and. i <= 23 ) then
     ch = ch_table(i)
  else
     ch = '?'
  end if

  return
end subroutine i4_to_amino_code
subroutine i4_to_base ( i4, base, s )

  !*****************************************************************************80
  !
  !! I4_TO_BASE represents an integer in any base up to 16.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !      Input            Output
  !    -------------    --------
  !    INTVAL   BASE           S
  !
  !         5    -1   '101010101'
  !         5     1       '11111'
  !         5     2         '101'
  !         5     3          '12'
  !         5     4          '11'
  !        -5     5         '-10'
  !         5     6           '5'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, the integer whose representation
  !    is desired.
  !
  !    Input, integer ( kind = 4 ) BASE, the base in which the representation is
  !    given.  BASE must be greater than 0 and no greater than 16.
  !
  !    Output, character ( len = * ) S, the string.
  !
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) icopy
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) jdig
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s = ' '
  s_length = len ( s )
  !
  !  Check the base.
  !
  if ( base < -1 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4_TO_BASE - Serious error!'
     write ( *, '(a)' ) '  The input base is less than -1!'
     return
  end if

  if ( base == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4_TO_BASE - Serious error!'
     write ( *, '(a)' ) '  The input base is zero.'
     return
  end if

  if ( 16 < base ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4_TO_BASE - Serious error!'
     write ( *, '(a)' ) '  The input base is greater than 16!'
     return
  end if
  !
  !  Special treatment for base 1 and -1.
  !
  if ( base == 1 ) then
     call i4_to_unary ( i4, s )
     return
  else if ( base == -1 ) then
     call i4_to_nunary ( i4, s )
     return
  end if
  !
  !  Do repeated mod's
  !
  jdig = 0
  icopy = abs ( i4 )

  do

     if ( ( 0 <= i4 .and. s_length <= jdig ) .or. &
          ( i4 < 0 .and. s_length - 1 <= jdig ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_BASE - Fatal error!'
        do i = 1, s_length
           s(i:i) = '*'
        end do
        return
     end if

     jdig = jdig + 1

     idig = mod ( icopy, base )
     icopy = ( icopy - idig ) / base

     call i4_to_hex_digit ( idig, s(s_length+1-jdig:s_length+1-jdig) )

     if ( icopy == 0 ) then
        exit
     end if

  end do
  !
  !  Take care of the minus sign.
  !
  if ( i4 < 0 ) then
     jdig = jdig + 1
     s(s_length+1-jdig:s_length+1-jdig) = '-'
  end if

  return
end subroutine i4_to_base
subroutine i4_to_binary ( i, s )

  !*****************************************************************************80
  !
  !! I4_TO_BINARY produces the binary representation of an I4.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !     I       S
  !    --  ------
  !     1      '1'
  !     2     '10'
  !     3     '11'
  !     4    '100'
  !    -9  '-1001'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, an integer to be represented.
  !
  !    Output, character ( len = * ) S, the binary representation.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_copy
  integer   ( kind = 4 ) j
  character ( len = * )  s

  i_copy = abs ( i )
  s = ' '
  j = len ( s )

  do while ( 0 < j )

     if ( mod ( i_copy, 2 ) == 1 ) then
        s(j:j) = '1'
     else
        s(j:j) = '0'
     end if

     i_copy = i_copy / 2

     if ( i_copy == 0 ) then
        exit
     end if

     j = j - 1

  end do

  if ( i_copy /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4_TO_BINARY - Serious error!'
     write ( *, '(a)' ) '  Not enough room in the string to represent the value.'
  end if

  if ( i < 0 ) then

     if ( 1 < j ) then
        j = j - 1
        s(j:j) = '-'
     else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_BINARY - Serious error!'
        write ( *, '(a)' ) '  No room to prefix minus sign!'
     end if

  end if

  return
end subroutine i4_to_binary
function i4_to_binhex ( i )

  !*****************************************************************************80
  !
  !! I4_TO_BINHEX returns the I-th character in the BINHEX encoding.
  !
  !  Example:
  !
  !    I  I4_TO_BINHEX
  !
  !    1  '!'
  !    2  '"'
  !    3  '#'
  !   ..
  !   64  'r'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the index of the character.
  !    1 <= I <= 64.
  !
  !    Output, character I4_TO_BINHEX, the requested character.
  !
  implicit none

  integer ( kind = 4 ) i
  character i4_to_binhex
  character ( len = 64 ), parameter :: string = &
       '!"#$%&''()*+,-012345689@ABCDEFGHIJKLMNPQRSTVWXYZ[`abcdefhijklmnpqr'

  if ( 1 <= i .and. i <= 64 ) then
     i4_to_binhex = string(i:i)
  else
     i4_to_binhex = ' '
  end if

  return
end function i4_to_binhex
subroutine i4_to_ch4 ( i4, ch4 )

  !*****************************************************************************80
  !
  !! I4_TO_CH4 converts an I4 to a 4 character string.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !    While most integers will be converted to unprintable strings,
  !    here are a few nice ones:
  !
  !    1097097581  Adam
  !    1114205292  Bill
  !    1131573111  Crow
  !    1147237989  Dave
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 May 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, the integer value.
  !
  !    Output, character ( len = 4 ) CH4, a corresponding character value.
  !
  implicit none

  character ( len = 4 ) ch4
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) j1
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j3
  integer   ( kind = 4 ) j4

  j1 = ibits ( i4,  0, 8 )
  j2 = ibits ( i4,  8, 8 )
  j3 = ibits ( i4, 16, 8 )
  j4 = ibits ( i4, 24, 8 )

  ch4 = achar ( j1 ) // achar ( j2 ) // achar ( j3 ) // achar ( j4 )

  return
end subroutine i4_to_ch4
subroutine i4_to_hex ( i4, s )

  !*****************************************************************************80
  !
  !! I4_TO_HEX produces the hexadecimal representation of an I4.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !        I4      S
  !       ---     ---
  !         0     '0'
  !         9     '9'
  !        10     'A'
  !        15     'F'
  !        16    '10'
  !       100    '64'
  !       -12    '-C'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, the integer to be represented.
  !
  !    Output, character ( len = * ) S, the hexadecimal representation.
  !
  implicit none

  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ), parameter :: i4_16 = 16
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) intcpy
  integer   ( kind = 4 ) isgn
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  intcpy = i4
  isgn = 1

  if ( intcpy < 0 ) then
     isgn = -1
     intcpy = -intcpy
  end if

  s = ' '
  !
  !  Point to the position just after the end of the string.
  !
  ichr = s_length + 1
  !
  !  Moving left, fill in the next digit of the string.
  !
  do

     ichr = ichr - 1

     if ( ichr <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_HEX - Serious error!'
        write ( *, '(a)' ) '  Ran out of room in the string!'
        return
     end if

     i1 = mod ( intcpy, i4_16 )
     intcpy = intcpy / 16

     call i4_to_hex_digit ( i1, s(ichr:ichr) )

     if ( intcpy == 0 ) then

        if ( isgn == -1 ) then

           if ( 1 < ichr ) then
              ichr = ichr - 1
              s(ichr:ichr) = '-'
           else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4_TO_HEX - Serious error!'
              write ( *, '(a)' ) '  No room to prefix minus sign!'
           end if

        end if

        return

     end if

  end do

  return
end subroutine i4_to_hex
subroutine i4_to_hex_digit ( i, ch )

  !*****************************************************************************80
  !
  !! I4_TO_HEX_DIGIT converts a (small) I4 to a hexadecimal digit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 August 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the integer, between 0 and 15.
  !
  !    Output, character CH, the hexadecimal digit corresponding to the integer.
  !
  implicit none

  character ch
  integer ( kind = 4 ) i

  if ( 0 <= i .and. i <= 9 ) then
     ch = achar ( i + 48 )
  else if ( 10 <= i .and. i <= 15 ) then
     ch = achar ( i + 55 )
  else
     ch = '*'
  end if

  return
end subroutine i4_to_hex_digit
function i4_to_isbn ( i )

  !*****************************************************************************80
  !
  !! I4_TO_ISBN converts an I4 to an ISBN digit.
  !
  !  Discussion:
  !
  !    Only the integers 0 through 10 can be input.  The representation
  !    of 10 is 'X'.
  !
  !    An I4 is an integer ( kind = 4 ) value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Book Industry Study Group,
  !    The Evolution in Product Identification:
  !    Sunrise 2005 and the ISBN-13,
  !    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, an integer between 0 and 10.
  !
  !    Output, character I4_TO_ISBN, the ISBN character code of the integer.
  !    If I is illegal, then I4_TO_ISBN is set to '?'.
  !
  implicit none

  integer   ( kind = 4 ) i
  character              i4_to_isbn

  if ( i == 0 ) then
     i4_to_isbn = '0'
  else if ( i == 1 ) then
     i4_to_isbn = '1'
  else if ( i == 2 ) then
     i4_to_isbn = '2'
  else if ( i == 3 ) then
     i4_to_isbn = '3'
  else if ( i == 4 ) then
     i4_to_isbn = '4'
  else if ( i == 5 ) then
     i4_to_isbn = '5'
  else if ( i == 6 ) then
     i4_to_isbn = '6'
  else if ( i == 7 ) then
     i4_to_isbn = '7'
  else if ( i == 8 ) then
     i4_to_isbn = '8'
  else if ( i == 9 ) then
     i4_to_isbn = '9'
  else if ( i == 10 ) then
     i4_to_isbn = 'X'
  else
     i4_to_isbn = '?'
  end if

  return
end function i4_to_isbn
subroutine i4_to_month_abb ( m, month_abb )

  !*****************************************************************************80
  !
  !! I4_TO_MONTH_ABB returns the 3 character abbreviation of a given month.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 June 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the index of the month, which should
  !    be between 1 and 12.
  !
  !    Output, character ( len = 3 ) MONTH_ABB, the month abbreviation
  !
  implicit none

  character ( len = 3 ), parameter, dimension(12) :: abb = (/ &
       'Jan', 'Feb', 'Mar', 'Apr', &
       'May', 'Jun', 'Jul', 'Aug', &
       'Sep', 'Oct', 'Nov', 'Dec' /)
  integer   ( kind = 4 ) m
  character ( len = 3 )  month_abb

  if ( m < 1 .or. 12 < m ) then

     month_abb = '???'

  else

     month_abb = abb(m)

  end if

  return
end subroutine i4_to_month_abb
subroutine i4_to_month_name ( m, month_name )

  !*****************************************************************************80
  !
  !! I4_TO_MONTH_NAME returns the name of a given month.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 April 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the index of the month, which should
  !    be between 1 and 12.
  !
  !    Output, character ( len = * ) MONTH_NAME, a string containing as much
  !    of the month's name as will fit.  To get the typical 3-letter abbreviations
  !    for the months, simply declare
  !      character ( len = 3 ) MONTH_NAME
  !    or pass in MONTH_NAME(1:3).
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) m
  character ( len = * )  month_name
  character ( len = 9 ), parameter, dimension(12) :: name = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)

  if ( m < 1 .or. 12 < m ) then

     do i = 1, len ( month_name )
        month_name(i:i) = '?'
     end do

  else

     month_name = name(m)

  end if

  return
end subroutine i4_to_month_name
subroutine i4_to_nunary ( intval, s )

  !*****************************************************************************80
  !
  !! I4_TO_NUNARY produces the "base -1" representation of an I4.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) INTVAL, an integer to be represented.
  !
  !    Output, character ( len = * ) S, the negative unary representation.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) intval
  character ( len = * )  s

  s = ' '

  if ( intval < 0 ) then

     do i = 1, abs ( intval )
        s(2*i-1:2*i) = '10'
     end do

  else if ( intval == 0 ) then

     s = '0'

  else if ( 0 < intval ) then

     s(1:1) = '1'
     do i = 2, intval
        s(2*i-2:2*i-1) = '01'
     end do

  end if

  s = adjustr ( s )

  return
end subroutine i4_to_nunary
subroutine i4_to_oct ( i4, s )

  !*****************************************************************************80
  !
  !! I4_TO_OCT produces the octal representation of an integer.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !        I4      S
  !
  !         0     '0'
  !         9    '11'
  !        10    '12'
  !        15    '17'
  !        16    '20'
  !       100   '144'
  !       -12   '-14'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, the integer to be represented.
  !
  !    Output, character ( len = * ) S, the octal representation.
  !
  implicit none

  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) ichr
  integer   ( kind = 4 ) intcpy
  integer   ( kind = 4 ) isgn
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  intcpy = i4
  isgn = 1

  if ( intcpy < 0 ) then
     isgn = - 1
     intcpy = -intcpy
  end if

  s = ' '
  !
  !  Point to the position just after the end of the string.
  !
  ichr = s_length + 1
  !
  !  Moving left, fill in the next digit of the string.
  !
  do

     ichr = ichr - 1

     if ( ichr <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_OCT - Fatal error!'
        write ( *, '(a)' ) '  Ran out of room in the string!'
        stop
     end if

     i1 = mod ( intcpy, 8 )
     intcpy = intcpy / 8

     call digit_oct_to_ch ( i1, s(ichr:ichr) )

     if ( intcpy == 0 ) then

        if ( isgn == -1 ) then

           if ( 1 < ichr ) then
              ichr = ichr - 1
              s(ichr:ichr) = '-'
           else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4_TO_OCT - Fatal error!'
              write ( *, '(a)' ) '  No room to prefix minus sign!'
              stop
           end if

        end if

        return

     end if

  end do

  return
end subroutine i4_to_oct
subroutine i4_to_s_left ( i4, s )

  !*****************************************************************************80
  !
  !! I4_TO_S_LEFT converts an I4 to a left-justified string.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !    Assume that S is 6 characters long:
  !
  !        I4  S
  !
  !         1  1
  !        -1  -1
  !         0  0
  !      1952  1952
  !    123456  123456
  !   1234567  ******  <-- Not enough room!
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, an integer to be converted.
  !
  !    Output, character ( len = * ) S, the representation of the integer.
  !    The integer will be left-justified.  If there is not enough space,
  !    the string will be filled with stars.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) ipos
  integer   ( kind = 4 ) ival
  character ( len = * )  s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
     return
  end if
  !
  !  Make a copy of the integer.
  !
  ival = i4
  !
  !  Handle the negative sign.
  !
  if ( ival < 0 ) then

     if ( ihi <= 1 ) then
        s(1:1) = '*'
        return
     end if

     ival = -ival
     s(1:1) = '-'
     ilo = 2

  end if
  !
  !  The absolute value of the integer goes into S(ILO:IHI).
  !
  ipos = ihi
  !
  !  Find the last digit of IVAL, strip it off, and stick it into the string.
  !
  do

     idig = mod ( ival, 10 )
     ival = ival / 10

     if ( ipos < ilo ) then
        do i = 1, ihi
           s(i:i) = '*'
        end do
        return
     end if

     call digit_to_ch ( idig, c )

     s(ipos:ipos) = c
     ipos = ipos - 1

     if ( ival == 0 ) then
        exit
     end if

  end do
  !
  !  Shift the string to the left.
  !
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end subroutine i4_to_s_left
subroutine i4_to_s_right ( intval, s )

  !*****************************************************************************80
  !
  !! I4_TO_S_RIGHT converts an I4 to a right justified string.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !    Assume that S is 6 characters long:
  !
  !    INTVAL       S
  !
  !         1       1
  !        -1      -1
  !         0       0
  !      1952    1952
  !    123456  123456
  !   1234567  ******  <-- Not enough room!
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
  !
  !    Output, character ( len = * ) S, the representation of the integer.
  !    The integer will be right-justified.  If there is not enough space,
  !    the string will be filled with stars.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ), parameter :: i4_ten = 10
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) ipos
  integer   ( kind = 4 ) ival
  character ( len = * )  s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
     return
  end if
  !
  !  Make a copy of the integer.
  !
  ival = intval
  !
  !  Handle the negative sign.
  !
  if ( ival < 0 ) then

     if ( ihi <= 1 ) then
        s(1:1) = '*'
        return
     end if

     ival = -ival
     s(1:1) = '-'
     ilo = 2

  end if
  !
  !  The absolute value of the integer goes into S(ILO:IHI).
  !
  ipos = ihi
  !
  !  Find the last digit of IVAL, strip it off, and stick it into the string.
  !
  do

     idig = mod ( ival, i4_ten )
     ival = ival / 10

     if ( ipos < ilo ) then
        do i = 1, ihi
           s(i:i) = '*'
        end do
        return
     end if

     call digit_to_ch ( idig, c )
     s(ipos:ipos) = c
     ipos = ipos - 1

     if ( ival == 0 ) then
        exit
     end if

  end do
  !
  !  Shift the minus sign, if any.
  !
  if ( s(1:1) == '-' ) then
     if ( ipos /= 1 ) then
        s(1:1) = ' '
        s(ipos:ipos) = '-'
     end if
  end if

  return
end subroutine i4_to_s_right
subroutine i4_to_s_right_comma ( i4, s )

  !*****************************************************************************80
  !
  !! I4_TO_S_RIGHT_COMMA converts an I4 to a right justified string with commas.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !    Assume that S is 10 characters long:
  !
  !          I4           S
  !
  !           1           1
  !          -1          -1
  !           0           0
  !        1952       1,952
  !      123456     123,456
  !     1234567   1,234,567
  !    12345678  12,345,678
  !   123456789  **********  <-- Not enough room!
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 April 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, an integer to be converted.
  !
  !    Output, character ( len = * ) S, the representation of the integer.
  !    The integer will be right-justified.  Commas will be used to separate
  !    sets of three digits.  If there is not enough space, the string will
  !    be filled with stars.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) digit
  integer   ( kind = 4 ) digit_num
  integer   ( kind = 4 ) hi
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ), parameter :: i4_3 = 3
  integer   ( kind = 4 ) lo
  integer   ( kind = 4 ) pos
  character ( len = * )  s
  integer   ( kind = 4 ) value

  s = ' '

  lo = 1
  hi = len ( s )

  if ( hi <= 0 ) then
     return
  end if
  !
  !  Make a copy of the integer.
  !
  value = i4
  !
  !  Handle the negative sign.
  !
  if ( value < 0 ) then

     if ( hi <= 1 ) then
        s(1:1) = '*'
        return
     end if

     value = -value
     s(1:1) = '-'
     lo = 2

  end if
  !
  !  The absolute value of the integer goes into S(LO:HI).
  !
  pos = hi
  !
  !  Find the last digit of VALUE, strip it off, and stick it into the string.
  !
  digit_num = 0

  do

     digit = mod ( value, 10 )
     value = value / 10
     digit_num = digit_num + 1

     if ( pos < lo ) then
        do i = 1, hi
           s(i:i) = '*'
        end do
        return
     end if
     !
     !  Insert a comma?
     !
     if ( 1 < digit_num .and. mod ( digit_num, i4_3 ) == 1 ) then

        if ( pos < lo ) then
           do i = 1, hi
              s(i:i) = '*'
           end do
           return
        end if

        s(pos:pos) = ','
        pos = pos - 1
     end if

     call digit_to_ch ( digit, c )
     s(pos:pos) = c
     pos = pos - 1

     if ( value == 0 ) then
        exit
     end if

  end do
  !
  !  Shift the minus sign, if any.
  !
  if ( s(1:1) == '-' ) then
     if ( pos /= 1 ) then
        s(1:1) = ' '
        s(pos:pos) = '-'
     end if
  end if

  return
end subroutine i4_to_s_right_comma
subroutine i4_to_s_roman ( intval, s )

  !*****************************************************************************80
  !
  !! I4_TO_S_ROMAN converts an I4 to a string of Roman numerals.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !    To generate numbers greater than 4999, the numeral 'V' had a bar
  !    above it, representing a value of 5000, a barred 'X' represented
  !    10,000 and so on.
  !
  !    In the subtractive representation of 4 by 'IV', 9 by 'IX' and so on,
  !    'I' can only subtract from 'V' or 'X',
  !    'X' can only subtract from 'L' or 'C',
  !    'C' can only subtract from 'D' or 'M'.
  !    Under these rules, 1999 cannot be written IMM!
  !
  !  Example:
  !
  !    INTVAL  S
  !
  !        -2  -II <-- Not a Roman numeral
  !        -1  -I  <-- Not a Roman numeral
  !         0   0  <-- Not a Roman numeral
  !         1   I
  !         2   II
  !         3   III
  !         4   IV
  !         5   V
  !        10   X
  !        20   XX
  !        30   XXX
  !        40   XL
  !        50   L
  !        60   LX
  !        70   LXX
  !        80   LXXX
  !        90   XC
  !       100   C
  !       500   D
  !      1000   M
  !      4999   MMMMCMLXLIX
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.  If the
  !    integer has absolute value greater than 4999, the string '?' will be
  !    returned.  If the integer is 0, then the string '0' will be returned.  If
  !    the integer is negative, then a minus sign will precede it, even
  !    though this has nothing to do with Roman numerals.
  !
  !    Output, character ( len = * ) S, the representation of the integer
  !    as a Roman numeral.
  !
  implicit none

  integer   ( kind = 4 ) icopy
  integer   ( kind = 4 ) intval
  character ( len = * )  s

  s = ' '
  icopy = intval

  if ( 4999 < abs ( icopy ) ) then
     s = '?'
     return
  end if

  if ( icopy == 0 ) then
     s = '0'
     return
  end if

  if ( icopy <= 0 ) then
     s = '-'
     icopy = -icopy
  end if

  do while ( 0 < icopy )

     if ( 1000 <= icopy ) then
        call s_cat ( s, 'M', s )
        icopy = icopy - 1000
     else if ( 900 <= icopy ) then
        call s_cat ( s, 'CM', s )
        icopy = icopy - 900
     else if ( 500 <= icopy ) then
        call s_cat ( s, 'D', s )
        icopy = icopy - 500
     else if ( 400 <= icopy ) then
        call s_cat ( s, 'CD', s )
        icopy = icopy - 400
     else if ( 100 <= icopy ) then
        call s_cat ( s, 'C', s )
        icopy = icopy - 100
     else if ( 90 <= icopy ) then
        call s_cat ( s, 'XC', s )
        icopy = icopy - 90
     else if ( 50 <= icopy ) then
        call s_cat ( s, 'L', s )
        icopy = icopy - 50
     else if ( 40 <= icopy ) then
        call s_cat ( s, 'XL', s )
        icopy = icopy - 40
     else if ( 10 <= icopy ) then
        call s_cat ( s, 'X', s )
        icopy = icopy - 10
     else if ( 9 <= icopy ) then
        call s_cat ( s, 'IX', s )
        icopy = icopy - 9
     else if ( 5 <= icopy ) then
        call s_cat ( s, 'V', s )
        icopy = icopy - 5
     else if ( 4 <= icopy ) then
        call s_cat ( s, 'IV', s )
        icopy = icopy - 4
     else
        call s_cat ( s, 'I', s )
        icopy = icopy - 1
     end if

  end do

  return
end subroutine i4_to_s_roman
subroutine i4_to_s_zero ( intval, s )

  !*****************************************************************************80
  !
  !! I4_TO_S_ZERO converts an I4 to a string, with zero padding.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !  Example:
  !
  !    Assume that S is 6 characters long:
  !
  !    INTVAL  S
  !
  !         1  000001
  !        -1  -00001
  !         0  000000
  !      1952  001952
  !    123456  123456
  !   1234567  ******  <-- Not enough room!
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
  !
  !    Output, character ( len = * ) S, the representation of the integer.
  !    The integer will be right justified, and zero padded.
  !    If there is not enough space, the string will be filled with stars.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) ipos
  integer   ( kind = 4 ) ival
  character ( len = * )  s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
     return
  end if
  !
  !  Make a copy of the integer.
  !
  ival = intval
  !
  !  Handle the negative sign.
  !
  if ( ival < 0 ) then

     if ( ihi <= 1 ) then
        s(1:1) = '*'
        return
     end if

     ival = -ival
     s(1:1) = '-'
     ilo = 2

  end if
  !
  !  Working from right to left, strip off the digits of the integer
  !  and place them into S(ILO:IHI).
  !
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

     idig = mod ( ival, 10 )
     ival = ival / 10

     if ( ipos < ilo ) then
        do i = 1, ihi
           s(i:i) = '*'
        end do
        return
     end if

     call digit_to_ch ( idig, c )

     s(ipos:ipos) = c
     ipos = ipos - 1

  end do
  !
  !  Fill the empties with zeroes.
  !
  do i = ilo, ipos
     s(i:i) = '0'
  end do

  return
end subroutine i4_to_s_zero
function i4_to_s32 ( i4 )

  !*****************************************************************************80
  !
  !! I4_TO_S32 converts an I4 to an S32.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ).
  !
  !    An S32 is a character ( len = 32 ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, the integer to be coded.
  !
  !    Output, character ( len = 32 ) I4_TO_S32, the character variable that
  !    corresponds to the integer.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) i4_copy
  character ( len = 32 ) i4_to_s32
  character ( len = 32 ) s32

  i4_copy = abs ( i4 )
  !
  !  Binary digits:
  !
  do i = 32, 2, -1

     if ( mod ( i4_copy, 2 ) == 1 ) then
        s32(i:i) = '1'
     else
        s32(i:i) = '0'
     end if

     i4_copy = i4_copy / 2

  end do
  !
  !  Sign bit
  !
  s32(1:1) = '0'
  !
  !  If original number was negative, then reverse all bits.
  !
  if ( i4 < 0 ) then
     do i = 1, 32
        if ( s32(i:i) == '0' ) then
           s32(i:i) = '1'
        else
           s32(i:i) = '0'
        end if
     end do
  end if

  i4_to_s32 = s32

  return
end function i4_to_s32
subroutine i4_to_unary ( i4, s )

  !*****************************************************************************80
  !
  !! I4_TO_UNARY produces the "base 1" representation of an I4.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I4, an integer to be represented.
  !
  !    Output, character ( len = * ) S, the unary representation.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )
  s = ' '

  if ( i4 < 0 ) then

     if ( s_length < i4 + 1 ) then
        s = '?'
        return
     end if

     s(1:1) = '-'

     do i = 2, abs ( i4 ) + 1
        s(i:i) = '1'
     end do

  else if ( i4 == 0 ) then

     s = '0'

  else if ( 0 < i4 ) then

     if ( s_length < i4 ) then
        s = '?'
        return
     end if

     do i = 1, i4
        s(i:i) = '1'
     end do

  end if

  s = adjustr ( s )

  return
end subroutine i4_to_unary
function i4_to_uudecode ( i )

  !*****************************************************************************80
  !
  !! I4_TO_UUDECODE returns the I-th character in the UUDECODE encoding.
  !
  !  Example:
  !
  !    I  I4_TO_UUDECODE
  !
  !    1  '`'
  !    2  '!'
  !    3  '"'
  !   ..
  !   64  '_'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the index of the character.
  !    1 <= I <= 64.
  !
  !    Output, character I4_TO_UUDECODE, the requested character.
  !
  implicit none

  integer   ( kind = 4 ) i
  character              i4_to_uudecode
  character ( len = 64 ), parameter :: string = &
       '`!"#$%&''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_'

  if ( 1 <= i .and. i <= 64 ) then
     i4_to_uudecode = string(i:i)
  else
     i4_to_uudecode = ' '
  end if

  return
end function i4_to_uudecode
function i4_to_xxdecode ( i )

  !*****************************************************************************80
  !
  !! I4_TO_XXDECODE returns the I-th character in the XXDECODE encoding.
  !
  !  Example:
  !
  !    I  I4_TO_UUDECODE
  !
  !    1  '+'
  !    2  '-'
  !    3  '0'
  !   ..
  !   64  'z'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the index of the character.
  !    1 <= I <= 64.
  !
  !    Output, character I4_TO_XXDECODE, the requested character.
  !
  implicit none

  integer   ( kind = 4 ) i
  character              i4_to_xxdecode
  character ( len = 64 ), parameter :: string = &
       '+-0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

  if ( 1 <= i .and. i <= 64 ) then
     i4_to_xxdecode = string(i:i)
  else
     i4_to_xxdecode = ' '
  end if

  return
end function i4_to_xxdecode
function i4_uniform ( a, b, seed )

  !*****************************************************************************80
  !
  !! I4_UNIFORM returns a scaled pseudorandom I4.
  !
  !  Discussion:
  !
  !    An I4 is an integer value.
  !
  !    The pseudorandom number will be scaled to be uniformly distributed
  !    between A and B.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 November 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Pierre L'Ecuyer,
  !    Random Number Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley Interscience, page 95, 1998.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) A, B, the limits of the interval.
  !
  !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
  !    should NOT be 0.  On output, SEED has been updated.
  !
  !    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
  !
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
     write ( *, '(a)' ) '  Input value of SEED = 0.'
     stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
     seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
  !
  !  Scale R to lie between A-0.5 and B+0.5.
  !
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
       +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
  !
  !  Use rounding to convert R to an integer between A and B.
  !
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end function i4_uniform
subroutine i4vec_indicator ( n, a )

  !*****************************************************************************80
  !
  !! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of elements of A.
  !
  !    Output, integer ( kind = 4 ) A(N), the array to be initialized.
  !
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
     a(i) = i
  end do

  return
end subroutine i4vec_indicator
subroutine i4vec_print ( n, a, title )

  !*****************************************************************************80
  !
  !! I4VEC_PRINT prints an I4VEC.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of components of the vector.
  !
  !    Input, integer ( kind = 4 ) A(N), the vector to be printed.
  !
  !    Input, character ( len = * ) TITLE, a title to be printed first.
  !    TITLE may be blank.
  !
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) big
  integer   ( kind = 4 ) i
  character ( len = * )  title

  if ( title /= ' ' ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
     do i = 1, n
        write ( *, '(i8,1x,i4)' ) i, a(i)
     end do
  else if ( big < 1000000 ) then
     do i = 1, n
        write ( *, '(i8,1x,i7)' ) i, a(i)
     end do
  else
     do i = 1, n
        write ( *, '(i8,i11)' ) i, a(i)
     end do
  end if

  return
end subroutine i4vec_print
subroutine i4vec_to_ch4vec ( n, i4vec, s )

  !*****************************************************************************80
  !
  !! I4VEC_TO_CH4VEC converts an I4VEC into a string.
  !
  !  Discussion:
  !
  !    This routine can be useful when trying to read character data from an
  !    unformatted direct access file, for instance.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of integers.
  !
  !    Input, integer ( kind = 4 ) I4VEC(N), the integers.
  !
  !    Output, character ( len = * ) S, a string of 4 * N characters
  !    representing the integer information.
  !
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4vec(n)
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) len_s
  character ( len = * )  s

  len_s = len ( s )

  if ( len_s < 4 * n ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4VEC_TO_CH4VEC - Fatal error!'
     write ( *, '(a)' ) '  String S is too small for the data.'
     stop
  end if

  s(1:4*n) = ' '

  do i = 1, n
     j = 4 * ( i - 1 ) + 1
     call i4_to_ch4 ( i4vec(i), s(j:j+3) )
  end do

  return
end subroutine i4vec_to_ch4vec
function ic_to_ibraille ( ic )

  !*****************************************************************************80
  !
  !! IC_TO_IBRAILLE converts an ASCII integer code to a Braille code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IC, the integer code for the ASCII character.
  !
  !    Output, integer ( kind = 4 ) IC_TO_IBRAILLE, the integer code for
  !    the Braille character, or -1 if no corresponding code is available.
  !
  implicit none

  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ic_to_ibraille
  integer ( kind = 4 ), parameter, dimension ( 0:255 ) :: junk = (/ &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       1, 33, 35, -1, -1, -1, 28, 36, 34, 34, -1, -1, 29, 37, 32, -1, &
       11, 02, 04, 05, 06, 07, 08, 09, 10, -1, 31, 30, -1, -1, -1, 35, &
       -1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, &
       17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, -1, -1, -1, -1, -1, &
       -1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, &
       17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /)

  if ( 0 <= ic .and. ic <= 255 ) then
     ic_to_ibraille = junk(ic)
  else
     ic_to_ibraille = -1
  end if

  return
end function ic_to_ibraille
function ic_to_iebcdic ( ic )

  !*****************************************************************************80
  !
  !! IC_TO_IEBCDIC converts an ASCII character code to an EBCDIC code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IC, the integer code for the ASCII character.
  !
  !    Output, integer ( kind = 4 ) IC_TO_IEBCDIC, the integer code for the
  !    EBCDIC character, or -1 if no corresponding EBCDIC code is available.
  !
  implicit none

  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ic_to_iebcdic
  integer ( kind = 4 ), parameter, dimension ( 0:255 ) :: junk = (/ &
       0,  1,  2,  3, 56, 45, 46, 47, 22,  5, 37, 11, 12, 13, 60, 61, &
       16, 17, 18, -1, -1, -1, 50, 38, 24, 25, 63, 39, 28, 29, 30, 31, &
       64, 90,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96, -1, 97, &
       240,241,242,243,244,245,246,247,248,249,122, 94, 76,126,110,111, &
       124,193,194,195,196,197,198,199,200,201,209,210,211,212,213,214, &
       215,216,217,226,227,228,229,230,231,232,233, -1, -1, -1, -1,109, &
       -1,129,130,131,132,133,134,135,136,137,145,146,147,148,149,150, &
       151,152,153,162,163,164,165,166,167,168,169, -1, 79, -1, -1,  7, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /)

  if ( 0 <= ic .and. ic <= 255 ) then
     ic_to_iebcdic = junk(ic)
  else
     ic_to_iebcdic = -1
  end if

  return
end function ic_to_iebcdic
function ic_to_imorse ( ic )

  !*****************************************************************************80
  !
  !! IC_TO_IMORSE converts an ASCII integer code to a Morse integer code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IC, the integer code for the ASCII character.
  !
  !    Output, integer ( kind = 4 ) IC_TO_IMORSE, the integer code for the
  !    Morse character, or -1 if no corresponding Morse code is available.
  !
  implicit none

  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ic_to_imorse
  integer ( kind = 4 ), parameter, dimension ( 0:255 ) :: junk = (/ &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       1, -1, 45, -1, -1, -1, -1, 42, -1, -1, -1, -1, 39, 43, 38, 44, &
       37, 28, 29, 30, 31, 32, 33, 34, 35, 36, 40, -1, -1, -1, -1, 41, &
       -1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, &
       17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, -1, -1, -1, -1, -1, &
       -1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, &
       17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /)

  if ( 0 <= ic .and. ic <= 255 ) then
     ic_to_imorse = junk(ic)
  else
     ic_to_imorse = -1
  end if

  return
end function ic_to_imorse
function ic_to_isoundex ( ic )

  !*****************************************************************************80
  !
  !! IC_TO_ISOUNDEX converts an ASCII integer code to a Soundex integer code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IC, the integer code for the ASCII character.
  !
  !    Output, integer ( kind = 4 ) IC_TO_ISOUNDEX, the integer code for the
  !    Soundex character, or -1 if no corresponding Soundex code is available.
  !
  implicit none

  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ic_to_isoundex
  integer ( kind = 4 ), parameter, dimension ( 0:255 ) :: junk = (/ &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, 48, 49, 50, 51, 48, 49, 50, 48, 48, 50, 50, 52, 53, 53, 48, &
       49, 50, 54, 50, 51, 48, 49, 48, 50, 48, 50, -1, -1, -1, -1, -1, &
       -1, 48, 49, 50, 51, 48, 49, 50, 48, 48, 50, 50, 52, 53, 53, 48, &
       49, 50, 54, 50, 51, 48, 49, 48, 50, 48, 50, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /)

  if ( 0 <= ic .and. ic <= 255 ) then
     ic_to_isoundex = junk(ic)
  else
     ic_to_isoundex = -1
  end if

  return
end function ic_to_isoundex
function iebcdic_to_ic ( iebcdic )

  !*****************************************************************************80
  !
  !! IEBCDIC_TO_IC converts an EBCDIC character code to ASCII.
  !
  !  Discussion:
  !
  !    What this actually means is the following:
  !
  !    If the letter "A" is entered into a file on an EBCDIC machine,
  !    it is coded internally as character 193.  Should this character
  !    be read on an ASCII machine, it will not be displayed as "A",
  !    but rather as an unprintable character!  But passing 193 in to
  !    IEBCDIC_TO_IC, out will come 65, the ASCII code for "A".  Thus, the
  !    correct character to display on an ASCII machine is
  !
  !      ACHAR ( IACHAR ( IEBCDIC_TO_IC ( EBCDIC Character ) ) ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IEBCDIC, the integer code for the EBCDIC
  !    character.
  !
  !    Output, integer ( kind = 4 ) IEBCDIC_TO_IC, the integer code for the
  !    ASCII character, or -1 if no corresponding ASCII code is available.
  !
  implicit none

  integer ( kind = 4 ) iebcdic
  integer ( kind = 4 ) iebcdic_to_ic
  integer ( kind = 4 ), parameter, dimension ( 0:255 ) :: junk = (/ &
       0,  1,  2,  3, -1,  9, -1,127, -1, -1, -1, 11, 12, 13, 14, 15, &
       16, 17, 18, -1, -1, -1,  8, -1, 24, 25, -1, -1, 28, 29, 30, 31, &
       -1, -1, -1, -1, -1, 10, 23, 27, -1, -1, -1, -1, -1,  5,  6,  7, &
       -1, -1, 22, -1, -1, -1, -1, -1,  4, -1, -1, -1, 14, 15, -1, 26, &
       32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 60, 40, 43,124, &
       38, -1, -1, -1, -1, -1, -1, -1, -1, -1, 33, 36, 42, 41, 59, -1, &
       45, 47, -1, -1, -1, -1, -1, -1, -1, -1, -1, 44, 37, 95, 62, 63, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 58, 35, 64, 39, 61, 34, &
       -1, 97, 98, 99,100,101,102,103,104,105, -1, -1, -1, -1, -1, -1, &
       -1,106,107,108,109,110,111,112,113,114, -1, -1, -1, -1, -1, -1, &
       -1, -1,115,116,117,118,119,120,121,122, -1, -1, -1, -1, -1, -1, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
       -1, 65, 66, 67, 68, 69, 70, 71, 72, 73, -1, -1, -1, -1, -1, -1, &
       -1, 74, 75, 76, 77, 78, 79, 80, 81, 82, -1, -1, -1, -1, -1, -1, &
       -1, -1, 83, 84, 85, 86, 87, 88, 89, 90, -1, -1, -1, -1, -1, -1, &
       48, 49, 50, 51, 52, 53, 54, 55, 56, 57, -1, -1, -1, -1, -1, -1 /)

  if ( 0 <= iebcdic .and. iebcdic <= 255 ) then
     iebcdic_to_ic = junk(iebcdic)
  else
     iebcdic_to_ic = -1
  end if

  return
end function iebcdic_to_ic
function isbn_to_i4 ( c )

  !*****************************************************************************80
  !
  !! ISBN_TO_I4 converts an ISBN character into an integer.
  !
  !  Discussion:
  !
  !    The characters '0' through '9' stand for themselves, but
  !    the character 'X' or 'x' stands for 10.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Book Industry Study Group,
  !    The Evolution in Product Identification:
  !    Sunrise 2005 and the ISBN-13,
  !    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
  !
  !  Parameters:
  !
  !    Input, character C, the ISBN character code to be converted.
  !
  !    Output, integer ( kind = 4 ) ISBN_TO_I4, the numeric value of the character
  !    code, between 0 and 10.  This value is returned as -1 if C is
  !    not a valid character code.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) isbn_to_i4

  if ( c == '0' ) then
     isbn_to_i4 = 0
  else if ( c == '1' ) then
     isbn_to_i4 = 1
  else if ( c == '2' ) then
     isbn_to_i4 = 2
  else if ( c == '3' ) then
     isbn_to_i4 = 3
  else if ( c == '4' ) then
     isbn_to_i4 = 4
  else if ( c == '5' ) then
     isbn_to_i4 = 5
  else if ( c == '6' ) then
     isbn_to_i4 = 6
  else if ( c == '7' ) then
     isbn_to_i4 = 7
  else if ( c == '8' ) then
     isbn_to_i4 = 8
  else if ( c == '9' ) then
     isbn_to_i4 = 9
  else if ( c == 'X' .or. c == 'x' ) then
     isbn_to_i4 = 10
  else
     isbn_to_i4 = -1
  end if

  return
end function isbn_to_i4
function istrcmp ( s1, s2 )

  !*****************************************************************************80
  !
  !! ISTRCMP compares two strings, returning +1, 0, or -1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to be compared.
  !
  !    Output, integer ( kind = 4 ) ISTRCMP:
  !    -1 if S1 < S2,
  !     0 if S1 = S2
  !    +1 if S2 < S1.
  !
  implicit none

  integer   ( kind = 4 ) istrcmp
  integer   ( kind = 4 ) nchar1
  integer   ( kind = 4 ) nchar2
  integer   ( kind = 4 ) s_length
  character ( len = * )  s1
  character ( len = * )  s2

  nchar1 = len ( s1 )
  nchar2 = len ( s2 )
  s_length = min ( nchar1, nchar2 )

  if ( llt ( s1(1:s_length), s2(1:s_length) ) ) then
     istrcmp = -1
  else if ( llt ( s2(1:s_length), s1(1:s_length) ) ) then
     istrcmp = 1
  else if ( s1(1:s_length) == s2(1:s_length) ) then

     if ( nchar1 == nchar2 ) then
        istrcmp = 0
     else if ( nchar1 < nchar2 ) then
        istrcmp = -1
     else
        istrcmp = 1
     end if

  end if

  return
end function istrcmp
function istrncmp ( s1, s2, nchar )

  !*****************************************************************************80
  !
  !! ISTRNCMP compares the start of two strings, returning +1, 0, or -1.
  !
  !  Discussion:
  !
  !    If either string is shorter than NCHAR characters, then it is
  !    treated as though it were padded with blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to be compared.
  !
  !    Input, integer ( kind = 4 ) NCHAR, the number of characters to be compared.
  !
  !    Output, integer ( kind = 4 ) ISTRNCMP:
  !    +1 if S1(1:NCHAR) is lexically greater than S2(1:NCHAR),
  !     0 if they are equal, and
  !    -1 if S1(1:NCHAR) is lexically less than S2(1:NCHAR).
  !
  implicit none

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) istrncmp
  integer   ( kind = 4 ) nchar
  integer   ( kind = 4 ) nchar1
  integer   ( kind = 4 ) nchar2
  character ( len = * )  s1
  character ( len = * )  s2
  !
  !  Figure out the maximum number of characters we will examine,
  !  which is the minimum of NCHAR and the lengths of the two
  !  strings.
  !
  istrncmp = 0

  nchar1 = len ( s1 )
  nchar2 = len ( s2 )

  do i = 1, nchar

     if ( i <= nchar1 ) then
        c1 = s1(i:i)
     else
        c1 = ' '
     end if

     if ( i <= nchar2 ) then
        c2 = s2(i:i)
     else
        c2 = ' '
     end if

     if ( llt ( c1, c2 ) ) then
        istrncmp = -1
        return
     else if ( lgt ( c1, c2 ) ) then
        istrncmp = 1
        return
     end if

  end do

  return
end function istrncmp
function len_nonnull ( s )

  !*****************************************************************************80
  !
  !! LEN_NONNULL returns the length of a string up to the last non-null character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to measure.
  !
  !    Output, integer ( kind = 4 ) LEN_NONNULL, the length of the string,
  !    up to the last non-null character.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) len_nonnull
  integer   ( kind = 4 ) len_s
  character, parameter :: NULL = achar ( 0 )
  character ( len = * ) s

  len_s = len ( s )

  do i = len_s, 1, -1
     if ( s(i:i) /= NULL ) then
        len_nonnull = i
        return
     end if
  end do

  len_nonnull = 0

  return
end function len_nonnull
function malphnum2 ( s )

  !*****************************************************************************80
  !
  !! MALPHNUM2 is TRUE if a string contains only alphanumerics and underscores.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    Alphanumeric characters are 'A' through 'Z', 'a' through 'z',
  !    '0' through '9' and the underscore character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, logical MALPHNUM2, is TRUE if the string contains only
  !    alphabetic characters, numerals, and underscores.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) itemp
  logical                malphnum2
  character ( len = * )  s

  malphnum2 = .false.

  do i = 1, len ( s )

     if ( s(i:i) /= '_' ) then

        itemp = iachar ( s(i:i) )

        if ( .not. ( 65 <= itemp .and. itemp <= 90 ) ) then
           if ( .not. ( 97 <= itemp .and. itemp <= 122 ) ) then
              if ( .not. ( 48 <= itemp .and. itemp <= 57 ) ) then
                 return
              end if
           end if
        end if
     end if

  end do

  malphnum2 = .true.

  return
end function malphnum2
subroutine military_to_ch ( military, ch )

  !*****************************************************************************80
  !
  !! MILITARY_TO_CH converts a Military code word to an ASCII character.
  !
  !  Example:
  !
  !    'Alpha'   'A'
  !    'Bravo'   'B'
  !    'Zulu'    'Z'
  !    'alpha'   'a'
  !    '7'       '7'
  !    '%'       '%'
  !    'Adam'    'A'
  !    'Anthrax' 'A'
  !    'amoeba'  'a'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = 8 ) MILITARY, the military code word.
  !
  !    Output, character CH, the ASCII character.  If MILITARY was not
  !    a recognized military code word, then CH is set to MILITARY(1:1).
  !
  implicit none

  integer ( kind = 4 ) a_to_i4
  character ch
  character ( len = 8 ), dimension ( 26 ) :: code = (/ &
       'alpha   ', 'bravo   ', 'charlie ', 'delta   ', 'echo    ', &
       'foxtrot ', 'golf    ', 'hotel   ', 'india   ', 'juliet  ', &
       'kilo    ', 'lima    ', 'mike    ', 'november', 'oscar   ', &
       'papa    ', 'quebec  ', 'romeo   ', 'sierra  ', 'tango   ', &
       'uniform ', 'victor  ', 'whiskey ', 'x-ray   ', 'yankee  ', &
       'zulu    ' /)
  integer   ( kind = 4 ) i
  character ( len = * )  military
  logical                s_eqi

  ch = military(1:1)

  i = a_to_i4 ( ch )

  if ( 1 <= i .and. i <= 26 ) then
     if ( s_eqi ( military, code(i) ) ) then
        ch = military(1:1)
     end if
  else if ( 27 <= i .and. i <= 52 ) then
     if ( s_eqi ( military, code(i-26) ) ) then
        ch = military(1:1)
     end if
  end if

  return
end subroutine military_to_ch
subroutine month_name_to_i4 ( month_name, month )

  !*****************************************************************************80
  !
  !! MONTH_NAME_TO_I4 returns the month number of a given month
  !
  !  Discussion:
  !
  !    Capitalization is ignored.  The month name has to match up to
  !    the unique beginning of a month name, and the rest is ignored.
  !    Here are the limits:
  !
  !      JAnuary
  !      February
  !      MARch
  !      APril
  !      MAY
  !      JUNe
  !      JULy
  !      AUgust
  !      September
  !      October
  !      November
  !      December
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) MONTH_NAME, a string containing a month
  !    name or abbreviation.
  !
  !    Output, integer ( kind = 4 ) MONTH, the number of the month,
  !    or -1 if the name could not be recognized.
  !
  implicit none

  integer   ( kind = 4 ) month
  character ( len = * )  month_name
  character ( len = 3 )  string

  string = month_name
  call s_cap ( string )

  if ( string(1:2) == 'JA' ) then
     month = 1
  else if ( string(1:1) == 'F' ) then
     month = 2
  else if ( string(1:3) == 'MAR' ) then
     month = 3
  else if ( string(1:2) == 'AP' ) then
     month = 4
  else if ( string(1:3) == 'MAY' ) then
     month = 5
  else if ( string(1:3) == 'JUN' ) then
     month = 6
  else if ( string(1:3) == 'JUL' ) then
     month = 7
  else if ( string(1:2) == 'AU' ) then
     month = 8
  else if ( string(1:1) == 'S' ) then
     month = 9
  else if ( string(1:1) == 'O' ) then
     month = 10
  else if ( string(1:1) == 'N' ) then
     month = 11
  else if ( string(1:1) == 'D' ) then
     month = 12
  else
     month = -1
  end if

  return
end subroutine month_name_to_i4
subroutine namefl ( s )

  !*****************************************************************************80
  !
  !! NAMEFL replaces "lastname, firstname" by "firstname lastname".
  !
  !  Discussion:
  !
  !    As part of the process, all commas and double blanks are
  !    removed, and the first character of the output string is
  !    never a blank.
  !
  !    A one-word name is left unchanged.
  !
  !  Example:
  !
  !      Input                      Output
  !
  !    Brown, Charlie             Charlie Brown
  !    Cher                       Cher
  !    Howell, James Thurston     James Thurston Howell
  !    Shakespeare Joe Bob        Joe Bob Shakespeare
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.
  !
  !    On input, a series of words separated by spaces.
  !
  !    On output, if S contained a single word, it is
  !    unchanged.  Otherwise, the first word has been moved
  !    to the end of S, and any trailing comma removed.
  !
  !    As part of this process, all double blanks are removed
  !    from S, and the output S never begins with
  !    a blank (unless the input S was entirely blank).
  !
  !    Any commas in the input string are deleted.
  !
  !    This routine cannot handle more than 255 characters in S.
  !
  implicit none

  character ( len = 255 ) s2
  integer   ( kind = 4 )  i
  character ( len = * )   s
  integer   ( kind = 4 )  s_length

  s2 = ' '
  !
  !  Remove all commas.
  !
  s_length = len_trim ( s )

  do i = 1, s_length
     if ( s(i:i) == ',') then
        s(i:i) = ' '
     end if
  end do
  !
  !  Remove double blanks.
  !  This also guarantees the string is flush left.
  !
  call s_blanks_delete ( s )
  !
  !  Get length of string.
  !
  s_length = len_trim ( s )

  if ( s_length <= 2 ) then
     return
  end if

  if ( 255 < s_length ) then
     s_length = len_trim ( s(1:255) )
  end if
  !
  !  Find the first blank in the string.
  !
  do i = 2, s_length - 1

     if ( s(i:i) == ' ' ) then
        s2(1:s_length-i) = s(i+1:s_length)
        s2(s_length-i+1:s_length-i+1) = ' '
        s2(s_length-i+2:s_length) = s(1:i-1)
        s = s2(1:s_length)
     end if

  end do

  return
end subroutine namefl
subroutine namelf ( s )

  !*****************************************************************************80
  !
  !! NAMELF replaces "firstname lastname" by "lastname, firstname".
  !
  !  Discussion:
  !
  !    A one-word name is left unchanged.
  !
  !  Example:
  !
  !      Input:                     Output:
  !
  !    Charlie Brown              Brown, Charlie
  !    Cher                       Cher
  !    James Thurston Howell      Howell, James Thurston
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.
  !
  !    On input, S contains a series of words separated by spaces.
  !
  !    On output, if S contained a single word, it is
  !    unchanged.  Otherwise, the last word has been moved
  !    to the beginning of S, and followed by a comma.
  !
  !    As part of this process, all double blanks are removed
  !    from S, and the output S never begins with
  !    a blank (unless the input S was entirely blank).
  !
  !    Moreover, any commas in the input string are deleted.
  !
  !    This routine cannot handle more than 255 characters
  !    in S.  If S is longer than that, only the
  !    first 255 characters will be considered.
  !
  implicit none

  character ( len = 255 ) s2
  integer   ( kind = 4 )  i
  character ( len = * )   s
  integer   ( kind = 4 )  s_length

  s2 = ' '
  !
  !  Remove all commas.
  !
  s_length = len_trim ( s )
  do i = 1, s_length
     if ( s(i:i) == ',' ) then
        s(i:i) = ' '
     end if
  end do
  !
  !  Remove all double blanks, and make string flush left.
  !
  call s_blanks_delete ( s )
  !
  !  Get length of string.
  !
  s_length = len_trim ( s )

  if ( s_length <= 2 ) then
     return
  end if

  if ( 255 < s_length ) then
     s_length = len_trim ( s(1:255) )
  end if
  !
  !  Find the last blank in the string.
  !
  do i = s_length, 2, -1
     if ( s(i:i) == ' ' ) then
        s2(1:s_length-i) = s(i+1:s_length)
        s2(s_length-i+1:s_length-i+2) = ', '
        s2(s_length-i+3:s_length+1) = s(1:i-1)
        s = s2(1:s_length+1)
     end if
  end do

  return
end subroutine namelf
subroutine namels ( name, ierror, rhs, value )

  !*****************************************************************************80
  !
  !! NAMELS reads a NAMELIST line, returning the variable name and value.
  !
  !  Discussion:
  !
  !    NAMELS is a simple program, and can only handle simple input.
  !    In particular, it cannot handle:
  !
  !      multiple assignments on one line,
  !      a single assignment extended over multiple lines,
  !      assignments to character or complex variables,
  !      assignments to arrays.
  !
  !    Typical input would be of the form:
  !
  !      name = value
  !
  !    including, for instance:
  !
  !      a = 1.0
  !      n = -17
  !      scale = +5.3E-2
  !
  !    Spaces are ignored, and case is not important.  Integer values
  !    will be returned as real, but this is never a
  !    problem as long as the integers are "small".
  !
  !    If a line begins with the character "#", it is assumed to be
  !    a comment, and is ignored.  IERROR is returned as 6.
  !
  !    If a line begins with the characters "end-of-input", it is
  !    assumed to be an "end-of-input" marker, and IERROR is returned
  !    as 7.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, character ( len = * ) NAME.
  !    NAME contains the left hand side of the assignment statement.
  !    Normally, this will be the name of a variable.
  !    If the input line was blank, then NAME will equal ' '.
  !    If an error occurred while trying to process the
  !    input line, NAME will contain the text of the line.
  !    If the line began with "#", then NAME will contain the
  !    text of the line.
  !    If the line equals "end-of-input", then NAME will contain
  !    the text of the line.
  !
  !    Output, integer ( kind = 4 ) IERROR.
  !    0, no errors were detected.
  !    1, the line was blank.
  !    2, the line did not contain an "=" sign.
  !    3, the line did not contain a variable name to the
  !       left of the "=" sign.
  !    4, the right hand side of the assignment did not make
  !       sense.
  !    5, end of input.
  !    6, the line began with "#", signifying a comment.
  !       The text of the line is returned in NAME.
  !    7, the line began with "end-of-input".
  !
  !    Output, character ( len = * ) RHS.
  !    RHS contains the right hand side of the assignment statement.
  !
  !    Output, real ( kind = 4 ) VALUE.
  !    VALUE contains the right hand side of the assignment statement.
  !    Normally, this will be a real value.
  !    But if the input line was blank, or if an error occurred
  !    while trying to process the input line, or if input
  !    terminated, then VALUE will simply be set to 0.
  !
  implicit none

  integer   ( kind = 4 ) iequal
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) length
  character ( len = 80 ) line
  character ( len = * )  name
  integer   ( kind = 4 ) pos
  character ( len = * )  rhs
  logical                s_eqi
  real      ( kind = 4 ) value
  !
  !  Set default values
  !
  ierror = 0
  name = ' '
  rhs = ' '
  value = 0
  !
  !  Read a line
  !
  read ( *, '(a)', iostat = ios ) line

  if ( ios /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'NAMELS - Reached end of input.'
     ierror = 5
     return
  end if
  !
  !  Empty lines are OK
  !
  if ( len_trim ( line ) <= 0 ) then
     ierror = 1
     return
  end if
  !
  !  Check for comment.
  !
  if ( line(1:1) == '#' ) then
     ierror = 6
     name = line
     return
  end if
  !
  !  Check for "end-of-line".
  !
  if ( s_eqi ( line, 'END-OF-INPUT' ) ) then
     ierror = 7
     name = line
     return
  end if
  !
  !  Does the line contain an = sign?
  !
  if ( index ( line, '=' ) <= 0 ) then
     ierror = 2
     value = 0
     name = line
     return
  end if
  !
  !  Find the name of the variable to be assigned.
  !
  iequal = index ( name, '=' )

  if ( 0 < iequal ) then
     rhs = line(iequal+1:)
  else
     rhs = line
  end if

  call s_before_ss_copy ( line, '=', name )
  call s_blank_delete ( name )

  if ( len_trim ( name ) <= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'NAMELS - Warning!'
     write ( *, '(a)' ) '  The following input line was ignored, because'
     write ( *, '(a)' ) '  there was no variable name on the left hand'
     write ( *, '(a)' ) '  side of the assignment statement:'
     write ( *, '(a)' ) line
     write ( *, '(a)' ) ' '
     ierror = 3
     return
  end if
  !
  !  Read the value, as a real number.
  !
  pos = index ( line, '=' )
  call s_to_r4 ( line(pos+1:), value, ierror, length )

  if ( ierror /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'NAMELS - Warning!'
     write ( *, '(a)' ) '  The following input line was ignored, because'
     write ( *, '(a)' ) '  the right hand side of the assignment '
     write ( *, '(a)' ) '  statement did not seem to make sense:'
     write ( *, '(a)' ) line
     write ( *, '(a)' ) ' '
     ierror = 4
  end if

  return
end subroutine namels
subroutine nexchr ( s, i, c )

  !*****************************************************************************80
  !
  !! NEXCHR returns the next nonblank character from a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, integer ( kind = 4 ) I.  If I is 0, then there were no
  !    nonblank characters in the string.  Otherwise I is
  !    the index of the first nonblank character in the string.
  !
  !    Output, character C, the first nonblank character in the string.
  !
  implicit none

  character              c
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_first_nonblank

  i = s_first_nonblank ( s )

  if ( 0 < i ) then
     c = s(i:i)
  else
     c = ' '
  end if

  return
end subroutine nexchr
subroutine nexstr ( s, nsub, isub, sub )

  !*****************************************************************************80
  !
  !! NEXSTR returns the next nonblank characters from a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Input, integer ( kind = 4 ) NSUB, the number of nonblank characters
  !    desired.
  !
  !    Output, integer ( kind = 4 ) ISUB, the index of the NSUB-th nonblank
  !    character.  However, if ISUB is 0, there were NO nonblank
  !    characters.  And if there are less than NSUB nonblank characters
  !    ISUB is the location of the last one of them.
  !
  !    Output, character ( len = NSUB ) SUB, the first NSUB nonblanks.
  !
  implicit none

  integer   ( kind = 4 )   nsub

  integer   ( kind = 4 )   i
  integer   ( kind = 4 )   isub
  integer   ( kind = 4 )   jsub
  integer   ( kind = 4 )   s_first_nonblank
  character ( len = * )    s
  character ( len = nsub ) sub

  sub = ' '
  isub = 0

  do i = 1, nsub

     jsub = s_first_nonblank ( s(isub+1:) )

     if ( jsub <= 0 ) then
        return
     end if

     isub = isub + jsub

     sub(i:i) = s(isub:isub)

  end do

  return
end subroutine nexstr
subroutine number_inc ( s )

  !*****************************************************************************80
  !
  !! NUMBER_INC increments the integer represented by a string.
  !
  !  Example:
  !
  !    Input      Output
  !    -----      ------
  !    '17'       '18'
  !    'cat3'     'cat4'
  !    '2for9'    '3for0'
  !    '99thump'  '00thump'
  !
  !  Discussion:
  !
  !    If the string contains characters that are not digits, they will
  !    simply be ignored.  If the integer is all 9's on input, then
  !    the output will be all 0's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 January 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a string representing an integer.
  !
  implicit none

  logical                ch_is_digit
  integer   ( kind = 4 ) i
  character ( len = * )  s

  do i = len ( s ), 1, -1

     if ( ch_is_digit ( s(i:i) ) ) then

        call digit_inc ( s(i:i) )

        if ( s(i:i) /= '0' ) then
           return
        end if

     end if

  end do

  return
end subroutine number_inc
subroutine oct_to_i4 ( s, intval )

  !*****************************************************************************80
  !
  !! OCT_TO_I4 converts an octal string to its integer value.
  !
  !  Warning:
  !
  !    If too many digits are strung together, the computation will overflow.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string of digits.
  !
  !    Output, integer ( kind = 4 ) INTVAL, the corresponding integer value.
  !
  implicit none

  integer   ( kind = 4 ) first
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) idig
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) isgn
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )
  !
  !  Determine if there is a plus or minus sign.
  !
  isgn = 1
  first = s_length

  do i = 1, s_length - 1

     if ( s(i:i) == '-' ) then
        isgn = -1
     else if ( s(i:i) == '+' ) then
        isgn = + 1
     else if ( s(i:i) /= ' ' ) then
        first = i
        exit
     end if

  end do
  !
  !  Read the numeric portion of the string.
  !
  intval = 0

  do i = first, s_length
     call ch_to_digit_oct ( s(i:i), idig )
     intval = intval * 8 + idig
  end do

  intval = isgn * intval

  return
end subroutine oct_to_i4
subroutine perm_check ( n, p, ierror )

  !*****************************************************************************80
  !
  !! PERM_CHECK checks that a vector represents a permutation.
  !
  !  Discussion:
  !
  !    The routine verifies that each of the integers from 1
  !    to N occurs among the N entries of the permutation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 August 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries.
  !
  !    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, the array does represent a permutation.
  !    nonzero, the array does not represent a permutation.  The smallest
  !    missing value is equal to IERROR.
  !
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

     ierror = iseek

     do ifind = 1, n
        if ( p(ifind) == iseek ) then
           ierror = 0
           exit
        end if
     end do

     if ( ierror /= 0 ) then
        return
     end if

  end do

  return
end subroutine perm_check
subroutine perm_inverse3 ( n, perm, perm_inv )

  !*****************************************************************************80
  !
  !! PERM_INVERSE3 produces the inverse of a given permutation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 October 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of items permuted.
  !
  !    Input, integer ( kind = 4 ) PERM(N), a permutation.
  !
  !    Output, integer ( kind = 4 ) PERM_INV(N), the inverse permutation.
  !
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)

  do i = 1, n
     perm_inv(perm(i)) = i
  end do

  return
end subroutine perm_inverse3
subroutine perm_uniform ( n, seed, p )

  !*****************************************************************************80
  !
  !! PERM_UNIFORM selects a random permutation of N objects.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 April 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Albert Nijenhuis, Herbert Wilf,
  !    Combinatorial Algorithms for Computers and Calculators,
  !    Academic Press, 1978, second edition,
  !    ISBN 0-12-519260-6,
  !    LC: QA164.N54.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
  !
  !    Input/output, integer ( kind = 4 ) SEED, a seed for the random
  !    number generator.
  !
  !    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
  !    location of the object originally at I.
  !
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call i4vec_indicator ( n, p )

  do i = 1, n
     j = i4_uniform ( i, n, seed )
     call i4_swap ( p(i), p(j) )
  end do

  return
end subroutine perm_uniform
subroutine r4_to_b4_ieee ( r, word )

  !*****************************************************************************80
  !
  !! R4_TO_B4_IEEE converts an R4 to a 4 byte IEEE word.
  !
  !  Discussion:
  !
  !    This routine does not seem to working reliably for unnormalized data.
  !
  !  Example:
  !
  !    0 00000000 00000000000000000000000 =  0
  !    1 00000000 00000000000000000000000 = -0
  !
  !    0 11111111 00000000000000000000000 =  Infinity
  !    1 11111111 00000000000000000000000 = -Infinity
  !
  !    0 11111111 00000100000000000000000 = NaN
  !    1 11111111 00100010001001010101010 = NaN
  !
  !    0 01111110 00000000000000000000000 = +1 * 2**(126-127) * 1.0   = 0.5
  !    0 01111111 00000000000000000000000 = +1 * 2**(127-127) * 1.0   = 1
  !    0 10000000 00000000000000000000000 = +1 * 2**(128-127) * 1.0   = 2
  !    0 10000001 00000000000000000000000 = +1 * 2**(129-127) * 1.0   = 4
  !
  !    0 10000001 10100000000000000000000 = +1 * 2**(129-127) * 1.101 =  6.5
  !    1 10000001 10100000000000000000000 = -1 * 2**(129-127) * 1.101 = -6.5
  !
  !    0 00000001 00000000000000000000000 = +1 * 2**(  1-127) * 1.0 = 2**(-126)
  !    0 00000000 10000000000000000000000 = +1 * 2**(  0-126) * 0.1 = 2**(-127)
  !    0 00000000 00000000000000000000001 = +1 * 2**(  0-126) *
  !                                          0.00000000000000000000001 =
  !                                          2**(-149)  (Smallest positive value)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 November 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    IEEE Standards Committee 754,
  !    IEEE Standard for Binary Floating Point Arithmetic,
  !    ANSI/IEEE Standard 754-1985,
  !    SIGPLAN Notices,
  !    Volume 22, Number 2, 1987, pages 9-25.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R, the real number to be converted.
  !
  !    Output, integer ( kind = 4 ) WORD, the IEEE representation of the number.
  !
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  real    ( kind = 4 ) r
  real    ( kind = 4 ) r_copy
  integer ( kind = 4 ) s
  integer ( kind = 4 ) word

  r_copy = r
  !
  !  Determine S, the sign bit.
  !
  if ( 0.0E+00 <= r_copy ) then
     s = 0
  else
     s = 1
     r_copy = -r_copy
  end if
  !
  !  Determine E, the exponent.
  !  (FOR NOW, IGNORE UNNORMALIZED NUMBERS)
  !
  e = 0

  if ( r == 0.0E+00 ) then

  else

     do while ( 2.0E+00 <= r_copy )
        e = e + 1
        r_copy = r_copy / 2.0E+00
     end do

     do while ( r_copy < 1.0E+00 .and. -127 < e )
        e = e - 1
        r_copy = r_copy * 2.0E+00
     end do

     e = e + 127

  end if
  !
  !  Determine F, the fraction.
  !
  if ( r == 0.0E+00 ) then

     f = 0

  else if ( 0 < e) then

     r_copy = r_copy - 1.0E+00
     f = int ( r_copy * 2.0E+00**23 )

  else if ( e == 0 ) then

     f = int ( r_copy * 2.0E+00**23 )

  end if
  !
  !  Set the bits corresponding to S, E, F.
  !
  call mvbits ( s, 0,  1, word, 31 )
  call mvbits ( e, 0,  8, word, 23 )
  call mvbits ( f, 0, 23, word,  0 )

  return
end subroutine r4_to_b4_ieee
subroutine r4_to_binary ( r, s )

  !*****************************************************************************80
  !
  !! R4_TO_BINARY represents an R4 as a string of binary digits.
  !
  !  Discussion:
  !
  !    No check is made to ensure that the string is long enough.
  !
  !    The binary digits are a faithful representation of the real
  !    number in base 2.
  !
  !  Example:
  !
  !      R           S
  !
  !    -10.75000    -1010.11
  !      0.4218750  0.011011
  !      0.3333333  0.01010101010101010101010
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 August 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R, the real number.
  !
  !    Output, character ( len = * ) S, the binary representation.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iexp
  integer   ( kind = 4 ) ilo
  real      ( kind = 4 ) r
  real      ( kind = 4 ) rcopy
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  if ( s_length < 1 ) then
     return
  end if

  rcopy = r
  s = ' '

  if ( rcopy == 0.0E+00 ) then
     s = '0'
     return
  end if

  ilo = 0

  if ( rcopy < 0.0E+00 ) then
     ilo = 1
     s(ilo:ilo) = '-'
     rcopy = -rcopy
  end if
  !
  !  Figure out the divisor.
  !
  iexp = 0

  do while ( 1.0E+00 <= rcopy )
     rcopy = rcopy / 2.0E+00
     iexp = iexp + 1
  end do

  do while ( rcopy < 0.5E+00 )
     rcopy = rcopy * 2.0E+00
     iexp = iexp - 1
  end do
  !
  !  Now 0.5 <= RCOPY < 1.
  !
  !  If IEXP < 0, print leading zeroes.
  !
  if ( iexp == 0 ) then
     ilo = ilo + 1
     s(ilo:ilo) = '0'
  else if ( iexp < 0 ) then
     ilo = ilo + 1
     s(ilo:ilo) = '0'
     ilo = ilo + 1
     s(ilo:ilo) = '.'
     do i = 1, -iexp
        ilo = ilo + 1
        s(ilo:ilo) = '0'
     end do

  end if
  !
  !  Now repeatedly double RCOPY.
  !  Every time you exceed 1, that's a '1' digit.
  !
  iexp = iexp + 1

  do

     rcopy = 2.0E+00 * rcopy

     iexp = iexp - 1

     if ( iexp == 0 ) then
        ilo = ilo + 1
        s(ilo:ilo) = '.'
        if ( s_length <= ilo ) then
           return
        end if
     end if

     ilo = ilo + 1

     if ( 1.0E+00 <= rcopy ) then
        rcopy = rcopy - 1.0E+00
        s(ilo:ilo) = '1'
     else
        s(ilo:ilo) = '0'
     end if

     if ( s_length <= ilo ) then
        return
     end if

     if ( rcopy == 0.0E+00 ) then
        exit
     end if

  end do

  return
end subroutine r4_to_binary
subroutine r4_to_ch4 ( r4, ch4 )

  !*****************************************************************************80
  !
  !! R4_TO_CH4 converts an R4 to a 4 character string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R4, the real value.
  !
  !    Output, character ( len = 4 ) CH4, a corresponding character value.
  !
  implicit none

  character ( len = 4 )  ch4
  integer   ( kind = 4 ) i4
  integer   ( kind = 4 ) j1
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j3
  integer   ( kind = 4 ) j4
  real      ( kind = 4 ) r4

  i4 = transfer ( r4, i4 )

  j1 = ibits ( i4,  0, 8 )
  j2 = ibits ( i4,  8, 8 )
  j3 = ibits ( i4, 16, 8 )
  j4 = ibits ( i4, 24, 8 )

  ch4 = achar ( j1 ) // achar ( j2 ) // achar ( j3 ) // achar ( j4 )

  return
end subroutine r4_to_ch4
subroutine r4_to_flt ( r4, isgn, mant, iexp, ndig )

  !*****************************************************************************80
  !
  !! R4_TO_FLT computes the scientific representation of an R4.
  !
  !  Discussion:
  !
  !    The routine is given a real number R and computes a sign ISGN,
  !    an integer mantissa MANT and an integer exponent IEXP so
  !    that
  !
  !      R4 = ISGN * MANT * 10 ** IEXP
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R4, the real number whose scientific
  !    representation is desired.
  !
  !    Output, integer ( kind = 4 ) ISGN, the sign of the number:
  !    -1, if R4 is negative.
  !     0, if R4 is zero.
  !     1, if R4 is positive.
  !
  !    Output, integer ( kind = 4 ) MANT, the mantissa of the representation.
  !    This is an integer between 0 and 10**NDIG, that is,
  !    0 <= MANT < 10**NDIG.
  !
  !    Output, integer ( kind = 4 ) IEXP, the exponent of 10 that multiplies MULT.
  !
  !    Input, integer ( kind = 4 ) NDIG, the number of decimal digits.
  !
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) mant
  integer ( kind = 4 ) ndig
  real    ( kind = 4 ) rmant
  real    ( kind = 4 ) r4

  mant = 0
  iexp = 0
  isgn = 0
  !
  !  Find the first digit.
  !  That is, write the value in the form RMANT * 10**IEXP
  !  where 1/10 < RMANT <= 1.
  !
  if ( r4 == 0.0E+00 ) then
     return
  else if ( r4 < 0.0E+00 ) then
     isgn = -1
     rmant = abs ( r4 )
  else
     isgn = 1
     rmant = r4
  end if

  do while ( 1.0E+00 < rmant )
     rmant = rmant / 10.0E+00
     iexp = iexp + 1
  end do

  do while ( rmant <= 0.1E+00 )
     rmant = rmant * 10.0E+00
     iexp = iexp - 1
  end do
  !
  !  Now read off NDIG digits of RMANT.
  !
  do i = 1, ndig
     rmant = rmant * 10.0E+00
     idig = int ( rmant )
     rmant = rmant - idig
     mant = 10 * mant + idig
     iexp = iexp - 1
  end do
  !
  !  Now do rounding.
  !
  idig = int ( rmant * 10.0E+00 )
  mant = 10 * mant + idig
  mant =  nint ( real ( mant, kind = 4 ) / 10.0E+00 )
  !
  !  Now chop off trailing zeroes.
  !
  do while ( mod ( mant, 10 ) == 0 )
     mant = mant / 10
     iexp = iexp + 1
  end do

  return
end subroutine r4_to_flt
subroutine r4_to_s_left ( r4, s )

  !*****************************************************************************80
  !
  !! R4_TO_S_LEFT writes an R4 into a left justified character string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 December 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R4, the real number to be written into the string.
  !
  !    Output, character ( len = * ) S, the string into which
  !    the real number is to be written.
  !
  implicit none

  real      ( kind = 4 ) r4
  character ( len = * )  s
  character ( len = 14 ) s2

  if ( real ( int ( r4 ), kind = 4 ) == r4 ) then
     write ( s2, '(i14)' ) int ( r4 )
  else if ( abs ( r4 ) < 999999.5E+00 ) then
     write ( s2, '(f14.6)' ) r4
  else
     write ( s2, '(g14.6)' ) r4
  end if

  s = adjustl ( s2 )

  return
end subroutine r4_to_s_left
subroutine r4_to_s_right ( r4, s )

  !*****************************************************************************80
  !
  !! R4_TO_S_RIGHT writes an R4 into a right justified character string.
  !
  !  Discussion:
  !
  !    Thanks to Bill Richmond for pointing out a programming error
  !    that stored the data in S2, and then failed to copy it to the
  !    output quantity S.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R4, the real number to be written into the string.
  !
  !    Output, character ( len = * ) S, the string into which
  !    the real number is to be written.
  !
  implicit none

  real      ( kind = 4 ) r4
  character ( len = * )  s
  character ( len = 14 ) s2

  if ( real ( int ( r4 ), kind = 4 ) == r4 ) then
     write ( s2, '(i14)' ) int ( r4 )
  else if ( abs ( r4 ) < 999999.5E+00 ) then
     write ( s2, '(f14.6)' ) r4
  else
     write ( s2, '(g14.6)' ) r4
  end if

  s = ' '
  s(1:14) = s2(1:14)

  call s_adjustr ( s )

  return
end subroutine r4_to_s_right
function r4_to_s32 ( r4 )

  !*****************************************************************************80
  !
  !! R4_TO_S32 encodes an R4 as 32 characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R4, the real number to be coded.
  !
  !    Output, character ( len = 32 ) R4_TO_S32, the character variable that
  !    corresponds to the real number.
  !
  implicit none

  character ( len = 32 ) chr32
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iexp
  integer   ( kind = 4 ) ii
  integer   ( kind = 4 ) j
  real      ( kind = 4 ) r4
  character ( len = 32 ) r4_to_s32
  real      ( kind = 4 ) rcopy

  rcopy = r4
  !
  !  Sign bit
  !
  if ( rcopy < 0.0E+00 ) then
     chr32(1:1) = '1'
  else
     chr32(1:1) = '0'
  end if

  rcopy = abs ( rcopy )
  !
  !  Exponent: 'excess 128' format, legal values of IEXP are 1 to 255.
  !
  if ( rcopy == 0.0E+00 ) then

     iexp = 0

  else

     iexp = 128

     if ( rcopy < 1.0E+00 ) then

        do while ( 1 < iexp )
           rcopy = 2.0E+00 * rcopy
           iexp = iexp - 1
        end do

     else if ( 2.0E+00 <= rcopy ) then

        do while ( iexp < 255 )
           rcopy = 0.5E+00 * rcopy
           iexp = iexp + 1
        end do

     end if

  end if
  !
  !  Write characters 2 through 9 that represent exponent.
  !
  do i = 1, 8
     ii = 10 - i
     j = mod ( iexp, 2 )
     iexp = iexp / 2
     if ( j == 0 ) then
        chr32(ii:ii) = '0'
     else
        chr32(ii:ii) = '1'
     end if
  end do
  !
  !  Write mantissa in positions 10 through 32.
  !  Note that, unless exponent equals 0, the most significant bit is
  !  assumed to be 1 and hence is not stored.
  !
  if ( rcopy /= 0.0E+00 ) then
     rcopy = rcopy - 1.0E+00
  end if

  do i = 10, 32
     rcopy = 2.0E+00 * rcopy
     if ( 1.0E+00 <= rcopy ) then
        chr32(i:i) = '1'
        rcopy = rcopy - 1.0E+00
     else
        chr32(i:i) = '0'
     end if
  end do

  r4_to_s32 = chr32

  return
end function r4_to_s32
subroutine r4_to_sef ( r4, s, e, f )

  !*****************************************************************************80
  !
  !! R4_TO_SEF represents an R4 as R = S * 2**E * F.
  !
  !  Discussion:
  !
  !    Assuming no arithmetic problems, in fact, this equality should be
  !    exact, that is, S, E and F should exactly express the value
  !    as stored on the computer.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 November 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 4 ) R4, the real number.
  !
  !    Output, integer ( kind = 4 ) S, the sign bit:
  !    0, if R is nonnegative;
  !    1, if R is negative.
  !
  !    Output, integer ( kind = 4 ) E, the exponent base 2.
  !
  !    Output, integer ( kind = 4 ) F, the mantissa.
  !
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  real    ( kind = 4 ) r4
  real    ( kind = 4 ) r4_copy
  integer ( kind = 4 ) s

  if ( r4 == 0.0E+00 ) then
     s = 0
     e = 0
     f = 0
     return
  end if

  r4_copy = r4
  !
  !  Set S.
  !
  if ( 0.0E+00 <= r4_copy ) then
     s = 0
  else
     s = 1
     r4_copy = -r4_copy
  end if
  !
  !  Extracting the exponent leaves 0.5 <= R4_COPY < 1.0.
  !
  e = 0

  do while ( r4_copy < 0.5E+00 )
     r4_copy = r4_copy * 2.0E+00
     e = e - 1
  end do

  do while ( 1.0E+00 <= r4_copy )
     r4_copy = r4_copy / 2.0E+00
     e = e + 1
  end do
  !
  !  Get the binary mantissa, adjusting the exponent as you go.
  !
  f = 0
  e = e + 1

  do

     f = 2 * f
     e = e - 1

     if ( 1.0E+00 <= r4_copy ) then
        f = f + 1
        r4_copy = r4_copy - 1.0E+00
     end if

     if ( r4_copy == 0.0E+00 ) then
        exit
     end if

     r4_copy = 2.0E+00 * r4_copy

  end do

  return
end subroutine r4_to_sef
function r4_uniform_01 ( seed )

  !*****************************************************************************80
  !
  !! R4_UNIFORM_01 returns a unit pseudorandom R4.
  !
  !  Discussion:
  !
  !    An R4 is a real ( kind = 4 ) value.
  !
  !    This routine implements the recursion
  !
  !      seed = 16807 * seed mod ( 2**31 - 1 )
  !      r4_uniform_01 = seed / ( 2**31 - 1 )
  !
  !    The integer arithmetic never requires more than 32 bits,
  !    including a sign bit.
  !
  !    If the initial seed is 12345, then the first three computations are
  !
  !      Input     Output      R4_UNIFORM_01
  !      SEED      SEED
  !
  !         12345   207482415  0.096616
  !     207482415  1790989824  0.833995
  !    1790989824  2035175616  0.947702
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 May 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Pierre L'Ecuyer,
  !    Random Number Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley Interscience, page 95, 1998.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
  !    should NOT be 0.  On output, SEED has been updated.
  !
  !    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
  !    strictly between 0 and 1.
  !
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
     write ( *, '(a)' ) '  Input value of SEED = 0.'
     stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
     seed = seed + i4_huge ( )
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end function r4_uniform_01
subroutine r8_extract ( s, r8, ierror )

  !*****************************************************************************80
  !
  !! R8_EXTRACT "extracts" an R8 from the beginning of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S; on input, a string from
  !    whose beginning a real is to be extracted.  On output,
  !    the real, if found, has been removed.
  !
  !    Output, real ( kind = 8 ) R8.  If IERROR is 0, then R4 contains the
  !    next real read from the string; otherwise R4 is 0.
  !
  !    Output, integer ( kind = 4 ) IERROR.
  !    0, no error.
  !    nonzero, a real could not be extracted from the beginning of the
  !    string.  R4 is 0.0 and S is unchanged.
  !
  implicit none

  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) length
  real      ( kind = 8 ) r8
  character ( len = * )  s

  r8 = 0.0D+00

  call s_to_r8 ( s, r8, ierror, length )

  if ( ierror /= 0 .or. length == 0 ) then
     ierror = 1
     r8 = 0.0D+00
  else
     call s_shift_left ( s, length )
  end if

  return
end subroutine r8_extract
subroutine r8_input ( string, value, ierror )

  !*****************************************************************************80
  !
  !! R8_INPUT prints a prompt string and reads an R8 from the user.
  !
  !  Discussion:
  !
  !    An R8 is a real ( kind = 8 ) value.
  !
  !    If the input line starts with a comment character ('#') or is blank,
  !    the routine ignores that line, and tries to read the next one.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 March 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) STRING, the prompt string.
  !
  !    Output, real ( kind = 8 ) VALUE, the value input by the user.
  !
  !    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero
  !    if no error occurred.
  !
  implicit none

  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) length
  character ( len = 80 ) line
  character ( len = * )  string
  real      ( kind = 8 ) value

  ierror = 0
  value = huge ( value )
  !
  !  Write the prompt.
  !
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

     read ( *, '(a)', iostat = ierror ) line

     if ( ierror /= 0 ) then
        return
     end if
     !
     !  If the line begins with a comment character, go back and read the next line.
     !
     if ( line(1:1) == '#' ) then
        cycle
     end if

     if ( len_trim ( line ) == 0 ) then
        cycle
     end if
     !
     !  Extract integer information from the string.
     !
     call s_to_r8 ( line, value, ierror, length )

     if ( ierror /= 0 ) then
        value = huge ( value )
        return
     end if

     exit

  end do

  return
end subroutine r8_input
subroutine r8_next ( s, r, done )

  !*****************************************************************************80
  !
  !! R8_NEXT "reads" R8's from a string, one at a time.
  !
  !  Discussion:
  !
  !    An R8 is a real ( kind = 8 ) value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string, presumably containing real
  !    numbers.  These may be separated by spaces or commas.
  !
  !    Output, real ( kind = 8 ) R.  If DONE is FALSE, then R contains the
  !    "next" real value read from the string.  If DONE is TRUE, then
  !    R is zero.
  !
  !    Input/output, logical DONE.
  !    On input with a fresh string, the user should set DONE to TRUE.
  !    On output, the routine sets DONE to FALSE if another real
  !    value was read, or TRUE if no more reals could be read.
  !
  implicit none

  logical                done
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ), save :: next = 1
  real      ( kind = 8 ) r
  character ( len = * )  s

  r = 0.0D+00

  if ( done ) then
     next = 1
     done = .false.
  end if

  if ( len ( s ) < next ) then
     done = .true.
     return
  end if

  call s_to_r8 ( s(next:), r, ierror, length )

  if ( ierror /= 0 ) then
     done = .true.
     next = 1
  else if ( length == 0 ) then
     done = .true.
     next = 1
  else
     done = .false.
     next = next + length
  end if

  return
end subroutine r8_next
subroutine r8_to_binary ( r, s )

  !*****************************************************************************80
  !
  !! R8_TO_BINARY represents an R8 as a string of binary digits.
  !
  !  Discussion:
  !
  !    No check is made to ensure that the string is long enough.
  !
  !    The binary digits are a faithful representation of the real
  !    number in base 2.
  !
  !  Example:
  !
  !      R           S
  !
  !    -10.75000    -1010.11
  !      0.4218750  0.011011
  !      0.3333333  0.01010101010101010101010
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 June 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) R, the real number.
  !
  !    Output, character ( len = * ) S, the binary representation.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iexp
  integer   ( kind = 4 ) ilo
  real      ( kind = 8 ) r
  real      ( kind = 8 ) rcopy
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  if ( s_length < 1 ) then
     return
  end if

  rcopy = r
  s = ' '

  if ( rcopy == 0.0D+00 ) then
     s = '0'
     return
  end if

  ilo = 0

  if ( rcopy < 0.0D+00 ) then
     ilo = 1
     s(ilo:ilo) = '-'
     rcopy = -rcopy
  end if
  !
  !  Figure out the divisor.
  !
  iexp = 0

  do while ( 1.0D+00 <= rcopy )
     rcopy = rcopy / 2.0D+00
     iexp = iexp + 1
  end do

  do while ( rcopy < 0.5D+00 )
     rcopy = rcopy * 2.0D+00
     iexp = iexp - 1
  end do
  !
  !  Now 0.5 <= RCOPY < 1.
  !
  !  If IEXP < 0, print leading zeroes.
  !
  if ( iexp == 0 ) then
     ilo = ilo + 1
     s(ilo:ilo) = '0'
  else if ( iexp < 0 ) then
     ilo = ilo + 1
     s(ilo:ilo) = '0'
     ilo = ilo + 1
     s(ilo:ilo) = '.'
     do i = 1, -iexp
        ilo = ilo + 1
        s(ilo:ilo) = '0'
     end do

  end if
  !
  !  Now repeatedly double RCOPY.
  !  Every time you exceed 1, that's a '1' digit.
  !
  iexp = iexp + 1

  do

     rcopy = 2.0D+00 * rcopy

     iexp = iexp - 1

     if ( iexp == 0 ) then
        ilo = ilo + 1
        s(ilo:ilo) = '.'
        if ( s_length <= ilo ) then
           return
        end if
     end if

     ilo = ilo + 1

     if ( 1.0D+00 <= rcopy ) then
        rcopy = rcopy - 1.0D+00
        s(ilo:ilo) = '1'
     else
        s(ilo:ilo) = '0'
     end if

     if ( s_length <= ilo ) then
        return
     end if

     if ( rcopy == 0.0D+00 ) then
        exit
     end if

  end do

  return
end subroutine r8_to_binary
subroutine r8_to_s_left ( r8, s )

  !*****************************************************************************80
  !
  !! R8_TO_S_LEFT writes an R8 into a left justified string.
  !
  !  Discussion:
  !
  !    An R8 is a real ( kind = 8 ) value.
  !
  !    A 'G14.6' format is used with a WRITE statement.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) R8, the number to be written into the string.
  !
  !    Output, character ( len = * ) S, the string into which
  !    the real number is to be written.  If the string is less than 14
  !    characters long, it will will be returned as a series of asterisks.
  !
  implicit none

  integer   ( kind = 4 ) i
  real      ( kind = 8 ) r8
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character ( len = 14 ) s2

  s_length = len ( s )

  if ( s_length < 14 ) then

     do i = 1, s_length
        s(i:i) = '*'
     end do

  else if ( r8 == 0.0D+00 ) then
     s(1:14) = '     0.0      '
  else
     write ( s2, '(g14.6)' ) r8
     s(1:14) = s2
  end if
  !
  !  Shift the string left.
  !
  s = adjustl ( s )

  return
end subroutine r8_to_s_left
subroutine r8_to_s_right ( d, s )

  !*****************************************************************************80
  !
  !! R8_TO_S_LEFT writes an R8 into a right justified string.
  !
  !  Discussion:
  !
  !    An R8 is a real ( kind = 8 ) value.
  !
  !    A 'G14.6' format is used with a WRITE statement.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) D, the number to be written into the string.
  !
  !    Output, character ( len = * ) S, the string into which
  !    the real number is to be written.  If the string is less than 14
  !    characters long, it will will be returned as a series of asterisks.
  !
  implicit none

  real      ( kind = 8 ) d
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character ( len = 14 ) s2

  s_length = len ( s )

  if ( s_length < 14 ) then

     do i = 1, s_length
        s(i:i) = '*'
     end do

  else if ( d == 0.0D+00 ) then
     s(1:14) = '     0.0      '
  else
     write ( s2, '(g14.6)' ) d
     s(1:14) = s2
  end if
  !
  !  Shift the string right.
  !
  call s_adjustr ( s )

  return
end subroutine r8_to_s_right
function r8_uniform_01 ( seed )

  !*****************************************************************************80
  !
  !! R8_UNIFORM_01 returns a unit pseudorandom R8.
  !
  !  Discussion:
  !
  !    An R8 is a real ( kind = 8 ) value.
  !
  !    This routine implements the recursion
  !
  !      seed = 16807 * seed mod ( 2**31 - 1 )
  !      r8_uniform_01 = seed / ( 2**31 - 1 )
  !
  !    The integer arithmetic never requires more than 32 bits,
  !    including a sign bit.
  !
  !    If the initial seed is 12345, then the first three computations are
  !
  !      Input     Output      R8_UNIFORM_01
  !      SEED      SEED
  !
  !         12345   207482415  0.096616
  !     207482415  1790989824  0.833995
  !    1790989824  2035175616  0.947702
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 July 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Pierre L'Ecuyer,
  !    Random Number Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley Interscience, page 95, 1998.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
  !    NOT be 0. On output, SEED has been updated.
  !
  !    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
  !    strictly between 0 and 1.
  !
  implicit none

  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
     write ( *, '(a)' ) '  Input value of SEED = 0.'
     stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
     seed = seed + 2147483647
  end if
  !
  !  Although SEED can be represented exactly as a 32 bit integer,
  !  it generally cannot be represented exactly as a 32 bit real number!
  !
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01
subroutine r8vec_to_s ( n, x, s )

  !*****************************************************************************80
  !
  !! R8VEC_TO_S "writes" an R8VEC into a string.
  !
  !  Discussion:
  !
  !    An R8VEC is a vector of "real ( kind = 8 )" values.
  !
  !    The values will be separated by commas and a single space.
  !    If the string is too short, then data will be lost.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the dimension of X.
  !
  !    Input, real ( kind = 8 ) X(N), a vector to be written to a string.
  !
  !    Output, character ( len = * ) S, a string to which the real vector
  !    has been written.
  !
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  character ( len = * )  s
  character ( len = 14 ) s2
  real      ( kind = 8 ) x(n)

  do i = 1, n

     if ( x(i) == 0.0D+00 ) then
        s2 = '0'
     else if ( 1.0D+10 <= abs ( x(i) ) ) then
        write ( s2, '(g14.6)' ) x(i)
        call s_trim_zeros ( s2 )
     else if ( real ( int ( x(i) ), kind = 8 ) == x(i) ) then
        write ( s2, '(i12)' ) int ( x(i) )
     else
        write ( s2, '(g14.6)' ) x(i)
        call s_trim_zeros ( s2 )
     end if

     if ( i == 1 ) then
        s = adjustl ( s2 )
     else
        s = trim ( s ) // ', ' // adjustl ( s2 )
     end if

  end do

  return
end subroutine r8vec_to_s
subroutine ranger ( s, maxval, nval, ival )

  !*****************************************************************************80
  !
  !! RANGER "understands" a range defined by a string like '4:8'.
  !
  !  Discussion:
  !
  !    The range can be much more complicated, as in
  !
  !      '4:8 10 2 14:20'
  !
  !    or (commas are optional)
  !
  !      '4:8,10, 2 , 14:20'
  !
  !    RANGER will return the values
  !
  !      4, 5, 6, 7, 8, 10, 2, 14, 15, 16, 17, 18, 19, 20
  !
  !    0 and negative integers are acceptable.  So are pairs
  !    of values that are equal, as in '4:4', which just represents
  !    4, and pairs that represent descending sequences, as
  !    in '4:-2'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, contains a string of integers,
  !    representing themselves, and pairs of integers representing
  !    themselves and all integers between them.
  !
  !    Input, integer ( kind = 4 ) MAXVAL, the dimension of the IVAL vector,
  !    which represents the maximum number of integers that may
  !    be read from the string.
  !
  !    Output, integer ( kind = 4 ) NVAL, the number of integers read from
  !    the string.
  !
  !    Output, integer ( kind = 4 ) IVAL(MAXVAL).  The first NVAL entries of
  !    IVAL contain the integers read from the string.
  !
  implicit none

  integer   ( kind = 4 ) maxval

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) intval
  integer   ( kind = 4 ) ival(maxval)
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) lens
  integer   ( kind = 4 ) next
  integer   ( kind = 4 ) nval
  character ( len = * )  s

  nval = 0
  !
  !  Replace all commas by blanks.
  !
  call s_ch_blank ( s, ',' )
  !
  !  Replace multiple consecutive blanks by one blank.
  !
  call s_blanks_delete ( s )
  !
  !  Get the length of the string to the last nonblank.
  !
  lens = len_trim ( s )
  !
  !  Set a pointer to the next location to be examined.
  !
  next = 1

  do while ( next <= lens )
     !
     !  Find the next integer in the string.
     !
     call s_to_i4 ( s(next:), intval, ierror, length )

     if ( ierror /= 0 ) then
        return
     end if
     !
     !  Move the pointer.
     !
     next = next + length
     !
     !  If there's room, add the value to the list.
     !
     if ( maxval <= nval ) then
        return
     end if

     nval = nval + 1
     ival(nval) = intval
     !
     !  Have we reached the end of the string?
     !
     if ( lens < next ) then
        return
     end if
     !
     !  Skip past the next character if it is a space.
     !
     if ( s(next:next) == ' ' ) then
        next = next + 1
        if ( lens < next ) then
           return
        end if
     end if
     !
     !  Is the next character a colon?
     !
     if ( s(next:next) /= ':' ) then
        cycle
     end if
     !
     !  Increase the pointer past the colon.
     !
     next = next + 1

     if ( lens < next ) then
        return
     end if
     !
     !  Find the next integer in the string.
     !
     call s_to_i4 ( s(next:), intval, ierror, length )

     if ( ierror /= 0 ) then
        return
     end if
     !
     !  Move the pointer.
     !
     next = next + length
     !
     !  Generate integers between the two values.
     !
     ilo = ival(nval)

     if ( ilo <= intval ) then
        inc = + 1
     else
        inc = -1
     end if

     do i = ilo+inc, intval, inc

        if ( maxval <= nval ) then
           return
        end if

        nval = nval + 1
        ival(nval) = i

     end do

  end do

  return
end subroutine ranger
subroutine rat_to_s_left ( ival, jval, s )

  !*****************************************************************************80
  !
  !! RAT_TO_S_LEFT returns a left-justified representation of IVAL/JVAL.
  !
  !  Discussion:
  !
  !    If the ratio is negative, a minus sign precedes IVAL.
  !    A slash separates IVAL and JVAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 February 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IVAL, JVAL, the integers whose ratio
  !    IVAL/JVAL is to be represented.
  !
  !    If IVAL is nonzero and JVAL is 0, the string will be returned as "Inf"
  !    or "-Inf" (Infinity), and if both IVAL and JVAL are zero, the string
  !    will be returned as "NaN" (Not-a-Number).
  !
  !    Output, character ( len = * ) S, a left-justified string
  !    containing the representation of IVAL/JVAL.
  !
  implicit none

  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) ival2
  integer   ( kind = 4 ) jval
  integer   ( kind = 4 ) jval2
  character ( len = * )  s
  character ( len = 22 ) s2
  !
  !  Take care of simple cases right away.
  !
  if ( ival == 0 ) then

     if ( jval /= 0 ) then
        s2 = '0'
     else
        s2 = 'NaN'
     end if

  else if ( jval == 0 ) then

     if ( 0 < ival ) then
        s2 = 'Inf'
     else
        s2 = '-Inf'
     end if
     !
     !  Make copies of IVAL and JVAL.
     !
  else

     ival2 = ival
     jval2 = jval

     if ( jval2 == 1 ) then
        write ( s2, '(i11)' ) ival2
     else
        write ( s2, '(i11, ''/'', i10)' ) ival2, jval2
     end if

     call s_blank_delete ( s2 )

  end if

  s = s2

  return
end subroutine rat_to_s_left
subroutine rat_to_s_right ( ival, jval, s )

  !*****************************************************************************80
  !
  !! RAT_TO_S_RIGHT returns a right-justified representation of IVAL/JVAL.
  !
  !  Discussion:
  !
  !    If the ratio is negative, a minus sign precedes IVAL.
  !    A slash separates IVAL and JVAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 February 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IVAL, JVAL, the two integers whose
  !    ratio IVAL/JVAL is to be represented.
  !
  !    Note that if IVAL is nonzero and JVAL is 0, the string will
  !    be returned as "Inf" or "-Inf" (Infinity), and if both
  !    IVAL and JVAL are zero, the string will be returned as "NaN"
  !    (Not-a-Number).
  !
  !    Output, character ( len = * ) S, a right-justified string
  !    containing the representation of IVAL/JVAL.
  !
  implicit none

  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) jval
  character ( len = * )  s

  call rat_to_s_left ( ival, jval, s )
  call s_adjustr ( s )

  return
end subroutine rat_to_s_right
subroutine s_adjustl ( s )

  !*****************************************************************************80
  !
  !! S_ADJUSTL flushes a string left.
  !
  !  Discussion:
  !
  !    Both blanks and tabs are treated as "white space".
  !
  !    This routine is similar to the FORTRAN90 ADJUSTL routine.
  !
  !  Example:
  !
  !    Input             Output
  !
  !    '     Hello'      'Hello     '
  !    ' Hi there!  '    'Hi there!   '
  !    'Fred  '          'Fred  '
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.
  !    On input, S is a string of characters.
  !    On output, any initial blank or tab characters have been cut.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) nonb
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character, parameter :: TAB = achar ( 9 )
  !
  !  Check the length of the string to the last nonblank.
  !  If nonpositive, return.
  !
  s_length = len_trim ( s )

  if ( s_length <= 0 ) then
     return
  end if
  !
  !  Find NONB, the location of the first nonblank, nontab.
  !
  nonb = 0

  do i = 1, s_length

     if ( s(i:i) /= ' ' .and. s(i:i) /= TAB ) then
        nonb = i
        exit
     end if

  end do

  if ( nonb == 0 ) then
     s = ' '
     return
  end if
  !
  !  Shift the string left.
  !
  if ( 1 < nonb ) then
     do i = 1, s_length + 1 - nonb
        s(i:i) = s(i+nonb-1:i+nonb-1)
     end do
  end if
  !
  !  Blank out the end of the string.
  !
  s(s_length+2-nonb:s_length) = ' '

  return
end subroutine s_adjustl
subroutine s_adjustr ( s )

  !*****************************************************************************80
  !
  !! S_ADJUSTR flushes a string right.
  !
  !  Example:
  !
  !    Input             Output
  !    'Hello     '      '     Hello'
  !    ' Hi there!  '    '   Hi there!'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, on output, trailing blank
  !    characters have been cut, and pasted back onto the front.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) nonb
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  !
  !  Check the full length of the string.
  !
  s_length = len ( s )
  !
  !  Find the occurrence of the last nonblank.
  !
  nonb = len_trim ( s )
  !
  !  Shift the string right.
  !
  do i = s_length, s_length + 1 - nonb, -1
     s(i:i) = s(i-s_length+nonb:i-s_length+nonb)
  end do
  !
  !  Blank out the beginning of the string.
  !
  s(1:s_length-nonb) = ' '

  return
end subroutine s_adjustr
subroutine s_after_ss_copy ( s1, ss, s2 )

  !*****************************************************************************80
  !
  !! S_AFTER_SS_COPY copies a string after a given substring.
  !
  !  Discussion:
  !
  !    S1 and S2 can be the same object, in which case the string is
  !    overwritten by a copy of itself after the substring.
  !
  !  Example:
  !
  !    Input:
  !
  !      S1 = 'ABCDEFGH'
  !      SS = 'EF'
  !
  !    Output:
  !
  !      S2 = 'GH'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be copied.
  !
  !    Input, character ( len = * ) SS, the substring after which the copy begins.
  !
  !    Output, character ( len = * ) S2, the copied portion of S.
  !
  implicit none

  integer   ( kind = 4 ) first
  integer   ( kind = 4 ) last
  integer   ( kind = 4 ) last_s2
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  character ( len = * )  ss
  !
  !  Find the first occurrence of the substring.
  !
  first = index ( s1, ss )
  !
  !  If the substring doesn't occur at all, then S2 is blank.
  !
  if ( first == 0 ) then
     s2 = ' '
     return
  end if
  !
  !  Redefine FIRST to point to the first character to copy after
  !  the substring.
  !
  first = first + len ( ss )
  !
  !  Measure the two strings.
  !
  s1_length = len ( s1 )
  last_s2 = len ( s2 )
  !
  !  Adjust effective length of S if S2 is short.
  !
  last = min ( s1_length, last_s2 + first - 1 )
  !
  !  Copy the string.
  !
  s2(1:s1_length+1-first) = s1(first:s1_length)
  !
  !  Clear out the rest of the copy.
  !
  s2(s1_length+2-first:last_s2) = ' '

  return
end subroutine s_after_ss_copy
subroutine s_alpha_last ( s, iloc )

  !*****************************************************************************80
  !
  !! S_ALPHA_LAST returns the location of the last alphabetic character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 May 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Output, integer ( kind = 4 ) ILOC, the location of the last alphabetic
  !    character in the string.  If there are no alphabetic
  !    characters, ILOC is returned as 0.
  !
  implicit none

  logical                ch_is_alpha
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iloc
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = s_length, 1, -1
     if ( ch_is_alpha ( s(i:i) ) ) then
        iloc = i
        return
     end if
  end do

  iloc = 0

  return
end subroutine s_alpha_last
function s_any_alpha ( s )

  !*****************************************************************************80
  !
  !! S_ANY_ALPHA is TRUE if a string contains any alphabetic character.
  !
  !  Example:
  !
  !    Input         Output
  !
  !    Riding Hood   TRUE
  !    123 + 34      FALSE
  !    Seven Eleven  TRUE
  !    1.0E+11       TRUE
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string to be checked.
  !
  !    Output, logical S_ANY_ALPHA is TRUE if any character in string
  !    is an alphabetic character.
  !
  implicit none

  logical                ch_is_alpha
  integer   ( kind = 4 ) i
  character ( len = * )  s
  logical                s_any_alpha
  integer   ( kind = 4 ) s_length

  s_any_alpha = .true.
  s_length = len_trim ( s )

  do i = 1, s_length
     if ( ch_is_alpha ( s(i:i) ) ) then
        return
     end if
  end do

  s_any_alpha = .false.

  return
end function s_any_alpha
function s_any_control ( s )

  !*****************************************************************************80
  !
  !! S_ANY_CONTROL is TRUE if a string contains any control characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, is the string to check.
  !
  !    Output, logical S_ANY_CONTROL, is TRUE if any character is a control
  !    character.
  !
  implicit none

  logical                ch_is_control
  integer   ( kind = 4 ) i
  character ( len = * )  s
  logical                s_any_control
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
     if ( ch_is_control ( s(i:i) ) ) then
        s_any_control = .true.
        return
     end if
  end do

  s_any_control = .false.

  return
end function s_any_control
subroutine s_b2u ( s )

  !*****************************************************************************80
  !
  !! S_B2U replaces interword blanks by underscores.
  !
  !  Discussion:
  !
  !    Initial blanks are deleted by shifting the string to be
  !    flush left.
  !
  !    This routine is useful for making a multiword name look
  !    like a single blank-delimited string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 December 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be
  !    transformed.
  !
  implicit none

  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s = adjustl ( s )
  s_length = len_trim ( s )

  do i = 1, s_length
     if ( s(i:i) == ' ' ) then
        s(i:i) = '_'
     end if
  end do

  return
end subroutine s_b2u
subroutine s_before_ss_copy ( s, ss, s2 )

  !*****************************************************************************80
  !
  !! S_BEFORE_SS_COPY copies a string up to a given substring.
  !
  !  Discussion:
  !
  !    S and S2 can be the same object, in which case the string is
  !    overwritten by a copy of itself up to the substring, followed
  !    by blanks.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = 'ABCDEFGH'
  !      SS = 'EF'
  !
  !    Output:
  !
  !      S2 = 'ABCD'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be copied.
  !
  !    Input, character ( len = * ) SS, the substring before which the copy stops.
  !
  !    Output, character ( len = * ) S2, the copied portion of S.
  !
  implicit none

  integer   ( kind = 4 ) last
  integer   ( kind = 4 ) last_s2
  character ( len = * )  s
  character ( len = * )  s2
  character ( len = * )  ss
  !
  !  Find the first occurrence of the substring.
  !
  last = index ( s, ss )
  !
  !  If the substring doesn't occur at all, behave as though it begins
  !  just after the string terminates.
  !
  !  Now redefine LAST to point to the last character to copy before
  !  the substring begins.
  !
  if ( last == 0 ) then
     last = len ( s )
  else
     last = last - 1
  end if
  !
  !  Now adjust again in case the copy holder is "short".
  !
  last_s2 = len ( s2 )

  last = min ( last, last_s2 )
  !
  !  Copy the beginning of the string.
  !  Presumably, compilers now understand that if LAST is 0, we don't
  !  copy anything.
  !  Clear out the rest of the copy.
  !
  s2(1:last) = s(1:last)
  s2(last+1:last_s2) = ' '

  return
end subroutine s_before_ss_copy
function s_begin ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_BEGIN is TRUE if one string matches the beginning of the other.
  !
  !  Discussion:
  !
  !    The strings are compared, ignoring blanks, spaces and capitalization.
  !
  !  Example:
  !
  !     S1              S2      S_BEGIN
  !
  !    'Bob'          'BOB'     TRUE
  !    '  B  o b '    ' bo b'   TRUE
  !    'Bob'          'Bobby'   TRUE
  !    'Bobo'         'Bobb'    FALSE
  !    ' '            'Bob'     FALSE    (Do not allow a blank to match
  !                                       anything but another blank string.)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 January 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to be compared.
  !
  !    Output, logical S_BEGIN, is TRUE if the strings match up to
  !    the end of the shorter string, ignoring case.
  !
  implicit none

  logical                ch_eqi
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) i2
  logical                s_begin
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  integer   ( kind = 4 ) s2_length

  s1_length = len_trim ( s1 )
  s2_length = len_trim ( s2 )
  !
  !  If either string is blank, then both must be blank to match.
  !  Otherwise, a blank string matches anything, which is not
  !  what most people want.
  !
  if ( s1_length == 0 .or. s2_length == 0 ) then

     if ( s1_length == 0 .and. s2_length == 0 ) then
        s_begin = .true.
     else
        s_begin = .false.
     end if

     return

  end if

  i1 = 0
  i2 = 0
  !
  !  Find the next nonblank in S1.
  !
  do

     do

        i1 = i1 + 1

        if ( s1_length < i1 ) then
           s_begin = .true.
           return
        end if

        if ( s1(i1:i1) /= ' ' ) then
           exit
        end if

     end do
     !
     !  Find the next nonblank in S2.
     !
     do

        i2 = i2 + 1

        if ( s2_length < i2 ) then
           s_begin = .true.
           return
        end if

        if ( s2(i2:i2) /= ' ' ) then
           exit
        end if

     end do
     !
     !  If the characters match, get the next pair.
     !
     if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
        exit
     end if

  end do

  s_begin = .false.

  return
end function s_begin
subroutine s_behead_substring ( s, sub )

  !*****************************************************************************80
  !
  !! S_BEHEAD_SUBSTRING "beheads" a string, removing a given substring.
  !
  !  Discussion:
  !
  !    Initial blanks in the string are removed first.
  !
  !    Then, if the initial part of the string matches the substring,
  !    that part is removed and the remainder shifted left.
  !
  !    Initial blanks in the substring are NOT ignored.
  !
  !    Capitalization is ignored.
  !
  !    If the substring is equal to the string, then the resultant
  !    string is returned as a single blank.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 January 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, character ( len = * ) SUB, the substring to be removed from
  !    the beginning of the string.
  !
  implicit none

  character ( len = * )  s
  logical                s_eqi
  integer   ( kind = 4 ) s_length
  character ( len = * )  sub
  integer   ( kind = 4 ) sub_length
  !
  !  Remove leading blanks from the string.
  !
  s = adjustl ( s )
  !
  !  Get lengths.
  !
  s_length = len_trim ( s )
  sub_length = len_trim ( sub )

  if ( s_length < sub_length ) then
     return
  end if
  !
  !  If the string begins with the substring, chop it off.
  !
  if ( s_eqi ( s(1:sub_length), sub(1:sub_length) ) ) then

     if ( sub_length < s_length ) then
        s = s(sub_length+1:s_length)
        s = adjustl ( s )
     else
        s = ' '
     end if

  end if

  return
end subroutine s_behead_substring
subroutine s_blank_delete ( s )

  !*****************************************************************************80
  !
  !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
  !
  !  Discussion:
  !
  !    All TAB characters are also removed.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character, parameter :: tab = achar ( 9 )

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

     ch = s(get:get)

     if ( ch /= ' ' .and. ch /= tab ) then
        put = put + 1
        s(put:put) = ch
     end if

  end do

  s(put+1:s_length) = ' '

  return
end subroutine s_blank_delete
subroutine s_blanks_delete ( s )

  !*****************************************************************************80
  !
  !! S_BLANKS_DELETE replaces consecutive blanks by one blank.
  !
  !  Discussion:
  !
  !    Thanks to Bill Richmond for pointing out a programming flaw which
  !    meant that, as characters were slid to the left through multiple
  !    blanks, their original images were not blanked out.  This problem
  !    is easiest resolved by using a copy of the string.
  !
  !    The remaining characters are left justified and right padded with blanks.
  !    TAB characters are converted to spaces.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  character              newchr
  character              oldchr
  character ( len = * )  s
  character ( len = len ( s ) ) s_copy
  integer   ( kind = 4 ) s_length
  character, parameter :: TAB = achar ( 9 )

  s_length = len ( s )

  j = 0
  s_copy(1:s_length) = s(1:s_length)
  s(1:s_length) = ' '

  newchr = ' '

  do i = 1, s_length

     oldchr = newchr
     newchr = s_copy(i:i)

     if ( newchr == TAB ) then
        newchr = ' '
     end if

     if ( oldchr /= ' ' .or. newchr /= ' ' ) then
        j = j + 1
        s(j:j) = newchr
     end if

  end do

  return
end subroutine s_blanks_delete
subroutine s_blanks_insert ( s, ilo, ihi )

  !*****************************************************************************80
  !
  !! S_BLANKS_INSERT inserts blanks into a string, sliding old characters over.
  !
  !  Discussion:
  !
  !    Characters at the end of the string "drop off" and are lost.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, integer ( kind = 4 ) ILO, the location where the first blank
  !    is to be inserted.
  !
  !    Input, integer ( kind = 4 ) IHI, the location where the last blank
  !    is to be inserted.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) imax
  integer   ( kind = 4 ) imin
  integer   ( kind = 4 ) put
  integer   ( kind = 4 ) nmove
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  if ( ihi < ilo .or. s_length < ilo ) then
     return
  end if

  if ( ihi <= s_length ) then
     imax = ihi
  else
     imax = s_length
  end if

  if ( 1 <= ilo ) then
     imin = ilo
  else
     imin = 1
  end if

  nmove = s_length - imax

  do i = 1, nmove
     put = s_length + 1 - i
     get = s_length - imax + imin - i
     ch = s(get:get)
     s(put:put) = ch
  end do

  do i = imin, imax
     s(i:i) = ' '
  end do

  return
end subroutine s_blanks_insert

subroutine s_cap ( s )

  !*****************************************************************************80
  !
  !! S_CAP replaces any lowercase letters by uppercase ones in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length

     ch = s(i:i)
     call ch_cap ( ch )
     s(i:i) = ch

  end do

  return
end subroutine s_cap
subroutine s_cat ( s1, s2, s3 )

  !*****************************************************************************80
  !
  !! S_CAT concatenates two strings to make a third string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 September 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the "prefix" string.
  !
  !    Input, character ( len = * ) S2, the "postfix" string.
  !
  !    Output, character ( len = * ) S3, the string made by
  !    concatenating S1 and S2, ignoring any trailing blanks.
  !
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
     s3 = ' '
  else if ( s1 == ' ' ) then
     s3 = s2
  else if ( s2 == ' ' ) then
     s3 = s1
  else
     s3 = trim ( s1 ) // trim ( s2 )
  end if

  return
end subroutine s_cat
subroutine s_cat1 ( s1, s2, s3 )

  !*****************************************************************************80
  !
  !! S_CAT1 concatenates two strings, with a single blank separator.
  !
  !  Example:
  !
  !    S1       S2       S
  !
  !    'cat'    'dog'    'cat dog'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 September 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the "prefix" string.
  !
  !    Input, character ( len = * ) S2, the "postfix" string.
  !
  !    Output, character ( len = * ) S3, the string made by concatenating
  !    S1, a blank, and S2, ignoring any trailing blanks.
  !
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
     s3 = ' '
  else if ( s1 == ' ' ) then
     s3 = s2
  else if ( s2 == ' ' ) then
     s3 = s1
  else
     s3 = trim ( s1 ) // ' ' // trim ( s2 )
  end if

  return
end subroutine s_cat1
subroutine s_center ( s )

  !*****************************************************************************80
  !
  !! S_CENTER centers the non-blank portion of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 October 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.  On input, a string to be
  !    centered.  On output, the centered string.
  !
  implicit none

  integer   ( kind = 4 ) l1
  integer   ( kind = 4 ) l2
  integer   ( kind = 4 ) n1
  integer   ( kind = 4 ) n2
  integer   ( kind = 4 ) n3
  character ( len = * )  s
  !
  !  How much space is in the string?
  !
  n1 = len ( s )
  !
  !  Shift the string flush left and find the last nonblank.
  !
  s = adjustl ( s )
  n2 = len_trim ( s )

  if ( n2 <= 0 ) then
     return
  end if

  if ( n2 == n1 .or. n2 == n1 - 1 ) then
     return
  end if

  n3 = n1 - n2
  l1 = n3 / 2
  l2 = l1 + n2 + 1

  s(l1+1:l2-1) = s(1:n2)

  s(1:l1) = ' '
  s(l2:n1) = ' '

  return
end subroutine s_center
subroutine s_center_insert ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_CENTER_INSERT inserts one string into the center of another.
  !
  !  Discussion:
  !
  !    The receiving string is not blanked out first.  Therefore, if there is
  !    already information in it, some of it may still be around
  !    after the insertion.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, a string to be inserted into S2.
  !
  !    Output, character ( len = * ) S2, the string to receive S1.
  !
  implicit none

  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  integer   ( kind = 4 ) m
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  integer   ( kind = 4 ) s2_length

  s1_length = len_trim ( s1 )
  s2_length = len ( s2 )

  if ( s1_length < s2_length ) then

     m = s2_length - s1_length
     ilo = 1
     ihi = s1_length
     jlo = ( m / 2 ) + 1
     jhi = jlo + s1_length - 1

  else if ( s2_length < s1_length ) then

     m = s1_length - s2_length
     ilo = ( m / 2 ) + 1
     ihi = ilo + s2_length - 1
     jlo = 1
     jhi = s2_length

  else

     ilo = 1
     ihi = s1_length
     jlo = 1
     jhi = s2_length

  end if

  s2(jlo:jhi) = s1(ilo:ihi)

  return
end subroutine s_center_insert
subroutine s_ch_blank ( s, ch )

  !*****************************************************************************80
  !
  !! S_CH_BLANK replaces each occurrence of a particular character by a blank.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, character CH, the character to be removed.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length

     if ( s(i:i) == ch ) then
        s(i:i) = ' '
     end if

  end do

  return
end subroutine s_ch_blank
subroutine s_ch_count ( s, ch, ch_count )

  !*****************************************************************************80
  !
  !! S_CH_COUNT counts occurrences of a particular character in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    13 January 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string.
  !
  !    Input, character CH, the character to be counted.
  !
  !    Output, integer ( kind = 4 ) CH_COUNT, the number of occurrences.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) ch_count
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  ch_count = 0

  s_length = len ( s )

  do i = 1, s_length

     if ( s(i:i) == ch ) then
        ch_count = ch_count + 1
     end if

  end do

  return
end subroutine s_ch_count
subroutine s_ch_delete ( s, ch )

  !*****************************************************************************80
  !
  !! S_CH_DELETE removes all occurrences of a character from a string.
  !
  !  Discussion:
  !
  !    Each time the given character is found in the string, the characters
  !    to the right of the string are shifted over one position.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, character CH, the character to be removed.
  !
  implicit none

  character              ch
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  put = 1

  do get = 1, s_length

     if ( s(get:get) == ch ) then

     else if ( put == get ) then
        put = put + 1
     else
        s(put:put) = s(get:get)
        put = put + 1
     end if

  end do

  s(put:s_length) = ' '

  return
end subroutine s_ch_delete
function s_ch_last ( s )

  !*****************************************************************************80
  !
  !! S_CH_LAST returns the last nonblank character in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, character S_CH_LAST, the last nonblank character in S,
  !    or ' ' if S is all blank.
  !
  implicit none

  character ( len = * )  s
  character              s_ch_last
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  if ( 0 < s_length ) then
     s_ch_last = s(s_length:s_length)
  else
     s_ch_last = ' '
  end if

  return
end function s_ch_last
subroutine s_chop ( s, ilo, ihi )

  !*****************************************************************************80
  !
  !! S_CHOP "chops out" a portion of a string, and closes up the hole.
  !
  !  Example:
  !
  !    S = 'Fred is not a jerk!'
  !
  !    call s_chop ( S, 9, 12 )
  !
  !    S = 'Fred is a jerk!    '
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
  !    characters to be removed.
  !
  implicit none

  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ihi2
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) ilo2
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  ilo2 = max ( ilo, 1 )
  ihi2 = min ( ihi, s_length )

  if ( ihi2 < ilo2 ) then
     return
  end if

  s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
  s(s_length+ilo2-ihi2:s_length) = ' '

  return
end subroutine s_chop

subroutine s_compare ( s1, s2, order )
  !*****************************************************************************80
  !
  !! S_COMPARE compares two strings.
  !
  !  Discussion:
  !
  !    The FORTRAN function LLT ( S1, S2 ) returns TRUE if S1 is lexically
  !    strictly less than S2, and FALSE otherwise.
  !
  !    There are related functions LLE, LGT, LGE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 December 2010
  !
  !  Author:
  !
  !   John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, integer ( kind = 4 ) ORDER:
  !    -1, S1 < S2.
  !     0, S1 = S2
  !    +1, S1 > S2
  !
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  integer ( kind = 4 ) s_len_trim
  character ( len = * ) s1
  integer ( kind = 4 ) s1_len
  character ( len = * ) s2
  integer ( kind = 4 ) s2_len

  s1_len = len_trim ( s1 )
  s2_len = len_trim ( s2 )

  order = 0

  do i = 1, min ( s1_len, s2_len )
     if ( s1(i:i) < s2(i:i) ) then
        order = -1
        return
     else if ( s2(i:i) < s1(i:i) ) then
        order = +1
        return
     end if

  end do
  !
  !  If one string is actually longer than the other, and nonblank,
  !  it must come after the other.
  !
  if ( s1_len < s2_len ) then
     order = -1
     return
  else if ( s2_len < s1_len ) then
     order = +1
     return
  end if

  return
end subroutine s_compare
subroutine s_control_blank ( s )

  !*****************************************************************************80
  !
  !! S_CONTROL_BLANK replaces control characters with blanks.
  !
  !  Discussion:
  !
  !    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  implicit none

  logical                ch_is_control
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
     if ( ch_is_control ( s(i:i) ) ) then
        s(i:i) = ' '
     end if
  end do

  return
end subroutine s_control_blank
subroutine s_control_count ( s, ifound )

  !*****************************************************************************80
  !
  !! S_CONTROL_COUNT returns the number of control characters in a string.
  !
  !  Discussion:
  !
  !    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Output, integer ( kind = 4 ) IFOUND, the number of control characters.
  !
  implicit none

  logical                ch_is_control
  integer   ( kind = 4 ) ifound
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  ifound = 0
  s_length = len_trim ( s )

  do i = 1, s_length

     if ( ch_is_control ( s(i:i) ) ) then
        ifound = ifound + 1
     end if

  end do

  return
end subroutine s_control_count
subroutine s_control_delete ( s )

  !*****************************************************************************80
  !
  !! S_CONTROL_DELETE removes all control characters from a string.
  !
  !  Discussion:
  !
  !    The string is collapsed to the left, and padded on the right with
  !    blanks to replace the removed characters.
  !
  !    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, is the string to be transformed.
  !
  implicit none

  logical                ch_is_control
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

     if ( .not. ch_is_control ( s(get:get) ) ) then
        put = put + 1
        s(put:put) = s(get:get)
     end if

  end do
  !
  !  Pad the end of the string with blanks
  !
  s(put+1:) = ' '

  return
end subroutine s_control_delete
subroutine s_copy ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_COPY copies one string into another.
  !
  !  Discussion:
  !
  !    If S1 is shorter than S2, the rest of S2 is blank.
  !    If S1 is longer than S2, then the excess information is lost.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 January 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be copied.
  !
  !    Output, character ( len = * ) S2, the copy.
  !
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2

  s2(1:min(len(s1),len(s2))) = s1(1:min(len(s1),len(s2)))
  s2(len(s1)+1:len(s2)) = ' '

  return
end subroutine s_copy
subroutine s_detag ( s )

  !*****************************************************************************80
  !
  !! S_DETAG removes from a string all substrings marked by angle brackets.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = '<I>This is Italic</I> whereas this is <B>boldly</B> not!'
  !
  !    Output:
  !
  !      S = ' whereas this is  not!'
  !
  !  Discussion:
  !
  !    This routine was written to help extract some data that was hidden
  !    inside an elaborate HTML table.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 September 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  implicit none

  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i3
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  do

     s_length = len_trim ( s )

     if ( len_trim ( s ) == 0 ) then
        exit
     end if

     i1 = index ( s, '<' )

     if ( i1 <= 0 .or. s_length <= i1 ) then
        exit
     end if

     i2 = index ( s(i1+1:), '>' )

     if ( i2 == 0 ) then
        exit
     end if

     i2 = i2 + i1
     !
     !  Shift.
     !
     s(i1:s_length+i1-i2-1) = s(i2+1:s_length)
     !
     !  Pad.
     !
     s(s_length+i1-i2:) = ' '

  end do

  return
end subroutine s_detag
subroutine s_detroff ( s )

  !*****************************************************************************80
  !
  !! S_DETROFF removes obnoxious "character" + backspace pairs from a string.
  !
  !  Discussion:
  !
  !    Given the string of characters:
  !      'AB#C#D#E'
  !    where we are using "#" to represent a backspace, the returned string
  !    will be
  !      'AE'.
  !
  !    This function was written for use in "cleaning up" UNICOS MAN pages.
  !    These MAN pages were formatted in the Byzantine TROFF printing format.
  !    Although the files were text, and would seem to "print" correctly to
  !    the screen, an unholy mess would emerge if the same file was sent
  !    to the printer.  This is because the screen handled the backspaces
  !    by backspacing, but most printers don't know anymore how to handle
  !    TROFF's backspaces, and so they just print them as blobs, instead of,
  !    say, spacing back.
  !
  !    In particular:
  !
  !      Passages which are to be underlined are written so:
  !      "_#T_#e_#x_#t" when what is meant is that "Text" is to be
  !      underlined if possible.  Note that the seemingly equivalent
  !      "T#_e#_x#_t#_" is NOT used.  This is because, in the olden
  !      days, certain screen terminals could backspace, but would only
  !      display the new character, obliterating rather than
  !      overwriting the old one.  This convention allows us to know
  !      that we want to delete "character" + Backspace, rather than
  !      Backspace + "character".
  !
  !      Passages which are meant to be in BOLDFACE are written so:
  !      "U#U#U#Ug#g#g#gl#l#l#ly#y#y#y", when what is meant is that
  !      "Ugly" is to be printed as boldly as possible.  These boldface
  !      passages may also be cleaned up using the same rule of
  !      removing all occurrences of "character" + Backspace.
  !
  !    It is truly a fright to look at the text of one of these MAN
  !    pages with all the ugly Backspace's, which display on VMS as ^H.
  !    These files print or type properly, but look awful in an editor.
  !    Moreoever, the lavish use of boldface means that text that is
  !    meant to fit in 80 columns can sometimes require 7 times as much
  !    space to describe.  This can cause a VMS editor to abort, or to
  !    skip the line, since 255 characters is the maximum for EDT.
  !
  !    A FORTRAN program that tries to read a long line like that will
  !    also fail if not careful, since a formatted sequential file
  !    on VMS has a default maximum record length of something like
  !    133 characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 January 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the line of text to
  !    be de-TROFF'ed.
  !
  implicit none

  character, parameter :: BS = achar ( 8 )
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )
  i = 1

  do while ( i <= s_length )

     if ( s(i:i) == BS ) then

        if ( i == 1 ) then
           s(1:s_length-1) = s(2:s_length)
           s(s_length:s_length) = ' '
           s_length = s_length - 1
           i = i - 1
        else
           s(i-1:s_length-2) = s(i+1:s_length)
           s(s_length-1:s_length) = ' '
           s_length = s_length - 2
           i = i - 2
        end if

     end if

     i = i + 1

  end do

  return
end subroutine s_detroff
function s_eqi ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_EQI is a case insensitive comparison of two strings for equality.
  !
  !  Discussion:
  !
  !    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_EQI, the result of the comparison.
  !
  implicit none

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  logical                s_eqi
  character ( len = *  ) s1
  integer   ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer   ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

  s_eqi = .false.

  do i = 1, lenc

     c1 = s1(i:i)
     c2 = s2(i:i)
     call ch_cap ( c1 )
     call ch_cap ( c2 )

     if ( c1 /= c2 ) then
        return
     end if

  end do

  do i = lenc + 1, s1_length
     if ( s1(i:i) /= ' ' ) then
        return
     end if
  end do

  do i = lenc + 1, s2_length
     if ( s2(i:i) /= ' ' ) then
        return
     end if
  end do

  s_eqi = .true.

  return
end function s_eqi
function s_eqidb ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_EQIDB compares two strings, ignoring case and blanks.
  !
  !  Example:
  !
  !    S_EQIDB ( 'Nor Way', 'NORway' ) is TRUE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Modified:
  !
  !    19 July 1998
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_EQIDB, the result of the comparison.
  !
  implicit none

  character              c1
  character              c2
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) len2
  logical                s_eqidb
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  !
  !  Get the length of each string to the last nonblank.
  !
  s1_length = len_trim ( s1 )
  len2 = len_trim ( s2 )
  !
  !  Assume we're going to fail.
  !
  s_eqidb = .false.
  !
  !  Initialize the pointers to characters in each string.
  !
  i1 = 0
  i2 = 0

  do
     !
     !  If we've matched all the nonblank characters in both strings,
     !  then return with S_EQIDB = TRUE.
     !
     if ( i1 == s1_length .and. i2 == len2 ) then
        s_eqidb = .true.
        return
     end if
     !
     !  Get the next nonblank character in the first string.
     !
     do

        i1 = i1 + 1

        if ( s1_length < i1 ) then
           return
        end if

        if ( s1(i1:i1) /= ' ' ) then
           exit
        end if

     end do

     c1 = s1(i1:i1)
     call ch_cap ( c1 )
     !
     !  Get the next nonblank character in the second string.
     !
     do

        i2 = i2 + 1
        if ( len2 < i2 ) then
           return
        end if

        c2 = s2(i2:i2)

        if ( c2 /= ' ' ) then
           exit
        end if

     end do

     call ch_cap ( c2 )

     if ( c1 /= c2 ) then
        exit
     end if

  end do

  return
end function s_eqidb

subroutine s_escape_tex ( s1, s2 )
 !*****************************************************************************80
  !
  !! S_ESCAPE_TEX de-escapes TeX escape sequences.
  !
  !  Discussion:
  !
  !    In particular, every occurrence of the characters '\', '_',
  !    '^', '{' and '}' will be replaced by '\\', '\_', '\^',
  !    '\{' and '\}'.  A TeX interpreter, on seeing these character
  !    strings, is then likely to return the original characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 January 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be de-escaped.
  !
  !    Output, character ( len = * ) S2, a copy of the string,
  !    modified to avoid TeX escapes.
  !
  implicit none

  character              ch
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  integer   ( kind = 4 ) s1_pos
  character ( len = * )  s2
  integer   ( kind = 4 ) s2_pos

  s1_length = len_trim ( s1 )

  s1_pos = 0
  s2_pos = 0
  s2 = ' '

  do while ( s1_pos < s1_length )
     s1_pos = s1_pos + 1
     ch = s1(s1_pos:s1_pos)
     if(ch == '\' .or. ch == '_' .or. &
          ch == '^' .or. &
          ch == '{' .or. &
          ch == '}' ) then
     s2_pos = s2_pos + 1
     s2(s2_pos:s2_pos) = '\'
     end if
     s2_pos = s2_pos + 1
    s2(s2_pos:s2_pos) = ch
   end do
  return
end 

subroutine s_fill ( s, ch )

  !*****************************************************************************80
  !
  !! S_FILL overwrites every character of a string by a given character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, character ( len = * ) S, the string to be overwritten.
  !
  !    Input, character CH, the overwriting character.
  !
implicit none

character              ch
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len ( s )

do i = 1, s_length
  s(i:i) = ch
end do

return
end subroutine s_fill
function s_first_nonblank ( s )

  !*****************************************************************************80
  !
  !! S_FIRST_NONBLANK returns the location of the first nonblank.
  !
  !  Discussion:
  !
  !    If all characters are blanks, a 0 is returned.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, integer ( kind = 4 ) S_FIRST_NONBLANK, the location of the first
  !    nonblank character in the string, or 0 if all are blank.
  !
implicit none

integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_first_nonblank
integer   ( kind = 4 ) s_length

s_length = len ( s )

do i = 1, s_length

  if ( s(i:i) /= ' ' ) then
     s_first_nonblank = i
     return
  end if

end do

s_first_nonblank = 0

return
end function s_first_nonblank
function s_gei ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_GEI = ( S1 is lexically greater than or equal to S2 ).
  !
  !  Discussion:
  !
  !    The comparison is done in a case-insensitive way.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_GEI, the result of the comparison.
  !
implicit none

character              c1
character              c2
integer   ( kind = 4 ) i
integer   ( kind = 4 ) lenc
logical                s_gei
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2
integer   ( kind = 4 ) s2_length

s1_length = len_trim ( s1 )
s2_length = len_trim ( s2 )
lenc = min ( s1_length, s2_length )

do i = 1, lenc

  c1 = s1(i:i)
  c2 = s2(i:i)
  call ch_cap ( c1 )
  call ch_cap ( c2 )

  if ( lgt ( c1, c2 ) ) then
     s_gei = .true.
     return
  else if ( llt ( c1, c2 ) ) then
     s_gei = .false.
     return
  end if

end do

if ( s1_length < s2_length ) then
  s_gei = .false.
else
  s_gei = .true.
end if

return
end function s_gei
function s_gti ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_GTI = S1 is lexically greater than S2.
  !
  !  Discussion:
  !
  !    The comparison is done in a case-insensitive way.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_GTI, the result of the comparison.
  !
implicit none

character              c1
character              c2
integer   ( kind = 4 ) i
integer   ( kind = 4 ) lenc
logical                s_gti
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2
integer   ( kind = 4 ) s2_length

s1_length = len ( s1 )
s2_length = len ( s2 )
lenc = min ( s1_length, s2_length )

do i = 1, lenc

  c1 = s1(i:i)
  c2 = s2(i:i)
  call ch_cap ( c1 )
  call ch_cap ( c2 )

  if ( lgt ( c1, c2 ) ) then
     s_gti = .true.
     return
  else if ( llt ( s1, s2 ) ) then
     s_gti = .false.
     return
  end if

end do

if ( s1_length <= s2_length ) then
  s_gti = .false.
else
  s_gti = .true.
end if

return
end function s_gti
function s_index ( s, sub )

  !*****************************************************************************80
  !
  !! S_INDEX seeks the first occurrence of a substring.
  !
  !  Discussion:
  !
  !    The function returns the location in the string at which the
  !    substring SUB is first found, or 0 if the substring does not
  !    occur at all.
  !
  !    The routine is trailing blank insensitive.  This is very
  !    important for those cases where you have stored information in
  !    larger variables.  If S is of length 80, and SUB is of
  !    length 80, then if S = 'FRED' and SUB = 'RED', a match would
  !    not be reported by the standard FORTRAN INDEX, because it treats
  !    both variables as being 80 characters long!  This routine assumes that
  !    trailing blanks represent garbage!
  !
  !    Because of the suppression of trailing blanks, this routine cannot be
  !    used to find, say, the first occurrence of the two-character
  !    string 'A '.  However, this routine treats as a special case the
  !    occurrence where S or SUB is entirely blank.  Thus you can
  !    use this routine to search for occurrences of double or triple blanks
  !    in a string, for example.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character ( len = * ) SUB, the substring to search for.
  !
  !    Output, integer ( kind = 4 ) S_INDEX.  0 if SUB does not occur in
  !    the string.  Otherwise S(S_INDEX:S_INDEX+LENS-1) = SUB,
  !    where LENS is the length of SUB, and is the first place
  !    this happens.
  !
implicit none

integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_index
integer   ( kind = 4 ) s_length
character ( len = * )  sub
integer   ( kind = 4 ) sub_length

s_index = 0

s_length = len_trim ( s )
sub_length = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
if ( s_length == 0 ) then
  s_length = len ( s )
end if

if ( sub_length == 0 ) then
  sub_length = len ( sub )
end if

if ( s_length < sub_length ) then
  return
end if

do i = 1, s_length + 1 - sub_length

  if ( s(i:i+sub_length-1) == sub(1:sub_length) ) then
     s_index = i
     return
  end if

end do

return
end function s_index
function s_index_set ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_INDEX_SET searches a string for any of a set of characters.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be examined.
  !
  !    Input, character ( len = * ) S2, the characters to search for.
  !
  !    Output, integer ( kind = 4 ) S_INDEX_SET, the first location of a
  !    character from S2 in S1, or 0 if no character from S2 occurs in S1.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) j
integer   ( kind = 4 ) k
character ( len = *  ) s1
integer   ( kind = 4 ) s1_length
character ( len = *  ) s2
integer   ( kind = 4 ) s_index_set

s1_length = len ( s1 )

j = s1_length + 1

do i = 1, len ( s2 )
  k = index ( s1, s2(i:i) )
  if ( k /= 0 ) then
     j = min ( j, k )
  end if
end do

if ( j == s1_length + 1 ) then
  j = 0
end if

s_index_set = j

return
end function s_index_set
function s_indexi ( s, sub )

  !*****************************************************************************80
  !
  !! S_INDEXI is a case-insensitive INDEX function.
  !
  !  Discussion:
  !
  !    The function returns the location in the string at which the
  !    substring SUB is first found, or 0 if the substring does not
  !    occur at all.
  !
  !    The routine is also trailing blank insensitive.  This is very
  !    important for those cases where you have stored information in
  !    larger variables.  If S is of length 80, and SUB is of
  !    length 80, then if S = 'FRED' and SUB = 'RED', a match would
  !    not be reported by the standard FORTRAN INDEX, because it treats
  !    both variables as being 80 characters long!  This routine assumes that
  !    trailing blanks represent garbage!
  !
  !    Because of the suppression of trailing blanks, this routine cannot be
  !    used to find, say, the first occurrence of the two-character
  !    string 'A '.  However, this routine treats as a special case the
  !    occurrence where S or SUB is entirely blank.  Thus you can
  !    use this routine to search for occurrences of double or triple blanks
  !    in a string, for example, although INDEX itself would be just as
  !    suitable for that problem.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character ( len = * ) SUB, the substring to search for.
  !
  !    Output, integer ( kind = 4 ) S_INDEXI.  0 if SUB does not occur in
  !    the string.  Otherwise S(S_INDEXI:S_INDEXI+LENS-1) = SUB,
  !    where LENS is the length of SUB, and is the first place
  !    this happens.  However, note that this routine ignores case,
  !    unlike the standard FORTRAN INDEX function.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) llen2
character ( len = * )  s
logical                s_eqi
integer   ( kind = 4 ) s_indexi
integer   ( kind = 4 ) s_length
character ( len = * )  sub

s_indexi = 0

s_length = len_trim ( s )
llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
if ( s_length == 0 ) then
  s_length = len ( s )
end if

if ( llen2 == 0 ) then
  llen2 = len ( sub )
end if

if ( s_length < llen2 ) then
  return
end if

do i = 1, s_length + 1 - llen2

  if ( s_eqi ( s(i:i+llen2-1), sub ) ) then
     s_indexi = i
     return
  end if

end do

return
end function s_indexi
function s_index_last ( s, sub )

  !*****************************************************************************80
  !
  !! S_INDEX_LAST finds the LAST occurrence of a given substring.
  !
  !  Discussion:
  !
  !    It returns the location in the string at which the substring SUB is
  !    first found, or 0 if the substring does not occur at all.
  !
  !    The routine is also trailing blank insensitive.  This is very
  !    important for those cases where you have stored information in
  !    larger variables.  If S is of length 80, and SUB is of
  !    length 80, then if S = 'FRED' and SUB = 'RED', a match would
  !    not be reported by the standard FORTRAN INDEX, because it treats
  !    both variables as being 80 characters long!  This routine assumes that
  !    trailing blanks represent garbage!
  !
  !    This means that this routine cannot be used to find, say, the last
  !    occurrence of a substring 'A ', since it assumes the blank space
  !    was not specified by the user, but is, rather, padding by the
  !    system.  However, as a special case, this routine can properly handle
  !    the case where either S or SUB is all blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character ( len = * ) SUB, the substring to search for.
  !
  !    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
  !    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
  !    where LENS is the length of SUB, and is the last place
  !    this happens.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) j
integer   ( kind = 4 ) llen2
character ( len = * )  s
integer   ( kind = 4 ) s_index_last
integer   ( kind = 4 ) s_length
character ( len = * )  sub

s_index_last = 0

s_length = len_trim ( s )
llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
if ( s_length == 0 ) then
  s_length = len ( s )
end if

if ( llen2 == 0 ) then
  llen2 = len ( sub )
end if

if ( s_length < llen2 ) then
  return
end if

do j = 1, s_length + 1 - llen2

  i = s_length + 2 - llen2 - j

  if ( s(i:i+llen2-1) == sub ) then
     s_index_last = i
     return
  end if

end do

return
end function s_index_last
function s_index_last_c ( s, c )

  !*****************************************************************************80
  !
  !! S_INDEX_LAST_C finds the LAST occurrence of a given character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 December 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character C, the character to search for.
  !
  !    Output, integer ( kind = 4 ) S_INDEX_LAST_C, the index in S where C occurs
  !    last, or -1 if it does not occur.
  !
implicit none

character              c
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length
integer   ( kind = 4 ) s_index_last_c

if ( c == ' ' ) then
  s_length = len ( s )
else
  s_length = len_trim ( s )
end if

do i = s_length, 1, -1

  if ( s(i:i) == c ) then
     s_index_last_c = i
     return
  end if

end do

s_index_last_c = -1

return
end function s_index_last_c
subroutine s_i_append ( s, i, done )

  !*****************************************************************************80
  !
  !! S_I_APPEND appends an integer to a string.
  !
  !  Discussion:
  !
  !    A blank space will separate the integer from the text already
  !    in the line.
  !
  !    The routine warns the user if the integer will not fit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 December 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a line of text.
  !    On input, the current string.  On output, the current string
  !    with the integer value appended.
  !
  !    Input, integer ( kind = 4 ) I, an integer to be appended to the line.
  !
  !    Output, logical DONE, is FALSE if there was not enough room
  !    to append the integer.
  !
implicit none

logical                done
integer   ( kind = 4 ) i
integer   ( kind = 4 ) lens
integer   ( kind = 4 ) lents
integer   ( kind = 4 ) lenw
integer   ( kind = 4 ) next
character ( len = * )  s
integer   ( kind = 4 ) s_length
character ( len = 13 ) w

done = .false.

s_length = len ( s )
lents = len_trim ( s )

call i4_to_s_left ( i, w )

lenw = len_trim ( w )

if ( lents == 0 ) then
  if ( s_length < lenw ) then
     done = .true.
     return
  end if
else
  if ( s_length < lents + 1 + lenw ) then
     done = .true.
     return
  end if
end if

if ( lents == 0 ) then
  next = 1
else
  next = lents + 1
  s(next:next) = ' '
  next = next + 1
end if

s(next:next+lenw-1) = w(1:lenw)

return
end subroutine s_i_append
subroutine s_inc_c ( s )

  !*****************************************************************************80
  !
  !! S_INC_C "increments" the characters in a string.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    The routine tries to produce the next string, in dictionary order,
  !    following the input value of a string.  Digits, spaces, and other
  !    nonalphabetic characters are ignored.  Case is respected; in other
  !    words, the case of every alphabetic character on input will be the
  !    same on output.
  !
  !    The following error conditions can occur:
  !
  !      There are no alphabetic characters in the string.  No
  !      incrementing is possible.
  !
  !      All alphabetic characters are equal to 'Z' or 'z'.  In this
  !      the string is also "wrapped around" so that all alphabetic
  !      characters are "A" or "a".
  !
  !    If the word "Tax" were input, the successive outputs would be
  !    "Tay", "Taz", "Tba", "Tbb", ...  If the input word "January 4, 1989"
  !    were input, the output would be "Januarz 4, 1989".
  !
  !    This routine could be useful when trying to create a unique file
  !    name or variable name at run time.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 April 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string whose
  !    alphabetic successor is desired.  On output, S has been replaced
  !    by its "successor".
  !
implicit none

integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) iloc
character ( len = * )  s

ilo = 1
ihi = len ( s )
!
!  Find the last alphabetic character in the string.
!
do

  call s_alpha_last ( s(ilo:ihi), iloc )
  !
  !  If there is no alphabetic character, we can't help.
  !
  if ( iloc == 0 ) then
     return
  end if

  if ( s(iloc:iloc) == achar ( 122 ) ) then

     s(iloc:iloc) = achar ( 97 )
     ihi = iloc - 1

     if ( ihi <= 0 ) then
        exit
     end if

  else if ( s(iloc:iloc) == achar ( 90 ) ) then

     s(iloc:iloc) = achar ( 65 )
     ihi = iloc - 1

     if ( ihi <= 0 ) then
        return
     end if

  else

     s(iloc:iloc) = achar ( iachar ( s(iloc:iloc) ) + 1 )
     exit

  end if

end do

return
end subroutine s_inc_c
subroutine s_inc_n ( s )

  !*****************************************************************************80
  !
  !! S_INC_N increments the digits in a string.
  !
  !  Discussion:
  !
  !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, which
  !    guarantees the ASCII collating sequence.
  !
  !    It is assumed that the digits in the name, whether scattered or
  !    connected, represent a number that is to be increased by 1 on
  !    each call.  If this number is all 9's on input, the output number
  !    is all 0's.  Non-numeric letters of the name are unaffected.
  !
  !    If the name is empty, then the routine stops.
  !
  !    If the name contains no digits, the empty string is returned.
  !
  !  Example:
  !
  !      Input            Output
  !      -----            ------
  !      'a7to11.txt'     'a7to12.txt'
  !      'a7to99.txt'     'a8to00.txt'
  !      'a9to99.txt'     'a0to00.txt'
  !      'cat.txt'        ' '
  !      ' '              STOP!
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 September 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.
  !    On input, a character string to be incremented.
  !    On output, the incremented string.
  !
implicit none

character              c
integer   ( kind = 4 ) change
integer   ( kind = 4 ) digit
integer   ( kind = 4 ) i
integer   ( kind = 4 ) lens
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

if ( s_length <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_INC_N - Fatal error!'
  write ( *, '(a)' ) '  The input string is empty.'
  stop
end if

change = 0

do i = s_length, 1, -1

  c = s(i:i)

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

     change = change + 1

     digit = iachar ( c ) - 48
     digit = digit + 1

     if ( digit == 10 ) then
        digit = 0
     end if

     c = achar ( digit + 48 )

     s(i:i) = c

     if ( c /= '0' ) then
        return
     end if

  end if

end do

if ( change == 0 ) then
  s = ' '
  return
end if

return
end subroutine s_inc_n
subroutine s_input ( string, value, ierror )

  !*****************************************************************************80
  !
  !! S_INPUT prints a prompt string and reads a string from the user.
  !
  !  Discussion:
  !
  !    If the input line starts with a comment character ('#'), or is blank,
  !    the routine ignores that line, and tries to read the next one.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    13 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) STRING, the prompt string.
  !
  !    Output, character ( len = * ) VALUE, the value input by the user.
  !
  !    Output, integer ( kind = 4 ) IERROR, an error flag, which is 0
  !    if no error occurred.
  !
implicit none

integer   ( kind = 4 ) ierror
character ( len = * )  string
character ( len = * )  value

ierror = 0
value = ' '
!
!  Write the prompt.
!
write ( *, '(a)' ) ' '
write ( *, '(a)' ) trim ( string )

do

  read ( *, '(a)', iostat = ierror ) value

  if ( ierror /= 0 ) then
     value = 'S_INPUT: Input error!'
     return
  end if
  !
  !  If the line begins with a comment character, go back and read the next line.
  !
  if ( value(1:1) == '#' ) then
     cycle
  end if

  if ( len_trim ( value ) == 0 ) then
     cycle
  end if

  exit

end do

return
end subroutine s_input
function s_is_alpha ( s )

  !*****************************************************************************80
  !
  !! S_IS_ALPHA returns TRUE if the string contains only alphabetic characters.
  !
  !  Discussion:
  !
  !    Here, alphabetic characters are 'A' through 'Z' and 'a' through 'z'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, logical S_IS_ALPHA, is TRUE if the string contains only
  !    alphabetic characters.
  !
implicit none

logical                ch_is_alpha
integer   ( kind = 4 ) i
character ( len = * )  s
logical                s_is_alpha
integer   ( kind = 4 ) s_length

s_is_alpha = .false.
s_length = len_trim ( s )

do i = 1, s_length

  if ( .not. ch_is_alpha ( s(i:i) ) ) then
     return
  end if

end do

s_is_alpha = .true.

return
end function s_is_alpha
function s_is_alphanumeric ( s )

  !*****************************************************************************80
  !
  !! S_IS_ALPHANUMERIC = string contains only alphanumeric characters.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    Alphanumeric characters are 'A' through 'Z', 'a' through 'z' and
  !    '0' through '9'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, logical S_IS_ALPHANUMERIC, is TRUE if the string contains only
  !    alphabetic characters and numerals.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) itemp
character ( len = * )  s
logical                s_is_alphanumeric
integer   ( kind = 4 ) s_length

s_is_alphanumeric = .false.
s_length = len_trim ( s )

do i = 1, s_length

  itemp = iachar ( s(i:i) )

  if ( .not. ( 65 <= itemp .and. itemp <= 90 ) ) then
     if ( .not. ( 97 <= itemp .and. itemp <= 122 ) ) then
        if ( .not. ( 48 <= itemp .and. itemp <= 57 ) ) then
           return
        end if
     end if
  end if

end do

s_is_alphanumeric = .true.

return
end function s_is_alphanumeric
function s_is_digit ( s )

  !*****************************************************************************80
  !
  !! S_IS_DIGIT returns TRUE if a string contains only decimal digits.
  !
  !  Discussion:
  !
  !    This is a strict comparison.
  !    The check is made from the first character to the last nonblank.
  !    Each character in between must be one of '0', '1', ..., '9'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, logical S_IS_DIGIT, is TRUE if S contains only digits.
  !
implicit none

character              c
integer   ( kind = 4 ) i
integer   ( kind = 4 ) lenc
character ( len = * )  s
logical                s_is_digit
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

s_is_digit = .false.

do i = 1, s_length
  c = s(i:i)
  if ( llt ( c, '0' ) .or. lgt ( c, '9' ) ) then
     return
  end if
end do

s_is_digit = .true.

return
end function s_is_digit
function s_is_f77_name ( s )

  !*****************************************************************************80
  !
  !! S_IS_F77_NAME = input string represent a legal FORTRAN77 identifier.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, logical S_IS_F77_NAME, is TRUE if the string is a legal FORTRAN77
  !    identifier.  That is, the string must begin with an alphabetic
  !    character, and all subsequent characters must be alphanumeric.
  !    The string may terminate with blanks.  No underscores are allowed.
  !
implicit none

logical                ch_is_alpha
integer   ( kind = 4 ) lenc
character ( len = * )  s
logical                s_is_alphanumeric
logical                s_is_f77_name
integer   ( kind = 4 ) s_length

s_is_f77_name = .false.

s_length = len_trim ( s )

if ( s_length <= 0 ) then
  return
end if

if ( .not. ch_is_alpha ( s(1:1) ) ) then
  return
end if

if ( .not. s_is_alphanumeric ( s(2:s_length) ) ) then
  return
end if

s_is_f77_name = .true.

return
end function s_is_f77_name
function s_is_f90_name ( s )

  !*****************************************************************************80
  !
  !! S_IS_F90_NAME = input string represent a legal FORTRAN90 identifier.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, logical S_IS_F90_NAME, is TRUE if the string is a legal
  !    FORTRAN90 identifier.  That is, the string must begin with an alphabetic
  !    character, and all subsequent characters must be alphanumeric
  !    or underscores.  The string may terminate with blanks.
  !
implicit none

logical                ch_is_alpha
integer   ( kind = 4 ) i
logical                malphnum2
character ( len = * )  s
logical                s_is_f90_name
integer   ( kind = 4 ) s_length

s_is_f90_name = .false.

s_length = len_trim ( s )

if ( s_length <= 0 ) then
  return
end if

if ( .not. ch_is_alpha ( s(1:1) ) ) then
  return
end if

do i = 2, s_length
  if ( .not. malphnum2 ( s(i:i) ) ) then
     return
  end if
end do

s_is_f90_name = .true.

return
end function s_is_f90_name
function s_is_i ( s, i )

  !*****************************************************************************80
  !
  !! S_IS_I is TRUE if a string represents an integer.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, integer ( kind = 4 ) I.  If the string represents an integer,
  !    I is the integer represented.  Otherwise I is 0.
  !
  !    Output, logical S_IS_I, is TRUE if the string represents an integer.
  !
implicit none

integer ( kind = 4 )  i
integer ( kind = 4 )  ierror
integer ( kind = 4 )  length
character ( len = * ) s
logical               s_is_i
integer ( kind = 4 )  s_length

s_length = len_trim ( s )

call s_to_i4 ( s, i, ierror, length )

if ( ierror == 0 .and. s_length <= length ) then
  s_is_i = .true.
else
  s_is_i = .false.
  i = 0
end if

return
end function s_is_i
subroutine s_is_r ( s, r, lval )

  !*****************************************************************************80
  !
  !! S_IS_R is TRUE if a string represents a real number.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, real ( kind = 4 ) R.  If the string represents a real number,
  !    then R is the real number represented.  Otherwise R is 0.
  !
  !    Output, logical LVAL, is TRUE if the string represents a real number.
  !
implicit none

integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) length
logical                lval
real      ( kind = 4 ) r
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

call s_to_r4 ( s, r, ierror, length )

if ( ierror == 0 .and. s_length <= length ) then
  lval = .true.
else
  lval = .false.
  r = 0.0E+00
end if

return
end subroutine s_is_r
subroutine s_left_insert ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_LEFT_INSERT inserts one string flush left into another.
  !
  !  Discussion:
  !
  !    S2 is not blanked out first.  Therefore, if there is
  !    already information in S2, some of it may still be around
  !    after S1 is written into S2.  Users may want to first
  !    assign S2 = ' ' if this is not the effect desired.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, a string to be inserted into S2.  Only
  !    the portion of S1 up to the last nonblank will be used.
  !
  !    Output, character ( len = * ) S2, a string which will contain,
  !    on output, a left flush copy of S1.
  !
implicit none

integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) jhi
integer   ( kind = 4 ) jlo
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2
integer   ( kind = 4 ) s2_length

s1_length = len_trim ( s1 )
s2_length = len ( s2 )

if ( s1_length < s2_length ) then
  ilo = 1
  ihi = s1_length
  jlo = 1
  jhi = s1_length
else if ( s2_length < s1_length ) then
  ilo = 1
  ihi = s2_length
  jlo = 1
  jhi = s2_length
else
  ilo = 1
  ihi = s1_length
  jlo = 1
  jhi = s2_length
end if

s2(jlo:jhi) = s1(ilo:ihi)

return
end subroutine s_left_insert
function s_lei ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_LEI = ( S1 is lexically less than or equal to S2 ).
  !
  !  Discussion:
  !
  !    The comparison is done in a case-insensitive way.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_LEI, the result of the comparison.
  !
implicit none

character              c1
character              c2
integer   ( kind = 4 ) i
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) lenc
logical                s_lei
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2

s1_length = len ( s1 )
len2 = len ( s2 )
lenc = min ( s1_length, len2 )

do i = 1, lenc

  c1 = s1(i:i)
  c2 = s2(i:i)
  call ch_cap ( c1 )
  call ch_cap ( c2 )

  if ( llt ( c1, c2 ) ) then
     s_lei = .true.
     return
  else if ( lgt ( c1, c2 ) ) then
     s_lei = .false.
     return
  end if

end do

if ( s1_length <= len2 ) then
  s_lei = .true.
else
  s_lei = .false.
end if

return
end function s_lei
subroutine s_low ( s )

  !*****************************************************************************80
  !
  !! S_LOW replaces all uppercase letters by lowercase ones.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be
  !    transformed.  On output, the string is all lowercase.
  !
implicit none

integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

do i = 1, s_length
  call ch_low ( s(i:i) )
end do

return
end subroutine s_low
function s_lti ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_LTI = ( S1 is lexically less than S2 ).
  !
  !  Discussion:
  !
  !    The comparison is done in a case-insensitive way.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_LTI, the result of the comparison.
  !
implicit none

character              c1
character              c2
integer   ( kind = 4 ) i
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) lenc
logical                s_lti
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2

s1_length = len ( s1 )
len2 = len ( s2 )
lenc = min ( s1_length, len2 )

do i = 1, lenc

  c1 = s1(i:i)
  c2 = s2(i:i)
  call ch_cap ( c1 )
  call ch_cap ( c2 )

  if ( llt ( c1, c2 ) ) then
     s_lti = .true.
     return
  else if ( lgt ( c1, c2 ) ) then
     s_lti = .false.
     return
  end if

end do

if ( s1_length < len2 ) then
  s_lti = .true.
else
  s_lti = .false.
end if

return
end function s_lti
function s_neqi ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_NEQI compares two strings for non-equality, ignoring case.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 November 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_NEQI, the result of the comparison.
  !
implicit none

character              c1
character              c2
integer   ( kind = 4 ) i
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) lenc
logical                s_neqi
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2

s1_length = len ( s1 )
len2 = len ( s2 )
lenc = min ( s1_length, len2 )

s_neqi = .true.

do i = 1, lenc

  c1 = s1(i:i)
  c2 = s2(i:i)
  call ch_cap ( c1 )
  call ch_cap ( c2 )

  if ( c1 /= c2 ) then
     return
  end if

end do

do i = lenc + 1, s1_length
  if ( s1(i:i) /= ' ' ) then
     return
  end if
end do

do i = lenc + 1, len2
  if ( s2(i:i) /= ' ' ) then
     return
  end if
end do

s_neqi = .false.

return
end function s_neqi
function s_no_control ( s )

  !*****************************************************************************80
  !
  !! S_NO_CONTROL = string contains no control characters.
  !
  !  Discussion:
  !
  !    Non-control characters are ASCII codes 32 through 127 inclusive.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 January 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, is the string to be checked.
  !
  !    Output, logical S_NO_CONTROL, is TRUE if S contains only printable
  !    characters, FALSE otherwise.
  !
implicit none

logical                ch_is_control
integer   ( kind = 4 ) i
logical                s_no_control
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_no_control = .false.
s_length = len_trim ( s )

do i = 1, s_length
  if ( ch_is_control ( s(i:i) ) ) then
     return
  end if
end do

s_no_control = .true.

return
end function s_no_control
subroutine s_nonalpha_delete ( s )

  !*****************************************************************************80
  !
  !! S_NONALPHA_DELETE removes nonalphabetic characters from a string.
  !
  !  Discussion:
  !
  !    The remaining characters are left justified and blank padded.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 August 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
implicit none

character              ch
integer   ( kind = 4 ) get
integer   ( kind = 4 ) put
character ( len = * )  s
integer   ( kind = 4 ) s_length

put = 0
s_length = len_trim ( s )

do get = 1, s_length

  ch = s(get:get)

  if ( ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) .or. &
       ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) ) then
     put = put + 1
     s(put:put) = ch
  end if

end do

s(put+1:s_length) = ' '

return
end subroutine s_nonalpha_delete
function s_of_i4 ( i )

  !*****************************************************************************80
  !
  !! S_OF_I4 converts an integer to a left-justified string.
  !
  !  Example:
  !
  !         I  S
  !
  !         1  1
  !        -1  -1
  !         0  0
  !      1952  1952
  !    123456  123456
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    13 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, an integer to be converted.
  !
  !    Output, character ( len = 11 ) S_OF_I4, the representation of the
  !    integer ( kind = 4 ).  The integer will be left-justified.
  !
implicit none

character              c
integer   ( kind = 4 ) i
integer   ( kind = 4 ) idig
integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) ipos
integer   ( kind = 4 ) ival
integer   ( kind = 4 ) j
character ( len = 11 ) s
character ( len = 11 ) s_of_i4

s = ' '

ilo = 1
ihi = 11
!
!  Make a copy of the integer.
!
ival = i
!
!  Handle the negative sign.
!
if ( ival < 0 ) then

  if ( ihi <= 1 ) then
     s(1:1) = '*'
     return
  end if

  ival = -ival
  s(1:1) = '-'
  ilo = 2

end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
ipos = ihi
!
!  Find the last digit, strip it off, and stick it into the string.
!
do

  idig = mod ( ival, 10 )
  ival = ival / 10

  if ( ipos < ilo ) then
     do j = 1, ihi
        s(j:j) = '*'
     end do
     return
  end if

  call digit_to_ch ( idig, c )

  s(ipos:ipos) = c
  ipos = ipos - 1

  if ( ival == 0 ) then
     exit
  end if

end do
!
!  Shift the string to the left.
!
s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
s(ilo+ihi-ipos:ihi) = ' '

s_of_i4 = s

return
end function s_of_i4
function s_only_alphab ( s )

  !*****************************************************************************80
  !
  !! S_ONLY_ALPHAB checks if a string is only alphabetic and blanks.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    Acceptable characters are 'A' through 'Z' and 'a' through 'z' and blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, logical S_ONLY_ALPHAB, is TRUE if the string contains only
  !    alphabetic characters and blanks.
  !
implicit none

character              c
integer   ( kind = 4 ) i
integer   ( kind = 4 ) itemp
character ( len = * )  s
integer   ( kind = 4 ) s_length
logical                s_only_alphab

s_only_alphab = .false.
s_length = len_trim ( s )

do i = 1, s_length

  c = s(i:i)

  if ( c /= ' ' ) then

     itemp = iachar ( c )

     if ( .not. ( 65 <= itemp .and. itemp <= 90 ) ) then
        if ( .not. ( 97 <= itemp .and. itemp <= 122 ) ) then
           return
        end if
     end if

  end if

end do

s_only_alphab = .true.

return
end function s_only_alphab
function s_only_digitb ( s )

  !*****************************************************************************80
  !
  !! S_ONLY_DIGITB returns TRUE if the string contains only digits or blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be checked.
  !
  !    Output, logical S_ONLY_DIGITB, is TRUE if the string contains only digits
  !    and blanks.
  !
implicit none

character              c
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length
logical                s_only_digitb

s_only_digitb = .false.
s_length = len_trim ( s )

do i = 1, s_length

  c = s(i:i)

  if ( c /= ' ' ) then
     if ( llt ( c, '0' ) .or. lgt ( c, '9' ) ) then
        return
     end if
  end if

end do

s_only_digitb = .true.

return
end function s_only_digitb
subroutine s_overlap ( s1, s2, overlap )

  !*****************************************************************************80
  !
  !! S_OVERLAP determines the overlap between two strings.
  !
  !  Discussion:
  !
  !    To determine the overlap, write the first word followed immediately
  !    by the second word.  Find the longest substring S which is both
  !    a suffix of S1 and a prefix of S2.  The length of this substring
  !    is the overlap.
  !
  !  Example:
  !
  !    S1              S2        OVERLAP
  !
  !    'timber'        'beret'   3
  !    'timber'        'timber'  6
  !    'beret'         'timber'  1
  !    'beret'         'berets'  5
  !    'beret'         'berth'   0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to be checked.
  !
  !    Output, integer ( kind = 4 ) OVERLAP, the length of the overlap.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) len3
integer   ( kind = 4 ) overlap
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2

overlap = 0

s1_length = len_trim ( s1 )
len2 = len_trim ( s2 )
len3 = min ( s1_length, len2 )

do i = 1, len3
  if ( s1(s1_length+1-i:s1_length) == s2(1:i) ) then
     overlap = i
  end if
end do

return
end subroutine s_overlap
function s_paren_check ( s )

  !*****************************************************************************80
  !
  !! S_PAREN_CHECK checks the parentheses in a string.
  !
  !  Discussion:
  !
  !    Blanks are removed from the string, and then the following checks
  !    are made:
  !
  !    1) as we read the string left to right, there must never be more
  !       right parentheses than left ones;
  !    2) there must be an equal number of left and right parentheses;
  !    3) there must be no occurrences of the abutting packages '...)(...'.
  !    4) there must be no occurrences of the empty package '()'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 November 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to check.
  !
  !    Output, logical S_PAREN_CHECK is TRUE if the string passed the checks.
  !
implicit none

integer   ( kind = 4 )  i
integer   ( kind = 4 )  isum
character ( len = * )   s
character ( len = 255 ) s_copy
integer   ( kind = 4 )  s_length
logical                 s_paren_check

s_copy = s
call s_blank_delete ( s_copy)

s_length = len_trim ( s_copy )
!
!  1) Letting '(' = +1 and ')' = -1, check that the running parentheses sum
!  is always nonnegative.
!
isum = 0
do i = 1, s_length

  if ( s_copy(i:i) == '(' ) then
     isum = isum + 1
  end if

  if ( s_copy(i:i) == ')' ) then

     isum = isum - 1

     if ( isum < 0 ) then
        s_paren_check = .false.
        return
     end if

  end if

end do
!
!  2) Check that the final parentheses sum is zero.
!
if ( isum /= 0 ) then
  s_paren_check = .false.
  return
end if
!
!  3) Check that there are no ")(" pairs.
!
do i = 2, s_length
  if ( s_copy(i-1:i) == ')(' ) then
     s_paren_check = .false.
     return
  end if
end do
!
!  4) Check that there are no "()" pairs.
!
do i = 2, s_length

  if ( s_copy(i-1:i) == '()' ) then
     s_paren_check = .false.
     return
  end if

end do
!
!  The checks were passed.
!
s_paren_check = .true.

return
end function s_paren_check
subroutine s_r_append ( s, r, done )

  !*****************************************************************************80
  !
  !! S_R_APPEND appends a real number to a string.
  !
  !  Discussion:
  !
  !    A blank space will separate the value from the text already
  !    in the line.
  !
  !    The routine warns the user if the value will not fit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 December 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a line of text.
  !    On input, the current string.  On output, the current string
  !    with the integer value appended.
  !
  !    Input, real ( kind = 4 ) R, the real number to be appended to the line.
  !
  !    Output, logical DONE, is FALSE if there was not enough room
  !    to append the data.
  !
implicit none

logical                done
integer   ( kind = 4 ) lens
integer   ( kind = 4 ) lents
integer   ( kind = 4 ) lenw
integer   ( kind = 4 ) next
real      ( kind = 4 ) r
character ( len = * )  s
character ( len = 14 ) w

done = .false.

lens = len ( s )
lents = len_trim ( s )

call r4_to_s_left ( r, w )

lenw = len_trim ( w )

if ( lents == 0 ) then
  if ( lens < lenw ) then
     done = .true.
     return
  end if
else
  if ( lens < lents + 1 + lenw ) then
     done = .true.
     return
  end if
end if

if ( lents == 0 ) then
  next = 1
else
  next = lents + 1
  s(next:next) = ' '
  next = next + 1
end if

s(next:next+lenw-1) = w(1:lenw)

return
end subroutine s_r_append
subroutine s_replace_ch ( s, c1, c2 )

  !*****************************************************************************80
  !
  !! S_REPLACE_CH replaces all occurrences of one character by another.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 March 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string.
  !
  !    Input, character C1, C2, the character to be replaced, and the
  !    replacement character.
  !
implicit none

character              c1
character              c2
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

do i = 1, s_length
  if ( s(i:i) == c1 ) then
     s(i:i) = c2
  end if
end do

return
end subroutine s_replace_ch
subroutine s_replace_one ( s1, sub1, sub2, s2 )

  !*****************************************************************************80
  !
  !! S_REPLACE_ONE replaces the first occurrence of SUB1 with SUB2.
  !
  !  Discussion:
  !
  !    The input and output strings may coincide.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the initial string.
  !
  !    Input, character ( len = * ) SUB1, the string to be replaced.
  !
  !    Input, character ( len = * ) SUB2, the replacement string.
  !
  !    Output, character ( len = * ) S2, the final string.
  !
implicit none

integer   ( kind = 4 )  i1
integer   ( kind = 4 )  i2
integer   ( kind = 4 )  i3
integer   ( kind = 4 )  i4
character ( len = * )   s1
character ( len = * )   s2
character ( len = 255 ) s3
character ( len = * )   sub1
character ( len = * )   sub2

s3 = ' '

i1 = index ( s1, sub1 )

if ( i1 == 0 ) then

  s3 = s1

else

  s3(1:i1-1) = s1(1:i1-1)

  i2 = len_trim ( sub2 )
  s3(i1:i1+i2-1) = sub2(1:i2)

  i3 = i1 + len_trim ( sub1 )
  i4 = len_trim ( s1 )

  s3(i1+i2:i1+i2+1+i4-i3) = s1(i3:i4)

end if

s2 = s3

return
end subroutine s_replace_one
subroutine s_replace_rec ( s, sub1, sub2, irep )

  !*****************************************************************************80
  !
  !! S_REPLACE_REC is a recursive replacement of one string by another.
  !
  !  Discussion:
  !
  !    All occurrences of SUB1 should be replaced by SUB2.
  !    This is not always true if SUB2 is longer than SUB1.
  !    The replacement is recursive.  In other words, replacing all
  !    occurrences of "ab" by "a" in "abbbbb" will return "a" rather
  !    than "abbbb".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.  On input,
  !    the string in which occurrences are to be replaced.  On
  !    output, the revised string.
  !
  !    Input, character ( len = * ) SUB1, the string which is to be replaced.
  !
  !    Input, character ( len = * ) SUB2, the replacement string.
  !
  !    Output, integer ( kind = 4 ) IREP, the number of replacements made.
  !    If IREP is negative, then its absolute value is the
  !    number of replacements made, and SUB2 is longer than
  !    SUB1, and at least one substring SUB1 could not be
  !    replaced by SUB2 because there was no more space in
  !    S.  (If S = 'aab' and SUB1 = 'a' and SUB2 = 'cc'
  !    then the result would be S = 'cca'.  The first 'a'
  !    was replaced, the 'b' fell off the end, the second 'a'
  !    was not replaced because the replacement 'cc' would
  !    have fallen off the end)
  !
implicit none

integer   ( kind = 4 ) irep
integer   ( kind = 4 ) len1
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) loc
character ( len = * )  s
integer   ( kind = 4 ) s_length
character ( len = * )  sub1
character ( len = * )  sub2

irep = 0
s_length = len ( s )

if ( s_length <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_REPLACE_REC - Serious error!'
  write ( *, '(a)' ) '  Null string not allowed!'
  return
end if

len1 = len ( sub1 )

if ( len1 <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_REPLACE_REC - Serious error!'
  write ( *, '(a)' ) '  Null SUB1 not allowed!'
  return
end if

len2 = len ( sub2 )

if ( len2 == len1 ) then

  if ( sub1 == sub2 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'S_REPLACE_REC - Warning!'
     write ( *, '(a)' ) '  Replacement = original!'
     return
  end if

  do

     loc = index ( s, sub1 )

     if ( loc == 0 ) then
        exit
     end if

     irep = irep + 1
     s(loc:loc+len1-1) = sub2

  end do

else if ( len2 < len1 ) then

  do

     loc = index ( s, sub1 )

     if ( loc == 0 ) then
        exit
     end if

     irep = irep + 1
     s(loc:loc+len2-1) = sub2
     call s_chop ( s, loc+len2, loc+len1-1 )

  end do

else

  do

     loc = index ( s, sub1 )

     if ( loc == 0 ) then
        exit
     end if

     irep = irep + 1

     if ( s_length < loc + len2 - 1 ) then
        irep = -irep
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_REPLACE_REC - Warning!'
        write ( *, '(a)' ) '  Some replaceable elements remain!'
        return
     end if

     call s_blanks_insert ( s, loc, loc+len2-1-len1 )
     s(loc:loc+len2-1) = sub2

  end do

end if

return
end subroutine s_replace_rec
subroutine s_replace ( s, sub1, sub2, irep )

  !*****************************************************************************80
  !
  !! S_REPLACE replaces all occurrences of SUB1 by SUB2 in a string.
  !
  !  Discussion:
  !
  !    This is not always true if SUB2 is longer than SUB1.  The
  !    replacement is NOT recursive.  In other words, replacing all
  !    occurrences of "ab" by "a" in "abbbbb" will return "abbbb"
  !    rather than "a".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.  On input,
  !    the string in which occurrences are to be replaced.
  !    On output, the revised string.
  !
  !    Input, character ( len = * ) SUB1, the string which is to be replaced.
  !    Trailing blank characters are ignored.  The routine is case sensitive.
  !
  !    Input, character ( len = * ) SUB2, the replacement string.
  !
  !    Output, integer ( kind = 4 ) IREP, the number of replacements made.
  !    If IREP is negative, then its absolute value is the
  !    number of replacements made, and SUB2 is longer than
  !    SUB1, and at least one substring SUB1 could not be
  !    replaced by SUB2 because there was no more space.
  !    (If S = 'aab' and SUB1 = 'a' and SUB2 = 'cc'
  !    then the result would be S = 'cca'.  The first 'a'
  !    was replaced, the 'b' fell off the end, the second 'a'
  !    was not replaced because the replacement 'cc' would have
  !    fallen off the end)
  !
implicit none

integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) irep
integer   ( kind = 4 ) len1
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) lens
integer   ( kind = 4 ) loc
character ( len = * )  s
character ( len = * )  sub1
character ( len = * )  sub2

irep = 0
lens = len ( s )

if ( lens <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_REPLACE - Serious error!'
  write ( *, '(a)' ) '  Null string not allowed!'
  return
end if

len1 = len_trim ( sub1 )

if ( len1 <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_REPLACE - Serious error!'
  write ( *, '(a)' ) '  Null SUB1 not allowed!'
  return
end if

len2 = len_trim ( sub2 )

if ( len2 == len1 ) then

  if ( sub1(1:len1) == sub2(1:len2) ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'S_REPLACE - Warning!'
     write ( *, '(a)' ) '  Replacement = original!'
     return
  end if

  ilo = 1

  do

     loc = index ( s(ilo:lens), sub1(1:len1) )

     if ( loc == 0 ) then
        exit
     end if

     loc = loc + ilo - 1
     irep = irep + 1
     s(loc:loc+len1-1) = sub2(1:len2)
     ilo = loc + len1

     if ( lens < ilo ) then
        exit
     end if

  end do

else if ( len2 < len1 ) then

  ilo = 1

  do

     loc = index ( s(ilo:lens), sub1(1:len1) )

     if ( loc == 0 ) then
        exit
     end if

     irep = irep + 1
     loc = loc + ilo - 1
     s(loc:loc+len2-1) = sub2(1:len2)
     call s_chop ( s, loc+len2, loc+len1-1 )
     ilo = loc + len2

     if ( lens < ilo ) then
        exit
     end if

  end do

else

  ilo = 1

  do

     loc = index ( s(ilo:lens), sub1(1:len1) )

     if ( loc == 0 ) then
        exit
     end if

     loc = loc + ilo - 1
     irep = irep + 1

     if ( lens < loc + len2 - 1 ) then
        irep = -irep
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_REPLACE - Warning!'
        write ( *, '(a)' ) '  Some replaceable elements remain!'
        exit
     end if

     call s_blanks_insert ( s, loc, loc+len2-1-len1 )

     s(loc:loc+len2-1) = sub2(1:len2)
     ilo = loc + len2

  end do

end if

return
end subroutine s_replace
subroutine s_replace_i ( s, sub1, sub2 )

  !*****************************************************************************80
  !
  !! S_REPLACE_I replaces all occurrences of SUB1 by SUB2 in a string.
  !
  !  Discussion:
  !
  !    Matches are made without regard to case.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.  On input,
  !    the string in which occurrences are to be replaced.
  !    On output, the revised string.
  !
  !    Input, character ( len = * ) SUB1, the string which is to be replaced.
  !
  !    Input, character ( len = * ) SUB2, the replacement string.
  !
  !    Output, integer ( kind = 4 ) IREP, the number of replacements made.
  !    If IREP is negative, then its absolute value is the
  !    number of replacements made, and SUB2 is longer than
  !    SUB1, and at least one substring SUB1 could not be
  !    replaced by SUB2 because there was no more space.
  !    (If S = 'aab' and SUB1 = 'a' and SUB2 = 'cc'
  !    then the result would be S = 'cca'.  The first 'a'
  !    was replaced, the 'b' fell off the end, the second 'a'
  !    was not replaced because the replacement 'cc' would have
  !    fallen off the end)
  !
implicit none

integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) len1
integer   ( kind = 4 ) len2
integer   ( kind = 4 ) lens
integer   ( kind = 4 ) s_indexi
character ( len = * )  s
character ( len = * )  sub1
character ( len = * )  sub2

lens = len ( s )

if ( lens <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_REPLACE_I - Serious error!'
  write ( *, '(a)' ) '  Null string not allowed!'
  return
end if

len1 = len ( sub1 )

if ( len1 <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_REPLACE_I - Serious error!'
  write ( *, '(a)' ) '  Null SUB1 not allowed!'
  return
end if

len2 = len ( sub2 )

ilo = s_indexi ( s, sub1 )
!
!  If the match string has been found, then insert the replacement.
!
if ( ilo /= 0 ) then
  s(ilo+len2:lens+len2-len1) = s(ilo+len1:lens)
  s(ilo:ilo+len2-1) = sub2(1:len2)
end if

return
end subroutine s_replace_i
subroutine s_reverse ( s )

  !*****************************************************************************80
  !
  !! S_REVERSE reverses the characters in a string.
  !
  !  Example:
  !
  !    Input        Output
  !
  !    ' Cat'       'taC '
  !    'Goo gol  '  'log ooG  '
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 November 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to reverse.
  !    Trailing blanks are ignored.
  !
implicit none

character              ch
integer   ( kind = 4 ) i
integer   ( kind = 4 ) j
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

do i = 1, s_length / 2
  j = s_length + 1 - i
  ch     = s(i:i)
  s(i:i) = s(j:j)
  s(j:j) = ch
end do

return
end subroutine s_reverse
subroutine s_right_insert ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_RIGHT_INSERT inserts a string flush right into another.
  !
  !  Discussion:
  !
  !    S2 is not blanked out first.  If there is already information in S2,
  !    some of it may still be around after S1 is written into S2.  Users may
  !    want to first assign S2 = ' ' if this is not the effect desired.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, a string which is to be
  !    inserted into S2.  Only the portion of S1 up to the last
  !    nonblank will be used.
  !
  !    Output, character ( len = * ) S2, a string whose length
  !    will be determined by a call to LEN, and which will
  !    contain, on output, a right flush copy of S1.
  !
implicit none

integer ( kind = 4 ) ihi
integer ( kind = 4 ) ilo
integer ( kind = 4 ) jhi
integer ( kind = 4 ) jlo
integer ( kind = 4 ) len1
integer ( kind = 4 ) len2
character ( len = * ) s1
character ( len = * ) s2

len1 = len_trim ( s1 )
len2 = len ( s2 )

if ( len1 < len2 ) then
  ilo = 1
  ihi = len1
  jlo = len2 + 1 - len1
  jhi = len2
else if ( len2 < len1 ) then
  ilo = len1 + 1 - len2
  ihi = len1
  jlo = 1
  jhi = len2
else
  ilo = 1
  ihi = len1
  jlo = 1
  jhi = len2
end if

s2(jlo:jhi) = s1(ilo:ihi)

return
end subroutine s_right_insert
subroutine s_roman_to_i4 ( s, i )

  !*****************************************************************************80
  !
  !! S_ROMAN_TO_I4 converts a Roman numeral to an integer.
  !
  !  Example:
  !
  !    S      I
  !
  !    X      10
  !    XIX    19
  !    MI     1001
  !    CXC    190
  !
  !  Discussion:
  !
  !    The subroutine does not check carefully as to whether the Roman numeral
  !    is properly formed.  In particular, it will accept a string like 'IM'
  !    and return 999, even though this is not a well formed Roman numeral.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string containing a Roman numeral.
  !
  !    Output, integer ( kind = 4 ) I, the corresponding value.
  !
implicit none

integer ( kind = 4 ) ch_roman_to_i4
character c1
character c2
logical done
integer ( kind = 4 ) i
integer ( kind = 4 ) i1
integer ( kind = 4 ) i2
character ( len = * ) s

i = 0
done = .true.

do

  call ch_next ( s, c2, done )

  if ( done ) then
     return
  end if

  i2 = ch_roman_to_i4 ( c2 )

  if ( i2 == 0 .and. c2 /= ' ' ) then
     return
  end if

  do

     c1 = c2
     i1 = i2

     call ch_next ( s, c2, done )

     if ( done ) then
        i = i + i1
        return
     end if

     i2 = ch_roman_to_i4 ( c2 )

     if ( i2 == 0 .and. c2 /= ' ' ) then
        i = i + i1
        return
     end if

     if ( i1 < i2 ) then
        i = i + i2 - i1
        c1 = ' '
        c2 = ' '
        exit
     end if

     i = i + i1

  end do

end do

return
end subroutine s_roman_to_i4
subroutine s_s_delete ( s, sub, irep )

  !*****************************************************************************80
  !
  !! S_S_DELETE removes all occurrences of a substring from a string.
  !
  !  Discussion:
  !
  !    The remainder is left justified and padded with blanks.
  !
  !    The deletion is not recursive.  Removing all occurrences of "ab" from
  !    "aaaaabbbbbQ" results in "aaaabbbbQ" rather than "Q".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, character ( len = * ) SUB1, the substring to be removed.
  !
  !    Output, integer ( kind = 4 ) IREP, the number of occurrences of SUB1
  !    which were found.
  !
implicit none

integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) irep
integer   ( kind = 4 ) loc
integer   ( kind = 4 ) nsub
character ( len = * )  s
character ( len = * )  sub

nsub = len_trim ( sub )

irep = 0
ilo = 1
ihi = len_trim ( s )

do while ( ilo <= ihi )

  loc = index ( s(ilo:ihi), sub )

  if ( loc == 0 ) then
     return
  end if

  irep = irep + 1
  loc = loc + ilo - 1
  call s_chop ( s, loc, loc+nsub-1 )
  ilo = loc
  ihi = ihi - nsub

end do

return
end subroutine s_s_delete
subroutine s_s_delete2 ( s, sub, irep )

  !*****************************************************************************80
  !
  !! S_S_DELETE2 recursively removes a substring from a string.
  !
  !  Discussion:
  !
  !    The remainder is left justified and padded with blanks.
  !
  !    The substitution is recursive, so
  !    that, for example, removing all occurrences of "ab" from
  !    "aaaaabbbbbQ" results in "Q".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
  !    Input, character ( len = * ) SUB, the substring to be removed.
  !
  !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
  !    the substring.
  !
implicit none

integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) irep
integer   ( kind = 4 ) loc
integer   ( kind = 4 ) nsub
character ( len = * )  s
integer   ( kind = 4 ) s_length
character ( len = * )  sub

s_length = len ( s )
nsub = len ( sub )

irep = 0
ihi = s_length

do while ( 0 < ihi )

  loc = index ( s(1:ihi), sub )

  if ( loc == 0 ) then
     return
  end if

  irep = irep + 1
  call s_chop ( s, loc, loc+nsub-1 )
  ihi = ihi - nsub

end do

return
end subroutine s_s_delete2
subroutine s_s_insert ( s1, ipos, s2 )

  !*****************************************************************************80
  !
  !! S_S_INSERT inserts a substring into a string.
  !
  !  Discussion:
  !
  !    Characters in the string are moved to the right to make room, and
  !    hence the trailing characters, if any, are lost.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S1, the string into which
  !    the second string is to be inserted.
  !
  !    Input, integer ( kind = 4 ) IPOS, the position in S at which S2 is
  !    to be inserted.
  !
  !    Input, character ( len = * ) S2, the string to be inserted.
  !
implicit none

integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ipos
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2
integer   ( kind = 4 ) s2_length

s1_length = len ( s1 )
s2_length = len_trim ( s2 )

ihi = min ( s1_length, ipos+s2_length-1 )

call s_blanks_insert ( s1, ipos, ihi )

s1(ipos:ihi) = s2

return
end subroutine s_s_insert
function s_s_subanagram ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_S_SUBANAGRAM determines if S2 is a "subanagram" of S1.
  !
  !  Discussion:
  !
  !    S2 is an anagram of S1 if S2 can be formed by permuting the letters
  !    of S1
  !
  !    S2 is an subanagram of S1 if S2 can be formed by selecting SOME of
  !    the letters of S1 and permuting them.
  !
  !    Blanks (trailing or otherwise), punctuation, and capitalization
  !    are all significant, so be sure to input exactly the information
  !    you want to check.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 June 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the master string.
  !
  !    Input, character ( len = * ) S2, the second string.
  !
  !    Output, logical S_S_SUBANAGRAM is TRUE if S2 is a subanagram of S1.
  !
implicit none

integer   ( kind = 4 ) i1
integer   ( kind = 4 ) i2
logical                s_s_subanagram
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2
integer   ( kind = 4 ) s2_length

s_s_subanagram = .false.
!
!  Sort both.
!
call s_sort_a ( s1 )
call s_sort_a ( s2 )

s1_length = len ( s1 )
s2_length = len ( s2 )

i1 = 0

do i2 = 1, s2_length

  do

     i1 = i1 + 1
     !
     !  Ran out of S1 before finishing.  No match is possible.
     !
     if ( s1_length < i1 ) then
        return
     end if
     !
     !  The current character in S1 is already greater than the character in S2.
     !  No match is possible.
     !
     if ( llt ( s2(i2:i2), s1(i1:i1) ) ) then
        return
     end if
     !
     !  Found an exact match for current character.  Keep going.
     !
     if ( s1(i1:i1) == s2(i2:i2) ) then
        exit
     end if
     !
     !  Didn't find a match, but one might be possible if we increase I1.
     !

  end do

end do
!
!  We matched every character of S2 with something in S1.
!
s_s_subanagram = .true.

return
end function s_s_subanagram
function s_s_subanagram_sorted ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_S_SUBANAGRAM_SORTED determines if S2 is a "subanagram" of S1.
  !
  !  Discussion:
  !
  !    This routine assumes that S1 and S2 have already been sorted.
  !
  !    S2 is an anagram of S1 if S2 can be formed by permuting the letters
  !    of S1
  !
  !    S2 is an subanagram of S1 if S2 can be formed by selecting SOME of
  !    the letters of S1 and permuting them.
  !
  !    Blanks (trailing or otherwise), punctuation, and capitalization
  !    are all significant, so be sure to input exactly the information
  !    you want to check.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 June 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the master string.
  !
  !    Input, character ( len = * ) S2, the second string.
  !
  !    Output, logical S_S_SUBANAGRAM_SORTED is TRUE if S2 is a subanagram of S1.
  !
implicit none

integer   ( kind = 4 ) i1
integer   ( kind = 4 ) i2
logical                s_s_subanagram_sorted
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2
integer   ( kind = 4 ) s2_length

s_s_subanagram_sorted = .false.

s1_length = len ( s1 )
s2_length = len ( s2 )

i1 = 0

do i2 = 1, s2_length

  do

     i1 = i1 + 1
     !
     !  Ran out of S1 before finishing.  No match is possible.
     !
     if ( s1_length < i1 ) then
        return
     end if
     !
     !  The current character in S1 is already greater than the character in S2.
     !  No match is possible.
     !
     if ( llt ( s2(i2:i2), s1(i1:i1) ) ) then
        return
     end if
     !
     !  Found an exact match for current character.  Keep going.
     !
     if ( s1(i1:i1) == s2(i2:i2) ) then
        exit
     end if
     !
     !  Didn't find a match, but one might be possible if we increase I1.
     !
  end do

end do
!
!  We matched every character of S2 with something in S1.
!
s_s_subanagram_sorted = .true.

return
end function s_s_subanagram_sorted
subroutine s_set_delete ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_SET_DELETE removes any characters in one string from another string.
  !
  !  Discussion:
  !
  !    When an element is removed, the rest of the string is shifted to the
  !    left, and padded with blanks on the right.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be examined.
  !
  !    Input, character ( len = * ) S2, the characters to be removed.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) j
integer   ( kind = 4 ) nset
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2

s1_length = len ( s1 )
nset = len ( s2 )

i = 0

do while ( i < s1_length )

  i = i + 1

  do j = 1, nset
     if ( s1(i:i) == s2(j:j) ) then
        call s_chop ( s1, i, i )
        s1_length = s1_length - 1
        i = i - 1
        exit
     end if
  end do

end do

return
end subroutine s_set_delete
subroutine s_shift_circular ( s, ishft )

  !*****************************************************************************80
  !
  !! S_SHIFT_CIRCULAR circular shifts the characters in a string to the right.
  !
  !  Discussion:
  !
  !    Thus, a shift of -1 would change "Violin" to "iolinV", and a shift
  !    of 1 would change it to "nVioli".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be shifted.
  !
  !    Input, integer ( kind = 4 ) ISHFT, the number of positions to the
  !    right to shift the characters.
  !
implicit none

character              chrin
character              chrout
integer   ( kind = 4 ) icycle
integer   ( kind = 4 ) idid
integer   ( kind = 4 ) igoto
integer   ( kind = 4 ) imove
integer   ( kind = 4 ) ishft
integer   ( kind = 4 ) jshft
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len ( s )

if ( s_length <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_SHIFT_CIRCULAR - Serious error!'
  write ( *, '(a)' ) '  String has nonpositive length!'
  return
end if
!
!  Force the shift to be positive and between 0 and S_LENGTH.
!
jshft = ishft

do while ( jshft < 0 )
  jshft = jshft + s_length
end do

do while ( s_length < jshft )
  jshft = jshft - s_length
end do

if ( jshft == 0 ) then
  return
end if
!
!  Shift the first character.  Shift the character that got
!  displaced by the first character...Repeat until you've shifted
!  all, or have "cycled" back to the first character early.
!
!  If you've cycled, start again at the second character, and
!  so on.
!
icycle = 0
idid = 0
imove = 0

do while ( idid < s_length )

  if ( imove == icycle ) then
     imove = imove + 1
     icycle = icycle + 1
     chrin = s(imove:imove)
  end if

  idid = idid + 1
  igoto = imove + jshft

  if ( s_length < igoto ) then
     igoto = igoto - s_length
  end if

  chrout = s(igoto:igoto)
  s(igoto:igoto) = chrin
  chrin = chrout

  imove = igoto

end do

return
end subroutine s_shift_circular
subroutine s_shift_left ( s, ishft )

  !*****************************************************************************80
  !
  !! S_SHIFT_LEFT shifts the characters in a string to the left and blank pads.
  !
  !  Discussion:
  !
  !    A shift of 2 would change "Violin" to "olin  ".
  !    A shift of -2 would change "Violin" to "  Violin".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be shifted.
  !
  !    Input, integer ( kind = 4 ) ISHFT, the number of positions to the
  !    left to shift the characters.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ishft
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len ( s )

if ( 0 < ishft ) then

  do i = 1, s_length - ishft
     s(i:i) = s(i+ishft:i+ishft)
  end do

  do i = s_length - ishft + 1, s_length
     s(i:i) = ' '
  end do

else if ( ishft < 0 ) then

  do i = s_length, - ishft + 1, - 1
     s(i:i) = s(i+ishft:i+ishft)
  end do

  do i = -ishft, 1, -1
     s(i:i) = ' '
  end do

end if

return
end subroutine s_shift_left
subroutine s_shift_right ( s, ishft )

  !*****************************************************************************80
  !
  !! S_SHIFT_RIGHT shifts the characters in a string to the right and blank pads.
  !
  !  Discussion:
  !
  !    A shift of 2 would change "Violin" to "  Viol".
  !    A shift of -2 would change "Violin" to "olin  ".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be shifted.
  !
  !    Input, integer ( kind = 4 ) ISHFT, the number of positions to the
  !    right to shift the characters.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ishft
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len ( s )

if ( s_length <= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_SHIFT_RIGHT - Serious error!'
  write ( *, '(a)' ) '  String has nonpositive length!'
  return
end if

if ( 0 < ishft ) then

  do i = s_length, ishft + 1, - 1
     s(i:i) = s(i-ishft:i-ishft)
  end do

  do i = ishft, 1, -1
     s(i:i) = ' '
  end do

else if ( ishft < 0 ) then

  do i = 1, s_length + ishft
     s(i:i) = s(i-ishft:i-ishft)
  end do

  do i = s_length + ishft + 1, s_length
     s(i:i) = ' '
  end do

end if

end subroutine s_shift_right
function s_skip_set ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_SKIP_SET finds the first entry of a string that is NOT in a set.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the string to be examined.
  !
  !    Input, character ( len = * ) S2, the characters to skip.
  !
  !    Output, integer ( kind = 4 ) S_SKIP_SET, the location of the first
  !    character in S1 that is not in S2, or 0 if no such index was found.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) s_skip_set
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
character ( len = * )  s2

s1_length = len_trim ( s1 )

do i = 1, s1_length

  if ( index ( s2, s1(i:i) ) == 0 ) then
     s_skip_set = i
     return
  end if

end do

s_skip_set = 0

return
end function s_skip_set
subroutine s_sort_a ( s )

  !*****************************************************************************80
  !
  !! S_SORT_A sorts a string into ascending order.
  !
  !  Discussion:
  !
  !    The string is assumed to be short, and so a simple bubble sort is used.
  !
  !    ALL the characters are sorted, including blanks and punctuation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 June 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be sorted.
  !
implicit none

character              c
integer   ( kind = 4 ) i
integer   ( kind = 4 ) j
integer   ( kind = 4 ) k
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len ( s )

do i = 1, s_length - 1

  c = s(i:i)
  j = i

  do k = i + 1, s_length
     if ( iachar ( s(k:k) ) < iachar ( s(j:j) ) ) then
        j = k
     end if
  end do

  if ( i /= j ) then
     s(i:i) = s(j:j)
     s(j:j) = c
  end if

end do

return
end subroutine s_sort_a
subroutine s_split ( s, sub, s1, s2, s3 )

  !*****************************************************************************80
  !
  !! S_SPLIT divides a string into three parts, given the middle.
  !
  !  Discussion:
  !
  !    This version of the routine is case-insensitive.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = 'aBCdEfgh'
  !      S2 = 'eF'
  !
  !    Output:
  !
  !      S1 = 'aBCd'
  !      S2 =  'gh'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be analyzed.
  !
  !    Input, character ( len = * ) SUB, the substring used to "split" S.
  !    Trailing blanks in SUB are ignored.
  !
  !    Output, character ( len = * ) S1, the entries in the string, up
  !    to, but not including, the first occurrence, if any,
  !    of SUB.  If SUB occurs immediately, then S1 = ' '.
  !    If SUB is not long enough, trailing entries will be lost.
  !
  !    Output, character ( len = * ) S2, the part of the string that matched SUB.
  !    If S2 is ' ', then there wasn't a match.
  !
  !    Output, character ( len = * ) S3, the part of the string after the match.
  !    If there was no match, then S3 is blank.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) lenm
character ( len = * )  s
integer   ( kind = 4 ) s_indexi
integer   ( kind = 4 ) s_length
character ( len = * )  s1
character ( len = * )  s2
character ( len = * )  s3
character ( len = * )  sub

s_length = len_trim ( s )

lenm = len_trim ( sub )

if ( lenm == 0 ) then
  lenm = 1
end if

i = s_indexi ( s, sub )
!
!  The substring did not occur.
!
if ( i == 0 ) then
  s1 = s
  s2 = ' '
  s3 = ' '
  !
  !  The substring begins immediately.
  !
else if ( i == 1 ) then
  s1 = ' '
  s2 = s(1:lenm)
  s3 = s(lenm+1:)
  !
  !  What am I checking here?
  !
else if ( s_length < i + lenm ) then
  s1 = s
  s2 = ' '
  s3 = ' '
  !
  !  The substring occurs in the middle.
  !
else
  s1 = s(1:i-1)
  s2 = s(i:i+lenm-1)
  s3 = s(i+lenm: )
end if
!
!  Drop leading blanks.
!
s1 = adjustl ( s1 )
s2 = adjustl ( s2 )
s3 = adjustl ( s3 )

return
end subroutine s_split
subroutine s_swap ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_SWAP swaps two strings.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S1, S2.  On output, the values of S1
  !    and S2 have been interchanged.
  !
implicit none

character ( len = * )   s1
character ( len = * )   s2
character ( len = 255 ) s3

s3 = s1
s1 = s2
s2 = s3

return
end subroutine s_swap
subroutine s_tab_blank ( s )

  !*****************************************************************************80
  !
  !! S_TAB_BLANK replaces each TAB character by one space.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
implicit none

integer ( kind = 4 ) i
character ( len = * ) s
integer ( kind = 4 ) s_length
character, parameter :: tab = achar ( 9 )

s_length = len_trim ( s )

do i = 1, s_length

  if ( s(i:i) == tab ) then
     s(i:i) = ' '
  end if

end do

return
end subroutine s_tab_blank
subroutine s_tab_blanks ( s )

  !*****************************************************************************80
  !
  !! S_TAB_BLANKS replaces TAB characters by 6 spaces.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be modified.  On
  !    output, some significant characters at the end of S may have
  !    been lost.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) get
integer   ( kind = 4 ) put
integer   ( kind = 4 ) lenc
integer   ( kind = 4 ) lens
integer   ( kind = 4 ) ntab
character ( len = * )  s
character, parameter :: tab = achar ( 9 )
!
!  If no TAB's occur in the line, there is nothing to do.
!
if ( index ( s, tab ) == 0 ) then
  return
end if
!
!  Otherwise, find out how long the string is.
!
lenc = len_trim ( s )
lens = len ( s )
!
!  Count the TAB's.
!
ntab = 0
do i = 1, lenc
  if ( s(i:i) == tab ) then
     ntab = ntab + 1
  end if
end do
!
!  Now copy the string onto itself, going backwards.
!  As soon as we've processed the first TAB, we're done.
!
put = lenc + 5 * ntab

do get = lenc, 1, - 1

  if ( s(get:get) /= tab ) then

     if ( put <= lens ) then
        s(put:put) = s(get:get)
     end if

     put = put - 1

  else

     do i = put, put - 5, -1
        if ( i <= lens ) then
           s(i:i) = ' '
        end if
     end do

     put = put - 6
     ntab = ntab - 1

     if ( ntab == 0 ) then
        return
     end if

  end if

end do

return
end subroutine s_tab_blanks
subroutine s_to_c4 ( s, cval, ierror, length )

  !*****************************************************************************80
  !
  !! S_TO_C4 reads a complex number from a string.
  !
  !  Discussion:
  !
  !    A C4 is simply a complex number to be stored as a
  !    "complex ( kind = 4 )" value.
  !
  !    This routine will read as many characters as possible until it reaches
  !    the end of the string, or encounters a character which cannot be
  !    part of the number.
  !
  !  Legal input is:
  !
  !     1 blanks,
  !     2 '+' or '-' sign,
  !     3 integer part,
  !     4 decimal point,
  !     5 fraction part,
  !     6 'E' or 'e' or 'D' or 'd', exponent marker,
  !     7 exponent sign,
  !     8 exponent integer part,
  !     9 exponent decimal point,
  !    10 exponent fraction part,
  !    11 blanks,
  !    12 '+' or '-' sign,
  !    13 integer part,
  !    14 decimal point,
  !    15 fraction part,
  !    16 'E' or 'e' or 'D' or 'd', exponent marker,
  !    17 exponent sign,
  !    18 exponent integer part,
  !    19 exponent decimal point,
  !    20 exponent fraction part,
  !    21 blanks,
  !    22 "*"
  !    23 spaces
  !    24 I
  !    25 comma or semicolon
  !
  !    with most quantities optional.
  !
  !  Example:
  !
  !    S               CVAL      IERROR     LENGTH
  !
  !    '1'               1         0          1
  !    '1+I'             1 + 1 i   0          3
  !    '1+1 i'           1 + 1 i   0          5
  !    '1+1*i'           1 + 1 i   0          5
  !    'i'                   1 i   0          1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate at the end of the string, or when no more
  !    characters can be read to form a legal real.  Blanks,
  !    commas, or other nonnumeric data will, in particular,
  !    cause the conversion to halt.
  !
  !    Output, complex ( kind = 4 ) CVAL, the value that was read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    1, the string was empty.
  !    2, could not read A correctly.
  !    3, could not read B correctly.
  !    4, could not read I correctly.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
  !    the string to form the number, including any terminating
  !    characters such as a trailing comma or blanks.
  !
implicit none

real      ( kind = 4 ) aval
real      ( kind = 4 ) bval
character              c
character              c2
logical                ch_eqi
complex   ( kind = 4 ) cval
integer   ( kind = 4 ) ichr
integer   ( kind = 4 ) ichr2
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) length
logical                s_neqi
character ( len = * )  s
!
!  Initialize the return arguments.
!
ierror = 0
aval = 0.0E+00
bval = 0.0E+00
cval = cmplx ( aval, bval, kind = 4 )
length = 0
!
!  Get the length of the line, and if it's zero, return.
!
if ( len_trim ( s ) <= 0 ) then
  ierror = 1
  return
end if

call nexchr ( s, ichr, c )
!
!  If the next character is "I", then this number is 0+I.
!
if ( ch_eqi ( c, 'I' ) ) then
  aval = 0.0E+00
  bval = 1.0E+00
  length = length + ichr
  cval = cmplx ( aval, bval, kind = 4 )
  return
end if
!
!  OK, the next string has to be a number!
!
call s_to_r4 ( s, aval, ierror, ichr )

if ( ierror /= 0 ) then
  ierror = 2
  length = 0
  return
end if

length = length + ichr
!
!  See if this is a pure real number, because:
!
!    1) There's no more input left.
!
if ( len_trim ( s(length+1:) ) == 0 ) then
  cval = cmplx ( aval, bval, kind = 4 )
  return
end if
!
!    2) The last character read was a comma.
!
if ( s(length:length) == ',' .or. s(length:length) == ';' ) then
  cval = cmplx ( aval, bval, kind = 4 )
  return
end if
!
!  If the very next character is "I", then this is a pure
!  imaginary number!
!
call nexchr ( s(length+1:), ichr, c )

if ( ch_eqi ( c, 'I' ) ) then
  bval = aval
  aval = 0.0E+00
  length = length + ichr
  cval = cmplx ( aval, bval, kind = 4 )
  return
end if
!
!  If the very next character is "*" and the one after that is
!  "I", then this is a pure imaginary number!
!
if ( c == '*' ) then

  call nexchr ( s(length+ichr+1:), ichr2, c2 )

  if ( ch_eqi ( c2, 'I' ) ) then
     bval = aval
     aval = 0.0E+00
     length = length + ichr + ichr2
  end if

  cval = cmplx ( aval, bval, kind = 4 )
  return

end if
!
!  OK, now we've got A.  We have to be careful because the next
!  thing we see MIGHT be "+ I" or "- I" which we can't let CHRCTR
!  see, because it will have fits.  So let's check these two
!  possibilities.
!
call nexchr ( s(length+1:), ichr, c )
call nexchr ( s(length+1+ichr:), ichr2, c2 )

if ( ch_eqi ( c2, 'I' ) ) then

  if ( c == '+' ) then
     bval = 1
     length = length + ichr + ichr2
     cval = cmplx ( aval, bval, kind = 4 )
     return
  else if ( c == '-' ) then
     bval = -1
     length = length + ichr + ichr2
     cval = cmplx ( aval, bval, kind = 4 )
     return
  end if

end if
!
!  Read the next real number.
!
call s_to_r4 ( s(length+1:), bval, ierror, ichr )

if ( ierror /= 0 ) then
  ierror = 3
  length = 0
  return
end if

length = length + ichr
!
!  If the next character is a "*", that's OK, advance past it.
!
call nexchr ( s(length+1:), ichr, c )

if ( c == '*' ) then
  length = length + ichr
end if
!
!  Now we really do want the next character to be "I".
!
call nexchr ( s(length+1:), ichr, c )

if ( s_neqi ( c, 'I' ) ) then
  ierror = 4
  length = 0
  return
end if
!
!  Form the complex number.
!
cval = cmplx ( aval, bval, kind = 4 )

return
end subroutine s_to_c4
subroutine s_to_chvec ( s, n, chvec )

  !*****************************************************************************80
  !
  !! S_TO_CHVEC converts a string to a character vector.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string of characters.
  !
  !    Input/output, integer ( kind = 4 ) N.
  !    if N is -1, extract characters from 1 to len(S);
  !    if N is 0, extract characters up to the last nonblank;
  !    if N is positive, extract characters from 1 to N.
  !
  !    On output, N is the number of characters successfully extracted.
  !
  !    Output, character CHVEC(N), the characters extracted from S.
  !
implicit none

character              chvec(*)
integer   ( kind = 4 ) i
integer   ( kind = 4 ) n
character ( len = * )  s

if ( n <= - 1 ) then
  n = len ( s )
else if ( n == 0 ) then
  n = len_trim ( s )
else
  n = min ( n, len ( s ) )
end if

do i = 1, n
  chvec(i) = s(i:i)
end do

return
end subroutine s_to_chvec
subroutine s_to_date ( s1, s2 )

  !*****************************************************************************80
  !
  !! S_TO_DATE converts the F90 date string to a more usual format.
  !
  !  Example:
  !
  !    S1        S2
  !    --------  ----------------
  !    20010204  4 February 2001
  !    17760704  4 July 1776
  !    19520310  10 March 1952
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 February 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = 8 ) S1, the F90 date string returned by
  !    the routine DATE_AND_TIME.
  !
  !    Output, character ( len = * ) S2, a more usual format for the date.
  !    Allowing 16 characters for S2 should be sufficient for the
  !    forseeable future.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) m
character ( len = 8 )  month
character ( len = * )  s1
character ( len = * )  s2

if ( s1(7:7) == '0' ) then
  s2(1:1) = s1(8:8)
  i = 1
else
  s2(1:2) = s1(7:8)
  i = 2
end if

i = i + 1
s2(i:i) = ' '

read ( s1(5:6), '(i2)' ) m
call i4_to_month_name ( m, month )

s2(i+1:) = month
i = i + len_trim ( month )

i = i + 1
s2(i:i) = ' '

s2(i+1:i+4) = s1(1:4)
i = i + 4

return
end subroutine s_to_date
subroutine s_to_dec ( s, itop, ibot, length )

  !*****************************************************************************80
  !
  !! S_TO_DEC reads a number from a string, returning a decimal result.
  !
  !  Discussion:
  !
  !    The integer may be in real format, for example '2.25'.  It
  !    returns ITOP and IBOT.  If the input number is an integer, ITOP
  !    equals that integer, and IBOT is 1.  But in the case of 2.25,
  !    the program would return ITOP = 225, IBOT = 100.
  !
  !    Legal input is
  !
  !          blanks,
  !       2  initial sign,
  !          blanks,
  !       3  whole number,
  !       4  decimal point,
  !       5  fraction,
  !       6  'E' or 'e' or 'D' or 'd', exponent marker,
  !       7  exponent sign,
  !       8  exponent,
  !          blanks
  !       9  comma or semicolon
  !      10end of information
  !
  !  Example:
  !
  !    S                 ITOP      IBOT     Length  Meaning
  !
  !    '1'                  1         0          1        1
  !    '     1   '          1         0          6        1
  !    '1A'                 1         0          1        1
  !    '12,34,56'          12         0          3       12
  !    '  34 7'            34         0          4       34
  !    '-1E2ABCD'          -1         2          4     -100
  !    '-1X2ABCD'          -1         0          2       -1
  !    ' 2E-1'              2        -1          5      0.2
  !    '23.45'           2345        -2          5    23.45
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 February 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading begins at position 1 and
  !    terminate when no more characters
  !    can be read to form a legal integer.  Blanks, commas,
  !    or other nonnumeric data will, in particular, cause
  !    the conversion to halt.
  !
  !    Output, integer ( kind = 4 ) ITOP, the integer read from the string,
  !    assuming that no negative exponents or fractional parts
  !    were used.  Otherwise, the 'integer' is ITOP/IBOT.
  !
  !    Output, integer ( kind = 4 ) IBOT, the integer divisor required to
  !    represent numbers which are in real format or have a
  !    negative exponent.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters used.
  !
implicit none

character              c
logical                ch_is_digit
integer   ( kind = 4 ) digit
integer   ( kind = 4 ) exponent
integer   ( kind = 4 ) exponent_sign
integer   ( kind = 4 ) ibot
integer   ( kind = 4 ) ihave
integer   ( kind = 4 ) iterm
integer   ( kind = 4 ) itop
integer   ( kind = 4 ) length
integer   ( kind = 4 ) mantissa_sign
character ( len = * )  s
logical                s_eqi

itop = 0
ibot = 0

if ( len ( s ) <= 0 ) then
  length = 0
  return
end if

length = -1
exponent_sign = 0
mantissa_sign = 1
exponent = 0
ihave = 1
iterm = 0
!
!  Consider the next character in the string.
!
do

  length = length + 1
  c = s(length+1:length+1)
  !
  !  Blank.
  !
  if ( c == ' ' ) then

     if ( ihave == 1 ) then

     else if ( ihave == 2 ) then

     else
        iterm = 1
     end if
     !
     !  Comma or semicolon.
     !
  else if ( c == ',' .or. c == ';' ) then

     if ( ihave /= 1 ) then
        iterm = 1
        ihave = 9
        length = length + 1
     end if
     !
     !  Minus sign.
     !
  else if ( c == '-' ) then

     if ( ihave == 1 ) then
        ihave = 2
        mantissa_sign = -1
     else if ( ihave == 6 ) then
        ihave = 7
        exponent_sign = -1
     else
        iterm = 1
     end if
     !
     !  Plus sign.
     !
  else if ( c == '+' ) then

     if ( ihave == 1 ) then
        ihave = 2
        mantissa_sign = +1
     else if ( ihave == 6 ) then
        ihave = 7
        exponent_sign = +1
     else
        iterm = 1
     end if
     !
     !  Decimal point.
     !
  else if ( c == '.' ) then

     if ( ihave < 4 ) then
        ihave = 4
     else
        iterm = 1
     end if
     !
     !  Exponent marker.
     !
  else if ( s_eqi ( c, 'E' ) .or. s_eqi ( c, 'D' ) ) then

     if ( ihave < 6 ) then
        ihave = 6
     else
        iterm = 1
     end if
     !
     !  Digit.
     !
  else if ( ch_is_digit ( c ) ) then

     if ( ihave <= 3 ) then
        ihave = 3
        call ch_to_digit ( c, digit )
        itop = 10 * itop + digit
     else if ( ihave <= 5 ) then
        ihave = 5
        call ch_to_digit ( c, digit )
        itop = 10 * itop + digit
        ibot = ibot - 1
     else if ( ihave <= 8 ) then
        ihave = 8
        call ch_to_digit ( c, digit )
        exponent = 10 * exponent + digit
     else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_DEC: Fatal error!'
        write ( *, '(a,i8)' ) '  IHAVE = ', ihave
        stop
     end if
     !
     !  Anything else is regarded as a terminator.
     !
  else
     iterm = 1
  end if

  if ( iterm == 1 ) then
     exit
  end if

  if ( len ( s ) <= length + 1 ) then
     length = len ( s )
     exit
  end if

end do
!
!  Number seems to have terminated.
!  Have we got a legal number?
!
if ( ihave == 1 ) then
  return
else if ( ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_TO_DEC - Serious error!'
  write ( *, '(a)' ) '  Illegal or nonnumeric input:'
  write ( *, '(a)' ) trim ( s )
  return
end if
!
!  Normalize.
!
if ( 0 < itop ) then

  do while ( mod ( itop, 10 ) == 0 )
     itop = itop / 10
     ibot = ibot + 1
  end do

end if
!
!  Consolidate the number in the form ITOP * 10**IBOT.
!
ibot = ibot + exponent_sign * exponent
itop = mantissa_sign * itop

if ( itop == 0 ) then
  ibot = 0
end if

return
end subroutine s_to_dec
subroutine s_to_ebcdic ( s )

  !*****************************************************************************80
  !
  !! S_TO_EBCDIC converts a character string from ASCII to EBCDIC.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S.  On input, the ASCII
  !    string, on output, the EBCDIC string.
  !
implicit none

character              ch_to_ebcdic
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len ( s )

do i = 1, s_length
  s(i:i) = ch_to_ebcdic ( s(i:i) )
end do

return
end subroutine s_to_ebcdic
subroutine s_to_format ( s, r, code, w, m )

  !*****************************************************************************80
  !
  !! S_TO_FORMAT reads a FORTRAN format from a string.
  !
  !  Discussion:
  !
  !    This routine will read as many characters as possible until it reaches
  !    the end of the string, or encounters a character which cannot be
  !    part of the format.  This routine is limited in its ability to
  !    recognize FORTRAN formats.  In particular, we are only expecting
  !    a single format specification, and cannot handle extra features
  !    such as 'ES' and 'EN' codes, '5X' spacing, and so on.
  !
  !    Legal input is:
  !
  !       0 nothing
  !       1 blanks
  !       2 optional '('
  !       3 blanks
  !       4 optional repeat factor R
  !       5 blanks
  !       6 CODE ( 'A', 'B', 'E', 'F', 'G', 'I', 'L', 'O', 'Z', '*' )
  !       7 blanks
  !       8 width W
  !       9 optional decimal point
  !      10 optional mantissa M
  !      11 blanks
  !      12 optional ')'
  !      13 blanks
  !
  !  Example:
  !
  !    S                 R   CODE   W    M
  !
  !    'I12              1   I      12   0
  !    'E8.0'            1   E       8   0
  !    'F10.5'           1   F      10   5
  !    '2G14.6'          2   G      14   6
  !    '*'               1   *      -1  -1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 November 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate at the end of the string, or when no more
  !    characters can be read.
  !
  !    Output, integer ( kind = 4 ) R, the repetition factor, which defaults to 1.
  !
  !    Output, character CODE, the format code.
  !
  !    Output, integer ( kind = 4 ) W, the field width.
  !
  !    Output, integer ( kind = 4 ) M, the mantissa width.
  !
implicit none

character              c
logical                ch_is_digit
logical                ch_is_format_code
character              code
integer   ( kind = 4 ) d
logical, parameter ::  debug = .true.
integer   ( kind = 4 ), parameter :: LEFT = 1
integer   ( kind = 4 ) m
integer   ( kind = 4 ) paren_sum
integer   ( kind = 4 ) pos
integer   ( kind = 4 ) r
integer   ( kind = 4 ), parameter :: RIGHT = -1
character ( len = * )  s
integer   ( kind = 4 ) s_length
integer   ( kind = 4 ) state
integer   ( kind = 4 ) w

state = 0
paren_sum = 0
pos = 0
s_length = len_trim ( s )

r = 0
w = 0
code = '?'
m = 0

do while ( pos < s_length )

  pos = pos + 1
  c = s(pos:pos)
  !
  !  BLANK character:
  !
  if ( c == ' ' ) then

     if ( state == 4 ) then
        state = 5
     else if ( state == 6 ) then
        state = 7
     else if ( state == 10 ) then
        state = 11
     else if ( state == 12 ) then
        state = 13
     end if
     !
     !  LEFT PAREN
     !
  else if ( c == '(' ) then

     if ( state < 2 ) then
        paren_sum = paren_sum + LEFT
     else
        if ( debug ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
           write ( *, '(a,i8)' ) '  Current state = ', state
           write ( *, '(a)' ) '  Input character = "' // c // '".'
        end if
        state = -1
        exit
     end if
     !
     !  DIGIT (R, F, or W)
     !
  else if ( ch_is_digit ( c ) ) then

     if ( state <= 3 ) then
        state = 4
        call ch_to_digit ( c, r )
     else if ( state == 4 ) then
        call ch_to_digit ( c, d )
        r = 10 * r + d
     else if ( state == 6 .or. state == 7 ) then
        if ( code == '*' ) then
           if ( debug ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
              write ( *, '(a,i8)' ) '  Current state = ', state
              write ( *, '(a,i8)' ) '  Current code = "' // code // '".'
              write ( *, '(a)' ) '  Input character = "' // c // '".'
           end if
           state = -1
           exit
        end if
        state = 8
        call ch_to_digit ( c, w )
     else if ( state == 8 ) then
        call ch_to_digit ( c, d )
        w = 10 * w + d
     else if ( state == 9 ) then
        state = 10
        call ch_to_digit ( c, m )
     else if ( state == 10 ) then
        call ch_to_digit ( c, d )
        m = 10 * m + d
     else
        if ( debug ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
           write ( *, '(a,i8)' ) '  Current state = ', state
           write ( *, '(a)' ) '  Input character = "' // c // '".'
        end if
        state = -1
        exit
     end if
     !
     !  DECIMAL POINT
     !
  else if ( c == '.' ) then

     if ( state == 8 ) then
        state = 9
     else
        if ( debug ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
           write ( *, '(a,i8)' ) '  Current state = ', state
           write ( *, '(a)' ) '  Input character = "' // c // '".'
        end if
        state = -1
        exit
     end if
     !
     !  RIGHT PAREN
     !
  else if ( c == ')' ) then

     paren_sum = paren_sum + RIGHT

     if ( paren_sum /= 0 ) then
        if ( debug ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
           write ( *, '(a,i8)' ) '  Current paren sum = ', paren_sum
           write ( *, '(a)' ) '  Input character = "' // c // '".'
        end if
        state = -1
        exit
     end if

     if ( state == 6 .and. code == '*' ) then
        state = 12
     else if ( 6 <= state ) then
        state = 12
     else
        if ( debug ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
           write ( *, '(a,i8)' ) '  Current state = ', state
           write ( *, '(a)' ) '  Input character = "' // c // '".'
        end if
        state = -1
        exit
     end if
     !
     !  Code
     !
  else if ( ch_is_format_code ( c ) ) then

     if ( state < 6 ) then
        state = 6
        code = c
     else
        if ( debug ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
           write ( *, '(a,i8)' ) '  Current state = ', state
           write ( *, '(a)' ) '  Input character = "' // c // '".'
        end if
        state = -1
        exit
     end if
     !
     !  Unexpected character
     !
  else

     if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
        write ( *, '(a,i8)' ) '  Current state = ', state
        write ( *, '(a)' ) '  Input character = "' // c // '".'
     end if
     state = -1
     exit

  end if

end do

if ( paren_sum /= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
  write ( *, '(a)' ) '  Parentheses mismatch.'
  stop
end if

if ( state < 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_TO_FORMAT - Fatal error!'
  write ( *, '(a)' ) '  Parsing error.'
  stop
end if

if ( r == 0 ) then
  r = 1
end if

return
end subroutine s_to_format
subroutine s_to_hex ( s, hex )

  !*****************************************************************************80
  !
  !! S_TO_HEX replaces a character string by a hexadecimal representation.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !    The string 'ABC' causes the hexadecimal string '414243' to be returned.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string of characters.
  !
  !    Output, character ( len = * ) HEX, the string of hex values.
  !
implicit none

character ( len = * )  hex
integer   ( kind = 4 ) i
integer   ( kind = 4 ) intval
integer   ( kind = 4 ) j
integer   ( kind = 4 ) ndo
integer   ( kind = 4 ) nhex
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )
nhex = len ( hex )

ndo = min ( nhex / 2, s_length )
hex = ' '

do i = 1, ndo

  j = 2 * i - 1
  intval = iachar ( s(i:i) )

  call i4_to_hex ( intval, hex(j:j+1) )

end do

return
end subroutine s_to_hex
subroutine s_to_i4 ( s, value, ierror, length )

  !*****************************************************************************80
  !
  !! S_TO_I4 reads an integer value from a string.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string to be examined.
  !
  !    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
  !    If the string is blank, then VALUE will be returned 0.
  !
  !    Output, integer ( kind = 4 ) IERROR, an error flag.
  !    0, no error.
  !    1, an error occurred.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters
  !    of S used to make the integer.
  !
implicit none

character              c
integer   ( kind = 4 ) i
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) isgn
integer   ( kind = 4 ) length
character ( len = * )  s
integer   ( kind = 4 ) state
character              :: TAB = achar ( 9 )
integer   ( kind = 4 ) value

value = 0
ierror = 0
length = 0

state = 0
isgn = 1

do i = 1, len_trim ( s )

  c = s(i:i)
  !
  !  STATE = 0, haven't read anything.
  !
  if ( state == 0 ) then

     if ( c == ' ' .or. c == TAB ) then

     else if ( c == '-' ) then
        state = 1
        isgn = -1
     else if ( c == '+' ) then
        state = 1
        isgn = +1
     else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
     else
        ierror = 1
        return
     end if
     !
     !  STATE = 1, have read the sign, expecting digits or spaces.
     !
  else if ( state == 1 ) then

     if ( c == ' ' .or. c == TAB ) then

     else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
     else
        ierror = 1
        return
     end if
     !
     !  STATE = 2, have read at least one digit, expecting more.
     !
  else if ( state == 2 ) then

     if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

     else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

     end if

  end if

end do
!
!  If we read all the characters in the string, see if we're OK.
!
if ( state == 2 ) then

  value = isgn * value
  ierror = 0
  length = len_trim ( s )

else

  value = 0
  ierror = 1
  length = 0

end if

return
end subroutine s_to_i4
subroutine s_to_i4vec ( s, n, i4vec, ierror )

  !*****************************************************************************80
  !
  !! S_TO_I4VEC reads an integer vector from a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 October 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be read.
  !
  !    Input, integer ( kind = 4 ) N, the number of values expected.
  !
  !    Output, integer ( kind = 4 ) I4VEC(N), the values read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    -K, could not read data for entries -K through N.
  !
implicit none

integer   ( kind = 4 ) n

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) i4vec(n)
integer   ( kind = 4 ) length
character ( len = * )  s

i = 0
ierror = 0
ilo = 1

do while ( i < n )

  i = i + 1

  call s_to_i4 ( s(ilo:), i4vec(i), ierror, length )

  if ( ierror /= 0 ) then
     ierror = -i
     exit
  end if

  ilo = ilo + length

end do

return
end subroutine s_to_i4vec
function s_to_l ( s )

  !*****************************************************************************80
  !
  !! S_TO_L reads a logical value from a string.
  !
  !  Discussion:
  !
  !    There are several ways of representing logical data that this routine
  !    recognizes:
  !
  !      False   True
  !      -----   ----
  !
  !      0       1
  !      F       T
  !      f       t
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 December 2010
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be read.
  !
  !    Output, logical S_TO_L, the logical value read from the string.
  !
implicit none

integer ( kind = 4 )  i
character ( len = * ) s
integer ( kind = 4 )  s_length
logical               s_to_l

s_length = len_trim ( s )

do i = 1, s_length
  if ( s(i:i) == '0' .or. s(i:i) == 'F' .or. s(i:i) == 'f' ) then
     s_to_l = .false.
     return
  else if ( s(i:i) == '1' .or. s(i:i) == 'T' .or. s(i:i) == 't' ) then
     s_to_l = .true.
     return
  end if
end do

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'S_TO_L - Fatal error!'
write ( *, '(a)' ) '  Input text did not contain logical data.'

stop
end function s_to_l
subroutine s_to_r4 ( s, r, ierror, length )

  !*****************************************************************************80
  !
  !! S_TO_R4 reads an R4 value from a string.
  !
  !  Discussion:
  !
  !    An "R4" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 4 )".
  !
  !    This routine will read as many characters as possible until it reaches
  !    the end of the string, or encounters a character which cannot be
  !    part of the real number.
  !
  !    Legal input is:
  !
  !       1 blanks,
  !       2 '+' or '-' sign,
  !       2.5 spaces
  !       3 integer part,
  !       4 decimal point,
  !       5 fraction part,
  !       6 'E' or 'e' or 'D' or 'd', exponent marker,
  !       7 exponent sign,
  !       8 exponent integer part,
  !       9 exponent decimal point,
  !      10 exponent fraction part,
  !      11 blanks,
  !      12 final comma or semicolon.
  !
  !    with most quantities optional.
  !
  !  Example:
  !
  !    S                 R
  !
  !    '1'               1.0
  !    '     1   '       1.0
  !    '1A'              1.0
  !    '12,34,56'        12.0
  !    '  34 7'          34.0
  !    '-1E2ABCD'        -100.0
  !    '-1X2ABCD'        -1.0
  !    ' 2E-1'           0.2
  !    '23.45'           23.45
  !    '-4.2E+2'         -420.0
  !    '17d2'            1700.0
  !    '-14e-2'         -0.14
  !    'e2'              100.0
  !    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 February 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate at the end of the string, or when no more
  !    characters can be read to form a legal real.  Blanks,
  !    commas, or other nonnumeric data will, in particular,
  !    cause the conversion to halt.
  !
  !    Output, real ( kind = 4 ) R, the real value that was read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    1, 2, 6 or 7, the input number was garbled.  The
  !    value of IERROR is the last type of input successfully
  !    read.  For instance, 1 means initial blanks, 2 means
  !    a plus or minus sign, and so on.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
  !    the string to form the number, including any terminating
  !    characters such as a trailing comma or blanks.
  !
implicit none

character              c
logical                ch_eqi
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) ihave
integer   ( kind = 4 ) isgn
integer   ( kind = 4 ) iterm
integer   ( kind = 4 ) jbot
integer   ( kind = 4 ) jsgn
integer   ( kind = 4 ) jtop
integer   ( kind = 4 ) length
integer   ( kind = 4 ) ndig
real      ( kind = 4 ) r
real      ( kind = 4 ) rbot
real      ( kind = 4 ) rexp
real      ( kind = 4 ) rtop
character ( len = * )  s
integer   ( kind = 4 ) s_length
character, parameter :: TAB = achar ( 9 )

s_length = len_trim ( s )
ierror = 0
r = 0.0E+00
length = -1
isgn = 1
rtop = 0.0E+00
rbot = 1.0E+00
jsgn = 1
jtop = 0
jbot = 1
ihave = 1
iterm = 0

do

  length = length + 1
  c = s(length+1:length+1)
  !
  !  Blank or TAB character.
  !
  if ( c == ' ' .or. c == TAB ) then

     if ( ihave == 2 ) then

     else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
     else if ( 1 < ihave ) then
        ihave = 11
     end if
     !
     !  Comma.
     !
  else if ( c == ',' .or. c == ';' ) then

     if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
     end if
     !
     !  Minus sign.
     !
  else if ( c == '-' ) then

     if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
     else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
     else
        iterm = 1
     end if
     !
     !  Plus sign.
     !
  else if ( c == '+' ) then

     if ( ihave == 1 ) then
        ihave = 2
     else if ( ihave == 6 ) then
        ihave = 7
     else
        iterm = 1
     end if
     !
     !  Decimal point.
     !
  else if ( c == '.' ) then

     if ( ihave < 4 ) then
        ihave = 4
     else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
     else
        iterm = 1
     end if
     !
     !  Exponent marker.
     !
  else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

     if ( ihave < 6 ) then
        ihave = 6
     else
        iterm = 1
     end if
     !
     !  Digit.
     !
  else if ( ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

     if ( ihave <= 2 ) then
        ihave = 3
     else if ( ihave == 4 ) then
        ihave = 5
     else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
     else if ( ihave == 9 ) then
        ihave = 10
     end if

     call ch_to_digit ( c, ndig )

     if ( ihave == 3 ) then
        rtop = 10.0E+00 * rtop + real ( ndig, kind = 4 )
     else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig, kind = 4 )
        rbot = 10.0E+00 * rbot
     else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
     else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
     end if
     !
     !  Anything else is regarded as a terminator.
     !
  else
     iterm = 1
  end if
  !
  !  If we haven't seen a terminator, and we haven't examined the
  !  entire string, go get the next character.
  !
  if ( iterm == 1 .or. s_length <= length + 1 ) then
     exit
  end if

end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
if ( iterm /= 1 .and. length + 1 == s_length ) then
  length = s_length
end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

  ierror = ihave

  return
end if
!
!  Number seems OK.  Form it.
!
if ( jtop == 0 ) then
  rexp = 1.0E+00
else

  if ( jbot == 1 ) then
     rexp = 10.0E+00**( jsgn * jtop )
  else
     rexp = jsgn * jtop
     rexp = rexp / jbot
     rexp = 10.0E+00**rexp
  end if

end if

r = isgn * rexp * rtop / rbot

return
end subroutine s_to_r4
subroutine s_to_r4vec ( s, n, r4vec, ierror )

  !*****************************************************************************80
  !
  !! S_TO_R4VEC reads an R4VEC from a string.
  !
  !  Discussion:
  !
  !    An R4VEC is a vector of real values, of type "real ( kind = 4 )".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 February 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be read.
  !
  !    Input, integer ( kind = 4 ) N, the number of values expected.
  !
  !    Output, real ( kind = 4 ) R4VEC(N), the values read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    -K, could not read data for entries -K through N.
  !
implicit none

integer   ( kind = 4 ) n

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) length
real      ( kind = 4 ) r4vec(n)
character ( len = * )  s

i = 0
ierror = 0
ilo = 1

do while ( i < n )

  i = i + 1

  call s_to_r4 ( s(ilo:), r4vec(i), ierror, length )

  if ( ierror /= 0 ) then
     ierror = -i
     exit
  end if

  ilo = ilo + length

end do

return
end subroutine s_to_r4vec
subroutine s_to_r8 ( s, dval, ierror, length )

  !*****************************************************************************80
  !
  !! S_TO_R8 reads an R8 value from a string.
  !
  !  Discussion:
  !
  !    An "R8" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 8 )".
  !
  !    The routine will read as many characters as possible until it reaches
  !    the end of the string, or encounters a character which cannot be
  !    part of the number.
  !
  !    Legal input is:
  !
  !       1 blanks,
  !       2 '+' or '-' sign,
  !       2.5 blanks
  !       3 integer part,
  !       4 decimal point,
  !       5 fraction part,
  !       6 'E' or 'e' or 'D' or 'd', exponent marker,
  !       7 exponent sign,
  !       8 exponent integer part,
  !       9 exponent decimal point,
  !      10 exponent fraction part,
  !      11 blanks,
  !      12 final comma or semicolon,
  !
  !    with most quantities optional.
  !
  !  Example:
  !
  !    S                 DVAL
  !
  !    '1'               1.0
  !    '     1   '       1.0
  !    '1A'              1.0
  !    '12,34,56'        12.0
  !    '  34 7'          34.0
  !    '-1E2ABCD'        -100.0
  !    '-1X2ABCD'        -1.0
  !    ' 2E-1'           0.2
  !    '23.45'           23.45
  !    '-4.2E+2'         -420.0
  !    '17d2'            1700.0
  !    '-14e-2'         -0.14
  !    'e2'              100.0
  !    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate at the end of the string, or when no more
  !    characters can be read to form a legal real.  Blanks,
  !    commas, or other nonnumeric data will, in particular,
  !    cause the conversion to halt.
  !
  !    Output, real ( kind = 8 ) DVAL, the value read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    1, 2, 6 or 7, the input number was garbled.  The
  !    value of IERROR is the last type of input successfully
  !    read.  For instance, 1 means initial blanks, 2 means
  !    a plus or minus sign, and so on.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read
  !    to form the number, including any terminating
  !    characters such as a trailing comma or blanks.
  !
implicit none

character              c
logical                ch_eqi
real      ( kind = 8 ) dval
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) ihave
integer   ( kind = 4 ) isgn
integer   ( kind = 4 ) iterm
integer   ( kind = 4 ) jbot
integer   ( kind = 4 ) jsgn
integer   ( kind = 4 ) jtop
integer   ( kind = 4 ) length
integer   ( kind = 4 ) ndig
real      ( kind = 8 ) rbot
real      ( kind = 8 ) rexp
real      ( kind = 8 ) rtop
character ( len = * )  s
integer   ( kind = 4 ) s_length
character           :: TAB = achar ( 9 )

s_length = len_trim ( s )

ierror = 0
dval = 0.0D+00
length = -1
isgn = 1
rtop = 0
rbot = 1
jsgn = 1
jtop = 0
jbot = 1
ihave = 1
iterm = 0

do

  length = length + 1

  if ( s_length < length + 1 ) then
     exit
  end if

  c = s(length+1:length+1)
  !
  !  Blank character.
  !
  if ( c == ' ' .or. c == TAB ) then

     if ( ihave == 2 ) then

     else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
     else if ( 1 < ihave ) then
        ihave = 11
     end if
     !
     !  Comma.
     !
  else if ( c == ',' .or. c == ';' ) then

     if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
     end if
     !
     !  Minus sign.
     !
  else if ( c == '-' ) then

     if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
     else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
     else
        iterm = 1
     end if
     !
     !  Plus sign.
     !
  else if ( c == '+' ) then

     if ( ihave == 1 ) then
        ihave = 2
     else if ( ihave == 6 ) then
        ihave = 7
     else
        iterm = 1
     end if
     !
     !  Decimal point.
     !
  else if ( c == '.' ) then

     if ( ihave < 4 ) then
        ihave = 4
     else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
     else
        iterm = 1
     end if
     !
     !  Scientific notation exponent marker.
     !
  else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

     if ( ihave < 6 ) then
        ihave = 6
     else
        iterm = 1
     end if
     !
     !  Digit.
     !
  else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

     if ( ihave <= 2 ) then
        ihave = 3
     else if ( ihave == 4 ) then
        ihave = 5
     else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
     else if ( ihave == 9 ) then
        ihave = 10
     end if

     call ch_to_digit ( c, ndig )

     if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
     else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
     else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
     else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
     end if
     !
     !  Anything else is regarded as a terminator.
     !
  else
     iterm = 1
  end if
  !
  !  If we haven't seen a terminator, and we haven't examined the
  !  entire string, go get the next character.
  !
  if ( iterm == 1 ) then
     exit
  end if

end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
if ( iterm /= 1 .and. length + 1 == s_length ) then
  length = s_length
end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
  ierror = ihave
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
  write ( *, '(a)' ) '  Illegal or nonnumeric input:'
  write ( *, '(a)' ) '    ' // trim ( s )
  return
end if
!
!  Number seems OK.  Form it.
!
if ( jtop == 0 ) then
  rexp = 1.0D+00
else
  if ( jbot == 1 ) then
     rexp = 10.0D+00 ** ( jsgn * jtop )
  else
     rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
          / real ( jbot, kind = 8 ) )
  end if
end if

dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

return
end subroutine s_to_r8
subroutine s_to_r8vec ( s, n, r8vec, ierror )

  !*****************************************************************************80
  !
  !! S_TO_R8VEC reads an R8VEC from a string.
  !
  !  Discussion:
  !
  !    An R8VEC is a vector of real values, of type "real ( kind = 8 )".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 January 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be read.
  !
  !    Input, integer ( kind = 4 ) N, the number of values expected.
  !
  !    Output, real ( kind = 8 ) R8VEC(N), the values read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    -K, could not read data for entries -K through N.
  !
implicit none

integer   ( kind = 4 ) n

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) lchar
real      ( kind = 8 ) r8vec(n)
character ( len = * )  s

i = 0
ierror = 0
ilo = 1

do while ( i < n )

  i = i + 1

  call s_to_r8 ( s(ilo:), r8vec(i), ierror, lchar )

  if ( ierror /= 0 ) then
     ierror = -i
     exit
  end if

  ilo = ilo + lchar

end do

return
end subroutine s_to_r8vec
subroutine s_to_rot13 ( s )

  !*****************************************************************************80
  !
  !! S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.
  !
  !  Discussion:
  !
  !    Two applications of the routine will return the original string.
  !
  !  Example:
  !
  !    Input:                      Output:
  !
  !    abcdefghijklmnopqrstuvwxyz  nopqrstuvwxyzabcdefghijklm
  !    Cher                        Pure
  !    James Thurston Howell       Wnzrf Guhefgba Ubjryy
  !    0123456789                  5678901234
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a string to be "rotated".
  !
implicit none

character              ch_to_rot13
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

do i = 1, s_length
  s(i:i) = ch_to_rot13 ( s(i:i) )
end do

return
end subroutine s_to_rot13
subroutine s_to_soundex ( s, code )

  !*****************************************************************************80
  !
  !! S_TO_SOUNDEX computes the Soundex code of a string.
  !
  !  Example:
  !
  !    Input:                      Output:
  !
  !    Ellery                      E460
  !    Euler                       E460
  !    Gauss                       G200
  !    Ghosh                       G200
  !    Heilbronn                   H416
  !    Hilbert                     H416
  !    Kant                        K530
  !    Knuth                       K530
  !    Ladd                        L300
  !    Lloyd                       L300
  !    Lissajous                   L222
  !    Lukasiewicz                 L222
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 May 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Donald Knuth,
  !    The Art of Computer Programming,
  !    Volume 3, Sorting and Searching,
  !    Second Edition,
  !    Addison Wesley, 1998,
  !    ISBN: 0201896850,
  !    LC: QA76.6.K64.
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string to be converted.
  !
  !    Output, character ( len = 4 ) CODE, the Soundex code for the string.
  !
implicit none

character              c
character              c_put
logical                ch_is_alpha
character              ch_s
character              ch_s_old
character ( len = 4 )  code
integer   ( kind = 4 ) get
integer   ( kind = 4 ) nget
integer   ( kind = 4 ) put
character ( len = * )  s

ch_s = '0'
code = ' '
nget = len_trim ( s )
get = 0
!
!  Try to fill position PUT of the code.
!
do put = 1, 4

  do

     if ( nget <= get ) then
        c_put = '0'
        exit
     end if

     get = get + 1
     c = s(get:get)
     call ch_cap ( c )

     if ( .not. ch_is_alpha ( c ) ) then
        cycle
     end if

     ch_s_old = ch_s

     call ch_to_soundex ( c, ch_s )

     if ( put == 1 ) then
        c_put = c
        exit
     else if ( ch_s /= ch_s_old .and. ch_s /= '0' ) then
        c_put = ch_s
        exit
     end if

  end do

  code(put:put) = c_put

end do

return
end subroutine s_to_soundex
subroutine s_to_w ( s, w, ierror, last )

  !*****************************************************************************80
  !
  !! S_TO_W reads the next blank-delimited word from a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 November 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string to be examined.
  !
  !    Output, character ( len = * ) W, the word that was read.
  !
  !    Output, integer ( kind = 4 ) IERROR, an error flag.
  !    0, no error.
  !    1, an error occurred.
  !
  !    Output, integer ( kind = 4 ) LAST, the last character of S used to make W.
  !
implicit none

character              c
integer   ( kind = 4 ) first
integer   ( kind = 4 ) i
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) last
character ( len = * )  s
integer   ( kind = 4 ) state
character ( len = * )  w

w = ' '
ierror = 0
state = 0
first = 0
last = 0
i = 0

do

  i = i + 1

  if ( len_trim ( s ) < i ) then

     if ( state == 0 ) then
        ierror = 1
        last = 0
     else
        last = i-1
        w = s(first:last)
     end if

     exit

  end if

  c = s(i:i)

  if ( state == 0 ) then

     if ( c /= ' ' ) then
        first = i
        state = 1
     end if

  else if ( state == 1 ) then

     if ( c == ' ' ) then
        last = i - 1
        w = s(first:last)
        exit
     end if

  end if

end do

return
end subroutine s_to_w
subroutine s_token_equal ( s, set, nset, iset )

  !*****************************************************************************80
  !
  !! S_TOKEN_EQUAL checks whether a string is equal to any of a set of strings.
  !
  !  Discussion:
  !
  !    The comparison is case-insensitive.
  !
  !    Trailing blanks in S and the elements of SET are ignored.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to check.
  !
  !    Input, character ( len = * ) SET(NSET), the set of strings.
  !
  !    Input, integer ( kind = 4 ) NSET, the number of elements of SET.
  !
  !    Output, integer ( kind = 4 ) ISET, equals 0 if no element of SET
  !    equals S.  If ISET is nonzero, then SET(ISET) equals
  !    S, disregarding case.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) iset
integer   ( kind = 4 ) nset
character ( len = * )  s
logical                s_eqi
character ( len = * )  set(*)

iset = 0
do i = 1, nset

  if ( s_eqi ( s, set(i) ) ) then
     iset = i
     return
  end if

end do

return
end subroutine s_token_equal
subroutine s_token_match ( s, token_num, token, match )

  !*****************************************************************************80
  !
  !! S_TOKEN_MATCH matches the beginning of a string and a set of tokens.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = 'TOMMYGUN'
  !      TOKEN = 'TOM', 'ZEBRA', 'TOMMY', 'TOMMYKNOCKER'
  !
  !    Output:
  !
  !      MATCH = 3
  !
  !  Discussion:
  !
  !    The longest possible match is taken.
  !    Matching is done without regard to case or trailing blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Input, integer ( kind = 4 ) TOKEN_NUM, the number of tokens to be compared.
  !
  !    Input, character ( len = * ) TOKEN(TOKEN_NUM), the tokens.
  !
  !    Output, integer ( kind = 4 ) MATCH, the index of the (longest)
  !    token that matched the string, or 0 if no match was found.
  !
implicit none

integer   ( kind = 4 ) token_num

integer   ( kind = 4 ) match
integer   ( kind = 4 ) match_length
character ( len = * )  s
logical                s_eqi
integer   ( kind = 4 ) s_length
integer   ( kind = 4 ) token_i
integer   ( kind = 4 ) token_length
character ( len = * )  token(token_num)

match = 0
match_length = 0

s_length = len_trim ( s )

do token_i = 1, token_num

  token_length = len_trim ( token ( token_i ) )

  if ( match_length < token_length ) then

     if ( token_length <= s_length ) then

        if ( s_eqi ( s(1:token_length), token(token_i)(1:token_length) ) ) then
           match_length = token_length
           match = token_i
        end if

     end if

  end if

end do

return
end subroutine s_token_match
subroutine s_trim_zeros ( s )

  !*****************************************************************************80
  !
  !! S_TRIM_ZEROS removes trailing zeros from a string.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = '1401.072500'
  !
  !    Output:
  !
  !      S = '1401.0725'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be operated on.
  !
implicit none

character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

do while ( 0 < s_length .and. s(s_length:s_length) == '0' )
  s(s_length:s_length) = ' '
  s_length = s_length - 1
end do

return
end subroutine s_trim_zeros
subroutine s_u2b ( s )

  !*****************************************************************************80
  !
  !! S_U2B replaces underscores by blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 December 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be
  !    transformed.
  !
implicit none

integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

do i = 1, s_length
  if ( s(i:i) == '_' ) then
     s(i:i) = ' '
  end if
end do

return
end subroutine s_u2b
subroutine s_word_append ( s, w, done )

  !*****************************************************************************80
  !
  !! S_WORD_APPEND appends a word to a string.
  !
  !  Discussion:
  !
  !    A blank space will separate the word from the text already
  !    in the line.
  !
  !    The routine warns the user if the word will not fit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 December 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a line of text.
  !    On input, the current string.  On output, the current string
  !    with the integer value appended.
  !
  !    Input, character ( len = * ) W, a word to be appended.
  !    Trailing blanks in the word are ignored.
  !
  !    Output, logical DONE, is FALSE if there was not enough room
  !    to append the word.
  !
implicit none

logical                done
integer   ( kind = 4 ) lens
integer   ( kind = 4 ) lents
integer   ( kind = 4 ) lenw
integer   ( kind = 4 ) next
character ( len = * )  s
character ( len = * )  w

done = .false.
lens = len ( s )
lents = len_trim ( s )

lenw = len_trim ( w )

if ( lents == 0 ) then
  if ( lens < lenw ) then
     done = .true.
     return
  end if
else
  if ( lens < lents + 1 + lenw ) then
     done = .true.
     return
  end if
end if

if ( lents == 0 ) then
  next = 1
else
  next = lents + 1
  s(next:next) = ' '
  next = next + 1
end if

s(next:next+lenw-1) = w(1:lenw)

return
end subroutine s_word_append
subroutine s_word_cap ( s )

  !*****************************************************************************80
  !
  !! S_WORD_CAP capitalizes the first character of each word in a string.
  !
  !  Example:
  !
  !    Input:
  !
  !      S = 'it is time to turn the page.'
  !
  !    Output:
  !
  !      S = 'It Is Time To Turn The Page.'
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 August 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !
implicit none

logical                blank
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length

s_length = len_trim ( s )

blank = .true.

do i = 1, s_length

  if ( blank ) then
     call ch_cap ( s(i:i) )
  else
     call ch_low ( s(i:i) )
  end if

  blank = ( s(i:i) == ' ' )

end do

return
end subroutine s_word_cap
subroutine s_word_count ( s, word_num )

  !*****************************************************************************80
  !
  !! S_WORD_COUNT counts the number of "words" in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 February 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, integer ( kind = 4 ) WORD_NUM, the number of "words" in the
  !    string.  Words are presumed to be separated by one or more blanks.
  !
implicit none

logical                blank
integer   ( kind = 4 ) i
character ( len = * )  s
integer   ( kind = 4 ) s_length
character, parameter :: TAB = achar ( 9 )
integer   ( kind = 4 ) word_num

word_num = 0
s_length = len ( s )

if ( s_length <= 0 ) then
  return
end if

blank = .true.

do i = 1, s_length

  if ( s(i:i) == ' ' .or. s(i:i) == TAB ) then
     blank = .true.
  else if ( blank ) then
     word_num = word_num + 1
     blank = .false.
  end if

end do

return
end subroutine s_word_count
subroutine s_word_extract_first ( s, w )

  !*****************************************************************************80
  !
  !! S_WORD_EXTRACT_FIRST extracts the first word from a string.
  !
  !  Discussion:
  !
  !    A "word" is a string of characters terminated by a blank or
  !    the end of the string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 January 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string.  On output, the first
  !    word has been removed, and the remaining string has been shifted left.
  !
  !    Output, character ( len = * ) W, the leading word of the string.
  !
implicit none

integer   ( kind = 4 ) get1
integer   ( kind = 4 ) get2
character ( len = * )  s
integer   ( kind = 4 ) s_length
character ( len = * )  w

w = ' '

s_length = len_trim ( s )

if ( s_length < 1 ) then
  return
end if
!
!  Find the first nonblank.
!
get1 = 0

do

  get1 = get1 + 1

  if ( s_length < get1 ) then
     return
  end if

  if ( s(get1:get1) /= ' ' ) then
     exit
  end if

end do
!
!  Look for the last contiguous nonblank.
!
get2 = get1

do

  if ( s_length <= get2 ) then
     exit
  end if

  if ( s(get2+1:get2+1) == ' ' ) then
     exit
  end if

  get2 = get2 + 1

end do
!
!  Copy the word.
!
w = s(get1:get2)
!
!  Shift the string.
!
s(1:get2) = ' '
s = adjustl ( s )

return
end subroutine s_word_extract_first
subroutine s_word_find ( s, iword, word, nchar )

  !*****************************************************************************80
  !
  !! S_WORD_FIND finds the word of a given index in a string.
  !
  !  Discussion:
  !
  !    A "word" is any string of nonblank characters, separated from other
  !    words by one or more blanks or TABS.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 April 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, integer ( kind = 4 ) IWORD, the index of the word to be
  !    searched for.  If IWORD is positive, then the IWORD-th
  !    word is sought.  If IWORD is zero or negative, then
  !    assuming that the string has N words in it, the
  !    N+IWORD-th word will be sought.
  !
  !    Output, character ( len = * ) WORD, the IWORD-th word of the
  !    string, or ' ' if the WORD could not be found.
  !
  !    Output, integer ( kind = 4 ) NCHAR, the number of characters in WORD,
  !    or 0 if the word could not be found.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) iblank
integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) iword
integer   ( kind = 4 ) jhi
integer   ( kind = 4 ) jlo
integer   ( kind = 4 ) jword
integer   ( kind = 4 ) kword
integer   ( kind = 4 ) nchar
character ( len = * )  s
integer   ( kind = 4 ) s_len
character, parameter :: TAB = achar ( 9 )
character ( len = * ) word

ilo = 0
ihi = 0
s_len = len_trim ( s )

if ( s_len <= 0 ) then
  return
end if

if ( 0 < iword ) then

  if ( s(1:1) == ' ' .or. s(1:1) == TAB ) then
     iblank = 1
     jword = 0
     jlo = 0
     jhi = 0
  else
     iblank = 0
     jword = 1
     jlo = 1
     jhi = 1
  end if

  i = 1

  do

     i = i + 1

     if ( s_len < i ) then

        if ( jword == iword ) then
           ilo = jlo
           ihi = s_len
           nchar = s_len + 1 - jlo
           word = s(ilo:ihi)
        else
           ilo = 0
           ihi = 0
           nchar = 0
           word = ' '
        end if

        return

     end if

     if ( ( s(i:i) == ' ' .or. s(i:i) == TAB ) .and. iblank == 0 ) then

        jhi = i - 1
        iblank = 1
        if ( jword == iword ) then
           ilo = jlo
           ihi = jhi
           nchar = jhi + 1 - jlo
           word = s(ilo:ihi)
           return
        end if

     else if ( s(i:i) /= ' ' .and. s(i:i) /= TAB .and. iblank == 1 ) then

        jlo = i
        jword = jword + 1
        iblank = 0

     end if

  end do

else

  iblank = 0
  kword = 1 - iword
  jword = 1
  jlo = s_len
  jhi = s_len
  i = s_len

  do

     i = i - 1

     if ( i <= 0 ) then

        if ( jword == kword ) then
           ilo = 1
           ihi = jhi
           nchar = jhi
           word = s(ilo:ihi)
        else
           ilo = 0
           ihi = 0
           nchar = 0
           word = ' '
        end if

        return

     end if

     if ( ( s(i:i) == ' ' .or. s == TAB ) .and. iblank == 0 ) then

        jlo = i + 1
        iblank = 1

        if ( jword == kword ) then
           ilo = jlo
           ihi = jhi
           nchar = jhi + 1 - jlo
           word = s(ilo:ihi)
           return
        end if

     else if ( s(i:i) /= ' ' .and. s(i:i) /= TAB .and. iblank == 1 ) then

        jhi = i
        jword = jword + 1
        iblank = 0

     end if

  end do

end if

return
end subroutine s_word_find
subroutine s_word_index ( s, indx, ilo, ihi )

  !*****************************************************************************80
  !
  !! S_WORD_INDEX finds the word of a given index in a string.
  !
  !  Discussion:
  !
  !    The routine returns in ILO and IHI the beginning and end of the INDX-th
  !    word, or 0 and 0 if there is no INDX-th word.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S is the string of words to be analyzed.
  !
  !    Input, integer ( kind = 4 ) INDX is the index of the desired token.
  !
  !    Output, integer ( kind = 4 ) ILO is the index of the first character
  !    of the INDX-th word, or 0 if there was no INDX-th word.
  !
  !    Output, integer ( kind = 4 ) IHI is the index of the last character
  !    of the INDX-th word, or 0 if there was no INDX-th word.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ) indx
character ( len = * )  s

ihi = 0
ilo = 0

do i = 1, indx

  call word_next ( s, ilo, ihi )

  if ( ilo == 0 ) then
     return
  end if

end do

return
end subroutine s_word_index
subroutine s_word_next ( s, word, done )

  !*****************************************************************************80
  !
  !! S_WORD_NEXT "reads" words from a string, one at a time.
  !
  !  Special cases:
  !
  !    The following characters are considered to be a single word,
  !    whether surrounded by spaces or not:
  !
  !      " ( ) { } [ ]
  !
  !    Also, if there is a trailing comma on the word, it is stripped off.
  !    This is to facilitate the reading of lists.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 May 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string, presumably containing words
  !    separated by spaces.
  !
  !    Output, character ( len = * ) WORD.
  !    If DONE is FALSE, then WORD contains the "next" word read.
  !    If DONE is TRUE, then WORD is blank, because there was no more to read.
  !
  !    Input/output, logical DONE.
  !    On input with a fresh string, set DONE to TRUE.
  !    On output, the routine sets DONE:
  !      FALSE if another word was read,
  !      TRUE if no more words could be read.
  !
implicit none

logical                done
integer   ( kind = 4 ) ilo
integer   ( kind = 4 ), save :: next = 1
character ( len = * )  s
integer   ( kind = 4 ), save :: s_length = 0
character, parameter :: TAB = achar ( 9 )
character ( len = * )  word
!
!  We "remember" S_LENGTH and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
if ( done ) then

  next = 1
  done = .false.
  s_length = len_trim ( s )

  if ( s_length <= 0 ) then
     done = .true.
     word = ' '
     return
  end if

end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
do

  if ( s_length < ilo ) then
     word = ' '
     done = .true.
     next = s_length + 1
     return
  end if
  !
  !  If the current character is blank, skip to the next one.
  !
  if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
     exit
  end if

  ilo = ilo + 1

end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
if ( s(ilo:ilo) == '"' .or. &
    s(ilo:ilo) == '(' .or. &
    s(ilo:ilo) == ')' .or. &
    s(ilo:ilo) == '{' .or. &
    s(ilo:ilo) == '}' .or. &
    s(ilo:ilo) == '[' .or. &
    s(ilo:ilo) == ']' ) then

  word = s(ilo:ilo)
  next = ilo + 1
  return

end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
next = ilo + 1

do while ( next <= s_length )

  if ( s(next:next) == ' ' ) then
     exit
  else if ( s(next:next) == TAB ) then
     exit
  else if ( s(next:next) == '"' ) then
     exit
  else if ( s(next:next) == '(' ) then
     exit
  else if ( s(next:next) == ')' ) then
     exit
  else if ( s(next:next) == '{' ) then
     exit
  else if ( s(next:next) == '}' ) then
     exit
  else if ( s(next:next) == '[' ) then
     exit
  else if ( s(next:next) == ']' ) then
     exit
  end if

  next = next + 1

end do
!
!  Ignore a trailing comma.
!
if ( s(next-1:next-1) == ',' ) then
  word = s(ilo:next-2)
else
  word = s(ilo:next-1)
end if

return
end subroutine s_word_next
subroutine s_word_permute ( s1, n, perm, s2 )

  !*****************************************************************************80
  !
  !! S_WORD_PERMUTE permutes the words in a string.
  !
  !  Discussion:
  !
  !    A word is a blank-delimited sequence of characters.
  !
  !    The string is assumed to contain N "words".  If more words are
  !    in the string, their position is not affected.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 January 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, a line of text.
  !
  !    Input, integer ( kind = 4 ) N, the number of words to permute.
  !
  !    Input, integer ( kind = 4 ) PERM(N), the permutation.  PERM(1) is the new
  !    location of the item whose original location was 1.
  !
  !    Output, character ( len = * ) S2, a copy of S1 with the
  !    first N words permuted.
  !
implicit none

integer ( kind = 4 ) n

character              c1
character              c2
integer   ( kind = 4 ) index1
integer   ( kind = 4 ) index2
integer   ( kind = 4 ) perm(n)
integer   ( kind = 4 ) perm_inv(n)
character ( len = * )  s1
integer   ( kind = 4 ) s1_length
integer   ( kind = 4 ) s1_pos
integer   ( kind = 4 ) s1_word_index(n)
integer   ( kind = 4 ) s1_word_length(n)
character ( len = * )  s2
integer   ( kind = 4 ) s2_pos
integer   ( kind = 4 ) word_length
!
!  Set up word position and length vectors.
!
s1_length = len ( s1 )

s1_word_length(1:n) = 0
s1_word_index(1:n) = 0

index1 = 0
c2 = ' '

do s1_pos = 1, s1_length

  c1 = c2
  c2 = s1(s1_pos:s1_pos)

  if ( s1_pos == 1 .or. ( c1 /= ' ' .and. c2 == ' ' ) ) then

     if ( n <= index1 ) then
        exit
     end if

     index1 = index1 + 1

     s1_word_index(index1) = s1_pos

  end if

  s1_word_length(index1) = s1_word_length(index1) + 1

end do
!
!  Invert the permutation.
!
call perm_inverse3 ( n, perm, perm_inv )
!
!  Copy S1 into S2, so we get any trailing information.
!
call s_copy ( s1, s2 )
!
!  Copy the first N words of S1 into S2 in permuted order.
!
s2_pos = 1

do index2 = 1, n

  index1 = perm_inv(index2)

  s1_pos = s1_word_index(index1)

  word_length = s1_word_length(index1)

  s2(s2_pos:s2_pos+word_length-1) = s1(s1_pos:s1_pos+word_length-1)

  s2_pos = s2_pos + word_length

end do

return
end subroutine s_word_permute
function s32_to_i4 ( s32 )

  !*****************************************************************************80
  !
  !! S32_TO_I4 returns an I4 equivalent to a 32 character string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = 32 ) S32, the character value.
  !
  !    Output, integer ( kind = 4 ) S32_TO_I4, a corresponding integer value.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) intval
character ( len = 32 ) s32
integer   ( kind = 4 ) s32_to_i4
character ( len = 32 ) scopy

scopy = s32

if ( scopy(1:1) == '1' ) then

  do i = 2, 32

     if ( scopy(i:i) == '0' ) then
        scopy(i:i) = '1'
     else
        scopy(i:i) = '0'
     end if

  end do

end if

intval = 0

do i = 2, 32

  intval = 2 * intval

  if ( scopy(i:i) == '1' ) then
     intval = intval + 1
  end if

end do

if ( scopy(1:1) == '1' ) then
  intval = -intval
end if

s32_to_i4 = intval

return
end function s32_to_i4
function s32_to_r4 ( s32 )

  !*****************************************************************************80
  !
  !! S32_TO_R4 converts a 32-character variable into an R4.
  !
  !  Discussion:
  !
  !    An "R4" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 4 )".
  !
  !    The first bit is 1 for a negative real, or 0 for a
  !    positive real.  Bits 2 through 9 are the exponent.  Bits 10
  !    through 32 are used for a normalized representation of the
  !    mantissa. Since it is assumed that normalization means the first
  !    digit of the mantissa is 1, this 1 is in fact not stored.
  !
  !    The special case of 0 is represented by all 0 bits.
  !
  !    It is believed that this method corresponds to the format used
  !    in VMS FORTRAN for reals.
  !
  !    Because of the limits on the mantissa, many Cray numbers are not
  !    representable at all by this method.  These numbers are very big
  !    or very small in magnitude.  Other numbers will simply be
  !    represented with less accuracy.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = 32 ) S32, the character variable to be decoded.
  !
  !    Output, real ( kind = 4 ) RCHAR32, the corresponding real value.
  !
implicit none

integer   ( kind = 4 ) i
integer   ( kind = 4 ) iexp
integer   ( kind = 4 ) j
integer   ( kind = 4 ) mant
character ( len = 32 ) s32
real      ( kind = 4 ) s32_to_r4
real      ( kind = 4 ) sgn
!
!  Read sign bit.
!
if ( s32(1:1) == '1' ) then
  sgn = -1.0E+00
else
  sgn = 1.0E+00
end if
!
!  Construct exponent from bits 2 through 9, subtract 128.
!
iexp = 0

do i = 2, 9

  if ( s32(i:i) == '0' ) then
     j = 0
  else
     j = 1
  end if

  iexp = 2 * iexp + j

end do

if ( iexp == 0 ) then
  s32_to_r4 = 0.0E+00
  return
end if

iexp = iexp - 128
!
!  Read mantissa from positions 10 through 32.
!  Note that, unless exponent equals 0, the most significant bit is
!  assumed to be 1 and hence is not stored.
!
mant = 1

do i = 10, 32
  mant = 2 * mant
  if ( s32(i:i) == '1' ) then
     mant = mant + 1
  end if
end do

s32_to_r4 = sgn * mant * ( 2.0E+00 ** ( iexp - 23 ) )

return
end function s32_to_r4
subroutine sef_to_b4_ieee  ( s, e, f, word )

  !*****************************************************************************80
  !
  !! SEF_TO_B4_IEEE converts SEF information to a 4 byte IEEE real word.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 November 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) S, the sign bit:
  !    0, if R is nonnegative;
  !    1, if R is negative.
  !
  !    Input, integer ( kind = 4 ) E, the exponent, base 2.
  !    Normally, -127 < E <= 127.
  !    If E = 128, then the data is interpreted as NaN, Inf, or -Inf.
  !    If -127 < E <= 127, the data is a normalized value.
  !    If E < -127, then the data is a denormalized value.
  !
  !    Input, integer ( kind = 4 ) F, the mantissa.
  !
  !    Output, integer ( kind = 4 ) WORD, the real number stored in IEEE format.
  !
implicit none

integer ( kind = 4 ) e
integer ( kind = 4 ) e2
integer ( kind = 4 ) f
integer ( kind = 4 ), parameter :: f_max = 2**24
integer ( kind = 4 ), parameter :: f_min = 2**23
integer ( kind = 4 ) f2
integer ( kind = 4 ) s
integer ( kind = 4 ) s2
integer ( kind = 4 ) word

s2 = s
e2 = e
f2 = f
!
!  Handle +Inf and -Inf.
!
if ( f /= 0 .and. e == 128 ) then
  e2 = e2 + 127
  f2 = 2**23 - 1
  call mvbits ( s2, 0,  1, word, 31 )
  call mvbits ( e2, 0,  8, word, 23 )
  call mvbits ( f2, 0, 23, word,  0 )
  return
end if
!
!  Handle NaN.
!
if ( f == 0 .and. e == 128 ) then
  e2 = e2 + 127
  f2 = 0
  call mvbits ( s2, 0,  1, word, 31 )
  call mvbits ( e2, 0,  8, word, 23 )
  call mvbits ( f2, 0, 23, word,  0 )
  return
end if
!
!  Handle +0 and -0.
!
if ( f == 0 ) then
  e2 = 0
  call mvbits ( s2, 0,  1, word, 31 )
  call mvbits ( e2, 0,  8, word, 23 )
  call mvbits ( f2, 0, 23, word,  0 )
  return
end if
!
!  Normalize.
!
if ( f < 0 ) then
  s2 = 1 - s2
  f2 = -f2
end if

e2 = e2 + 127 + 23

do while ( f_max <= f2 )
  f2 = f2 / 2
  e2 = e2 + 1
end do

do while ( f2 < f_min )
  f2 = f2 * 2
  e2 = e2 - 1
end do
!
!  The biased exponent cannot be negative.
!  Shift it up to zero, and reduce F2.
!
do while ( e2 < 0 .and. f2 /= 0 )
  e2 = e2 + 1
  f2 = f2 / 2
end do
!
!  Normalized values drop the leading 1.
!
if ( 0 < e2 ) then
  call mvbits ( s2, 0,  1, word, 31 )
  call mvbits ( e2, 0,  8, word, 23 )
  f2 = f2 - f_min
  call mvbits ( f2, 0, 23, word,  0 )
  !
  !  Denormalized values have a biased exponent of 0.
  !
else
  call mvbits ( s2, 0,  1, word, 31 )
  call mvbits ( e2, 0,  8, word, 23 )
  call mvbits ( f2, 0, 23, word,  0 )
end if

return
end subroutine sef_to_b4_ieee
subroutine sef_to_r4 ( s, e, f, r )

  !*****************************************************************************80
  !
  !! SEF_TO_R4 converts SEF information to an R4 = S * 2.0**E * F.
  !
  !  Discussion:
  !
  !    An "R4" value is simply a real number to be stored as a
  !    variable of type "real ( kind = 4 )".
  !
  !    Assuming no arithmetic problems, in fact, this equality should be
  !    exact, that is, S, E and F should exactly express the value
  !    as stored on the computer.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 November 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) S, the sign bit:
  !    0, if R is nonnegative;
  !    1, if R is negative.
  !
  !    Input, integer ( kind = 4 ) E, the exponent, base 2.
  !
  !    Input, integer ( kind = 4 ) F, the mantissa.
  !
  !    Output, real ( kind = 4 ) R, the real number.
  !
implicit none

integer ( kind = 4 ) e
integer ( kind = 4 ) f
integer ( kind = 4 ) i
real    ( kind = 4 ) r
integer ( kind = 4 ) s

if ( f == 0 ) then
  r = 0.0E+00
  return
end if

if ( s == 0 ) then
  r = 1.0E+00
else if ( s == 1 ) then
  r = -1.0E+00
else
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SEF_TO_R4 - Fatal error!'
  write ( *, '(a,i8)' ) '  Illegal input value of S = ', s
  stop
end if

r = r * real ( f, kind = 4 )

if ( 0 < e ) then
  do i = 1, e
     r = r * 2.0E+00
  end do
else if ( e < 0 ) then
  do i = 1, -e
     r = r / 2.0E+00
  end do
end if

return
end subroutine sef_to_r4
subroutine sort_heap_external ( n, indx, i, j, isgn )

  !*****************************************************************************80
  !
  !! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
  !
  !  Discussion:
  !
  !    The actual list of data is not passed to the routine.  Hence this
  !    routine may be used to sort integers, reals, numbers, names,
  !    dates, shoe sizes, and so on.  After each call, the routine asks
  !    the user to compare or interchange two items, until a special
  !    return value signals that the sorting is completed.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 February 2004
  !
  !  Author:
  !
  !    Original version by Albert Nijenhuis, Herbert Wilf.
  !    FORTRAN90 version by John Burkardt
  !
  !  Reference:
  !
  !    Albert Nijenhuis, Herbert Wilf,
  !    Combinatorial Algorithms for Computers and Calculators,
  !    Academic Press, 1978, second edition,
  !    ISBN 0-12-519260-6,
  !    LC: QA164.N54.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of items to be sorted.
  !
  !    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
  !
  !    The user must set INDX to 0 before the first call.
  !    Thereafter, the user should not change the value of INDX until
  !    the sorting is done.
  !
  !    On return, if INDX is
  !
  !      greater than 0,
  !      * interchange items I and J;
  !      * call again.
  !
  !      less than 0,
  !      * compare items I and J;
  !      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
  !      * call again.
  !
  !      equal to 0, the sorting is done.
  !
  !    Output, integer ( kind = 4 ) I, J, the indices of two items.
  !    On return with INDX positive, elements I and J should be interchanged.
  !    On return with INDX negative, elements I and J should be compared, and
  !    the result reported in ISGN on the next call.
  !
  !    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
  !    I and J. (Used only when the previous call returned INDX less than 0).
  !    ISGN <= 0 means I is less than or equal to J;
  !    0 <= ISGN means I is greater than or equal to J.
  !
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ), save :: i_save = 0
integer ( kind = 4 ) indx
integer ( kind = 4 ) isgn
integer ( kind = 4 ) j
integer ( kind = 4 ), save :: j_save = 0
integer ( kind = 4 ), save :: k = 0
integer ( kind = 4 ), save :: k1 = 0
integer ( kind = 4 ) n
integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
if ( indx == 0 ) then

  i_save = 0
  j_save = 0
  k = n / 2
  k1 = k
  n1 = n
  !
  !  INDX < 0: The user is returning the results of a comparison.
  !
else if ( indx < 0 ) then

  if ( indx == -2 ) then

     if ( isgn < 0 ) then
        i_save = i_save + 1
     end if

     j_save = k1
     k1 = i_save
     indx = -1
     i = i_save
     j = j_save
     return

  end if

  if ( 0 < isgn ) then
     indx = 2
     i = i_save
     j = j_save
     return
  end if

  if ( k <= 1 ) then

     if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
     else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
     end if

     i = i_save
     j = j_save
     return

  end if

  k = k - 1
  k1 = k
  !
  !  0 < INDX, the user was asked to make an interchange.
  !
else if ( indx == 1 ) then

  k1 = k

end if

do

  i_save = 2 * k1

  if ( i_save == n1 ) then
     j_save = k1
     k1 = i_save
     indx = -1
     i = i_save
     j = j_save
     return
  else if ( i_save <= n1 ) then
     j_save = i_save + 1
     indx = -2
     i = i_save
     j = j_save
     return
  end if

  if ( k <= 1 ) then
     exit
  end if

  k = k - 1
  k1 = k

end do

if ( n1 == 1 ) then
  i_save = 0
  j_save = 0
  indx = 0
  i = i_save
  j = j_save
else
  i_save = n1
  n1 = n1 - 1
  j_save = 1
  indx = 1
  i = i_save
  j = j_save
end if

return
end subroutine sort_heap_external
function state_id ( state )

  !*****************************************************************************80
  !
  !! STATE_ID returns the 2 letter Postal Code for one of the 50 states.
  !
  !  Discussion:
  !
  !    The states are listed in order of their admission to the union.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 April 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) STATE, the index of a state.
  !
  !    Output, character ( len = 2 ) STATE_ID, the 2 letter code.
  !
implicit none

character ( len = 2 ), parameter, dimension ( 50 ) :: id = (/ &
    'DE', 'PA', 'NJ', 'GA', 'CT', &
    'MA', 'MD', 'SC', 'NH', 'VA', &
    'NY', 'NC', 'RI', 'VT', 'KY', &
    'TN', 'OH', 'LA', 'IN', 'MS', &
    'IL', 'AL', 'ME', 'MO', 'AR', &
    'MI', 'FL', 'TX', 'IA', 'WI', &
    'CA', 'MN', 'OR', 'KS', 'WV', &
    'NV', 'NE', 'CO', 'ND', 'SD', &
    'MT', 'WA', 'ID', 'WY', 'UT', &
    'OK', 'NM', 'AZ', 'AL', 'HI' /)
integer ( kind = 4 ) state
character ( len = 2 ) state_id

if ( state < 1 ) then
  state_id = '??'
else if ( state <= 50 ) then
  state_id = id(state)
else
  state_id = '??'
end if

return
end function state_id
function state_name ( state )

  !*****************************************************************************80
  !
  !! STATE_NAME returns the name of one of the 50 states.
  !
  !  Discussion:
  !
  !    The states are listed in order of their admission to the union.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 April 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) STATE, the index of a state.
  !
  !    Output, character ( len = 14 ) STATE_NAME, the name of the state.
  !
implicit none

character ( len = 14 ), parameter, dimension ( 50 ) :: name = (/ &
    'Delaware      ', &
    'Pennsylvania  ', &
    'New Jersey    ', &
    'Georgia       ', &
    'Connecticut   ', &
    'Massachusetts ', &
    'Maryland      ', &
    'South Carolina', &
    'New Hampshire ', &
    'Virginia      ', &
    'New York      ', &
    'North Carolina', &
    'Rhode Island  ', &
    'Vermont       ', &
    'Kentucky      ', &
    'Tennessee     ', &
    'Ohio          ', &
    'Louisiana     ', &
    'Indiana       ', &
    'Missippi      ', &
    'Illinois      ', &
    'Alabama       ', &
    'Maine         ', &
    'Missouri      ', &
    'Arkansas      ', &
    'Michigan      ', &
    'Florida       ', &
    'Texas         ', &
    'Iowa          ', &
    'Wisconsin     ', &
    'California    ', &
    'Minnesota     ', &
    'Oregon        ', &
    'Kansas        ', &
    'West Virginia ', &
    'Nevada        ', &
    'Nebraska      ', &
    'Colorado      ', &
    'North Dakota  ', &
    'South Dakota  ', &
    'Montana       ', &
    'Washington    ', &
    'Idaho         ', &
    'Wyoming       ', &
    'Utah          ', &
    'Oklahoma      ', &
    'New Mexico    ', &
    'Arizona       ', &
    'Alaska        ', &
    'Hawaii        ' /)
integer ( kind = 4 ) state
character ( len = 14 ) state_name

if ( state < 1 ) then
  state_name = '??????????????'
else if ( state <= 50 ) then
  state_name = name(state)
else
  state_name = '??????????????'
end if

return
end function state_name
subroutine svec_lab ( n, nuniq, svec, ident )

  !*****************************************************************************80
  !
  !! SVEC_LAB makes an index array for an array of (repeated) strings.
  !
  !  Discussion:
  !
  !    The routine is given an array of strings.  It assigns an integer
  !    to each unique string, and returns an equivalent array of
  !    these values.
  !
  !    Note that blank strings are treated specially.  Any blank
  !    string gets an identifier of 0.  Blank strings are not
  !    counted in the value of NUNIQ.
  !
  !  Example:
  !
  !    SVEC    IDENT
  !
  !    ALPHA       1
  !    ALPHA      -1
  !    BETA        2
  !    ALPHA      -1
  !    BETA       -2
  !    GAMMA       3
  !    ALPHA      -1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 February 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries.
  !
  !    Output, integer ( kind = 4 ) NUNIQ, the number of unique nonblank entries.
  !
  !    Input, character ( len = * ) SVEC(N), the list of strings.
  !
  !    Output, integer ( kind = 4 ) IDENT(N), the identifiers assigned to the
  !    strings.  If SVEC(I) is blank, then IDENT(I) is 0.
  !    Otherwise, if SVEC(I) is the first occurrence of a
  !    given string, then it is assigned a positive identifier.
  !    If SVEC(I) is a later occurrence of a string, then
  !    it is assigned a negative identifier, whose absolute
  !    value is the identifier of the first occurrence.
  !
implicit none

integer   ( kind = 4 ) n

integer   ( kind = 4 ) i
integer   ( kind = 4 ) ident(n)
integer   ( kind = 4 ) j
integer   ( kind = 4 ) match
integer   ( kind = 4 ) nuniq
character ( len = * )  svec(n)

nuniq = 0

do i = 1, n

  if ( svec(i) == ' ' ) then

     ident(i) = 0

  else

     match = 0

     do j = 1, i-1
        if ( 0 < ident(j) ) then
           if ( svec(j) == svec(i) ) then
              ident(i) = -ident(j)
              match = j
              exit
           end if
        end if
     end do

     if ( match == 0 ) then
        nuniq = nuniq + 1
        ident(i) = nuniq
     end if

  end if

end do

return
end subroutine svec_lab
subroutine svec_merge_a ( na, a, nb, b, nc, c )

  !*****************************************************************************80
  !
  !! SVEC_MERGE_A merges two ascending sorted string arrays.
  !
  !  Discussion:
  !
  !    The elements of A and B should be sorted in ascending order.
  !
  !    The elements in the output array C will be in ascending order, and unique.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NA, the dimension of A.
  !
  !    Input, character ( len = * ) A(NA), the first sorted array.
  !
  !    Input, integer ( kind = 4 ) NB, the dimension of B.
  !
  !    Input, character ( len = * ) B(NB), the second sorted array.
  !
  !    Output, integer ( kind = 4 ) NC, the number of elements in the output
  !    array.  Note that C should usually be dimensioned at least NA+NB in the
  !    calling routine.
  !
  !    Output, character ( len = * ) C(NC), the merged unique sorted array.
  !
implicit none

integer   ( kind = 4 ) na
integer   ( kind = 4 ) nb

character ( len = * )  a(na)
character ( len = * )  b(nb)
character ( len = * )  c(na+nb)
integer   ( kind = 4 ) j
integer   ( kind = 4 ) ja
integer   ( kind = 4 ) jb
integer   ( kind = 4 ) na2
integer   ( kind = 4 ) nb2
integer   ( kind = 4 ) nc

na2 = na
nb2 = nb

ja = 0
jb = 0
nc = 0

do
   !
   !  If we've used up all the entries of A, stick the rest of B on the end.
   !
  if ( na2 <= ja ) then

     do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 ) then
           nc = nc + 1
           c(nc) = b(jb)
        else if ( llt ( c(nc), b(jb) ) ) then
           nc = nc + 1
           c(nc) = b(jb)
        end if
     end do

     exit
     !
     !  If we've used up all the entries of B, stick the rest of A on the end.
     !
  else if ( nb2 <= jb ) then

     do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 ) then
           nc = nc + 1
           c(nc) = a(ja)
        else if ( llt ( c(nc), a(ja) ) ) then
           nc = nc + 1
           c(nc) = a(ja)
        end if
     end do

     exit
     !
     !  Otherwise, if the next entry of A is smaller, that's our candidate.
     !
  else if ( lle ( a(ja+1), b(jb+1) ) ) then

     ja = ja + 1
     if ( nc == 0 ) then
        nc = nc + 1
        c(nc) = a(ja)
     else if ( llt ( c(nc), a(ja) ) ) then
        nc = nc + 1
        c(nc) = a(ja)
     end if
     !
     !  ...or if the next entry of B is the smaller, consider that.
     !
  else

     jb = jb + 1
     if ( nc == 0 ) then
        nc = nc + 1
        c(nc) = b(jb)
     else if ( llt ( c(nc), b(jb) ) ) then
        nc = nc + 1
        c(nc) = b(jb)
     end if
  end if

end do

return
end subroutine svec_merge_a
subroutine svec_permute ( n, a, p )

  !*****************************************************************************80
  !
  !! SVEC_PERMUTE permutes a string vector in place.
  !
  !  Example:
  !
  !    Input:
  !
  !      N = 5
  !      P = (  3,     2,     4,       2,      1 )
  !      A = ( 'ONE', 'TWO', 'THREE', 'FOUR', 'FIVE' )
  !
  !    Output:
  !
  !      A    = ( 'FIVE', 'FOUR', 'ONE', 'THREE', 'TWO' ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of objects.
  !
  !    Input/output, character ( len = * ) A(N), the array to be permuted.
  !
  !    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
  !    that the I-th element of the output array should be the J-th
  !    element of the input array.  P must be a legal permutation
  !    of the integers from 1 to N, otherwise the algorithm will
  !    fail catastrophically.
  !
implicit none

integer ( kind = 4 ) n

character ( len = * )   a(n)
character ( len = 255 ) a_temp
integer   ( kind = 4 ) ierror
integer   ( kind = 4 ) iget
integer   ( kind = 4 ) iput
integer   ( kind = 4 ) istart
integer   ( kind = 4 ) p(n)

call perm_check ( n, p, ierror )

if ( ierror /= 0 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVEC_PERMUTE - Fatal error!'
  write ( *, '(a)' ) '  The input array does not represent'
  write ( *, '(a)' ) '  a proper permutation.  In particular, the'
  write ( *, '(a,i8)' ) '  array is missing the value ', ierror
  stop
end if
!
!  Search for the next element of the permutation that has not been used.
!
do istart = 1, n

  if ( p(istart) < 0 ) then

     cycle

  else if ( p(istart) == istart ) then

     p(istart) = -p(istart)
     cycle

  else

     a_temp = a(istart)
     iget = istart
     !
     !  Copy the new value into the vacated entry.
     !
     do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'SVEC_PERMUTE - Fatal error!'
           write ( *, '(a)' ) '  A permutation index is out of range.'
           write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
           stop
        end if

        if ( iget == istart ) then
           a(iput) = a_temp
           exit
        end if

        a(iput) = a(iget)

     end do

  end if

end do
!
!  Restore the signs of the entries.
!
p(1:n) = -p(1:n)

return
end subroutine svec_permute
subroutine svec_reverse ( n, a )

  !*****************************************************************************80
  !
  !! SVEC_REVERSE reverses the elements of a string vector.
  !
  !  Example:
  !
  !    Input:
  !
  !      N = 4,
  !      A = ( 'Bob', 'Carol', 'Ted', 'Alice' ).
  !
  !    Output:
  !
  !      A = ( 'Alice', 'Ted', 'Carol', 'Bob' ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in the array.
  !
  !    Input/output, character ( len = * ) A(N), the array to be reversed.
  !
implicit none

integer   ( kind = 4 )  n

character ( len = * )   a(n)
character ( len = 255 ) a_temp
integer   ( kind = 4 )  i

do i = 1, n/2
  a_temp   = a(i)
  a(i)     = a(n+1-i)
  a(n+1-i) = a_temp
end do

return
end subroutine svec_reverse
subroutine svec_search_binary_a ( n, a, b, indx )

  !*****************************************************************************80
  !
  !! SVEC_SEARCH_BINARY_A searches an ascending sorted string vector.
  !
  !  Discussion:
  !
  !    Binary search is used.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Donald Kreher, Douglas Simpson,
  !    Algorithm 1.9,
  !    Combinatorial Algorithms,
  !    CRC Press, 1998, page 26.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of elements in the vector.
  !
  !    Input, character ( len = * ) A(N), the array to be searched.  A must
  !    be sorted in increasing order.
  !
  !    Input, character ( len = * ) B, the value to be searched for.
  !
  !    Output, integer ( kind = 4 ) INDX, the result of the search.
  !    0, B does not occur in A.
  !    I, A(I) = B.
  !
implicit none

integer   ( kind = 4 ) n

character ( len = * )  a(n)
character ( len = * )  b
integer   ( kind = 4 ) high
integer   ( kind = 4 ) indx
integer   ( kind = 4 ) low
integer   ( kind = 4 ) mid

indx = 0

low = 1
high = n

do while ( low <= high )

  mid = ( low + high ) / 2

  if ( a(mid) == b ) then
     indx = mid
     exit
  else if ( llt ( a(mid), b ) ) then
     low = mid + 1
  else if ( lgt ( a(mid), b ) ) then
     high = mid - 1
  end if

end do

return
end subroutine svec_search_binary_a
subroutine svec_sort_heap_a ( n, a )

  !*****************************************************************************80
  !
  !! SVEC_SORT_HEAP_A ascending sorts an SVEC using heap sort.
  !
  !  Discussion:
  !
  !    The ASCII collating sequence is used.  This means
  !      A < B < C < .... < Y < Z < a < b < .... < z.
  !    Numbers and other symbols may also occur, and will be sorted according to
  !    the ASCII ordering.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of strings
  !
  !    Input/output, character ( len = * ) A(N);
  !    On input, an array of strings to be sorted.
  !    On output, the sorted array.
  !
implicit none

integer   ( kind = 4 ) n

character ( len = * )  a(n)
integer   ( kind = 4 ) i
integer   ( kind = 4 ) indx
integer   ( kind = 4 ) isgn
integer   ( kind = 4 ) j
!
!  Do the sorting using the external heap sort routine.
!
i = 0
indx = 0
isgn = 0
j = 0

do

  call sort_heap_external ( n, indx, i, j, isgn )

  if ( 0 < indx ) then

     call s_swap ( a(i), a(j) )

  else if ( indx < 0 ) then

     if ( lle ( a(i), a(j) ) ) then
        isgn = -1
     else
        isgn = +1
     end if

  else if ( indx == 0 ) then

     exit

  end if

end do

return
end subroutine svec_sort_heap_a
subroutine svec_sort_heap_a_index ( n, sarray, indx )

  !*****************************************************************************80
  !
  !! SVEC_SORT_HEAP_A_INDEX: case-sensitive indexed heap sort of an SVEC.
  !
  !  Discussion:
  !
  !    The sorting is not actually carried out.
  !    Rather an index array is created which defines the sorting.
  !    This array may be used to sort or index the array, or to sort or
  !    index related arrays keyed on the original array.
  !
  !    The ASCII collating sequence is used, and case is significant.
  !    This means
  !
  !      A < B < C < .... < Y < Z < a < b < .... < z.
  !
  !    Numbers and other symbols may also occur, and will be sorted according to
  !    the ASCII ordering.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 July 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in SARRAY.
  !
  !    Input, character ( len = * ) SARRAY(N), an array to be sorted.
  !
  !    Output, integer ( kind = 4 ) INDX(N), contains the sort index.  The
  !    I-th element of the sorted array is SARRAY ( INDX(I) ).
  !
implicit none

integer   ( kind = 4 )  n

integer   ( kind = 4 )  i
integer   ( kind = 4 )  indx(n)
integer   ( kind = 4 )  indxt
integer   ( kind = 4 )  ir
integer   ( kind = 4 )  j
integer   ( kind = 4 )  l
character ( len = * )   sarray(n)
character ( len = 255 ) string

do i = 1, n
  indx(i) = i
end do

l = n / 2 + 1
ir = n

do

  if ( 1 < l ) then

     l = l - 1
     indxt = indx(l)
     string = sarray(indxt)

  else

     indxt = indx(ir)
     string = sarray(indxt)
     indx(ir) = indx(1)
     ir = ir - 1

     if ( ir == 1 ) then
        indx(1) = indxt
        return
     end if

  end if

  i = l
  j = l + l

  do while ( j <= ir )

     if ( j < ir ) then
        if ( llt ( sarray ( indx(j) ), sarray ( indx(j+1) ) ) ) then
           j = j + 1
        end if
     end if

     if ( llt ( string, sarray ( indx(j) ) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
     else
        j = ir + 1
     end if

  end do

  indx(i) = indxt

end do

return
end subroutine svec_sort_heap_a_index
subroutine svec_sorted_unique ( n, a, unique_num )

  !*****************************************************************************80
  !
  !! SVEC_SORTED_UNIQUE: number of unique entries in a sorted SVEC.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of elements in the array.
  !
  !    Input/output, character ( len = * ) A(N).
  !    On input, the sorted list of strings.
  !    On output, the unique elements, in sorted order.
  !
  !    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
  !    in the array.
  !
implicit none

integer   ( kind = 4 ) n

character ( len = * )  a(n)
integer   ( kind = 4 ) itest
integer   ( kind = 4 ) unique_num

if ( n <= 0 ) then
  unique_num = 0
  return
end if

unique_num = 1

do itest = 2, n

  if ( a(itest) /= a(unique_num) ) then
     unique_num = unique_num + 1
     a(unique_num) = a(itest)
  end if

end do

return
end subroutine svec_sorted_unique
subroutine sveci_search_binary_a ( n, a, b, indx )

  !*****************************************************************************80
  !
  !! SVECI_SEARCH_BINARY_A: search ascending sorted implicitly capitalized SVEC
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Donald Kreher, Douglas Simpson,
  !    Algorithm 1.9,
  !    Combinatorial Algorithms,
  !    CRC Press, 1998, page 26.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of elements in the vector.
  !
  !    Input, character ( len = * ) A(N), the array to be searched.  A must
  !    be sorted in increasing order.
  !
  !    Input, character ( len = * ) B, the value to be searched for.
  !
  !    Output, integer ( kind = 4 ) INDX, the result of the search.
  !    0, B does not occur in A.
  !    I, A(I) = B, ignoring capitalization.
  !
implicit none

integer   ( kind = 4 ) n

character ( len = * )  a(n)
character ( len = * )  b
integer   ( kind = 4 ) high
integer   ( kind = 4 ) indx
integer   ( kind = 4 ) low
integer   ( kind = 4 ) mid
logical                s_eqi
logical                s_gti
logical                s_lti

indx = 0

low = 1
high = n

do while ( low <= high )

  mid = ( low + high ) / 2

  if ( s_eqi ( a(mid), b ) ) then
     indx = mid
     exit
  else if ( s_lti ( a(mid), b ) ) then
     low = mid + 1
  else if ( s_gti ( a(mid), b ) ) then
     high = mid - 1
  end if

end do

return
end subroutine sveci_search_binary_a
subroutine sveci_sort_heap_a ( n, sarray )

  !*****************************************************************************80
  !
  !! SVECI_SORT_HEAP_A heap sorts an SVEC of implicitly capitalized strings.
  !
  !  Discussion:
  !
  !    The characters in an implicitly capitalized string are treated as
  !    though they had been capitalized.  Thus, the letters 'a' and 'A'
  !    are considered equal, both 'a' and 'A' precede 'B', and
  !    'Fox' and 'fOx' are considered equal.
  !
  !    The ASCII collating sequence is used, except that all
  !    alphabetic characters are treated as though they were uppercase.
  !
  !    This means
  !
  !      A = a < B = b < C = c < .... < Y = y < Z = z.
  !
  !    Numbers and other symbols may also occur, and will be sorted
  !    according to the ASCII ordering.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in SARRAY.
  !
  !    Input/output, character ( len = * ) SARRAY(N), the array to be sorted.
  !
implicit none

integer   ( kind = 4 ) n

integer   ( kind = 4 ) l
integer   ( kind = 4 ) l1
integer   ( kind = 4 ) m
integer   ( kind = 4 ) n1
logical                s_gei
logical                s_lti
character ( len = * )  sarray(n)
character ( len = 255 ) s

n1 = n
l = n / 2
s = sarray(l)
l1 = l

do

  m = 2 * l1

  if ( m <= n1 ) then

     if ( m < n1 ) then
        if ( s_gei ( sarray(m+1), sarray(m) ) ) then
           m = m + 1
        end if
     end if

     if ( s_lti ( s, sarray(m) ) ) then
        sarray(l1) = sarray(m)
        l1 = m
        cycle
     end if

  end if

  sarray(l1) = s

  if ( 1 < l ) then
     l = l - 1
     s = sarray(l)
     l1 = l
     cycle
  end if

  if ( n1 < 2 ) then
     exit
  end if

  s = sarray(n1)
  sarray(n1) = sarray(1)

  n1 = n1 - 1
  l1 = l

end do

return
end subroutine sveci_sort_heap_a
subroutine sveci_sort_heap_a_index ( n, sarray, indx )

  !*****************************************************************************80
  !
  !! SVECI_SORT_HEAP_A_INDEX index heap sorts an SVECI.
  !
  !  Discussion:
  !
  !    The sorting is not actually carried out,
  !    but rather an index vector is returned, which defines the
  !    sorting.  This index vector may be used to sort the array, or
  !    to sort related arrays keyed on the first one.
  !
  !    The ASCII collating sequence is used, except that all
  !    alphabetic characters are treated as though they were uppercase.
  !
  !    This means
  !
  !      A = a < B = b < C = c < .... < Y = y < Z = z.
  !
  !    Numbers and other symbols may also occur, and will be sorted according to
  !    the ASCII ordering.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of entries in SARRAY.
  !
  !    Input, character ( len = * ) SARRAY(N), an array to be sorted.
  !
  !    Output, integer ( kind = 4 ) INDX(N), contains the sort index.  The
  !    I-th element of the sorted array is SARRAY ( INDX(I) ).
  !
implicit none

integer   ( kind = 4 )  n

integer   ( kind = 4 )  i
integer   ( kind = 4 )  indx(n)
integer   ( kind = 4 )  indxt
integer   ( kind = 4 )  ir
integer   ( kind = 4 )  j
integer   ( kind = 4 )  l
logical                 s_lti
character ( len = * )   sarray(n)
character ( len = 255 ) s

do i = 1, n
  indx(i) = i
end do

l = n / 2 + 1
ir = n

do

  if ( 1 < l ) then

     l = l - 1
     indxt = indx(l)
     s = sarray(indxt)

  else

     indxt = indx(ir)
     s = sarray(indxt)
     indx(ir) = indx(1)
     ir = ir - 1

     if ( ir == 1 ) then
        indx(1) = indxt
        return
     end if

  end if

  i = l
  j = l + l

  do while ( j <= ir )

     if ( j < ir ) then
        if ( s_lti ( sarray ( indx(j) ), sarray ( indx(j+1) ) ) ) then
           j = j + 1
        end if
     end if

     if ( s_lti ( s, sarray ( indx(j) ) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
     else
        j = ir + 1
     end if

  end do

  indx(i) = indxt

end do

return
end subroutine sveci_sort_heap_a_index
subroutine sym_to_ch ( sym, c, ihi )

  !*****************************************************************************80
  !
  !! SYM_TO_CH returns the character represented by a symbol.
  !
  !  Discussion:
  !
  !    Instead of ICHAR, we now use the IACHAR function, which
  !    guarantees the ASCII collating sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 April 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) SYM is a string containing printable symbols.
  !
  !    Output, character C, is the ASCII character represented by the
  !    first symbol in SYM.
  !
  !    Output, integer ( kind = 4 ) IHI, C is represented by SYM(1:IHI).
  !    IHI = 0 if there was a problem.
  !
implicit none

character              c
integer   ( kind = 4 ) ialt
integer   ( kind = 4 ) ichr
integer   ( kind = 4 ) ictl
integer   ( kind = 4 ) ihi
logical                s_eqi
character ( len = * )  sym
integer   ( kind = 4 ) sym_length

c = ' '
sym_length = len_trim ( sym )

if ( sym_length <= 0 ) then
  c = ' '
  ihi = 0
  return
end if

ialt = 0
ictl = 0
ihi = 1
!
!  Could it be an ALT character?
!
if ( sym(ihi:ihi) == '!' .and. ihi < sym_length ) then
  ialt = 1
  ihi = ihi + 1
end if
!
!  Could it be a control character?
!
if ( sym(ihi:ihi) == '^' .and. ihi < sym_length ) then
  ictl = 1
  ihi = ihi + 1
end if
!
!  Could it be a DEL character?
!
ichr = iachar ( sym(ihi:ihi) )

if ( ihi+2 <= sym_length ) then
  if ( s_eqi ( sym(ihi:ihi+2), 'DEL' ) ) then
     ichr = 127
     ihi = ihi + 2
  end if
end if
!
!  Could it be an SP character?
!
if ( ihi + 1 <= sym_length ) then
  if ( s_eqi ( sym(ihi:ihi+1), 'SP' ) ) then
     ichr = 32
     ihi = ihi + 1
  end if
end if
!
!  Interpret the character.
!
if ( ialt == 1 ) then
  ichr = ichr + 128
end if

if ( ictl == 1 ) then
  ichr = ichr - 64
end if

c = achar ( ichr )

return
end subroutine sym_to_ch


subroutine token_expand ( s, tokens )

  !*****************************************************************************80
  !
  !! TOKEN_EXPAND makes sure certain tokens have spaces surrounding them.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be examined.
  !
  !    Input, character ( len = * ) TOKENS, a string of characters.  Every
  !    occurrence of a character from TOKENS in S must be
  !    preceded and followed by a blank space, except if the occurrence
  !    is in the first or last positions of S, in which a
  !    preceding or trailing blank space is implicit.
  !
implicit none

character               c1
character               c2
character               c3
integer   ( kind = 4 )  i
integer   ( kind = 4 )  put
integer   ( kind = 4 )  j
integer   ( kind = 4 )  lenc
integer   ( kind = 4 )  lent
character ( len = * )   s
character ( len = 255 ) s2
character ( len = * )   tokens

lenc = len_trim ( s )
lent = len_trim ( tokens )
s2 = ' '
put = 0
c2 = ' '
c3 = s(1:1)

do i = 1, lenc

  c1 = c2
  c2 = c3

  if ( i < lenc ) then
     c3 = s(i+1:i+1)
  else
     c3 = ' '
  end if

  do j = 1, lent

     if ( c2 == tokens(j:j) ) then
        if ( c1 /= ' ' ) then
           put = put + 1
           if ( put <= 255 ) then
              s2(put:put) = ' '
           end if
        end if
     end if

  end do

  put = put + 1

  if ( put <= 255 ) then
     s2(put:put) = c2
  end if

  do j = 1, lent

     if ( c2 == tokens(j:j) ) then
        if ( c3 /= ' ' ) then
           put = put + 1
           if ( put <= 255 ) then
              s2(put:put) = ' '
           end if
        end if
     end if

  end do

end do

s = s2

return
end subroutine token_expand
subroutine token_extract ( s, token_num, token, match )

  !*****************************************************************************80
  !
  !! TOKEN_EXTRACT "extracts" a token from the beginning of a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 November 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S; on input, a string from
  !    whose beginning a token is to be extracted.  On output,
  !    the token, if found, has been removed.
  !
  !    Input, integer ( kind = 4 ) TOKEN_NUM, the number of tokens to be
  !    compared.
  !
  !    Input, character ( len = * ) TOKEN(TOKEN_NUM), the tokens.
  !
  !    Output, integer ( kind = 4 ) MATCH, the index of the (longest) token
  !    that matched the string, or 0 if no match was found.
  !
implicit none

integer   ( kind = 4 ) token_num

integer   ( kind = 4 ) left
integer   ( kind = 4 ) match
character ( len = * )  s
character ( len = * )  token(token_num)

call s_token_match ( s, token_num, token, match )

if ( match /= 0 ) then
  left = len_trim ( token(match) )
  call s_shift_left ( s, left )
end if

return
end subroutine token_extract
subroutine token_index ( s, indx, ilo, ihi )

  !*****************************************************************************80
  !
  !! TOKEN_INDEX finds the N-th FORTRAN variable name in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S is the string of words to be analyzed.
  !
  !    Input, integer ( kind = 4 ) INDX is the index of the desired token.
  !
  !    Output, integer ( kind = 4 ) ILO is the index of the first character
  !    of the INDX-th token, or 0 if there was no INDX-th token.
  !
  !    Output, integer ( kind = 4 ) IHI is the index of the last character
  !    of the INDX-th token, or 0 if there was no INDX-th token.
  !
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ) ihi
integer ( kind = 4 ) ilo
integer ( kind = 4 ) indx
character ( len = * ) s

ihi = 0
ilo = 0

do i = 1, indx

  call token_next ( s, ilo, ihi)

  if ( ilo == 0 ) then
     return
  end if

end do

return
end subroutine token_index
subroutine token_next ( s, ilo, ihi )

  !*****************************************************************************80
  !
  !! TOKEN_NEXT finds the next FORTRAN variable name in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S is the string of words to be analyzed.
  !
  !    Output, integer ( kind = 4 ) ILO is the location of the first character
  !    of the next word, or 0 if there was no next word.
  !
  !    Input/output, integer ( kind = 4 ) IHI.
  !    On input, IHI is taken to be the LAST character of the
  !    PREVIOUS word, or 0 if the first word is sought.
  !
  !    On output, IHI is the index of the last character of
  !    the next word, or 0 if there was no next word.
  !
implicit none

integer   ( kind = 4 ) ihi
integer   ( kind = 4 ) ilo
character ( len = * )  s
integer   ( kind = 4 ) s_len
logical                s_only_alphab
logical                s_only_digitb

s_len = len_trim ( s )

ilo = ihi

if ( ilo < 0 ) then
  ilo = 0
end if
!
!  Find ILO, the index of the next alphabetic character.
!
do

  ilo = ilo + 1

  if ( s_len < ilo ) then
     ilo = 0
     ihi = 0
     return
  end if

  if ( s_only_alphab ( s(ilo:ilo) ) ) then
     exit
  end if

end do
!
!  Find the index of the next character which is neither
!  alphabetic nor numeric.
!
ihi = ilo

do

  ihi = ihi + 1

  if ( s_len < ihi ) then
     ihi = s_len
     return
  end if

  if ( .not. ( s_only_alphab ( s(ihi:ihi) ) ) .and. &
       .not. ( s_only_digitb ( s(ihi:ihi) ) ) ) then
     exit
  end if

end do

ihi = ihi - 1

return
end subroutine token_next
subroutine word_bounds ( line, word_num, word_start, word_end )

  !*****************************************************************************80
  !
  !! WORD_BOUNDS returns the start and end of each word in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 October 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) LINE, a string containing words
  !    separated by spaces.
  !
  !    Input, integer ( kind = 4 ) WORD_NUM, the number of words in the line.
  !
  !    Output, integer ( kind = 4 ) WORD_START(WORD_NUM), WORD_END(WORD_NUM),
  !    the locations in LINE of the beginning and end of each word.
  !
implicit none

integer ( kind = 4 ) word_num

logical blank
character c
logical, parameter :: debug = .true.
integer ( kind = 4 ) i
character ( len = * ) line
integer ( kind = 4 ) line_len
integer ( kind = 4 ) w
integer ( kind = 4 ) word_end(word_num)
integer ( kind = 4 ) word_start(word_num)

i = 0
w = 0
blank = .true.

line_len = len_trim ( line )

do i = 1, line_len + 1

  if ( i <= line_len ) then
     c = line(i:i)
  else
     c = ' '
  end if

  if ( c == ' ' ) then

     if ( .not. blank ) then
        word_end(w) = i-1
        if ( w == word_num ) then
           exit
        end if
     end if

     blank = .true.

  else

     if ( blank ) then
        w = w + 1
        word_start(w) = i
     end if

     blank = .false.

  end if

end do

if ( w /= word_num ) then
  if ( debug ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'WORD_BOUNDS - Warning:'
     write ( *, '(a)' ) '  Found fewer words than requested.'
  end if
end if

return
end subroutine word_bounds
subroutine word_last_read ( s, word )

  !*****************************************************************************80
  !
  !! WORD_LAST_READ returns the last word from a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 April 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string containing words separated
  !    by spaces.
  !
  !    Output, character ( len = * ) WORD, the last word.
  !
implicit none

integer ( kind = 4 ) first
integer ( kind = 4 ) last
character ( len = * ) s
character ( len = * ) word

last = len_trim ( s )

if ( last <= 0 ) then
  word = ' '
  return
end if

first = last

do

  if ( first <= 1 ) then
     exit
  end if

  if ( s(first-1:first-1) == ' ' ) then
     exit
  end if

  first = first - 1

end do

word = s(first:last)

return
end subroutine word_last_read
subroutine word_next ( s, ilo, ihi )

  !*****************************************************************************80
  !
  !! WORD_NEXT finds the next (blank separated) word in a string.
  !
  !  Discussion:
  !
  !    This routine is usually used repetitively on a fixed string.  On each
  !    call, it accepts IHI, the index of the last character of the
  !    previous word extracted from the string.
  !
  !    It then computes ILO and IHI, the first and last characters of
  !    the next word in the string.
  !
  !    It is assumed that words are separated by one or more spaces.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 April 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string of words to be analyzed.
  !
  !    Output, integer ( kind = 4 ) ILO is the location of the first character
  !    of the next word, or 0 if there was no next word.
  !
  !    Input/output, integer ( kind = 4 ) IHI.
  !    On input, IHI is taken to be the LAST character of the
  !    PREVIOUS word, or 0 if the first word is sought.
  !    On output, IHI is the index of the last character of
  !    the next word, or 0 if there was no next word.
  !
implicit none

integer ( kind = 4 ) ihi
integer ( kind = 4 ) ilo
character ( len = * ) s
integer ( kind = 4 ) s_len

s_len = len_trim ( s )
!
!  Find ILO, the index of the first nonblank character after
!  (the old value of) IHI.
!
if ( ihi < 0 ) then
  ilo = 0
else
  ilo = ihi
end if

do

  ilo = ilo + 1

  if ( s_len < ilo ) then
     ilo = 0
     ihi = 0
     return
  end if

  if ( s(ilo:ilo) /= ' ') then
     exit
  end if

end do
!
!  Find IHI, the index of the next blank character, or end of line.
!
ihi = ilo

do

  ihi = ihi + 1

  if ( s_len <= ihi ) then
     ihi = s_len
     return
  end if

  if ( s(ihi:ihi) == ' ' ) then
     exit
  end if

end do
!
!  Decrement IHI to point to the previous, nonblank, character.
!
ihi = ihi - 1

return
end subroutine word_next
subroutine word_next_read ( s, word, done )

  !*****************************************************************************80
  !
  !! WORD_NEXT_READ "reads" words from a string, one at a time.
  !
  !  Discussion:
  !
  !    This routine was written to process tokens in a file.
  !    A token is considered to be an alphanumeric string delimited
  !    by whitespace, or any of various "brackets".
  !
  !    The following characters are considered to be a single word,
  !    whether surrounded by spaces or not:
  !
  !      " ( ) { } [ ]
  !
  !    Also, if there is a trailing comma on the word, it is stripped off.
  !    This is to facilitate the reading of lists.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 May 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string, presumably containing words
  !    separated by spaces.
  !
  !    Output, character ( len = * ) WORD.
  !    If DONE is FALSE, then WORD contains the "next" word read.
  !    If DONE is TRUE, then WORD is blank, because there was no more to read.
  !
  !    Input/output, logical DONE.
  !    On input with a fresh string, set DONE to TRUE.
  !    On output, the routine sets DONE:
  !      FALSE if another word was read,
  !      TRUE if no more words could be read.
  !
implicit none

logical done
integer ( kind = 4 ) ilo
integer ( kind = 4 ), save :: lenc = 0
integer ( kind = 4 ), save :: next = 1
character ( len = * ) s
character, parameter :: TAB = achar ( 9 )
character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
if ( done ) then

  next = 1
  done = .false.
  lenc = len_trim ( s )

  if ( lenc <= 0 ) then
     done = .true.
     word = ' '
     return
  end if

end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
do

  if ( lenc < ilo ) then
     word = ' '
     done = .true.
     next = lenc + 1
     return
  end if
  !
  !  If the current character is blank, skip to the next one.
  !
  if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
     exit
  end if

  ilo = ilo + 1

end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
if ( s(ilo:ilo) == '"' .or. &
    s(ilo:ilo) == '(' .or. &
    s(ilo:ilo) == ')' .or. &
    s(ilo:ilo) == '{' .or. &
    s(ilo:ilo) == '}' .or. &
    s(ilo:ilo) == '[' .or. &
    s(ilo:ilo) == ']' ) then

  word = s(ilo:ilo)
  next = ilo + 1
  return

end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
next = ilo + 1

do while ( next <= lenc )

  if ( s(next:next) == ' ' ) then
     exit
  else if ( s(next:next) == TAB ) then
     exit
  else if ( s(next:next) == '"' ) then
     exit
  else if ( s(next:next) == '(' ) then
     exit
  else if ( s(next:next) == ')' ) then
     exit
  else if ( s(next:next) == '{' ) then
     exit
  else if ( s(next:next) == '}' ) then
     exit
  else if ( s(next:next) == '[' ) then
     exit
  else if ( s(next:next) == ']' ) then
     exit
  end if

  next = next + 1

end do

if ( s(next-1:next-1) == ',' ) then
  word = s(ilo:next-2)
else
  word = s(ilo:next-1)
end if

return
end subroutine word_next_read
subroutine word_next2 ( s, first, last )

  !*****************************************************************************80
  !
  !! WORD_NEXT2 returns the first word in a string.
  !
  !  Discussion:
  !
  !    "Words" are any string of characters, separated by commas or blanks.
  !
  !    The routine returns:
  !    * FIRST, the first string of nonblank, noncomma characters;
  !    * LAST, the characters of the string that occur after FIRST and
  !      the commas and blanks.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string of words to be analyzed.
  !
  !    Output, character ( len = * ) FIRST, the next word in the string.
  !
  !    Output, character ( len = * ) LAST, the remaining string.
  !
implicit none

character c
character ( len = * ) first
integer ( kind = 4 ) i
integer ( kind = 4 ) ido
integer ( kind = 4 ) ifirst
integer ( kind = 4 ) ilast
character ( len = * ) last
integer ( kind = 4 ) lenf
integer ( kind = 4 ) lenl
integer ( kind = 4 ) lens
character ( len = * ) s

first = ' '
last = ' '

ifirst = 0
ilast = 0

lens = len_trim ( s )
lenf = len ( first )
lenl = len ( last )

ido = 0

do i = 1, lens

  c = s(i:i)

  if ( ido == 0 ) then
     if ( c /= ' ' .and. c /= ',' ) then
        ido = 1
     end if
  end if

  if ( ido == 1 ) then
     if ( c /= ' ' .and. c /= ',' ) then
        ifirst = ifirst + 1
        if ( ifirst <= lenf ) then
           first(ifirst:ifirst) = c
        end if
     else
        ido = 2
     end if
  end if

  if ( ido == 2 ) then
     if ( c /= ' ' .and. c /= ',' ) then
        ido = 3
     end if
  end if

  if ( ido == 3 ) then
     ilast = ilast + 1
     if ( ilast <= lenl ) then
        last(ilast:ilast) = c
     end if
  end if

end do

return
end subroutine word_next2
subroutine word_swap ( s, i1, i2 )

  !*****************************************************************************80
  !
  !! WORD_SWAP swaps two words in a given string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, a string of characters.
  !    "Words" in the string are presumed to be separated by blanks.
  !
  !    Input, integer ( kind = 4 ) I1, I2, the indices of the words to be swapped.
  !    If either I1 or I2 is nonpositive, or greater than the number of
  !    words in the string, then nothing is done to the string.  Otherwise,
  !    words I1 and I2 are swapped.
  !
implicit none

logical blank
integer ( kind = 4 ) i
integer ( kind = 4 ) i1
integer ( kind = 4 ) i2
integer ( kind = 4 ) j1
integer ( kind = 4 ) j1beg
integer ( kind = 4 ) j1end
integer ( kind = 4 ) j2
integer ( kind = 4 ) j2beg
integer ( kind = 4 ) j2end
integer ( kind = 4 ) lens
character ( len = * ) s
character ( len = 80 ) s2
integer ( kind = 4 ) word_num

lens = len_trim ( s )
if ( lens <= 0 ) then
  return
end if

if ( 80 < lens ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WORD_SWAP - Warning!'
  write ( *, '(a)' ) '  The internal temporary string is too short'
  write ( *, '(a)' ) '  to copy your string.  Errors may result!'
  stop
end if
!
!  We need to make a copy of the input arguments, because we
!  might alter them.  We want to ensure that J1 <= J2.
!
j1 = min ( i1, i2)
j2 = max ( i1, i2)

if ( j1 <= 0 ) then
  return
else if ( j2 <= 0 ) then
  return
else if ( j1 == j2 ) then
  return
end if

j1beg = 0
j1end = 0
j2beg = 0
j2end = 0
word_num = 0
blank = .true.

do i = 1, lens

  if ( s(i:i) == ' ' ) then

     if ( j1beg /= 0 .and. j1end == 0 ) then
        j1end = i - 1
     else if ( j2beg /= 0 .and. j2end == 0 ) then
        j2end = i - 1
     end if

     blank = .true.

  else if ( blank ) then

     word_num = word_num + 1

     if ( word_num == j1 ) then
        j1beg = i
     else if ( word_num == j2 ) then
        j2beg = i
     end if

     blank = .false.

  end if

end do

if ( j1beg /= 0 .and. j1end == 0 ) then
  j1end = lens
else if ( j2beg /= 0 .and. j2end == 0 ) then
  j2end = lens
end if

if ( word_num < j1 .or. word_num < j2 ) then
  return
end if
!
!  OK, we can swap words J1 and J2.
!
s2 = s
!
!  Copy word 2.
!
s( j1beg : j1beg + j2end - j2beg )  = s2 ( j2beg : j2end )
!
!  Copy (possibly null) string between word 1 and word 2.
!
s ( j1beg + j2end - j2beg + 1 : j1beg + j2end - j1end - 1 ) &
    = s2 ( j1end + 1 : j2beg - 1 )
!
!  Copy word 1.
!
s ( j1beg + j2end - j1end : j2end )  = s2 ( j1beg : j1end )

return
end subroutine word_swap

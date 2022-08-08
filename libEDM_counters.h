#pragma once

#include <libEDM_library.h>
#include <libEDM_matrix.h>

class Error_Counter {
public:

    Error_Counter(const long ignorefirst = 0L) : _bit_errors(0L), _block_errors(0L), _bits(0L), _blocks(0L), _ignorefirst(ignorefirst) {}

    void count(const bVector &in1, const bVector &in2);

    void   clear()        {_bit_errors = 0L; _block_errors = 0L; _bits = 0L; _blocks = 0L;}
    double BER()          {return divide(_bit_errors,   _bits);}
    double BLER()         {return divide(_block_errors, _blocks);}
    size_t bits()         {return _bits;}
    size_t blocks()       {return _blocks;}
    size_t bit_errors()   {return _bit_errors;}
    size_t block_errors() {return _block_errors;}

  private:
    size_t _ignorefirst;
    size_t _bit_errors, _block_errors;
    size_t _bits, _blocks;
};
#include <libEDM_counters.h>

void Error_Counter::count(const bVector &in1, const bVector &in2)
{
    assert( in1.size() == in2.size() );

    if (_blocks >= _ignorefirst)
    {
        _blocks++;
        _bits += in1.size();

        size_t previous_bit_errors = _bit_errors;
        for (size_t i=0; i<in1.size(); i++)
            if (in1[i] != in2[i])
                _bit_errors++;

        if (_bit_errors != previous_bit_errors)
            _block_errors++;
    }
}
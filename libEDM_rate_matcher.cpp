#include <libEDM_rate_matcher.h>

#include <assert.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include <libEDM_library.h>

#include <boost/math/common_factor.hpp>

using std::min;

using boost::math::gcd;

RateMatcher::RateMatchingBlock::ActionType RateMatcher::RateMatchingBlock::action()
{
	_e -= _a * abs(_delta_N);

	if (_e <= 0)
	{
		_e += _a * _N;

		if (_delta_N < 0)
			return Puncture;
		else
			return Repeat;
	}
	else
		return NoChange;
}

dVector RateMatcher::dematch (const dVector &input)
{
	return process(input, &RateMatcher::dematcher);
}

void RateMatcher::dematcher (dVector &output, const RateMatcher::RateMatchingBlock::ActionType action, const dMatrix &frames, const size_t frame, size_t &bit)
{
	switch (action)
	{
	case RateMatchingBlock::Puncture :
		// bit is missing, so output 0.0
		output.push_back(0.0);
		bit--;
		break;

	case RateMatchingBlock::NoChange :
		// copy bit to output
		output.push_back(frames[frame][bit]);
		break;

	case RateMatchingBlock::Repeat :
		// bit is repeated, so add next two bits together and output
		double value = frames[frame][bit];
		bit++;
		value += frames[frame][bit];
		output.push_back(value);
		break;
	}

}

size_t RateMatcher::interleave (const size_t column_id)
{
	if ((column_id < (_interleaver_depth / 2)) && (column_id % 2 == 1))
	{
		// odd column in first half of list
		return column_id + _interleaver_depth / 2;
	}
	else
	{
		if ((column_id >= (_interleaver_depth / 2)) && (column_id % 2 == 0))
			// even column in second half of list
			return column_id - _interleaver_depth / 2;
		else
			// column position is not changed
			return column_id;
	}
}

vector<RateMatcher::SequenceType> RateMatcher::initialise_initial_parity_table()
{
	vector<SequenceType> initial_parity_table(_interleaver_depth);

	switch (_interleaver_depth)
	{
	case 1:
		initial_parity_table[0] = X;
		break;

	case 2:
		initial_parity_table[0] = X;
		initial_parity_table[1] = Y;
		break;

	case 4:
		initial_parity_table[0] = X;
		initial_parity_table[1] = Y_dash;
		initial_parity_table[2] = Y;
		initial_parity_table[3] = X;
		break;

	case 8:
		initial_parity_table[0] = X;
		initial_parity_table[1] = Y;
		initial_parity_table[2] = Y_dash;
		initial_parity_table[3] = X;
		initial_parity_table[4] = Y;
		initial_parity_table[5] = Y_dash;
		initial_parity_table[6] = X;
		initial_parity_table[7] = Y;
		break;
	}

	return initial_parity_table;
}

map<RateMatcher::SequenceType, RateMatcher::SequenceType> RateMatcher::initialise_next_sequence_table()
{
	map<SequenceType, SequenceType> next_sequence_table;
	switch (_interleaver_depth)
	{
	case 1:
	case 4:
		next_sequence_table[X]      = Y;
		next_sequence_table[Y]      = Y_dash;
		next_sequence_table[Y_dash] = X;

	case 2:
	case 8:
		next_sequence_table[X]      = Y_dash;
		next_sequence_table[Y]      = X;
		next_sequence_table[Y_dash] = Y;
	}

	return next_sequence_table;
}

map<RateMatcher::SequenceType, vector<RateMatcher::RateMatchingBlock*> > RateMatcher::initialise_rate_matching_blocks()
{
	map<SequenceType, vector<RateMatchingBlock*> > rate_matching_blocks;

	if (_coder_type == Convolutional || _matched_block_size >= _dematched_block_size)
	{
		// one rate matching block per frame required

		// compute a
		const size_t a = 2;

		// compute e_ini
		const int N       = divide(_dematched_block_size, _interleaver_depth);
		const int delta_N = divide(_matched_block_size,   _interleaver_depth) - N;

		int q;
		if (delta_N == 0)
			q = -1;
		else
			q = N / abs(delta_N);	// rounding down here is intentional

		double q_dash = q;
		if (q % 2 == 0)
			// q is even
            q_dash -= gcd(static_cast<size_t>(q), _interleaver_depth) / _interleaver_depth;

		uVector S(_interleaver_depth, 0);
		for (size_t x = 0; x < _interleaver_depth; x++)
		{
			ldiv_t temp = div(static_cast<long>(ceil(x * q_dash)), _interleaver_depth);
			S[interleave(temp.rem)] = temp.quot;
		}

		vector<RateMatchingBlock*> frame_rate_matching_blocks;
		for (size_t frame = 0; frame < _interleaver_depth; frame++)
		{
			size_t e_ini = (a * S[frame] * abs(delta_N) + N) % (a * N);
			if (e_ini == 0)
				e_ini = a * N;

			// initialise rate matching block and add to frame_rate_matching_blocks
			frame_rate_matching_blocks.push_back(new RateMatchingBlock(N, delta_N, a, e_ini));
		}

		rate_matching_blocks[X] = frame_rate_matching_blocks;
	}
	else
	{
		// Turbo code requiring puncturing

		// compute number of bits per frame per sequence
		map<SequenceType, uVector > N;
		for (size_t frame = 0; frame < _interleaver_depth; frame++)
		{
			SequenceType sequence = _initial_parity_table[frame];

			// compute and store N for initial parity type
			const size_t sequence1_bits_per_frame = static_cast<size_t>(ceil(dematched_frame_size() / 3.0));
			N[sequence].push_back(sequence1_bits_per_frame);

			// compute and store N for next parity type
			sequence = _next_sequence_table[sequence];
			const size_t sequence2_bits_per_frame = static_cast<size_t>(ceil((dematched_frame_size() - sequence1_bits_per_frame) / 2.0));
			N[sequence].push_back(sequence2_bits_per_frame);

			// compute and store N for final parity type
			sequence = _next_sequence_table[sequence];
			const size_t sequence3_bits_per_frame = dematched_frame_size() - sequence1_bits_per_frame - sequence2_bits_per_frame;
			N[sequence].push_back(sequence3_bits_per_frame);
		}

		const int    sequence_frame_size = _matched_block_size / (3 * _interleaver_depth); // rounding down is intentional
		const double unrounded_delta_N   = (static_cast<int>(_matched_block_size) - static_cast<int>(_dematched_block_size)) / (2.0 * _interleaver_depth);

		for (size_t sequence_index = 0; sequence_index < 2; sequence_index++)
		{
			// sequence 0 = Y; sequence 1 = Y dash
			SequenceType sequence;
			if (sequence_index == 0)
				sequence = Y;
			else
				sequence = Y_dash;

			size_t a;
			         int delta_N;
			if (sequence == Y)
			{
				a       = 2;
				delta_N = static_cast<int>(floor(unrounded_delta_N));
			}
			else
			{
				a       = 1;
				delta_N = static_cast<int>(ceil(unrounded_delta_N));
			}

			vector<RateMatchingBlock*> frame_rate_matching_blocks;
			for (size_t frame = 0; frame < _interleaver_depth; frame++)
			{
				const size_t q = N[sequence][frame] / abs(delta_N);

				uVector S(_interleaver_depth, 0);
				if (q <= 2)
				{
					for (size_t x = 0; x < _interleaver_depth; x++)
						if (sequence == 0)
							S[interleave((3*x + 1) % _interleaver_depth)] = x % 2;
						else
							S[interleave((3*x + 2) % _interleaver_depth)] = x % 2;
				}
				else
				{
					double q_dash = q;
					if (q % 2 == 0)
						// q is even
						q_dash -= gcd(q, _interleaver_depth) / _interleaver_depth;

					for (size_t x = 0; x < _interleaver_depth; x++)
					{
						const size_t r = static_cast<size_t>(ceil(x * q_dash)) % _interleaver_depth;
						if (sequence == 0)
							S[interleave((3*r + 1) % _interleaver_depth)] = ceil(x * q_dash) / _interleaver_depth;
						else
							S[interleave((3*r + 2) % _interleaver_depth)] = ceil(x * q_dash) / _interleaver_depth;
					}
				}

				size_t e_ini = (a * S[frame] * abs(delta_N) + N[sequence][frame]) % (a * N[sequence][frame]);
				if (e_ini == 0)
					e_ini = a * N[sequence][frame];

				// initialise rate matching block and add to frame_rate_matching_blocks
				frame_rate_matching_blocks.push_back(new RateMatchingBlock(N[sequence][frame], delta_N, a, e_ini));
			}

			rate_matching_blocks[sequence] = frame_rate_matching_blocks;
		}
	}

	return rate_matching_blocks;
}

dVector RateMatcher::match (const dVector &input)
{
	return process(input, &RateMatcher::matcher);
}

void RateMatcher::matcher (dVector &output, const RateMatcher::RateMatchingBlock::ActionType action, const dMatrix &frames, const size_t frame, size_t &bit)
{
	switch (action)
	{
	case RateMatchingBlock::Puncture :
		// do nothing
		break;

	case RateMatchingBlock::NoChange :
		// copy bit to output
		output.push_back(frames[frame][bit]);
		break;

	case RateMatchingBlock::Repeat :
		// copy bit twice to output
		output.push_back(frames[frame][bit]);
		output.push_back(frames[frame][bit]);
		break;
	}
}

dVector RateMatcher::process (const dVector &input, RateMatcher::ProcessFunction process_function)
{
	reset();

	const dMatrix frames = sort_into_frames(input);

	dVector output;
	for (size_t frame = 0; frame < _interleaver_depth; frame++)
	{
		SequenceType sequence = _initial_parity_table[frame];

		for (size_t bit = 0; bit < frames[frame].size(); bit++)
		{
			RateMatchingBlock::ActionType action;

			if (_coder_type == Convolutional || _matched_block_size >= _dematched_block_size)
			{
				action = _rate_matching_blocks[X][frame]->action();
			}
			else
			{
				// turbo coding requiring puncturing
				switch (sequence)
				{
				case X:
					// no puncturing of X sequence
					action = RateMatchingBlock::NoChange;
					break;

				case Y:
				case Y_dash:
					action = _rate_matching_blocks[sequence][frame]->action();
					assert(action != RateMatchingBlock::Repeat);
					break;
				}

				// find sequence of next bit
				sequence = _next_sequence_table[sequence];
			}

			(this->*process_function)(output, action, frames, frame, bit);
		}

		assert(output.size() == (frame + 1) * matched_frame_size() || output.size() == (frame + 1) * dematched_frame_size());
	}

    return output;
}


void RateMatcher::reset()
{
	for (map<SequenceType, vector<RateMatchingBlock*> >::iterator rmb = _rate_matching_blocks.begin(); rmb != _rate_matching_blocks.end(); rmb++)
		for (size_t frame = 0; frame < _interleaver_depth; frame++)
			rmb->second[frame]->reset();
}


dMatrix RateMatcher::sort_into_frames (const dVector &input)
{
	dMatrix output(_interleaver_depth);

	const size_t bits_per_frame = input.size() / _interleaver_depth;

	for (size_t frame = 0, input_index = 0; frame < _interleaver_depth; frame++)
		for (size_t bit = 0; bit < bits_per_frame; bit++, input_index++)
			output[frame].push_back(input[input_index]);

	return output;
}
#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <vector>

#include <libEDM_matrix.h>

using std::map;
using std::vector;

class RateMatcher {
public:
	enum CoderType {Convolutional, Turbo};

private:
	class RateMatchingBlock {
	private:
		const size_t _N, _a, _e_ini;
		const int    _delta_N;
		      int    _e;

	public:
		typedef enum{Puncture, NoChange, Repeat} ActionType;

		// constructor
		RateMatchingBlock(const size_t N, const int delta_N, const size_t a, const size_t e_ini) : _N(N), _delta_N(delta_N), _a(a), _e_ini(e_ini), _e(e_ini) {}

		ActionType action();
		void       reset() {_e = _e_ini;}
	};

	enum SequenceType {X, Y, Y_dash};

	typedef void (RateMatcher::*ProcessFunction)(dVector &output, const RateMatchingBlock::ActionType action, const dMatrix &frames, const size_t frame, size_t &bit);

	const size_t    _interleaver_depth;
	const size_t    _matched_block_size;
	const size_t    _dematched_block_size;
	const CoderType _coder_type;

	vector<SequenceType>            _initial_parity_table;
	map<SequenceType, SequenceType> _next_sequence_table;

	map<SequenceType, vector<RateMatchingBlock*> > _rate_matching_blocks;

	size_t matched_frame_size   () const {return _matched_block_size / _interleaver_depth;}
	size_t dematched_frame_size () const {return _dematched_block_size / _interleaver_depth;}

	vector<SequenceType>                           initialise_initial_parity_table();
	map<SequenceType, vector<RateMatchingBlock*> > initialise_rate_matching_blocks();
	map<SequenceType, SequenceType>                initialise_next_sequence_table();

	void    dematcher        (dVector &output, const RateMatchingBlock::ActionType action, const dMatrix &frames, const size_t frame, size_t &bit);
	void    matcher          (dVector &output, const RateMatchingBlock::ActionType action, const dMatrix &frames, const size_t frame, size_t &bit);
    dVector process          (const dVector &input, ProcessFunction process_function);
	dMatrix sort_into_frames (const dVector &input);

	size_t interleave (const size_t column_id);

public:
	// constructors
	RateMatcher(const CoderType coder_type, const size_t interleaver_depth, const size_t matched_block_size, const size_t dematched_block_size)
				: _coder_type(coder_type), _interleaver_depth(interleaver_depth), _matched_block_size(matched_block_size), _dematched_block_size(dematched_block_size), _rate_matching_blocks(initialise_rate_matching_blocks()), _initial_parity_table(initialise_initial_parity_table()), _next_sequence_table(initialise_next_sequence_table()) {}

	dVector match   (const dVector &input);
	dVector dematch (const dVector &input);

	void reset();
};
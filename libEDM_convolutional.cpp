#include <libEDM_convolutional.h>

#include <cassert>
#include <limits>

#include <libEDM_library.h>

#define INFINITY numeric_limits<double>::infinity()

using std::abs;
using std::exp;
using std::log;
using std::make_pair;
using std::max;
using std::min;
using std::numeric_limits;

inline double com_logmap(const double x, const double y)
{
    double output = max(x,y);
    if ((x != -INFINITY) && (y != -INFINITY))
        output += log(1.0 + exp(-fabs(y-x)));

    return output;
}

inline double com_logmax(const double x, const double y)
{
    return max(x,y);
}

size_t ConvolutionalCodec::calc_state_transition (const size_t inputState, const bool inputBit, dynamic_bitset<> &parityBits) const
{
    // initialise state
    dynamic_bitset<> state(m+1, inputState);

    // compute feedback bit
    bool feedbackBit = (state & feedbackGenerator).count() % 2;

    // modify state by feedback bit
    state[m] = feedbackBit ^ inputBit;

    // compute parity bits
    parityBits.resize(n);
    for (size_t i=0; i<n; i++)
       parityBits[i] = (state & parityGenerators[i]).count() % 2;

    // shift state right and return
    state >>= 1;
    return state.to_ulong();
}

ConvolutionalCodec::ConvolutionalCodec(const size_t constraintLength, const uVector &parityGenerators, const size_t feedbackGenerator, const bool terminated, const bool systematic, const size_t decoderDelay)
    : terminated                 (terminated),
      systematic                 (systematic),
      decoderDelay               (decoderDelay),
      Lc                         (1.0),
      n                          (parityGenerators.size()),
      m                          (constraintLength-1),
      numStates                  (1 << m),
      incrementalMetricsTable    (1 << n),
      incrementalSOVAMetricsTable(1 << n),
      incrementalDistanceTable   (1 << n),
      outputParity               (Matrix<dynamic_bitset<> >(numStates, 2, dynamic_bitset<>(n))),
      stateTransition            (uMatrix(numStates, 2, 0))
{
    // store feedback generator as dynamic_bitset
    this->feedbackGenerator = dynamic_bitset<>(constraintLength, feedbackGenerator);

    // store parity generators as dynamic_bitsets
    for (size_t i=0; i<n; i++)
    {
        dynamic_bitset<> parityGenerator(constraintLength, static_cast<size_t>(parityGenerators[i]));
        this->parityGenerators.push_back(parityGenerator);
    }

    reverseStateTransition.resize(numStates);
    for (size_t bit=0; bit < 2; bit++)
        for (size_t state=0; state < numStates; state++)
        {
            const size_t nextState = calc_state_transition(state, bit, outputParity[state][bit]);
            stateTransition[state][bit] = nextState;
            reverseStateTransition[nextState].push_back(pair<size_t,bool>(state,bit));
        }

    // create and initialise the states of the Viterbi decoder
	oddStates .reserve(numStates);
	evenStates.reserve(numStates);
	for (size_t stateIndex = 0; stateIndex < numStates; stateIndex++)
	{
		oddStates .push_back(new State(stateIndex, this, &evenStates));
		evenStates.push_back(new State(stateIndex, this, &oddStates));
	}
}

dVector ConvolutionalCodec::decode(const dVector &input, const dVector &extrinsicInput, const Metric metric)
{
    switch ( metric )
    {
    case MAP :
        return decode_MAP(input, extrinsicInput);

    case LOGMAP :
    case LOGMAX :
        return decode_LOG(input, extrinsicInput, metric);

    case SOVA :
        return decode_SOVA(input, dVector(input.size(),1.0));

    default:
		error("ConvolutionalCode::decode metric not supported");
    }
}

bVector ConvolutionalCodec::decode(const dVector &input, const Metric metric)
{
    bVector output;

    const dVector llr = decode(input, dVector(0), metric);

    for (size_t i=0; i<llr.size()-m; i++)
        output.push_back(llr[i] < 0.0);

    return output;
}

bVector ConvolutionalCodec::encode(const bVector &input, bVector &tail, bMatrix &parityBits) const
{
    bVector output;

    size_t numTailBits = 0;
    if (terminated)
        numTailBits = m;

    size_t numParityBits = input.size() + numTailBits;

    parityBits.set_size(numParityBits, n);
    tail.resize(numTailBits);

    size_t encoderState = 0;
    for (size_t i=0; i<input.size(); i++)
    {
        if ( systematic )
            output.push_back(input[i]);

        for (size_t j=0; j<n; j++)
        {
	        const bool parityBit = outputParity[encoderState][input[i]][j];
            parityBits[i][j] = parityBit;
            output.push_back(parityBit);
        }

        encoderState = stateTransition[encoderState][input[i]];
    }
  
    for (size_t i=0; i<numTailBits; i++)
    {
        // find tail bit that will shift state to right
        const size_t targetState = encoderState >> 1;
        tail[i] = (stateTransition[encoderState][true] == targetState);

        if ( systematic )
            output.push_back(tail[i]);

        for (size_t j=0; j<n; j++)
        {
            const bool parityBit = outputParity[encoderState][tail[i]][j];
	        parityBits[input.size()+i][j] = parityBit;
            output.push_back(parityBit);
        }

        encoderState = targetState;
    }

    return output;
}

void ConvolutionalCodec::unsplice(const dVector &input, dVector &rxSystematic, dMatrix &rxParity) const
{
    const size_t blockLength = input.size() / (n + systematic);

    rxSystematic.clear();
    rxParity.set_size(n,0);

    size_t bitIndex = 0;
    for (size_t k=0; k<blockLength; k++)
    {
        if ( systematic )
            rxSystematic.push_back(input[bitIndex++]);
        else
            rxSystematic.push_back(0.0);

        for (size_t i=0; i<n; i++)
            rxParity[i].push_back(input[bitIndex++]);
    }

    assert( rxSystematic.size() == rxParity[0].size() );
}

dVector ConvolutionalCodec::decode_MAP(const dVector &input, const dVector &extrinsicInput)
{
    dVector output;

    // extract systematic and parity information
    dVector rxSystematic;
    dMatrix rxParity;
    unsplice(input, rxSystematic, rxParity);

    // compute block length, including tail bits
    const size_t blockLength = rxSystematic.size();

    // initialise gamma
    dCubrix gamma(numStates, blockLength, 2, 0.0);

    // merge parity information
    dCubrix parity(numStates, blockLength, 2, 0.0);
    for (size_t state=0; state < numStates; state++)
        for (size_t k=0; k < blockLength; k++)
            for (size_t bit=0; bit<2; bit++)
            {
                double temp = 0.0;
                for (size_t i=0; i<n; i++)
                    temp += (outputParity[state][bit][i] ? -rxParity[i][k] : rxParity[i][k]);

                parity[state][k][bit] = exp(0.5 * temp);
            }

    // calculate gamma 
    for (size_t k=0; k<blockLength; k++)
    {
        double ex = 0.5 * rxSystematic[k];
        if (k < extrinsicInput.size())
            ex += 0.5 * extrinsicInput[k];

        for (size_t state=0; state<numStates; state++)
        {
            gamma[state][k][false] = parity[state][k][false] * exp( ex);
            gamma[state][k][true ] = parity[state][k][true ] * exp(-ex);
        }
    }

    // initialise alpha
    dMatrix alpha(numStates, blockLength+1, 0.0);
    alpha[0][0] = 1.0;

    // calculate alpha going forward through the trellis
    for (size_t k=0; k<blockLength; k++)
    {
        double sum = 0.0;
        for (size_t state=0; state<numStates; state++)
            for (size_t bit = 0; bit < 2; bit++)
            {
	            const size_t nextState = stateTransition[state][bit];

                double temp = alpha[state][k] * gamma[state][k][bit];
                alpha[nextState][k+1] += temp;
                sum                   += temp;
            }

        // normalise alpha
        for (size_t state=0; state < numStates; state++)
            alpha[state][k+1] /= sum;
    }
    
    // initialise beta
    dMatrix beta(numStates, blockLength+1, 0.0);
    if (terminated)
        beta[0][blockLength] = 1.0;
    else
        for (size_t state=0; state<numStates; state++)
            beta[state][blockLength] = alpha[state][blockLength];
  
    // calculate beta going backward through the trellis
    for (size_t k=blockLength; k>0; k--)
    {
        double sum = 0.0;
        for (size_t state=0; state < numStates; state++)
            for (size_t bit = 0; bit < 2; bit++)
            {
                const size_t nextState = stateTransition[state][bit];

                const double temp = beta[nextState][k] * gamma[state][k-1][bit];
                beta[state][k-1] += temp;
                sum              += temp;
            }

        // normalise beta
        for (size_t state=0; state < numStates; state++)
            beta[state][k-1] /= sum;
    }
  
    // calculate extrinsic output for each bit
    for (size_t k=0; k<extrinsicInput.size(); k++)
    {
        dVector prob(2, 0.0);
        for (size_t state=0; state<numStates; state++)
            for (size_t bit = 0; bit < 2; bit++)
            {
   	            const size_t nextState = stateTransition[state][bit];
                prob[bit] += alpha[state][k] * parity[state][k][bit] * beta[nextState][k+1];
            }

        output.push_back(log(prob[false]) - log(prob[true]));
    }

    return output;
}

dVector ConvolutionalCodec::decode_LOG(const dVector &input, const dVector &extrinsicInput, const Metric metric)
{
    dVector output;

    dVector rxSystematic;
    dMatrix rxParity;
    unsplice(input, rxSystematic, rxParity);

    const size_t blockLength = rxSystematic.size();

    double (*com_log)(const double, const double);

    // set the internal metric
    if ( metric == ConvolutionalCodec::LOGMAX )
        com_log = com_logmax;
    else
        com_log = com_logmap;
 
    // initialise gamma
    dCubrix gamma(numStates, blockLength, 2, 0.0);

    // merge parity information
    dCubrix parity(numStates, blockLength, 2, 0.0);
    for (size_t state=0; state < numStates; state++)
        for (size_t k=0; k < blockLength; k++)
            for (size_t bit=0; bit<2; bit++)
                for (size_t i=0; i<n; i++)
                    parity[state][k][bit] += (outputParity[state][bit][i] ? -rxParity[i][k] : rxParity[i][k]);

    // calculate gamma
    for (size_t k=0; k<blockLength; k++)
    {
        double ex = rxSystematic[k];
        if (k < extrinsicInput.size())
            ex += extrinsicInput[k];

        for (size_t state=0; state<numStates; state++)
        {
	        gamma[state][k][false] = 0.5 * (parity[state][k][false] + ex);
	        gamma[state][k][true ] = 0.5 * (parity[state][k][true ] - ex);
        }
    }

    // initialise alpha
    dMatrix alpha(numStates, blockLength+1, -INFINITY);
    alpha[0][0] = 0.0;
  
    // calculate alpha going forward through the trellis
    dVector denom(blockLength+1, -INFINITY);
    for (size_t k=0; k<blockLength; k++)
    {
        for (size_t state=0; state<numStates; state++)
        {
            dVector temp(2);
            for (size_t bit = 0; bit < 2; bit++)
            {
	            temp[bit] = alpha[state][k] + gamma[state][k][bit];

                const size_t nextState = stateTransition[state][bit];
	            alpha[nextState][k+1] = com_log(alpha[nextState][k+1], temp[bit]);
                //  denom[k]    = com_log(alpha[s][k], denom[k]);
            }
        }

        // normalise alpha
        const double norm = alpha[0][k+1];
        for (size_t state=0; state<numStates; state++)
//          alpha[state][k] -= denom[k];
            alpha[state][k+1] -= norm;
    }
  
    // initialise beta
    dMatrix beta(numStates, blockLength+1, -INFINITY);
    if (terminated)
        beta[0][blockLength] = 0.0;
    else
        for (size_t s=0; s<numStates; s++)
            beta[s][blockLength] = alpha[s][blockLength];
  
    // calculate beta going backward through the trellis
    for (size_t k=blockLength; k>0; k--)
    { 
        for (size_t state=0; state<numStates; state++)
        {
            dVector temp(2);
            for (size_t bit = 0; bit < 2; bit++)
            {
	            const size_t nextState = stateTransition[state][bit];
                temp[bit] = beta[nextState][k] + gamma[state][k-1][bit];
            }
	        beta[state][k-1] = com_log(temp[false], temp[true]);
        }

        // normalise beta
        const double norm = beta[0][k-1];
        for (size_t state=0; state<numStates; state++)
//          beta[state][k] -= denom[k];
            beta[state][k-1] -= norm;
    }

    // calculate extrinsic output for each bit
    for (size_t k=0; k<extrinsicInput.size(); k++)
    {
        dVector prob(2, -INFINITY);
        for (size_t state=0; state<numStates; state++)
            for (size_t bit=0; bit<2; bit++)
            {
                const size_t nextState = stateTransition[state][bit];
                prob[bit] = com_log(prob[bit], alpha[state][k] + 0.5 * parity[state][k][bit] + beta[nextState][k+1]);
            }

        output.push_back(prob[false] - prob[true]);
    }

    return output;
}

double ConvolutionalCodec::State::process_transitions(const size_t sourceBitIndex)
{
	// update state information

	const double metricA = transitionA.cumulativeMetric();
	const double metricB = transitionB.cumulativeMetric();

	const double SOVAMetricA = transitionA.cumulativeSOVAMetric();
	const double SOVAMetricB = transitionB.cumulativeSOVAMetric();

	const State *acceptedState, *rejectedState;
	if (metricA < metricB)
	{
		// choose transitionA
		_cumulativeMetric                = metricA;
		_cumulativeSOVAMetric            = SOVAMetricA;
		_cumulativeDistance              = transitionA.cumulativeDistance();
		acceptedState                    = transitionA.previousState();
        sequenceEstimate[sourceBitIndex] = transitionA.inputBit;
		rejectedState                    = transitionB.previousState();
	}
	else
	{
		// choose transitionB
		_cumulativeMetric                = metricB;
		_cumulativeSOVAMetric            = SOVAMetricB;
		_cumulativeDistance              = transitionB.cumulativeDistance();
		acceptedState                    = transitionB.previousState();
        sequenceEstimate[sourceBitIndex] = transitionB.inputBit;
		rejectedState                    = transitionA.previousState();
	}

	const size_t returnIndex = sourceBitIndex - min(codec->decoderDelay, sourceBitIndex);

	// update sequence_estimate and sequence_prob_error
	const double probNewError = 1.0 / (1.0 + exp(fabs(SOVAMetricA - SOVAMetricB)));
	for (size_t bitIndex = returnIndex; bitIndex < sourceBitIndex; bitIndex++)
	{
		const bool currentBit = acceptedState->sequenceEstimate[bitIndex];
		sequenceEstimate[bitIndex] = currentBit;

		// update probability of bit error if rejected_sequence_estimate bit is different, otherwise copy value
		const double currentProbError = acceptedState->sequenceProbError[bitIndex];
		sequenceProbError[bitIndex] = currentProbError;
		if (currentBit != rejectedState->sequenceEstimate[bitIndex])
			// modify probability that bit is in error
			sequenceProbError[bitIndex] += probNewError - 2.0 * probNewError * currentProbError;
	}

    // return likelihood ratio
    if ( codec->decoderDelay <= sourceBitIndex )
        // need to return hard decision
        if ( sequenceEstimate[returnIndex] )
            // sequenceProbError is probability that bit is not true (i.e. is false)
            return sequenceProbError[returnIndex] / (1.0 - sequenceProbError[returnIndex]);
        else
            return (1.0 - sequenceProbError[returnIndex]) / sequenceProbError[returnIndex];
    else
        // hard decision not required at this stage
        return 0.0;
}

void ConvolutionalCodec::State::reset(const size_t maxSize)
{
	// reset _cumulative_metric
	// setting to infinity means that path will always lose out to a path that started at state 0
	if (index == 0)
		_cumulativeMetric = 0.0;
	else
		_cumulativeMetric = INFINITY;

	_cumulativeSOVAMetric = 0.0;
	_cumulativeDistance    = 0;

	// reset sequence_estimate and sequence_prob_error
	sequenceEstimate .assign(maxSize, false);
    sequenceProbError.assign(maxSize, 0.0);
}

void ConvolutionalCodec::compute_incremental_metrics_table(dVector::const_iterator receivedTransitionBits, dVector::const_iterator snrEstimates)
{
	// incremental metrics do not depend on state, so can be computed and stored in a LUT indexed by the transition output bits
	for (size_t transition = 0; transition < incrementalMetricsTable.size(); transition++)
	{
		double incrementalMetric = 0.0, incrementalSOVAMetric = 0.0;
		size_t incrementalDistance = 0;

		for (size_t bitIndex = 0; bitIndex < n; bitIndex++)
		{
			const double sentTransitionBit     = transition & (1<<bitIndex) ? -1.0 : 1.0;
			const double receivedTransitionBit = *(receivedTransitionBits + bitIndex);

			// metric is computed as described by Hagenauer
			const double increment = sqr(receivedTransitionBit - sentTransitionBit);
			const double EsNo      = *(snrEstimates + bitIndex);

			// Don't scale metric used for hard decisions by Es/No as this only screws things up
			incrementalMetric += increment;

            if ( (EsNo == INFINITY) || (EsNo == 1.0) )
				// ignore Es_No if infinite or unity
			    incrementalSOVAMetric += increment;
            else
                incrementalSOVAMetric += EsNo * increment;

			if ( sgn(sentTransitionBit) != sgn(receivedTransitionBit) )
				incrementalDistance++;
		}

		incrementalMetricsTable    [transition] = incrementalMetric;
		incrementalSOVAMetricsTable[transition] = incrementalSOVAMetric;
		incrementalDistanceTable   [transition] = incrementalDistance;
	}
}

dVector ConvolutionalCodec::decode_SOVA(const dVector &input, const dVector &snrEstimates)
{
	dVector output;

	size_t sourceBitsPerBlock  = input.size() / n;	        // includes (constraint_length - 1) tail bits
	size_t uncodedBitsPerBlock = sourceBitsPerBlock - m;	// includes CRC check

	// reset codec
	reset(sourceBitsPerBlock);

	dVector codedBits;
	codedBits.reserve(n);

	// iterate over all source bits
	vector<State*> *states;
    size_t numStates = this->numStates;
	for (size_t sourceBitIndex = 0; sourceBitIndex < sourceBitsPerBlock; sourceBitIndex++)
	{
		if (sourceBitIndex > uncodedBitsPerBlock)
			numStates /= 2;

		// set states to either even_states or odd_states, depending on source_bit_index
		states = sourceBitIndex % 2 ? &oddStates : &evenStates;

		// compute incremental metrics table
		size_t offset = sourceBitIndex * n;
		compute_incremental_metrics_table(input.begin() + offset, snrEstimates.begin() + offset);

		double minMetric = INFINITY;
		double bitLRDecision;
		for (size_t stateIndex = 0; stateIndex < numStates; stateIndex++)
		{		
			double bitLR = (*states)[stateIndex]->process_transitions(sourceBitIndex);

			// force a decision after traceback delay
			if (sourceBitIndex >= decoderDelay)
				if ((*states)[stateIndex]->cumulativeMetric() < minMetric)
				{
					minMetric = (*states)[stateIndex]->cumulativeMetric();
					bitLRDecision = bitLR;
			    }
		}

		if (sourceBitIndex >= decoderDelay)
			output.push_back(log(bitLRDecision));
	}


	// copy remaining LLRs (including tail bits) from state 0 to output vector
	const size_t startBit = max(0, static_cast<int>(sourceBitsPerBlock) - static_cast<int>(decoderDelay));
	for (size_t outputBitIndex = startBit; outputBitIndex < sourceBitsPerBlock; outputBitIndex++)
	{
        const double probError = (*states)[0]->sequenceProbError[outputBitIndex];
        if ( (*states)[0]->sequenceEstimate[outputBitIndex] )
            // probError is probability of bit not being true (i.e. is false)
		    output.push_back(log(probError / (1.0 - probError)));
        else
		    output.push_back(log((1.0 - probError) / probError));
	}

	// output.estimated_rawBER = static_cast<double>((*states)[0]->cumulative_distance()) / input.size();

	return output;
}

void ConvolutionalCodec::reset(const size_t maxSize)
{
	// reset each state
	for (size_t stateIndex = 0; stateIndex < numStates; stateIndex++)
	{
		oddStates [stateIndex]->reset(maxSize);
		evenStates[stateIndex]->reset(maxSize);
	}
}
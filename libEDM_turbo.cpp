#include <libEDM_turbo.h>

#include <boost/math/common_factor.hpp>

using boost::math::gcd;

bVector TurboCodec::encode (const bVector &input) const
{
    bVector output;

    // initializations
    size_t numBlocks = input.size() / numUncoded;

    // encode all code blocks
    for (size_t i=0; i<numBlocks; i++)
    {    
        bMatrix tail  (numCoders,0,0.0);
        bCubrix parity(numCoders,0,0,0.0);

        // encode block
        bVector systematic = input.mid(i*numUncoded, numUncoded);

        rscCodec.encode(systematic, tail[0], parity[0]);
        for (size_t coder = 1; coder < numCoders; coder++)
        {
            bVector interleaved = bInterleavers[coder-1].interleave(systematic); 
            rscCodec.encode(interleaved, tail[coder], parity[coder]);
        }

        // the data part
        for (size_t k=0; k<numUncoded; k++)
        {
            // systematic bits
	        output.push_back(systematic[k]);

            // parity bits
            for (size_t coder = 0; coder < numCoders; coder++)
	            for (size_t j=0; j<n; j++)
                    output.push_back(parity[coder][k][j]);
        }
    
        // tail bits
        for (size_t coder = 0; coder < numCoders; coder++)
            for (size_t k=0; k<m; k++)
            {
                // systematic tail bits
                output.push_back(tail[coder][k]);

                // parity tail bits
	            for (size_t j=0; j<n; j++)
                    output.push_back(parity[coder][numUncoded+k][j]);
            }
    }

    return output;
}

dVector TurboCodec::unsplice (dVector::const_iterator rxSignal, dMatrix &received) const
{
    received.set_size(numCoders,0);

    dMatrix rxParity(numCoders,0,0.0), rxSystematic(numCoders,0,0.0);

    // data part
    for (size_t k=0; k<numUncoded; k++)
    {
	    rxSystematic[0].push_back(*rxSignal++);

        for (size_t coder = 0; coder < numCoders; coder++)
	        rxParity[coder].push_back(*rxSignal++);
    }

    for (size_t coder=1; coder<numCoders; coder++)
        rxSystematic[coder] = dInterleavers[coder-1].interleave(rxSystematic[0]);

    for (size_t coder = 0; coder < numCoders; coder++)
        for (size_t k=0; k<numUncoded; k++)
        {
            received[coder].push_back(rxSystematic[coder][k]);
            received[coder].push_back(rxParity    [coder][k]);
        }

    // tail bits
    for (size_t coder = 0; coder < numCoders; coder++)
        for (size_t k=0; k<m; k++)
        {
	        received[coder].push_back(*rxSignal++);
	        received[coder].push_back(*rxSignal++);
        }

    // scale the input data
    if (Lc != 1.0)
        for (size_t coder = 0; coder < numCoders; coder++)
	        received[coder] *= Lc;

    return rxSystematic[0];
}

bVector TurboCodec::decode(const dVector &rxSignal, const bVector &trueBits)
{
    bVector output;
    const size_t numBlocks = rxSignal.size() / numCoded;

    // initialise rxSignal iterator
    dVector::const_iterator rxSignalIterator = rxSignal.begin();

    // split received code block into systematic and parity bits
    for (size_t block=0; block<numBlocks; block++)
    {
        dMatrix received;
        dVector rxSystematic = unsplice(rxSignalIterator, received);

        // decode the block
        dVector extrinsic(numUncoded);
        bVector lastDecodedBlock, decodedBlock(numUncoded);

        for (size_t iteration=0; iteration<numIterations; iteration++)
        {
            // decode first RSC code and store extrinsic information
	        extrinsic = rscCodec.decode(received[0], extrinsic, metric);

            // update LLR
            dVector LLR = rxSystematic + extrinsic;

            for (size_t coder = 1; coder < numCoders; coder++)
            {
                // interleave current extrinsic information
                dVector interleavedExtrinsic = dInterleavers[coder-1].interleave(extrinsic);

                // decode received coder bits, passing current extrinsic information (suitably) interleaved to decoder
                interleavedExtrinsic = rscCodec.decode(received[coder], interleavedExtrinsic, metric);

                // deinterleave updated extrinsic information
                extrinsic = dInterleavers[coder-1].deinterleave(interleavedExtrinsic);

                // update LLR
                LLR += extrinsic;
            }

            // threshold extrinsic information to make bit decisions
	        for (size_t k=0; k<numUncoded; k++)
                decodedBlock[k] = (LLR[k]<0.0);

            if ( adaptiveStop )
            {
	            if ( decodedBlock == lastDecodedBlock )
	                break;
                lastDecodedBlock = decodedBlock;
            }
            else
                if ( (trueBits.size() == numBlocks * numUncoded) && (decodedBlock == trueBits.mid(block*numUncoded, numUncoded)) )
	                break;
        }
      
        // copy final bit decisions to decoded bits vector
        output.ins(output.size(), decodedBlock);
    }
    return output;
}

uVector wcdma_turbo_interleaver_sequence(size_t interleaverSize)
{
    const size_t MAX_INTERLEAVER_SIZE = 5114;
    const size_t MIN_INTERLEAVER_SIZE = 40;

    assert( (interleaverSize >= MIN_INTERLEAVER_SIZE) && (interleaverSize <= MAX_INTERLEAVER_SIZE) );

    size_t K = interleaverSize;
    
    //Definitions of primes and associated primitive roots:
    size_t prime_array[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257};
    size_t root_array[]  = {0, 0, 0, 3, 2, 2, 3, 2, 5, 2, 3, 2, 6, 3, 5, 2, 2, 2, 2, 7, 5, 3, 2, 3, 5, 2, 5, 2, 6, 3, 3, 2, 3, 2, 2, 6, 5, 2, 5, 2, 2, 2, 19, 5, 2, 3, 2, 3, 2, 6, 3, 7, 7, 6, 3};
    uVector primes(prime_array, 55);
    uVector roots (root_array,  55);
    
    //Determine R
    size_t rows = 20;
    if ((K>=40) && (K<=159))
	    rows = 5;
    else
        if ( ((K>=160)&&(K<=200)) || ((K>=481)&&(K<=530)) )
	        rows = 10;
    
    //Determine cols
    size_t cols, p = 0, v = 0;
    if ((K>=481) && (K<=530))
    {
	    p = 53;
	    v = 2;
	    cols = p;
    }
    else
    {
	    //Find minimum prime p such that (p+1) - K/rows >= 0 ...
	    for (size_t i=0; i<primes.size(); i++)
	        //if ( (primes[i] + 1 - static_cast<double>(K)/rows) >= 0.0 )
	        if ( (primes[i] + 1) * rows >= K )
            {
	            p = primes[i];
	            v = roots[i];
	            break;
	        }

	    //... and determine cols such that
	    if ( p * rows >= K )
	        if ( (p - 1) * rows >= K )
	            cols = p-1;
            else
	            cols = p;
        else
	        cols = p+1;
    }
    
    //Construct the base sequences for intra-row permutaions
    uVector s(p-1);
    s[0] = 1;
    for (size_t i=1; i<=(p-2); i++)
	    s[i] = (v * s[i-1]) % p;
    
    //Let q(0) = 1 be the first prime integer in {q(j)}, and select the consecutive 
    //minimum prime integers {q(j)}, j = 1, 2, ..., (rows-1) such that gcd( q(j), p-1) == 1, q(j) > 6, and q(j) > q(j-1)
    uVector q(rows);
    q[0] = 1;
    for (size_t j=1; j<rows; j++)
	    for (size_t i=0; i<primes.size(); i++)
        {
	        size_t qj = primes[i];
	        if ( (qj>6) && (qj>q[j-1]) )
                if (gcd(qj, p-1) == 1)
                {
	                q[j] = qj;
	                break;
	            }
        }
    
    //Definitions of Pat1, Pat2, Pat3, and Pat4:
    size_t pattern1_array[] = {4, 3, 2, 1, 0};
    size_t pattern2_array[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    size_t pattern3_array[] = {19, 9, 14, 4, 0, 2, 5, 7, 12, 18, 16, 13, 17, 15, 3, 1, 6, 11, 8, 10};
    size_t pattern4_array[] = {19, 9, 14, 4, 0, 2, 5, 7, 12, 18, 10, 8, 13, 17, 3, 1, 16, 6, 15, 11};
    
    uVector pattern1(pattern1_array, 5);
    uVector pattern2(pattern2_array, 10);
    uVector pattern3(pattern3_array, 20);
    uVector pattern4(pattern4_array, 20);

    //T(j) is the inter-row permutation patters defined as one of the following four
    //kinds of patterns: Pat1, Pat2, Pat3, and Pat4 depending on the number of input bits K
    uVector T;
    switch(rows) {
        case 5:
            T = pattern1;
            break;
        case 10:
            T = pattern2;
            break;
        case 20:
            if ( ((K>2281) && (K<2481)) || ((K>3160) && (K<3211)) )
                T = pattern3;
            else
                T = pattern4;
    }

    //Permute {q(j)} to make {r(j)} such that r(T(j)) = q(j), j = 0, 1, ..., (rows-1),
    //where T(j) indicates the original row position of the j-th permuted row
    uVector r(rows);
    for (size_t i=0; i<rows; i++)
        r[i] = q[T[i]];
    
    //U(j,i) is the input bit position of i-th output after the permutation of j-th row
    //Perform the j-th (j=0, 1, 2, ..., (rows-1)) intra-row permutation as
    uMatrix U(rows, cols, 0);
    if ( cols == p )
	    for (size_t row=0; row<rows; row++)
        {
	        for (size_t col=0; col<(p-1); col++)
	            U[row][col] = s[(col * r[row]) % (p-1)];
	        U[row][p-1] = 0;
	    }
    else
        if ( cols == (p+1) )
        {
	        for (size_t row=0; row<rows; row++)
            {
	            for (size_t col=0; col<(p-1); col++)
	                U[row][col] = s[(col * r[row]) % (p-1)];

	            U[row][p-1] = 0;
	            U[row][p]   = p;
	        }
	        if ( K == (rows*cols) )
            {
	            size_t temp = U[rows-1][p];
	            U[rows-1][p] = U[rows-1][0];
	            U[rows-1][0] = temp;
	        }
        } 
        else
            if ( cols == (p-1) )
	            for (size_t row=0; row<rows; row++)
	                for (size_t col=0; col<cols; col++)
	                    U[row][col] = s[(col * r[row]) % cols] - 1;

    
    //Calculate the interleaver sequence:
    uVector I(K);
    size_t count = 0;
    for (size_t col=0; col<cols; col++)
	    for (size_t row=0; row<rows; row++)
        {
	        size_t index = T[row] * cols + U[T[row]][col];
	        if (index < K)
	            I[count++] = index;
	    }
    
    return I;    
}
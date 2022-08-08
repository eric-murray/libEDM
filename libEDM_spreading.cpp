#include <libEDM_spreading.h>

dMatrix wcdma_spreading_codes(const size_t SF)
{
    dMatrix codes(SF,SF);

    if (SF == 1)
        codes[0][0] = 1.0;
    else
	{
        dMatrix prevCodes = wcdma_spreading_codes(SF/2);
        for (size_t i=0; i<SF/2; i++)
		{
			codes[2*i]  .replace_mid(0,    prevCodes[i]);
			codes[2*i]  .replace_mid(SF/2, prevCodes[i]);
			codes[2*i+1].replace_mid(0,    prevCodes[i]);
			codes[2*i+1].replace_mid(SF/2, prevCodes[i] * -1.0);
        }
    }
    return codes;
}

#pragma once

#include <libEDM_types.h>

//
// class degrees
//

degrees::degrees (const radians &angle) : unit(bound(180.0 * angle() * M_1_PI)) {}

double degrees::bound (double value) const
{
   while ( value > 180.0 )
        value -= 360.0;
    while ( value < -180.0 )
        value += 360.0;
    return value;
}

//
// class radians
//

radians::radians (const degrees &angle) : unit(bound(M_PI * angle() / 180.0)) {}

double radians::bound (double value) const
{
    while ( value > M_PI )
        value -= 2.0 * M_PI;
    while ( value < -M_PI )
        value += 2.0 * M_PI;
    return value;
}

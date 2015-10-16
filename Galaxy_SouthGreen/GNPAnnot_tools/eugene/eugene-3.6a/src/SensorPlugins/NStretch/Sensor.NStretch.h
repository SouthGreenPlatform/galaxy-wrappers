// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id: Sensor.NStretch.h,v 1.2 2009-04-30 12:45:22 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.NStretch.h
// Contents: Sensor NStretch
// ------------------------------------------------------------------

#ifndef  SENSOR_NSTRETCH_INCLUDED
#define  SENSOR_NSTRETCH_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorNStretch                    **
 *************************************************************/
class SensorNStretch : public Sensor
{
private:
    unsigned char *insideNStretch;
    double stretchPenalty;
    int maxLengthWithoutPenalty; // Penalize region of N with a length upper than maxLengthWithoutPenalty

public:
    SensorNStretch  (int n, DNASeq *X);
    virtual ~SensorNStretch   ();
    virtual void Init       (DNASeq *);
    virtual void GiveInfo   (DNASeq *, int, DATA *);
    virtual void Plot       (DNASeq *);
    virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorNStretch * builder0(int n, DNASeq *X)
{
    return new SensorNStretch(n, X);
}

#endif

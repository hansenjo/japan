/**********************************************************\
* File: QwFakeHelicity.h                                  *
*                                                         *
* Author: B.Waidyawansa                                   *
* Time-stamp:03-06-2010                                   *
\**********************************************************/
/**
  The QwFakeHelicity class uses a pregenerated random seed to generate
  the fake helicity signal that then can be used to perfrom helicity
  related calculations.
*/


#ifndef QwFAKEHELICITY_H
#define QwFAKEHELICITY_H

#include "QwHelicity.h"

class QwFakeHelicity: public QwHelicity {
 public:
  explicit QwFakeHelicity(const TString& region_tmp)
    : VQwSubsystem(region_tmp)
    , QwHelicity(region_tmp)
    , fMinPatternPhase(1)

    {
      // using the constructor of the base class
    };

    virtual ~QwFakeHelicity() = default;

    void    ClearEventData();
    Bool_t  IsGoodHelicity();
    void    ProcessEvent();

    Bool_t CheckForBurpFail(const VQwSubsystem *subsys){
      return kFALSE;
    };

 protected:
    Int_t fMinPatternPhase;

    Bool_t CollectRandBits();
    UInt_t GetRandbit(UInt_t& ranseed);

};

#endif

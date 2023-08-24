/**********************************************************\
* File: QwBeamLine.C                                      *
*                                                         *
* Author:                                                 *
* Time-stamp:                                             *
\**********************************************************/

#include "QwBeamLine.h"

// System headers
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cassert>

// Qweak headers
#include "QwLog.h"

#ifdef __USE_DATABASE__
#define MYSQLPP_SSQLS_NO_STATICS
#include "QwParitySSQLS.h"
#include "QwParityDB.h"
#endif // __USE_DATABASE__

#include "QwPromptSummary.h"

// Forward declarations
#ifdef __USE_DATABASE__
class QwParityDB;
//class QwDBInterface;
#endif // __USE_DATABASE__

// Register this subsystem with the factory
RegisterSubsystemFactory(QwBeamLine);

// Helper functors
namespace {

#ifndef ALL
#define ALL( c ) (c).begin(), (c).end()
#endif

// Pretty much all the helper functors  below should be made member functions
// in a common base class
template<class InputIt, class InputIt2, class UnaryFunction>
constexpr UnaryFunction QwForEachInBoth( InputIt first, InputIt last,
                                         InputIt2 first2, InputIt2 last2,
                                         UnaryFunction f ) {
  assert(std::distance(first, last) == std::distance(first2, last2));
  for( ; first != last && first2 != last2; ++first, ++first2 )
    f(*first, *first2);

  return f;
}

struct QwSetCutMode {
  explicit QwSetCutMode( Int_t mode ) : mode_(mode) {}
  template<class T>
  void operator()( boost::shared_ptr<T>& item ) { item->SetEventCutMode(mode_); }
  template<class T>
  void operator()( T& item ) { item.SetEventCutMode(mode_); }
private:
  Int_t mode_;
};

struct QwApplySngEvCut {
  explicit QwApplySngEvCut( Bool_t& status, Bool_t debug, TString name )
    : status_{status}, name_{std::move(name)}, i_{0} {}

  template<class T>
  void PrintStatus( const T& item ) const {
    std::cout << "******* QwBeamLine::SingleEventCuts()->" << name_ <<
              "[ " << i_ << " , "
              << item.GetElementName() << " ] ******\n";
  }
  template<class T>
  void operator()( boost::shared_ptr<T>& item ) {
    status_ &= item->ApplySingleEventCuts();
    PrintStatus(*item.get());
    ++i_;
  }
  template<class T>
  void operator()( T& item ) {
    status_ &= item.ApplySingleEventCuts();
    PrintStatus(item);
    ++i_;
  }
private:
  Bool_t& status_;
  TString name_;
  size_t  i_;
};

struct QwScale {
  explicit QwScale( Double_t factor ) : factor_(factor) {}
  template<class T>
  void operator()( boost::shared_ptr<T>& item ) { item->Scale(factor_); }
  template<class T>
  void operator()( T& item ) { item.Scale(factor_); }
private:
  Double_t factor_;
};

struct QwCalcRunAvg {
  QwCalcRunAvg() = default;
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) { ptr->CalculateRunningAverage(); }
  template<class T>
  void operator()( T& item ) { item.CalculateRunningAverage(); }
};

struct QwPrintValue {
  QwPrintValue() = default;
  template<class T>
  void operator()( const boost::shared_ptr<T>& ptr ) { ptr->PrintValue(); }
  template<class T>
  void operator()( const T& item ) { item.PrintValue(); }
};

struct QwConstructHistograms {
  QwConstructHistograms( TDirectory* folder, TString prefix )
    : folder_{folder}, prefix_{std::move(prefix)} {}
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) {
    ptr->ConstructHistograms(folder_, prefix_);
  }
  template<class T>
  void operator()( T& a ) { a.ConstructHistograms(folder_, prefix_); }
private:
  TDirectory* folder_;
  TString prefix_;
};

struct QwFillHistograms {
  QwFillHistograms() = default;
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) { ptr->FillHistograms(); }
  template<class T>
  void operator()( T& item ) { item.FillHistograms(); }
};

struct QwConstructBranchAndVector {
  QwConstructBranchAndVector( TTree* tree, TString prefix, std::vector<Double_t>& values )
    : tree_{tree}, prefix_{std::move(prefix)}, values_{values} {}
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) {
    ptr->ConstructBranchAndVector(tree_, prefix_, values_);
  }
  template<class T>
  void operator()( T& a ) { a.ConstructBranchAndVector(tree_, prefix_, values_); }
private:
  TTree* tree_;
  TString prefix_;
  std::vector<Double_t>& values_;
};

struct QwConstructBranch {
  QwConstructBranch( TTree* tree, TString prefix )
    : tree_{tree}, prefix_{std::move(prefix)} {}
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) {
    ptr->ConstructBranch(tree_, prefix_);
  }
  template<class T>
  void operator()( T& a ) { a.ConstructBranch(tree_, prefix_); }
private:
  TTree* tree_;
  TString prefix_;
};

struct QwFillTreeVector {
  explicit QwFillTreeVector( std::vector<Double_t>& values )
    : values_{values} {}
  template<class T>
  void operator()( const boost::shared_ptr<T>& ptr ) {
    ptr->FillTreeVector(values_);
  }
  template<class T>
  void operator()( const T& a ) { a.FillTreeVector(values_); }
private:
  std::vector<Double_t>& values_;
};

struct QwClearEventData {
  QwClearEventData() = default;
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) { ptr->ClearEventData(); }
  template<class T>
  void operator()( T& item ) { item.ClearEventData(); }
};

struct QwProcessEvent {
  QwProcessEvent() = default;
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) { ptr->ProcessEvent(); }
  template<class T>
  void operator()( T& item ) { item.ProcessEvent(); }
};

struct QwProcessEventAndPrint {
  QwProcessEventAndPrint() = default;
  template<class T>
  void operator()( boost::shared_ptr<T>& ptr ) {
    ptr->ProcessEvent();
    ptr->PrintInfo();
  }
  template<class T>
  void operator()( T& item ) {
    item.ProcessEvent();
    item.PrintInfo();
  }
};

struct QwCheckForBurpFail {
  explicit QwCheckForBurpFail( Bool_t& status ) : status_{status} {}
  template<class T>
  void operator()( boost::shared_ptr<T>& pa, const boost::shared_ptr<T>& pb ) const {
    status_ |= pa->CheckForBurpFail(pb.get());
  }
  template<class T>
  void operator()( T& a, const T& b ) const { status_ |= a.CheckForBurpFail(&b); }
private:
  Bool_t& status_;
};

struct QwAccumulate {
  QwAccumulate( Int_t count, Int_t mask ) : count_{count}, mask_{mask} {}
  template<class T>
  void operator()( boost::shared_ptr<T>& pa, boost::shared_ptr<T>& pb ) {
    pa->AccumulateRunningSum(*pb, count_, mask_);
  }
  template<class T>
  void operator()( T& a, T& b ) { a.AccumulateRunningSum(b, count_, mask_); }
private:
  Int_t count_;
  Int_t mask_;
};

struct QwDeaccumulate {
  explicit QwDeaccumulate( Int_t mask ) : mask_{mask} {}
  template<class T>
  void operator()( boost::shared_ptr<T>& pa, boost::shared_ptr<T>& pb ) {
    pa->DeaccumulateRunningSum(*pb, mask_);
  }
  template<class T>
  void operator()( T& a, T& b ) { a.DeaccumulateRunningSum(b, mask_); }
private:
  Int_t mask_;
};

} // namespace

//*****************************************************************//
void QwBeamLine::ProcessOptions( QwOptions& options ) {
  //Handle command line options
}

//*****************************************************************//
template<>
Int_t QwBeamLine::AddToElementList<VQwClock_ptr>( std::vector<VQwClock_ptr>& elementlist,
                                                  QwBeamDetectorID& detector_id ) {
  if( detector_id.fTypeID == kQwClock ) {
    VQwClock_ptr element(
      VQwClock::Create(GetName(),
                       detector_id.fdetectorname,
                       detector_id.fmoduletype));
    elementlist.push_back(element);
  }
  detector_id.fIndex = elementlist.size() - 1;
  return detector_id.fIndex;
}
//*****************************************************************//
template<>
Int_t QwBeamLine::AddToElementList<VQwBCM_ptr>( std::vector<VQwBCM_ptr>& elementlist,
                                                QwBeamDetectorID& detector_id ) {
  if( detector_id.fTypeID == kQwBCM ) {
    VQwBCM_ptr element(
      VQwBCM::Create(GetName(),
                     detector_id.fdetectorname,
                     detector_id.fmoduletype));
    elementlist.push_back(element);
  }
  detector_id.fIndex = elementlist.size() - 1;
  return detector_id.fIndex;
}
//*****************************************************************//
template<>
Int_t QwBeamLine::AddToElementList<VQwBPM_ptr>( std::vector<VQwBPM_ptr>& elementlist,
                                                QwBeamDetectorID& detector_id ) {
  if( detector_id.fTypeID == kQwBPMStripline ) {
    VQwBPM_ptr element(
      VQwBPM::CreateStripline(GetName(),
                              detector_id.fdetectorname,
                              detector_id.fmoduletype));
    elementlist.push_back(element);
  }
  detector_id.fIndex = elementlist.size() - 1;
  return detector_id.fIndex;
}

//*****************************************************************//
template<typename TT>
Int_t QwBeamLine::AddToElementList( std::vector<TT>& elementlist,
                                    QwBeamDetectorID& detector_id ) {
  TT element(GetName(), detector_id.fdetectorname);
  elementlist.push_back(element);
  detector_id.fIndex = elementlist.size() - 1;
  return detector_id.fIndex;
}


//*****************************************************************//
Int_t QwBeamLine::LoadChannelMap( const TString& mapfile ) {
  Bool_t ldebug = kFALSE;

  TString varname, varvalue;

  TString combotype, comboname, dev_name;

  Int_t index = 0;
  Bool_t combolistdecoded;
  Bool_t deviceok;

  std::vector<TString> fDeviceName;
  std::vector<TString> fProperty;
  std::vector<TString> fType;
  std::vector<Double_t> fQWeight;
  std::vector<Double_t> fXWeight;
  std::vector<Double_t> fYWeight;
  Double_t sumQweights = 0.0;

  std::vector<QwBeamDetectorID> clock_needed_list;

  QwParameterFile mapstr(mapfile.Data());  //Open the file
  fDetectorMaps.insert(mapstr.GetParamFileNameContents());
  mapstr.EnableGreediness();
  mapstr.SetCommentChars("!");
  mapstr.AddBreakpointKeyword("begin");
  mapstr.AddBreakpointKeyword("end");

  Int_t buffer_offset = 0; /*  Allow some extra words at the start of a bank.
			    *  The buffer_offset value will be reset at the
			    *  start of each ROC or bank declaration, so should
			    *  be relisted for each bank.
			    */
  while( mapstr.ReadNextLine() && mapstr.SkipSection("PUBLISH") ) {
    RegisterRocBankMarker(mapstr);
    //  Remove the "vqwk_buffer_offset" and "scaler_buffer_offset"
    //  keywords from the parameter file's listing
    if( mapstr.ReturnValue("vqwk_buffer_offset", buffer_offset) ) {
      QwDebug << "QwBeamLine::LoadChannelMap: "
              << "ROC " << fCurrentROC_ID
              << ", Bank " << fCurrentBank_ID
              << "; vqwk_buffer_offset:  " << buffer_offset << QwLog::endl;
    }

    if( mapstr.HasVariablePair("=", varname, varvalue) ) { //  This is a declaration line.  Decode it.
      varname.ToLower();

      if( varname == "begin" ) {

        // Start to decode derived beamline devices combined+energy
        deviceok = kTRUE;
        combotype = varvalue;
        combolistdecoded = kFALSE;

        TString dettype;

        while( mapstr.ReadNextLine() && !combolistdecoded ) {
          if( mapstr.HasVariablePair("=", varname, varvalue) ) {
            varname.ToLower();
            if( varname == "end" ) {
              // calculate the total weights of the charge
              sumQweights = 0.0;
              for( size_t i = 0; i < fDeviceName.size(); i++ )
                sumQweights += abs(fQWeight[i]);
              combolistdecoded = kTRUE;
              break;
            }
          }

          if( mapstr.PopValue("name", varvalue) ) {
            comboname = varvalue;
          }

          dev_name = mapstr.GetTypedNextToken<TString>();
          dev_name.ToLower();
          dettype = mapstr.GetTypedNextToken<TString>();
          dettype.ToLower();

          // Check to see if the device being read is a valid physical device.
          // If not, discard the combination.
          index = GetDetectorIndex(GetQwBeamInstrumentType(dettype), dev_name);

          if( index == -1 ) {
            QwError << "QwBeamLine::LoadChannelMap:  Unknown device: "
                    << dev_name << " used in " << comboname
                    << ". This combination  will not be decoded!"
                    << QwLog::endl;
            deviceok = kFALSE;
            combolistdecoded = kTRUE;
          } else {
            // Found the device
            // Add to the array of names
            fDeviceName.push_back(dev_name);

            // Read in the weights.
            // For combined bpms and combined bcms these are charge weights.
            // For the energy calculator these are the ratios of the transport matrix elements.
            fQWeight.push_back(mapstr.GetTypedNextToken<Double_t>());

            // For combined BPMs,in addition, there are weights for the X & Y positions.
            if( combotype == "combinedbpm" ) {
              fXWeight.push_back(mapstr.GetTypedNextToken<Double_t>());
              fYWeight.push_back(mapstr.GetTypedNextToken<Double_t>());
            }

            // For the enrgy calculator there are device type and the specified beam parameters.
            if( combotype == "energycalculator" ) {
              fType.push_back(dettype);
              fProperty.push_back(mapstr.GetTypedNextToken<TString>());
            }
          }
        }

        // Now create the combined device
        QwBeamDetectorID localComboID(-1, -1, comboname, combotype,
                                      fBeamDetectorID.at(index).fmoduletype);

        localComboID.fdetectorname = comboname(0, comboname.Sizeof() - 1);
        localComboID.fIndex = GetDetectorIndex(localComboID.fTypeID, localComboID.fdetectorname);

        if( localComboID.fTypeID == kQwUnknownDeviceType ) {
          QwError << "QwBeamLine::LoadChannelMap:  Unknown detector type: "
                  << combotype << ", the detector " << comboname << " will not be decoded "
                  << QwLog::endl;
          deviceok = kFALSE;
          continue;
        }

        if( (localComboID.fIndex == -1) && deviceok ) {

          // Decoding combined BCM array
          if( localComboID.fTypeID == kQwCombinedBCM ) {

            VQwBCM_ptr localbcmcombo(
              VQwBCM::CreateCombo(GetName(),
                                  localComboID.fdetectorname, localComboID.fmoduletype));
            fBCMCombo.push_back(localbcmcombo);

            for( size_t i = 0; i < fDeviceName.size(); i++ ) {
              index = GetDetectorIndex(GetQwBeamInstrumentType(dettype), fDeviceName[i]);
              fBCMCombo[fBCMCombo.size() - 1].get()->SetBCMForCombo(fBCM.at(index).get(),
                                                                    fQWeight[i], sumQweights);
            }
            fDeviceName.clear();
            fQWeight.clear();
            localComboID.fIndex = fBCMCombo.size() - 1;
          }

          // Decoding combined BPM array.
          if( localComboID.fTypeID == kQwCombinedBPM ) {
            VQwBPM_ptr localbpmcombo(
              VQwBPM::CreateCombo(GetName(),
                                  localComboID.fdetectorname, localComboID.fmoduletype));
            fBPMCombo.push_back(localbpmcombo);

            for( size_t i = 0; i < fDeviceName.size(); i++ ) {
              index = GetDetectorIndex(GetQwBeamInstrumentType(dettype), fDeviceName[i]);
              fBPMCombo[fBPMCombo.size() - 1].get()->SetBPMForCombo(
                fStripline.at(index).get(), fQWeight[i], fXWeight[i],
                fYWeight[i], sumQweights);

            }
            fDeviceName.clear();
            fQWeight.clear();
            fXWeight.clear();
            fYWeight.clear();
            localComboID.fIndex = fBPMCombo.size() - 1;
          }

          // Decoding energy calculator.
          if( localComboID.fTypeID == kQwEnergyCalculator ) {

            QwEnergyCalculator localecalculator(GetName(), localComboID.fdetectorname);
            fECalculator.push_back(localecalculator);

            for( size_t i = 0; i < fDeviceName.size(); i++ ) {
              index = GetDetectorIndex(GetQwBeamInstrumentType(fType[i]), fDeviceName[i]);

              if( GetQwBeamInstrumentType(fType[i]) == kQwBPMStripline )
                fECalculator[fECalculator.size() - 1].Set(fStripline.at(index).get(), fType[i], fProperty[i],
                                                          fQWeight[i]);

              if( GetQwBeamInstrumentType(fType[i]) == kQwCombinedBPM )
                fECalculator[fECalculator.size() - 1].Set(fBPMCombo.at(index).get(), fType[i], fProperty[i],
                                                          fQWeight[i]);

            }

            fDeviceName.clear();
            fQWeight.clear();

            fProperty.clear();
            fType.clear();
            localComboID.fIndex = fECalculator.size() - 1;
          }
        }
        // Use only the combinations that are of known type and has known physical devices.
        if( deviceok )
          fBeamDetectorID.push_back(localComboID);
      }

      QwDebug << "At end of processing the combined device " << QwLog::endl;

    } else {
      // Start to decode the physical beamline devices
      QwBeamDetectorID localBeamDetectorID(GetSubbankIndex(), mapstr);
      Bool_t lineok = localBeamDetectorID.ReportInitErrors();
      if( !lineok ) continue;

      localBeamDetectorID.fIndex =
        GetDetectorIndex(localBeamDetectorID.fTypeID,
                         localBeamDetectorID.fdetectorname);

      if( localBeamDetectorID.fIndex == -1 ) {
        Int_t index;
        VQwDataElement* local_element = nullptr;

        switch( localBeamDetectorID.fTypeID ) {
          case kQwHaloMonitor:
            index = AddToElementList(fHaloMonitor, localBeamDetectorID);
            local_element = &(fHaloMonitor.at(index));
            break;
          case kQwBPMCavity:
            index = AddToElementList(fCavity, localBeamDetectorID);
            local_element = &(fCavity.at(index));
            break;
          case kQwBPMStripline:
            index = AddToElementList(fStripline, localBeamDetectorID);
            local_element = fStripline.at(index).get();
            break;
          case kQwBCM:
            index = AddToElementList(fBCM, localBeamDetectorID);
            local_element = fBCM.at(index).get();
            break;
          case kQwClock:
            index = AddToElementList(fClock, localBeamDetectorID);
            local_element = fClock.at(index).get();
            break;
          case kQwQPD:
            index = AddToElementList(fQPD, localBeamDetectorID);
            local_element = &(fQPD.at(index));
            break;
          case kQwLinearArray:
            index = AddToElementList(fLinearArray, localBeamDetectorID);
            local_element = &(fLinearArray.at(index));
            break;
          default:
            // TODO: can this happen?
            break;
        }
        if( !local_element )
          // TODO: report error?
          continue;

        local_element->LoadChannelParameters(mapstr);
        if( local_element->NeedsExternalClock() ) {
          QwDebug << "Try to push device "
                  << local_element->GetElementName()
                  << " (address=" << std::hex << local_element
                  << ") onto the clock_needed_list"
                  << QwLog::endl;
          clock_needed_list.push_back(localBeamDetectorID);
        }
      }

      fBeamDetectorID.push_back(localBeamDetectorID);
    }
  }

  // Now load the variables to publish
  mapstr.RewindToFileStart();
  QwParameterFile* section;
  std::vector<TString> publishinfo;
  while( (section = mapstr.ReadNextSection(varvalue)) ) {
    if( varvalue == "PUBLISH" ) {
      fPublishList.clear();
      while( section->ReadNextLine() ) {
        section->TrimComment(); // Remove everything after a comment character
        section->TrimWhitespace(); // Get rid of leading and trailing spaces
        for( int ii = 0; ii < 4; ii++ ) {
          varvalue = section->GetTypedNextToken<TString>();
          if( varvalue.Length() ) {
            publishinfo.push_back(varvalue);
          }
        }
        if( publishinfo.size() == 4 )
          fPublishList.push_back(publishinfo);
        publishinfo.clear();
      }
    }
    delete section;
  }
  // Print list of variables to publish
  if( !fPublishList.empty() ) {
    QwMessage << "Variables to publish:" << QwLog::endl;
    for( const auto& pubinfo: fPublishList ) {
      assert(pubinfo.size() == 4);
      QwMessage << pubinfo[0] << " " << pubinfo[1] << " "
                << pubinfo[2] << " " << pubinfo[3] << QwLog::endl;
    }
  }

  if( ldebug ) {
    std::cout << "QwBeamLine::Done with Load map channel \n";
    for( const auto& detectorId: fBeamDetectorID )
      detectorId.Print();
  }

  // Now propagate clock pointers to those channels that need it
  index = 0;
  VQwDataElement* local_element;
  std::string clockname;
  for( const auto& detectorId: clock_needed_list ) {
    local_element = GetElement(detectorId);
    if( !local_element ) continue;
    clockname = local_element->GetExternalClockName();
    if( clockname.empty() ) {
      QwWarning << "QwBeamLine::LoadChannelMap  "
                << "Device, " << local_element->GetElementName()
                << " needs a reference clock, but the reference clock name is empty"
                << QwLog::endl;
    } else {
      index = GetDetectorIndex(kQwClock, clockname);
      if( index >= 0 ) {
        QwMessage << "QwBeamLine::LoadChannelMap  "
                  << "Setting " << fClock.at(index).get()->GetElementName()
                  << " as the reference clock for channel "
                  << local_element->GetElementName()
                  << QwLog::endl;
        local_element->SetExternalClockPtr(fClock.at(index).get()->GetTime());
      } else {
        QwWarning << "QwBeamLine::LoadChannelMap  "
                  << "Cannot find clock, " << local_element->GetExternalClockName()
                  << ", needed by device, " << local_element->GetElementName()
                  << QwLog::endl;
      }
    }
  }
  ldebug = kFALSE;

  mapstr.Close(); // Close the file (ifstream)
  return 0;
}


//*****************************************************************//
Int_t QwBeamLine::LoadEventCuts( const TString& filename ) {
  Int_t eventcut_flag = 1;

  // Open the file
  QwParameterFile mapstr(filename.Data());
  fDetectorMaps.insert(mapstr.GetParamFileNameContents());
  while( mapstr.ReadNextLine() ) {
    mapstr.TrimComment('!');   // Remove everything after a '!' character.
    mapstr.TrimWhitespace();   // Get rid of leading and trailing spaces.
    if( mapstr.LineIsEmpty() ) { continue; }

    TString varname, varvalue;
    if( mapstr.HasVariablePair("=", varname, varvalue) ) {
      if( varname == "EVENTCUTS" ) {
        eventcut_flag = QwParameterFile::GetUInt(varvalue);
      }
    } else {
      auto device_type = mapstr.GetTypedNextToken<TString>();
      device_type.ToLower();
      auto device_name = mapstr.GetTypedNextToken<TString>();
      device_name.ToLower();

      Int_t det_index = GetDetectorIndex(GetQwBeamInstrumentType(device_type), device_name);
      if( det_index == -1 ) {
        QwWarning << " Device not found " << device_name << " of type " << device_type << QwLog::endl;
        continue;
      }

      if( device_type == GetQwBeamInstrumentTypeName(kQwBCM) ) {
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for BCM value
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for BCM value
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwBCM " << errorflag
                  << ", burp = " << burplevel << QwLog::endl;
        fBCM[det_index]->SetSingleEventCuts(errorflag, LLX, ULX,
                                            stabilitycut, burplevel);//(fBCMEventCuts);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwHaloMonitor) ) {
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for HaloMonitor value
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for HaloMonitor value
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwHaloMonitor " << errorflag << QwLog::endl;
        fHaloMonitor[det_index].SetSingleEventCuts(errorflag, LLX, ULX,
                                                   stabilitycut, burplevel);//(fBCMEventCuts);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwEnergyCalculator) ) {
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for energy
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for energy
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwEnergyCalculator " << errorflag << QwLog::endl;
        fECalculator[det_index].SetSingleEventCuts(errorflag, LLX, ULX,
                                                   stabilitycut, burplevel);//(fEnergyEventCuts);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwBPMStripline) ) {
        auto channel_name = mapstr.GetTypedNextToken<TString>();
        channel_name.ToLower();
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for BPMStripline X
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for BPMStripline X
        varvalue = mapstr.GetTypedNextToken<TString>();//global/local
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine:QwBPMStripline " << channel_name << " " << varvalue << " " << stabilitycut
                  << QwLog::endl;
        //QwMessage<<"QwBeamLine Error Code passing to QwBPMStripline "<<GetGlobalErrorFlag(varvalue,eventcut_flag,stabilitycut)<<" stability  "<<stabilitycut <<QwLog::endl;
        fStripline[det_index]->SetSingleEventCuts(channel_name, errorflag,
                                                  LLX, ULX, stabilitycut, burplevel);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwQPD) ) {
        auto channel_name = mapstr.GetTypedNextToken<TString>();
        channel_name.ToLower();
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for QPD X
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for QPD X
        varvalue = mapstr.GetTypedNextToken<TString>();//global/local
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwQPD " << errorflag << QwLog::endl;
        fQPD[det_index].SetSingleEventCuts(channel_name, errorflag,
                                           LLX, ULX, stabilitycut, burplevel);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwLinearArray) ) {
        auto channel_name = mapstr.GetTypedNextToken<TString>();
        channel_name.ToLower();
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for LinearArray X
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for LinearArray X
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwLinearArray " << errorflag << QwLog::endl;
        fLinearArray[det_index].SetSingleEventCuts(channel_name, errorflag,
                                                   LLX, ULX, stabilitycut, burplevel);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwBPMCavity) ) {
        auto channel_name = mapstr.GetTypedNextToken<TString>();
        channel_name.ToLower();
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for cavity bpm X
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for cavity bpm X
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwBPMCavity " << errorflag
                  << " " << det_index << QwLog::endl;
        fCavity[det_index].SetSingleEventCuts(channel_name, errorflag,
                                              LLX, ULX, stabilitycut, burplevel);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwCombinedBCM) ) {
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for BCM value
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for BCM value
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwCombinedBCM " << errorflag << QwLog::endl;
        fBCMCombo[det_index]->PrintInfo();
        fBCMCombo[det_index]->SetSingleEventCuts(errorflag, LLX, ULX,
                                                 stabilitycut, burplevel);
      } else if( device_type == GetQwBeamInstrumentTypeName(kQwCombinedBPM) ) {
        auto channel_name = mapstr.GetTypedNextToken<TString>();
        channel_name.ToLower();
        auto LLX = mapstr.GetTypedNextToken<Double_t>();  //lower limit for combined bpm X
        auto ULX = mapstr.GetTypedNextToken<Double_t>();  //upper limit for combined bpm X
        varvalue = mapstr.GetTypedNextToken<TString>();//global/loacal
        auto stabilitycut = mapstr.GetTypedNextToken<Double_t>();
        auto burplevel = mapstr.GetTypedNextToken<Double_t>();
        varvalue.ToLower();
        auto errorflag = GetGlobalErrorFlag(varvalue, eventcut_flag, stabilitycut);
        QwMessage << "QwBeamLine Error Code passing to QwCombinedBPM " << errorflag << QwLog::endl;
        fBPMCombo[det_index]->SetSingleEventCuts(channel_name, errorflag,
                                                 LLX, ULX, stabilitycut, burplevel);
      }
    }
  }


  //update the event cut ON/OFF for all the devices
  //std::cout<<"EVENT CUT FLAG"<<eventcut_flag<<std::endl;

  QwSetCutMode setmode(eventcut_flag);

  for_each(ALL(fStripline),   setmode);
  for_each(ALL(fQPD),         setmode);
  for_each(ALL(fLinearArray), setmode);
  for_each(ALL(fCavity),      setmode);
  for_each(ALL(fBCM),         setmode);
  for_each(ALL(fClock),       setmode);
  for_each(ALL(fHaloMonitor), setmode);
  for_each(ALL(fBCMCombo),    setmode);
  for_each(ALL(fBPMCombo),    setmode);
  for_each(ALL(fECalculator), setmode);

  fQwBeamLineErrorCount = 0; //set the error counter to zero

  mapstr.Close(); // Close the file (ifstream)
  return 0;
}


//*****************************************************************//
Int_t QwBeamLine::LoadGeometryDefinition( const TString& mapfile )
{
  Bool_t ldebug = kFALSE;
  TString varname, varvalue;
  Int_t lineread = 1;
  Int_t index = 0;
  TString devname, devtype;
  TString melement;
  Double_t devOffsetX = 0, devOffsetY = 0, devOffsetZ = 0;
  Double_t devSENfactor = 0, devAlphaX = 0, devAlphaY = 0;
  TString rotation_stat;
  VQwBPM* bpm = nullptr;

  if( ldebug )std::cout << "QwBeamLine::LoadGeometryParameters(" << mapfile << ")\n";

  QwParameterFile mapstr(mapfile.Data());  //Open the file
  fDetectorMaps.insert(mapstr.GetParamFileNameContents());
  while( mapstr.ReadNextLine() ) {
    lineread += 1;
    if( ldebug )std::cout << " line read so far =" << lineread << "\n";
    mapstr.TrimComment('!');
    mapstr.TrimWhitespace();

    if( mapstr.LineIsEmpty() ) continue;

    Bool_t notfound = kTRUE;

    devtype = mapstr.GetTypedNextToken<TString>();
    devtype.ToLower();
    devtype.Remove(TString::kBoth, ' ');
    devname = mapstr.GetTypedNextToken<TString>();
    devname.ToLower();
    devname.Remove(TString::kBoth, ' ');

    if( GetQwBeamInstrumentType(devtype) == kQwUnknownDeviceType ) {
      QwError << "Error! Unknown detector type '" << devtype
              << "' in Geometry file!" << QwLog::endl;
      /*If the device type is unknown there is no point in going through the rest of the specs for that device*/
      continue;
    }

    index = GetDetectorIndex(GetQwBeamInstrumentType(devtype), devname);
    if( index < 0 ) {
      /*Detector name isn't recognized. Ignore it!*/
      QwWarning << "Unrecognized detector name '" << devname
                << "' in Geometry file.  This may not be a problem, "
                << "if we're using a reduced channel map."
                << QwLog::endl;
    } else {
      devOffsetX = mapstr.GetTypedNextToken<Double_t>(); // X offset
      devOffsetY = mapstr.GetTypedNextToken<Double_t>(); // Y offset
      devOffsetZ = mapstr.GetTypedNextToken<Double_t>(); // Z offset
      devSENfactor = mapstr.GetTypedNextToken<Double_t>(); // sensivity scaling factor
      devAlphaX = mapstr.GetTypedNextToken<Double_t>(); // alpha X
      devAlphaY = mapstr.GetTypedNextToken<Double_t>(); // alpha Y


      /*If the device is a bpm stripline, assign the rotations and gains*/
      if( GetQwBeamInstrumentType(devtype) == kQwBPMStripline ) {
        bpm = fStripline.at(index).get();
        AssignGeometry(&mapstr, bpm);
      }


      if( ldebug == 1 ) {
        std::cout << "####################\n";
        std::cout << "! device type, device_name, Xoffset, Yoffset, Zoffset, BSEN scaling factor, AlpaX, AlpaY\n"
                  << std::endl;
        std::cout << GetQwBeamInstrumentType(devtype) << " / "
                  << devname << " / "
                  << devOffsetX << " / "
                  << devOffsetY << " / "
                  << devOffsetZ << " / "
                  << devSENfactor << " / "
                  << devAlphaX << " / "
                  << devAlphaY << " / "
                  << std::endl;
      }


      while( notfound ) {
        if( GetQwBeamInstrumentType(devtype) == kQwBPMStripline ) {
          //Load bpm offsets
          if( index == -1 ) {
            QwWarning << "QwBeamLine::LoadGeometryDefinition:  Unknown bpm in qweak_beamline_geometry.map: "
                      << devname
                      << QwLog::endl;
            notfound = kFALSE;
            continue;
          }

          TString localname = fStripline.at(index).get()->GetElementName();
          localname.ToLower();
          if( ldebug )
            std::cout << "element name ==" << localname
                      << "== to be compared to ==" << devname << "== \n";

          if( localname == devname ) {
            if( ldebug ) std::cout << " I found the bpm !\n";
            bpm->GetSurveyOffsets(devOffsetX, devOffsetY, devOffsetZ);
            bpm->GetElectronicFactors(devSENfactor, devAlphaX, devAlphaY);

            // If nothing is specified, a default rotation of 45 degrees is implied.
            notfound = kFALSE;
          }
        } else if( GetQwBeamInstrumentType(devtype) == kQwCombinedBPM ) {
          //Load combined bpm offsets which are, ofcourse, target position in the beamline
          if( index == -1 ) {
            QwError << "QwBeamLine::LoadGeometryDefinition:  Unknown combined bpm in qweak_beamline_geometry.map: "
                    << devname << " Check the combined bpm names!\n "
                    << QwLog::endl;
            notfound = kFALSE;
            continue;
          }

          TString localname = fBPMCombo.at(index).get()->GetElementName();
          localname.ToLower();
          if( ldebug )
            std::cout << "element name ==" << localname << "== to be compared to ==" << devname << "== \n";

          if( localname == devname ) {
            if( ldebug ) std::cout << " I found the combinedbpm !\n";
            fBPMCombo.at(index).get()->GetSurveyOffsets(devOffsetX, devOffsetY, devOffsetZ);
            notfound = kFALSE;
          }
        } else if( GetQwBeamInstrumentType(devtype) == kQwBPMCavity ) {
          //Load cavity bpm offsets
          if( index == -1 ) {
            QwError << "QwBeamLine::LoadGeometryDefinition:  Unknown bpm : "
                    << devname << " will not be asigned with geometry parameters. \n"
                    << QwLog::endl;
            notfound = kFALSE;
            continue;
          }
          TString localname = fCavity.at(index).GetElementName();
          localname.ToLower();
          if( ldebug )
            std::cout << "element name ==" << localname
                      << "== to be compared to ==" << devname << "== \n";

          if( localname == devname ) {
            if( ldebug ) std::cout << " I found the cavity bpm !\n";
            fCavity.at(index).GetSurveyOffsets(devOffsetX, devOffsetY, devOffsetZ);
            notfound = kFALSE;
          }
        } else if( GetQwBeamInstrumentType(devtype) == kQwQPD ) {
          //Load QPD calibration factors
          if( index == -1 ) {
            QwError << "QwBeamLine::LoadGeometryDefinition:  Unknown QPD : "
                    << devname << " will not be asigned with calibration factors. \n"
                    << QwLog::endl;
            notfound = kFALSE;
            continue;
          }
          TString localname = fQPD.at(index).GetElementName();
          localname.ToLower();
          if( ldebug )
            std::cout << "element name ==" << localname
                      << "== to be compared to ==" << devname << "== \n";

          if( localname == devname ) {
            if( ldebug ) std::cout << "I found the QPD !\n";
            fQPD.at(index).GetCalibrationFactors(devAlphaX, devAlphaY);
            notfound = kFALSE;
          }
        } else
          QwError << "QwBeamLine::LoadGeometryDefinition: Unknown device type : " << devtype <<
                  ". Are you sure we have this in the beamline? I am skipping this." << QwLog::endl;
      }
    }
  }

  if( ldebug ) std::cout << "line read in the geometry file = " << lineread << " \n";

  ldebug = kFALSE;
  mapstr.Close(); // Close the file (ifstream)
  return 0;

}

//--------------------------------------------------------------------------------------------------------
Int_t QwBeamLine::LoadMockDataParameters( const TString& mapfile )
{
  Bool_t ldebug = kFALSE;
  TString varname, varvalue;
  Int_t lineread = 1;
  Int_t index = 0;
  TString devname, devtype;
  TString melement;

  if( ldebug ) std::cout << "QwBeamLine::LoadMockDataParameters(" << mapfile << ") \n" << std::endl;

  QwParameterFile mapstr(mapfile.Data());  //Open the file
  fDetectorMaps.insert(mapstr.GetParamFileNameContents());

  while( mapstr.ReadNextLine() ) {
    lineread += 1;
    //    std::cerr << "Line:  " << mapstr.GetLine() << std::endl;
    if( ldebug ) std::cout << "Line read so far = " << lineread << "\n" << std::endl;
    mapstr.TrimComment('!');
    mapstr.TrimWhitespace();

    if( mapstr.LineIsEmpty() ) continue;

    devtype = mapstr.GetTypedNextToken<TString>();
    devtype.ToLower();
    devtype.Remove(TString::kBoth, ' ');
    devname = mapstr.GetTypedNextToken<TString>();
    devname.ToLower();
    devname.Remove(TString::kBoth, ' ');

    if( GetQwBeamInstrumentType(devtype) == kQwUnknownDeviceType ) {
      /*If the device type is unknown there is no point in going through the rest of the specs for that device*/
      QwError << "Error! Unknown detector type '" << devtype << "' in MockDataParameters file!" << QwLog::endl;
      continue;
    }
    index = GetDetectorIndex(GetQwBeamInstrumentType(devtype), devname);
    if( index < 0 ) {
      /*Detector name isn't recognized. Ignore it!*/
      QwWarning << "Unrecognized detector name '" << devname << "' in MockDataParameters file." << QwLog::endl;
      continue;
    }
    //  The device should process the reminder of this line.
    GetElement(GetQwBeamInstrumentType(devtype), index)->LoadMockDataParameters(mapstr);
  }
  return 0;
}

//--------------------------------------------------------------------------------------------------------

void QwBeamLine::AssignGeometry( QwParameterFile* mapstr, VQwBPM* bpm )
{
  Bool_t ldebug = kFALSE;

  TString token = "0";
  TString angle, xgain, ygain;
  Double_t rotation_angle = 0;

  while( token != "" ) {
    token = mapstr->GetTypedNextToken<TString>();
    token.Remove(TString::kBoth, '\0');

    if( token.Contains("unrotated") ) {
      if( ldebug ) std::cout << " unrotated " << std::endl;
      bpm->SetRotationOff();
    } else if( token.Contains("rotation") ) {
      // If the status is 'rotated'

      // If a specific rotation angle is given read that
      if( token.Contains("=") ) {
        angle = token.Remove(0, 9);
        rotation_angle = atof(angle);
        if( ldebug ) std::cout << "Rotation angle = " << rotation_angle << std::endl;
        bpm->SetRotation(rotation_angle);
      }
    }
    // If nothing is specified for rotation, a default rotation of 45 degrees is implied.

    if( token.Contains("xgain") ) {
      xgain = token.Remove(0, 6);
      if( ldebug ) std::cout << " xgain =" << xgain << std::endl;
      bpm->SetGains("X", atof(xgain));
    }

    if( token.Contains("ygain") ) {
      ygain = token.Remove(0, 6);
      if( ldebug ) std::cout << " ygain =" << ygain << std::endl;
      bpm->SetGains("Y", atof(ygain));
    }
  }
}

//*****************************************************************//
Int_t QwBeamLine::LoadInputParameters( const TString& pedestalfile )
{
  Bool_t ldebug = kFALSE;

  Int_t lineread = 1;

  if( ldebug )std::cout << "QwBeamLine::LoadInputParameters(" << pedestalfile << ")\n";

  QwParameterFile mapstr(pedestalfile.Data());  //Open the file
  fDetectorMaps.insert(mapstr.GetParamFileNameContents());

  while( mapstr.ReadNextLine() ) {
    lineread += 1;
    if( ldebug )std::cout << " line read so far =" << lineread << "\n";
    mapstr.TrimComment('!');   // Remove everything after a '!' character.
    mapstr.TrimWhitespace();   // Get rid of leading and trailing spaces.
    if( mapstr.LineIsEmpty() ) continue;
    else {
      auto varname = mapstr.GetTypedNextToken<TString>();  //name of the channel
      varname.ToLower();
      varname.Remove(TString::kBoth, ' ');
      auto varped = mapstr.GetTypedNextToken<Double_t>(); // value of the pedestal
      auto varcal = mapstr.GetTypedNextToken<Double_t>(); // value of the calibration factor
      /*Double_t varweight = */ mapstr.GetTypedNextToken<Double_t>(); // value of the statistical weight

      //if(ldebug) std::cout<<"inputs for channel "<<varname
      //	      <<": ped="<<varped<<": cal="<<varcal<<": weight="<<varweight<<"\n";
      Bool_t notfound = kTRUE;

      if( notfound ) {
        for( size_t i = 0; i < fStripline.size(); i++ ) {
          for( int j = 0; j < 4; j++ ) {
            TString localname = fStripline[i].get()->GetSubElementName(j);
            localname.ToLower();
            if( ldebug )
              std::cout << "Stripline element name ==" << localname
                        << "== to be compared to ==" << varname << "== \n";
            if( localname == varname ) {
              if( ldebug ) std::cout << " I found it !\n";
              fStripline[i].get()->SetSubElementPedestal(j, varped);
              fStripline[i].get()->SetSubElementCalibrationFactor(j, varcal);
              notfound = kFALSE;
              j = 5;
              i = fStripline.size() + 1;
            }
          }
        }
        for( size_t i = 0; i < fQPD.size(); i++ ) {
          for( int j = 0; j < 4; j++ ) {
            TString localname = fQPD[i].GetSubElementName(j);
            localname.ToLower();
            if( ldebug )
              std::cout << "QPD element name ==" << localname
                        << "== to be compared to ==" << varname << "== \n";
            if( localname == varname ) {
              if( ldebug ) std::cout << " I found it !\n";
              fQPD[i].SetSubElementPedestal(j, varped);
              fQPD[i].SetSubElementCalibrationFactor(j, varcal);
              notfound = kFALSE;
              j = 5;
              i = fQPD.size() + 1;
            }
          }
        }
        for( size_t i = 0; i < fLinearArray.size(); i++ ) {
          for( int j = 0; j < 8; j++ ) {
            TString localname = fLinearArray[i].GetSubElementName(j);
            localname.ToLower();
            if( ldebug )
              std::cout << "LinearArray element name ==" << localname
                        << "== to be compared to ==" << varname << "== \n";
            if( localname == varname ) {
              if( ldebug ) std::cout << " I found it !\n";
              fLinearArray[i].SetSubElementPedestal(j, varped);
              fLinearArray[i].SetSubElementCalibrationFactor(j, varcal);
              notfound = kFALSE;
              j = 9;
              i = fLinearArray.size() + 1;
            }
          }
        }
        for( size_t i = 0; i < fCavity.size(); i++ ) {
          for( size_t j = 0; j < QwBPMCavity::kNumElements; j++ ) {
            TString localname = fCavity[i].GetSubElementName(j);
            localname.ToLower();
            if( ldebug )
              std::cout << "Cavity element name ==" << localname
                        << "== to be compared to ==" << varname << "== \n";
            if( localname == varname ) {
              if( ldebug ) std::cout << " I found it !\n";
              fCavity[i].SetSubElementPedestal(j, varped);
              fCavity[i].SetSubElementCalibrationFactor(j, varcal);
              notfound = kFALSE;
              j = 3;
              i = fCavity.size() + 1;
            }
          }
        }
        for( size_t i = 0; i < fBCM.size(); i++ ) {
          if( fBCM[i].get()->GetElementName() == varname ) {
            fBCM[i].get()->SetPedestal(varped);
            fBCM[i].get()->SetCalibrationFactor(varcal);
            i = fBCM.size() + 1;
            notfound = kFALSE;
            i = fBCM.size() + 1;
          }
        }
        for( size_t i = 0; i < fClock.size(); i++ ) {
          if( fClock[i].get()->GetElementName() == varname ) {
            fClock[i].get()->SetPedestal(varped);
            fClock[i].get()->SetCalibrationFactor(varcal);
            i = fClock.size() + 1;
            notfound = kFALSE;
            i = fClock.size() + 1;
          }
        }
        for( size_t i = 0; i < fHaloMonitor.size(); i++ ) {
          if( fHaloMonitor[i].GetElementName() == varname ) {
            std::cout << varname << " I found it ! " << varcal << " ped. " << varped << "\n";
            fHaloMonitor[i].SetPedestal(varped);
            fHaloMonitor[i].SetCalibrationFactor(varcal);
            i = fHaloMonitor.size() + 1;
            notfound = kFALSE;
            i = fHaloMonitor.size() + 1;
          }
        }
      }
    }
  }
  if( ldebug ) std::cout << " line read in the pedestal + cal file =" << lineread << " \n";

  ldebug = kFALSE;
  mapstr.Close(); // Close the file (ifstream)
  return 0;
}


//*****************************************************************//
void QwBeamLine::RandomizeEventData( int helicity, double time ) {
  // Randomize all QwBPMStripline buffers
  for( const auto& stripline: fStripline ) {
    stripline->RandomizeEventData(helicity, time);
    //    stripline->PrintInfo();
  }

  for( auto& cavity: fCavity )
    cavity.RandomizeEventData(helicity, time);

  // Randomize all QwBCM buffers
  for( const auto& bcm: fBCM ) {
    bcm->RandomizeEventData(helicity, time);
    //    bcm->PrintInfo();
  }

  // Randomize all QwHaloMonitor buffers
  //for (size_t i = 0; i < fHaloMonitor.size(); i++)
  //fHaloMonitor[i].RandomizeEventData(helicity, time);

//-------------------------------------------------------------------------------------------
  for( const auto& bcmCombo: fBCMCombo ) {
    bcmCombo->RandomizeEventData(helicity, time);
    for( size_t j = 0; j < bcmCombo->GetNumberOfElements(); j++ ) {
      VQwBCM* bcm = GetBCM(bcmCombo->GetSubElementName(j));
      if( bcm ) {
        bcmCombo->GetProjectedCharge(bcm);
      }
    }
  }
//-------------------------------------------------------------------------------------------
  for( const auto& bpmCombo: fBPMCombo ) {
    bpmCombo->RandomizeEventData(helicity, time);
    //fBPMCombo[i].get()->PrintValue();
    for( size_t j = 0; j < bpmCombo->GetNumberOfElements(); j++ ) {
      VQwBPM* bpm = GetBPMStripline(bpmCombo->GetSubElementName(j));
      if( bpm ) {
        bpmCombo->GetProjectedPosition(bpm);
//  print the bpm device name, and get the x and y and z? values and print them
//      std::cout << "bpm " << bpm->GetElementName() << std::endl;
//      std::cout << "xpos= " << bpm->GetPosition(VQwBPM::kXAxis) << std::endl;
//      std::cout << "ypos= " << bpm->GetPosition(VQwBPM::kYAxis) << std::endl; // << "ypos= " <<  << "zpos= " <<  << std::endl;
        //bpm->PrintInfo();
        //  Call the new function in stripline class to fill all internal variables from the fAbsPos for the bpm object
      }
    }
  }
  //-------------------------------------------------------------------------------------------
  for( auto& calc: fECalculator ) {
    calc.RandomizeEventData(helicity, time);
    //calc.PrintValue();
    for( size_t j = 0; j < calc.GetNumberOfElements(); j++ ) {
      // std::cout << "In loop over ec subelements; device name = " << calc.GetSubElementName(j) << std::endl;
      VQwBPM* bpm = GetBPMStripline(calc.GetSubElementName(j));
      if( bpm ) {
        calc.GetProjectedPosition(bpm);
      }
    }
  }
}


//*****************************************************************//
void QwBeamLine::EncodeEventData( std::vector<UInt_t>& buffer ) {
  std::vector<UInt_t> elements;
  elements.clear();

  // Get all buffers in the order they are defined in the map file
  for( const auto& id: fBeamDetectorID ) {
    // This is a QwBCM
    if( id.fTypeID == kQwBCM ) {
      fBCM[id.fIndex]->EncodeEventData(elements);
      //std::cout << "" << fBCM[id.fIndex]->GetElementName() << std::endl;
    }
    // This is a QwBPMStripline (which has 4 entries, only process the first one)
    if( id.fTypeID == kQwBPMStripline
        && id.fSubelement == 0 ) {
      fStripline[id.fIndex]->EncodeEventData(elements);
      //  Print the HWsum absolute position values for the BPM
      //fStripline[id.fIndex]->PrintValue();
    }

    //  If this is a combined BPM, let's try to print the position and angle HWsum values
  }

  // If there is element data, generate the subbank header
  std::vector<UInt_t> subbankheader;
  std::vector<UInt_t> rocheader;
  if( !elements.empty() ) {

    // Form CODA subbank header
    subbankheader.clear();
    subbankheader.push_back(elements.size() + 1);  // subbank size
    subbankheader.push_back((fCurrentBank_ID << 16) | (0x01 << 8) | (1 & 0xff));
    // subbank tag | subbank type | event number

    // Form CODA bank/roc header
    rocheader.clear();
    rocheader.push_back(subbankheader.size() + elements.size() + 1);  // bank/roc size
    rocheader.push_back((fCurrentROC_ID << 16) | (0x10 << 8) | (1 & 0xff));
    // bank tag == ROC | bank type | event number

    // Add bank header, subbank header and element data to output buffer
    buffer.insert(buffer.end(), rocheader.begin(), rocheader.end());
    buffer.insert(buffer.end(), subbankheader.begin(), subbankheader.end());
    buffer.insert(buffer.end(), elements.begin(), elements.end());
  }
}

//*****************************************************************//
Int_t QwBeamLine::ProcessEvBuffer( const ROCID_t roc_id, const BankID_t bank_id, UInt_t* buffer, UInt_t num_words ) {
  Bool_t lkDEBUG = kFALSE;

  Int_t index = GetSubbankIndex(roc_id, bank_id);
  if( index >= 0 && num_words > 0 ) {
    //  We want to process this ROC.  Begin looping through the data.
    if( lkDEBUG )
      std::cout << "QwBeamLine::ProcessEvBuffer:  "
                << "Begin processing ROC" << roc_id
                << " and subbank " << bank_id
                << " number of words=" << num_words << std::endl;
    if( buffer[0] == 0xf0f0f0f0 && num_words % 2 == 1 ) {
      buffer++;
      if( lkDEBUG )
        std::cout << "QwBeamLine::ProcessEvBuffer:  "
                  << "Skipped padding word 0xf0f0f0f0 at beginning of buffer."
                  << std::endl;
    }

    for( auto& id: fBeamDetectorID ) {
      if( id.fSubbankIndex == index ) {

        if( id.fTypeID == kQwBPMStripline ) {
          if( lkDEBUG ) {
            std::cout << "found stripline data for " << id.fdetectorname << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fStripline[id.fIndex]->
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank,
                            id.fSubelement);
        }

        if( id.fTypeID == kQwQPD ) {
          if( lkDEBUG ) {
            std::cout << "found qpd data for " << id.fdetectorname << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fQPD[id.fIndex].
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank,
                            id.fSubelement);
        }

        if( id.fTypeID == kQwLinearArray ) {
          if( lkDEBUG ) {
            std::cout << "found linear array data for " << id.fdetectorname << id.fIndex << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fLinearArray[id.fIndex].
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank,
                            id.fSubelement);

        }

        if( id.fTypeID == kQwBPMCavity ) {
          if( lkDEBUG ) {
            std::cout << "found stripline data for " << id.fdetectorname << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fCavity[id.fIndex].
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank,
                            id.fSubelement);
        }

        if( id.fTypeID == kQwBCM ) {
          if( lkDEBUG ) {
            std::cout << "found bcm data for " << id.fdetectorname << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fBCM[id.fIndex]->
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank);
        }

        if( id.fTypeID == kQwClock ) {
          if( lkDEBUG ) {
            std::cout << "found clock data for " << id.fdetectorname << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fClock[id.fIndex]->
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank);
        }

        if( id.fTypeID == kQwHaloMonitor ) {
          if( lkDEBUG ) {
            std::cout << "found halo monitor data for " << id.fdetectorname << std::endl;
            std::cout << "word left to read in this buffer:" << num_words - id.fWordInSubbank << std::endl;
          }
          fHaloMonitor[id.fIndex].
            ProcessEvBuffer(&(buffer[id.fWordInSubbank]),
                            num_words - id.fWordInSubbank);
        }

      }
    }
  }

  return 0;
}


//*****************************************************************//
Bool_t QwBeamLine::ApplySingleEventCuts() {

  Bool_t status = kTRUE;

  for_each(ALL(fBCM),         QwApplySngEvCut(status, bDEBUG, "BCM"));
  for_each(ALL(fClock),       QwApplySngEvCut(status, bDEBUG, "Clock"));
  for_each(ALL(fHaloMonitor), QwApplySngEvCut(status, bDEBUG, "HaloMonitor"));
  for_each(ALL(fStripline),   QwApplySngEvCut(status, bDEBUG, "BPMStripline"));
  for_each(ALL(fQPD),         QwApplySngEvCut(status, bDEBUG, "QPD"));
  for_each(ALL(fLinearArray), QwApplySngEvCut(status, bDEBUG, "LinearArray"));
  for_each(ALL(fCavity),      QwApplySngEvCut(status, bDEBUG, "BPMCavity"));
  for_each(ALL(fBCMCombo),    QwApplySngEvCut(status, bDEBUG, "CombinedBCM"));
  for_each(ALL(fBPMCombo),    QwApplySngEvCut(status, bDEBUG, "CombinedBPM"));
  for_each(ALL(fECalculator), QwApplySngEvCut(status, bDEBUG, "EnergyCalculator"));

  //If at least one of the devices failed  event cuts, increment error counter for QwBeamLine
  if( !status )
    fQwBeamLineErrorCount++;

  return status;
}

//*****************************************************************//
Bool_t QwBeamLine::CheckForBurpFail( const VQwSubsystem* subsys ) {
  Bool_t burpstatus = kFALSE;
  auto* tmp = const_cast<VQwSubsystem*>(subsys);
  if( Compare(tmp) ) {
    const auto* input = dynamic_cast<const QwBeamLine*>(subsys);
    QwCheckForBurpFail check(burpstatus);

    QwForEachInBoth(ALL(fClock),       ALL(input->fClock),       check);
    QwForEachInBoth(ALL(fStripline),   ALL(input->fStripline),   check);
    QwForEachInBoth(ALL(fQPD),         ALL(input->fQPD),         check);
    QwForEachInBoth(ALL(fLinearArray), ALL(input->fLinearArray), check);
    QwForEachInBoth(ALL(fCavity),      ALL(input->fCavity),      check);
    QwForEachInBoth(ALL(fBCM),         ALL(input->fBCM),         check);
    QwForEachInBoth(ALL(fBCMCombo),    ALL(input->fBCMCombo),    check);
    QwForEachInBoth(ALL(fBPMCombo),    ALL(input->fBPMCombo),    check);
    QwForEachInBoth(ALL(fECalculator), ALL(input->fECalculator), check);
    QwForEachInBoth(ALL(fHaloMonitor), ALL(input->fHaloMonitor), check);
  }
  return burpstatus;
}


//*****************************************************************//
void
QwBeamLine::PrintErrorCounters() const {//inherited from the VQwSubsystemParity; this will display the error summary

  QwMessage << "*********QwBeamLine Error Summary****************" << QwLog::endl;
  QwVQWK_Channel::PrintErrorCounterHead();

  for(size_t i=0;i<fClock.size();i++){
    fClock[i].get()->PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fBCM.size();i++){
    fBCM[i].get()->PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fHaloMonitor.size();i++){
    fHaloMonitor[i].PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fStripline.size();i++){
    fStripline[i].get()->PrintErrorCounters();
    printf("\n");
  }
  for(size_t i=0;i<fQPD.size();i++){
    fQPD[i].PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fLinearArray.size();i++){
    fLinearArray[i].PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fCavity.size();i++){
    fCavity[i].PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fBCMCombo.size();i++){
    fBCMCombo[i].get()->PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fBPMCombo.size();i++){
    fBPMCombo[i].get()->PrintErrorCounters();
  }
  printf("\n");
  for(size_t i=0;i<fECalculator.size();i++){
    fECalculator[i].PrintErrorCounters();
  }
  QwVQWK_Channel::PrintErrorCounterTail();
}

//*****************************************************************//
void QwBeamLine::IncrementErrorCounters()
{
  for(size_t i=0;i<fClock.size();i++){
    fClock[i].get()->IncrementErrorCounters();
  }
  for(size_t i=0;i<fBCM.size();i++){
    fBCM[i].get()->IncrementErrorCounters();
  }
  for(size_t i=0;i<fHaloMonitor.size();i++){
    fHaloMonitor[i].IncrementErrorCounters();
  }
  for(size_t i=0;i<fStripline.size();i++){
    fStripline[i].get()->IncrementErrorCounters();
    }
  for(size_t i=0;i<fQPD.size();i++){
    fQPD[i].IncrementErrorCounters();
  }
  for(size_t i=0;i<fLinearArray.size();i++){
    fLinearArray[i].IncrementErrorCounters();
  }
  for(size_t i=0;i<fCavity.size();i++){
    fCavity[i].IncrementErrorCounters();
  }
  for(size_t i=0;i<fBCMCombo.size();i++){
    fBCMCombo[i].get()->IncrementErrorCounters();
  }
  for(size_t i=0;i<fBPMCombo.size();i++){
    fBPMCombo[i].get()->IncrementErrorCounters();
  }
  for(size_t i=0;i<fECalculator.size();i++){
    fECalculator[i].IncrementErrorCounters();
  }
}

//*****************************************************************//
UInt_t QwBeamLine::GetEventcutErrorFlag() {//return the error flag
  UInt_t ErrorFlag;
  UInt_t ErrorFlagtmp;
  ErrorFlag=0;
  for(size_t i=0;i<fBCM.size();i++){
    ErrorFlagtmp = fBCM[i].get()->GetEventcutErrorFlag();
    ErrorFlag |=ErrorFlagtmp;
  }
  for(size_t i=0;i<fStripline.size();i++){
    ErrorFlag |= fStripline[i].get()->GetEventcutErrorFlag();
  }
  for(size_t i=0;i<fQPD.size();i++){
    ErrorFlag |= fQPD[i].GetEventcutErrorFlag();
  }
  for(size_t i=0;i<fLinearArray.size();i++){
    ErrorFlag |= fLinearArray[i].GetEventcutErrorFlag();
  }
  for(size_t i=0;i<fCavity.size();i++){
    ErrorFlag |= fCavity[i].GetEventcutErrorFlag();
  }
  for(size_t i=0;i<fBCMCombo.size();i++){
    ErrorFlag |= fBCMCombo[i].get()->GetEventcutErrorFlag();
  }
  for(size_t i=0;i<fBPMCombo.size();i++){
    ErrorFlag |= fBPMCombo[i].get()->GetEventcutErrorFlag();
  }
  for(size_t i=0;i<fECalculator.size();i++){
    ErrorFlag |= fECalculator[i].GetEventcutErrorFlag();
  }

  return ErrorFlag;

}

//*****************************************************************//
UInt_t QwBeamLine::UpdateErrorFlag() {//return the error flag
  UInt_t ErrorFlag;
  UInt_t ErrorFlagtmp;
  ErrorFlag=0;
  for(size_t i=0;i<fBCM.size();i++){
    ErrorFlagtmp = fBCM[i].get()->UpdateErrorFlag();
    ErrorFlag |=ErrorFlagtmp;
  }
  for(size_t i=0;i<fStripline.size();i++){
    ErrorFlag |= fStripline[i].get()->UpdateErrorFlag();
  }
  for(size_t i=0;i<fQPD.size();i++){
    ErrorFlag |= fQPD[i].UpdateErrorFlag();
  }
  for(size_t i=0;i<fLinearArray.size();i++){
    ErrorFlag |= fLinearArray[i].UpdateErrorFlag();
  }
  for(size_t i=0;i<fCavity.size();i++){
    ErrorFlag |= fCavity[i].UpdateErrorFlag();
  }
  for(size_t i=0;i<fBCMCombo.size();i++){
    ErrorFlag |= fBCMCombo[i].get()->UpdateErrorFlag();
  }
  for(size_t i=0;i<fBPMCombo.size();i++){
    ErrorFlag |= fBPMCombo[i].get()->UpdateErrorFlag();
  }
  for(size_t i=0;i<fECalculator.size();i++){
    ErrorFlag |= fECalculator[i].UpdateErrorFlag();
  }

  return ErrorFlag;

}

//*****************************************************************//
void QwBeamLine::UpdateErrorFlag( const VQwSubsystem* ev_error )
{
  auto* tmp = const_cast<VQwSubsystem*>(ev_error);
  if( Compare(tmp) ) {

    const auto* input = dynamic_cast<const QwBeamLine*>(ev_error);

    /*
   for(size_t i=0;i<input->fClock.size();i++)
*(fClock[i].get())=*(input->fClock[i].get());
*/
    for( size_t i = 0; i < input->fStripline.size(); i++ )
      (fStripline[i].get())->UpdateErrorFlag(input->fStripline[i].get());
    for( size_t i = 0; i < input->fQPD.size(); i++ )
      (fQPD[i]).UpdateErrorFlag(&(input->fQPD[i]));
    for( size_t i = 0; i < input->fLinearArray.size(); i++ )
      (fLinearArray[i]).UpdateErrorFlag(&(input->fLinearArray[i]));
    for( size_t i = 0; i < input->fCavity.size(); i++ )
      (fCavity[i]).UpdateErrorFlag(&(input->fCavity[i]));
    for( size_t i = 0; i < input->fBCM.size(); i++ ) {
      (fBCM[i].get())->UpdateErrorFlag(input->fBCM[i].get());
    }
    for( size_t i = 0; i < input->fBCMCombo.size(); i++ )
      (fBCMCombo[i].get())->UpdateErrorFlag(
        input->fBCMCombo[i].get()); //*(fBCMCombo[i].get())=*(input->fBCMCombo[i].get());
    for( size_t i = 0; i < input->fBPMCombo.size(); i++ )
      (fBPMCombo[i].get())->UpdateErrorFlag(input->fBPMCombo[i].get()); //=*(input->fBPMCombo[i].get());
    for( size_t i = 0; i < input->fECalculator.size(); i++ )
      (fECalculator[i]).UpdateErrorFlag(&(input->fECalculator[i]));
    for( size_t i = 0; i < input->fHaloMonitor.size(); i++ )
      (fHaloMonitor[i]).UpdateErrorFlag(&(input->fHaloMonitor[i]));
  }
}


//*****************************************************************//
void QwBeamLine::ProcessEvent()
{
  QwProcessEvent         process;
  //QwProcessEventAndPrint process_and_print;

  // Make sure this one comes first! The clocks are needed by
  // other elements.
  for_each(ALL(fClock),       process);

  for_each(ALL(fStripline),   process);
  for_each(ALL(fCavity),      process);
  for_each(ALL(fBCM),         process);
  for_each(ALL(fQPD),         process);
  for_each(ALL(fLinearArray), process);
  for_each(ALL(fHaloMonitor), process);
  for_each(ALL(fBCMCombo),    process);
  for_each(ALL(fBPMCombo),    process);
  for_each(ALL(fECalculator), process);
}


//*****************************************************************//
Int_t QwBeamLine::ProcessConfigurationBuffer( const ROCID_t roc_id, const BankID_t bank_id, UInt_t* buffer,
                                              UInt_t num_words ) {

  return 0;
}

//*****************************************************************//
Bool_t QwBeamLine::PublishInternalValues() const {
  // Publish variables
  Bool_t status = kTRUE;

  // Publish variables through map file
  // This should work with bcm, bpmstripline, bpmcavity, combo bpm and combo bcm
  for( const auto& names: fPublishList ) {
    const TString& publish_name = names.at(0);
    TString device_type = names.at(1);
    const TString& device_name = names.at(2);
    TString device_prop = names.at(3);
    device_type.ToLower();
    device_prop.ToLower();

    const VQwHardwareChannel* tmp_channel = nullptr;

    EQwBeamInstrumentType type_id;
    if( device_type == "combobpm" )
      type_id = kQwCombinedBPM;
    else if( device_type == "combobcm" )
      type_id = kQwCombinedBCM;
    else if( device_type == "comboenergy" )
      type_id = kQwEnergyCalculator;
    else if( device_type == "scaler" )
      type_id = kQwHaloMonitor;
    else
      type_id = GetQwBeamInstrumentType(device_type);

    Int_t index = GetDetectorIndex(type_id, device_name);

    if( type_id != kQwUnknownDeviceType && index != -1 ) {
      tmp_channel = GetChannel(type_id, index, device_prop);
    } else
      QwError << "QwBeamLine::PublishInternalValues() error " << QwLog::endl;

    if( !tmp_channel ) {
      QwError << "QwBeamLine::PublishInternalValues(): " << publish_name << " not found" << QwLog::endl;
      status = kFALSE;
    } else {
      QwDebug << "QwBeamLine::PublishInternalValues(): " << publish_name << " found" << QwLog::endl;

      status = PublishInternalValue(publish_name, "published-value", tmp_channel);
    }
  }

  return status;
}

//*****************************************************************//
Bool_t QwBeamLine::PublishByRequest( TString device_name ) {
  Bool_t status = kFALSE;
  const VQwHardwareChannel* tmp_channel = nullptr;

  std::vector<TString> publishinfo(4, TString(""));
  publishinfo.at(0) = device_name;

  EQwBeamInstrumentType type_id;
  Int_t index = -1;

  TString name = device_name;
  TString device_prop = "value";
  if( device_name.EndsWith("WS") ) {
    name = device_name(0, device_name.Length() - 2);
    device_prop = "ef";
  } else if( device_name.EndsWith("Q") ) {
    name = device_name(0, device_name.Length() - 1);
    device_prop = "ef";
  } else if( device_name.EndsWith("XSlope") ) {
    name = device_name(0, device_name.Length() - 6);
    device_prop = "xp";
  } else if( device_name.EndsWith("YSlope") ) {
    name = device_name(0, device_name.Length() - 6);
    device_prop = "yp";
  } else if( device_name.EndsWith("X") ) {
    name = device_name(0, device_name.Length() - 1);
    device_prop = "x";
  } else if( device_name.EndsWith("Y") ) {
    name = device_name(0, device_name.Length() - 1);
    device_prop = "y";
  }

  for( auto& id: fBeamDetectorID ) {
    if( id.fdetectorname == name
        || id.fdetectorname == device_name ) {
      index = id.fIndex;
      type_id = id.fTypeID;

      publishinfo.at(1) = GetQwBeamInstrumentTypeName(type_id);
      publishinfo.at(2) = id.fdetectorname;
      publishinfo.at(3) = device_prop;
      break;
    }
  }

  if( index != -1 ) {
    tmp_channel = GetChannel(type_id, index, publishinfo.at(3));
    fPublishList.push_back(publishinfo);
    status = PublishInternalValue(publishinfo.at(0), "published-by-request",
                                  tmp_channel);
  }

  return status;
}


//*****************************************************************//
void QwBeamLine::ClearEventData()
{
  QwClearEventData clear;

  for_each(ALL(fClock),       clear);
  for_each(ALL(fStripline),   clear);
  for_each(ALL(fCavity),      clear);
  for_each(ALL(fBCM),         clear);
  for_each(ALL(fQPD),         clear);
  for_each(ALL(fLinearArray), clear);
  for_each(ALL(fHaloMonitor), clear);
  for_each(ALL(fBCMCombo),    clear);
  for_each(ALL(fBPMCombo),    clear);
  for_each(ALL(fECalculator), clear);
}

//*****************************************************************//
Int_t QwBeamLine::GetDetectorIndex( EQwBeamInstrumentType type_id, const TString& name ) const {
  Bool_t ldebug = kFALSE;
  Int_t result = -1;
  if( ldebug ) {
    std::cout << "QwBeamLine::GetDetectorIndex\n";
    std::cout << "type_id==" << type_id << " name=" << name << "\n";
    std::cout << fBeamDetectorID.size() << " already registered detector\n";
  }
  for(const auto & detectorId : fBeamDetectorID) {
    if( ldebug ) {
      std::cout << "testing against (" << detectorId.fTypeID
                << "," << detectorId.fdetectorname << ")=>" << result << "\n";
    }
    if( detectorId.fTypeID == type_id
        && detectorId.fdetectorname == name ) {
      result = detectorId.fIndex;
      break;
    }
  }
  return result;
}

VQwDataElement* QwBeamLine::GetElement( const QwBeamDetectorID& det_id ) {
  return GetElement(det_id.fTypeID, det_id.fIndex);
}

VQwDataElement* QwBeamLine::GetElement( EQwBeamInstrumentType TypeID, const TString& name ) {
  Int_t index = GetDetectorIndex(TypeID, name);
  return GetElement(TypeID, index);
}

VQwDataElement* QwBeamLine::GetElement( EQwBeamInstrumentType TypeID, Int_t index ) {
  VQwDataElement* tmp_ptr;
  switch( TypeID ) {
    case kQwBPMStripline:
      tmp_ptr = fStripline.at(index).get();
      break;
    case kQwQPD:
      tmp_ptr = &(fQPD.at(index));
      break;
    case kQwLinearArray:
      tmp_ptr = &(fLinearArray.at(index));
      break;
    case kQwBCM:
      tmp_ptr = fBCM.at(index).get();
      break;
    case kQwCombinedBCM:
      tmp_ptr = fBCMCombo.at(index).get();
      break;
    case kQwCombinedBPM:
      tmp_ptr = fBPMCombo.at(index).get();
      break;
    case kQwEnergyCalculator:
      tmp_ptr = &(fECalculator.at(index));
      break;
    case kQwHaloMonitor:
      tmp_ptr = &(fHaloMonitor.at(index));
      break;
    case kQwBPMCavity:
      tmp_ptr = &(fCavity.at(index));
      break;
    case kQwClock:
      tmp_ptr = fClock.at(index).get();
      break;
    default:
      TString loc = "QwBeamLine::GetElement called by "
                    + GetName() + " with invalid arguements: "
                    + GetQwBeamInstrumentTypeName(TypeID) + " "
                    + Form("%d", index);
      throw std::invalid_argument(loc.Data());
  }
  return tmp_ptr;
}

const VQwDataElement* QwBeamLine::GetElement( EQwBeamInstrumentType TypeID, Int_t index ) const {
  return const_cast<QwBeamLine*>(this)->GetElement(TypeID, index);
}

const VQwHardwareChannel*
QwBeamLine::GetChannel( EQwBeamInstrumentType TypeID, Int_t index, const TString& device_prop ) const {
  const VQwHardwareChannel* tmp_channel = nullptr;

  if( TypeID == kQwBPMStripline || TypeID == kQwBPMCavity || TypeID == kQwQPD ) {
    const auto* tmp_ptr = dynamic_cast<const VQwBPM*>(GetElement(TypeID, index));
    assert(tmp_ptr);
    if( tmp_ptr ) {
      if( device_prop == "x" )
        tmp_channel = tmp_ptr->GetPosition(VQwBPM::kXAxis);
      else if( device_prop == "y" )
        tmp_channel = tmp_ptr->GetPosition(VQwBPM::kYAxis);
      else if( device_prop == "ef" )
        tmp_channel = tmp_ptr->GetEffectiveCharge();
    }
  } else if( TypeID == kQwCombinedBPM ) {
    const auto* tmp_ptr = dynamic_cast<const VQwBPM*>(GetElement(TypeID, index));
    assert(tmp_ptr);
    if( tmp_ptr ) {
      if( device_prop == "x" )
        tmp_channel = tmp_ptr->GetPosition(VQwBPM::kXAxis);
      else if( device_prop == "y" )
        tmp_channel = tmp_ptr->GetPosition(VQwBPM::kYAxis);
      else if( device_prop == "ef" )
        tmp_channel = fBPMCombo.at(index)->GetEffectiveCharge();
      else if( device_prop == "xp" )
        tmp_channel = fBPMCombo.at(index)->GetAngleX();
      else if( device_prop == "yp" )
        tmp_channel = fBPMCombo.at(index)->GetAngleY();
    }
  } else if( TypeID == kQwLinearArray ) {
    /// TODO: QwBeamLine::GetChannel
    /// How do we access linear array channel outputs?
  } else if( TypeID == kQwBCM || TypeID == kQwCombinedBCM ) {
    tmp_channel = dynamic_cast<const VQwBCM*>(GetElement(TypeID, index))->GetCharge();
  } else if( TypeID == kQwEnergyCalculator ) {
    tmp_channel = fECalculator.at(index).GetEnergy();
  } else if( TypeID == kQwHaloMonitor ) {
    tmp_channel = fHaloMonitor.at(index).GetScaler();
  } else if( TypeID == kQwClock ) {
    tmp_channel = fClock.at(index)->GetTime();
  } else {
    TString loc = "QwBeamLine::GetChannel called by "
                  + GetName() + " with invalid arguements: "
                  + GetQwBeamInstrumentTypeName(TypeID) + " "
                  + Form("%d", index);
    throw std::invalid_argument(loc.Data());
  }
  return tmp_channel;
}

//*****************************************************************//
VQwBPM* QwBeamLine::GetBPMStripline( const TString& name ) {
  for( const auto& stripline: fStripline ) {
    if( stripline->GetElementName() == name ) {
      return stripline.get();
    }
  }
  return nullptr;
}

//*****************************************************************//

QwBPMCavity* QwBeamLine::GetBPMCavity( const TString& name ) {
  for( auto& cavity: fCavity ) {
    if( cavity.GetElementName() == name ) {
      return &cavity;
    }
  }
  return nullptr;
}


//*****************************************************************//
VQwBCM* QwBeamLine::GetBCM( const TString& name ) {
  //QwWarning << "QwBeamLine::GetBCM" << QwLog::endl;
  for( const auto& bcm: fBCM ) {
    if( bcm->GetElementName() == name ) {
      return bcm.get();
    }
  }
  //QwWarning << "BCM Not Found" << QwLog::endl;
  return nullptr;
}


//*****************************************************************//
VQwClock* QwBeamLine::GetClock( const TString& name ) {
  //QwWarning << "QwBeamLine::GetClock" << QwLog::endl;
  for( auto& clock: fClock ) {
    if( clock->GetElementName() == name ) {
      return clock.get();
    }
  }

  //QwWarning << "Clock Not Found" << QwLog::endl;
  return nullptr;
}

//*****************************************************************//
VQwBCM* QwBeamLine::GetCombinedBCM( const TString& name ) {
  //QwWarning << "QwBeamLine::GetCombinedBCM" << QwLog::endl;

  for( auto& cbcm: fBCMCombo ) {
    if( cbcm->GetElementName() == name ) {
      return cbcm.get();
    }
  }
  return nullptr;
}

//*****************************************************************//
VQwBPM* QwBeamLine::GetCombinedBPM( const TString& name ) {
  //QwWarning << "QwBeamLine::GetCombinedBPM" << QwLog::endl;
  for( auto& cbpm: fBPMCombo ) {
    if( cbpm->GetElementName() == name ) {
      return cbpm.get();
    }
  }
  return nullptr;
}

//*****************************************************************//
QwEnergyCalculator* QwBeamLine::GetEnergyCalculator( const TString& name ) {
  for( auto& ecal: fECalculator ) {
    if( ecal.GetElementName() == name ) {
      return &ecal;
    }
  }
  return nullptr;
}

//*****************************************************************//
QwHaloMonitor* QwBeamLine::GetScalerChannel( const TString& name ) {
  for( auto& halo: fHaloMonitor ) {
    if( halo.GetElementName() == name ) {
      return &halo;
    }
  }
  return nullptr;
}

//*****************************************************************//
const VQwBPM* QwBeamLine::GetBPMStripline( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetBPMStripline(name);
}

//*****************************************************************//
const QwBPMCavity* QwBeamLine::GetBPMCavity( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetBPMCavity(name);
}

//*****************************************************************//
const VQwBCM* QwBeamLine::GetBCM( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetBCM(name);
}

//*****************************************************************//
const VQwClock* QwBeamLine::GetClock( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetClock(name);
}

//*****************************************************************//
const VQwBCM* QwBeamLine::GetCombinedBCM( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetCombinedBCM(name);
}

//*****************************************************************//
const VQwBPM* QwBeamLine::GetCombinedBPM( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetCombinedBPM(name);
}

//*****************************************************************//
const QwEnergyCalculator* QwBeamLine::GetEnergyCalculator( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetEnergyCalculator(name);
}

//*****************************************************************//
const QwHaloMonitor* QwBeamLine::GetScalerChannel( const TString& name ) const {
  return const_cast<QwBeamLine*>(this)->GetScalerChannel(name);
}


//*****************************************************************//
VQwSubsystem& QwBeamLine::operator=( VQwSubsystem* value ) {
  //  std::cout<<" here in QwBeamLine::operator= \n";
  if( Compare(value) ) {

    auto* input = dynamic_cast<QwBeamLine*>(value);

    for( size_t i = 0; i < input->fClock.size(); i++ )
      *fClock[i] = *input->fClock[i];
    for( size_t i = 0; i < input->fStripline.size(); i++ )
      *fStripline[i] = *input->fStripline[i];
    for( size_t i = 0; i < input->fQPD.size(); i++ )
      fQPD[i] = input->fQPD[i];
    for( size_t i = 0; i < input->fLinearArray.size(); i++ )
      fLinearArray[i] = input->fLinearArray[i];
    for( size_t i = 0; i < input->fCavity.size(); i++ )
      fCavity[i] = input->fCavity[i];
    for( size_t i = 0; i < input->fBCM.size(); i++ )
      *fBCM[i] = *input->fBCM[i];
    for( size_t i = 0; i < input->fHaloMonitor.size(); i++ )
      fHaloMonitor[i] = input->fHaloMonitor[i];
    for( size_t i = 0; i < input->fBCMCombo.size(); i++ )
      *fBCMCombo[i] = *input->fBCMCombo[i];
    for( size_t i = 0; i < input->fBPMCombo.size(); i++ )
      *fBPMCombo[i] = *input->fBPMCombo[i];
    for( size_t i = 0; i < input->fECalculator.size(); i++ )
      fECalculator[i] = input->fECalculator[i];

    if( !input->fPublishList.empty() ) {
      fPublishList.resize(input->fPublishList.size());
      for( size_t i = 0; i < input->fPublishList.size(); i++ ) {
        fPublishList[i].resize(input->fPublishList[i].size());
        for( size_t j = 0; j < input->fPublishList[i].size(); j++ ) {
          fPublishList[i][j] = input->fPublishList[i][j];
        }
      }
    }

  }
  return *this;
}


//*****************************************************************//
VQwSubsystem& QwBeamLine::operator+=( VQwSubsystem* value ) {
  if( Compare(value) ) {
    //QwBeamLine* input= (QwBeamLine*)value ;
    auto* input = dynamic_cast<QwBeamLine*>(value);

    for( size_t i = 0; i < input->fClock.size(); i++ )
      *(fClock[i].get()) += *(input->fClock[i].get());
    for( size_t i = 0; i < input->fStripline.size(); i++ )
      *(fStripline[i].get()) += *(input->fStripline[i].get());
    for( size_t i = 0; i < input->fCavity.size(); i++ )
      fCavity[i] += input->fCavity[i];
    for( size_t i = 0; i < input->fQPD.size(); i++ )
      fQPD[i] += input->fQPD[i];
    for( size_t i = 0; i < input->fLinearArray.size(); i++ )
      fLinearArray[i] += input->fLinearArray[i];
    for( size_t i = 0; i < input->fBCM.size(); i++ )
      *(fBCM[i].get()) += *(input->fBCM[i].get());
    for( size_t i = 0; i < input->fHaloMonitor.size(); i++ )
      fHaloMonitor[i] += input->fHaloMonitor[i];
    for( size_t i = 0; i < input->fBCMCombo.size(); i++ )
      *(fBCMCombo[i].get()) += *(input->fBCMCombo[i].get());
    for( size_t i = 0; i < input->fBPMCombo.size(); i++ )
      *(fBPMCombo[i].get()) += *(input->fBPMCombo[i].get());
    for( size_t i = 0; i < input->fECalculator.size(); i++ )
      fECalculator[i] += input->fECalculator[i];

    if( !input->fPublishList.empty() ) {
      fPublishList.resize(input->fPublishList.size());
      for( size_t i = 0; i < input->fPublishList.size(); i++ ) {
        fPublishList.at(i).resize(input->fPublishList.at(i).size());
        for( size_t j = 0; j < input->fPublishList.at(i).size(); j++ ) {
          fPublishList.at(i).at(j) = input->fPublishList.at(i).at(j);
        }
      }
    }

  }
  return *this;
}

//*****************************************************************//
VQwSubsystem& QwBeamLine::operator-=( VQwSubsystem* value ) {

  if( Compare(value) ) {
    auto* input = dynamic_cast<QwBeamLine*>(value);

    for( size_t i = 0; i < input->fClock.size(); i++ )
      *(fClock[i].get()) -= *(input->fClock[i].get());
    for( size_t i = 0; i < input->fStripline.size(); i++ )
      *(fStripline[i].get()) -= *(input->fStripline[i].get());
    for( size_t i = 0; i < input->fCavity.size(); i++ )
      fCavity[i] -= input->fCavity[i];
    for( size_t i = 0; i < input->fQPD.size(); i++ )
      fQPD[i] -= input->fQPD[i];
    for( size_t i = 0; i < input->fLinearArray.size(); i++ )
      fLinearArray[i] -= input->fLinearArray[i];
    for( size_t i = 0; i < input->fBCM.size(); i++ )
      *(fBCM[i].get()) -= *(input->fBCM[i].get());
    for( size_t i = 0; i < input->fHaloMonitor.size(); i++ )
      fHaloMonitor[i] -= input->fHaloMonitor[i];
    for( size_t i = 0; i < input->fBCMCombo.size(); i++ )
      *(fBCMCombo[i].get()) -= *(input->fBCMCombo[i].get());
    for( size_t i = 0; i < input->fBPMCombo.size(); i++ )
      *(fBPMCombo[i].get()) -= *(input->fBPMCombo[i].get());
    for( size_t i = 0; i < input->fECalculator.size(); i++ )
      fECalculator[i] -= input->fECalculator[i];

    if( !input->fPublishList.empty() ) {
      fPublishList.resize(input->fPublishList.size());
      for( size_t i = 0; i < input->fPublishList.size(); i++ ) {
        fPublishList.at(i).resize(input->fPublishList.at(i).size());
        for( size_t j = 0; j < input->fPublishList.at(i).size(); j++ ) {
          fPublishList.at(i).at(j) = input->fPublishList.at(i).at(j);
        }
      }
    }

  }
  return *this;
}


//*****************************************************************//
void QwBeamLine::Sum( VQwSubsystem* value1, VQwSubsystem* value2 ) {
  if( Compare(value1) && Compare(value2) ) {
    *this = value1;
    *this += value2;
  }
}


//*****************************************************************//
void QwBeamLine::Difference( VQwSubsystem* value1, VQwSubsystem* value2 ) {
  if( Compare(value1) && Compare(value2) ) {
    *this = value1;
    *this -= value2;
  }
}


//*****************************************************************//
void QwBeamLine::Ratio( VQwSubsystem* numer, VQwSubsystem* denom ) {
  if( Compare(numer) && Compare(denom) ) {
    auto* innumer = dynamic_cast<QwBeamLine*>(numer);
    auto* indenom = dynamic_cast<QwBeamLine*>(denom);

    for( size_t i = 0; i < innumer->fClock.size(); i++ )
      fClock[i].get()->Ratio(*(innumer->fClock[i].get()),
                             *(indenom->fClock[i].get()));
    for( size_t i = 0; i < innumer->fStripline.size(); i++ )
      fStripline[i].get()->Ratio(*(innumer->fStripline[i].get()),
                                 *(indenom->fStripline[i].get()));
    for( size_t i = 0; i < innumer->fCavity.size(); i++ )
      fCavity[i].Ratio(innumer->fCavity[i], indenom->fCavity[i]);
    for( size_t i = 0; i < innumer->fQPD.size(); i++ )
      fQPD[i].Ratio(innumer->fQPD[i], indenom->fQPD[i]);
    for( size_t i = 0; i < innumer->fLinearArray.size(); i++ )
      fLinearArray[i].Ratio(innumer->fLinearArray[i], indenom->fLinearArray[i]);
    for( size_t i = 0; i < innumer->fBCM.size(); i++ )
      fBCM[i].get()->Ratio(*(innumer->fBCM[i].get()),
                           *(indenom->fBCM[i].get()));
    for( size_t i = 0; i < innumer->fHaloMonitor.size(); i++ )
      fHaloMonitor[i].Ratio(innumer->fHaloMonitor[i], indenom->fHaloMonitor[i]);
    for( size_t i = 0; i < innumer->fBCMCombo.size(); i++ )
      fBCMCombo[i].get()->Ratio(*(innumer->fBCMCombo[i].get()),
                                *(indenom->fBCMCombo[i]));
    for( size_t i = 0; i < innumer->fBPMCombo.size(); i++ )
      fBPMCombo[i].get()->Ratio(*(innumer->fBPMCombo[i].get()),
                                *(indenom->fBPMCombo[i].get()));
    for( size_t i = 0; i < innumer->fECalculator.size(); i++ )
      fECalculator[i].Ratio(innumer->fECalculator[i], indenom->fECalculator[i]);

    // For the combined bcm, maybe we might want to think about getting
    // the asymmetry using the asymmetries of the individual bcms with a
    // weight. But right now it is unclear if really need to have that
    // option.
  }
}

//*****************************************************************//
void QwBeamLine::Scale( Double_t factor ) {
  QwScale scale(factor);

  for_each(ALL(fClock),       scale);
  for_each(ALL(fStripline),   scale);
  for_each(ALL(fCavity),      scale);
  for_each(ALL(fQPD),         scale);
  for_each(ALL(fLinearArray), scale);
  for_each(ALL(fBCM),         scale);
  for_each(ALL(fHaloMonitor), scale);
  for_each(ALL(fBCMCombo),    scale);
  for_each(ALL(fBPMCombo),    scale);
  for_each(ALL(fECalculator), scale);
}

//*****************************************************************//
void QwBeamLine::CalculateRunningAverage() {
  QwCalcRunAvg avg;

  for_each(ALL(fClock),       avg);
  for_each(ALL(fStripline),   avg);
  for_each(ALL(fCavity),      avg);
  for_each(ALL(fQPD),         avg);
  for_each(ALL(fLinearArray), avg);
  for_each(ALL(fBCM),         avg);
  for_each(ALL(fHaloMonitor), avg);
  for_each(ALL(fBCMCombo),    avg);
  for_each(ALL(fBPMCombo),    avg);
  for_each(ALL(fECalculator), avg);
}

//*****************************************************************//
void QwBeamLine::PrintValue() const {
  QwPrintValue print;

  QwMessage << "=== QwBeamLine: " << GetName() << " ===" << QwLog::endl;
  QwMessage << "Clock"         << QwLog::endl;
  for_each(ALL(fClock),       print);
  QwMessage << "BPM stripline" << QwLog::endl;
  for_each(ALL(fStripline),   print);
  QwMessage << "QPD"           << QwLog::endl;
  for_each(ALL(fQPD),         print);
  QwMessage << "LinearArray"   << QwLog::endl;
  for_each(ALL(fLinearArray), print);
  QwMessage << "BPM cavity"    << QwLog::endl;
  for_each(ALL(fCavity),      print);
  QwMessage << "BCM"           << QwLog::endl;
  for_each(ALL(fBCM),         print);
  QwMessage << "HaloMonitor"   << QwLog::endl;
  for_each(ALL(fHaloMonitor), print);
  QwMessage << "BPM combo"     << QwLog::endl;
  for_each(ALL(fBCMCombo),    print);
  QwMessage << "BPM combo"     << QwLog::endl;
  for_each(ALL(fBPMCombo),    print);
  QwMessage << "Energy "       << QwLog::endl;
  for_each(ALL(fECalculator), print);
}

//*****************************************************************//
void QwBeamLine::AccumulateRunningSum( VQwSubsystem* value1, Int_t count, Int_t ErrorMask ) {
  if( Compare(value1) ) {
    auto* value = dynamic_cast<QwBeamLine*>(value1);

    QwAccumulate acc(count, ErrorMask);

    QwForEachInBoth(ALL(fClock), ALL(value->fClock), acc);
    QwForEachInBoth(ALL(fStripline), ALL(value->fStripline), acc);
    QwForEachInBoth(ALL(fCavity), ALL(value->fCavity), acc);
    QwForEachInBoth(ALL(fBCM), ALL(value->fBCM), acc);
    QwForEachInBoth(ALL(fBCMCombo), ALL(value->fBCMCombo), acc);
    QwForEachInBoth(ALL(fBPMCombo), ALL(value->fBPMCombo), acc);
    QwForEachInBoth(ALL(fECalculator), ALL(value->fECalculator), acc);
    QwForEachInBoth(ALL(fQPD), ALL(value->fQPD), acc);
    QwForEachInBoth(ALL(fLinearArray), ALL(value->fLinearArray), acc);
    QwForEachInBoth(ALL(fHaloMonitor), ALL(value->fHaloMonitor), acc);
  }
}

//*****************************************************************//
void QwBeamLine::DeaccumulateRunningSum( VQwSubsystem* value1, Int_t ErrorMask ) {
  if( Compare(value1) ) {
    auto* value = dynamic_cast<QwBeamLine*>(value1);
    QwDeaccumulate deacc(ErrorMask);

    QwForEachInBoth(ALL(fClock), ALL(value->fClock), deacc);
    QwForEachInBoth(ALL(fStripline), ALL(value->fStripline), deacc);
    QwForEachInBoth(ALL(fCavity), ALL(value->fCavity), deacc);
    QwForEachInBoth(ALL(fBCM), ALL(value->fBCM), deacc);
    QwForEachInBoth(ALL(fBCMCombo), ALL(value->fBCMCombo), deacc);
    QwForEachInBoth(ALL(fBPMCombo), ALL(value->fBPMCombo), deacc);
    QwForEachInBoth(ALL(fECalculator), ALL(value->fECalculator), deacc);
    QwForEachInBoth(ALL(fQPD), ALL(value->fQPD), deacc);
    QwForEachInBoth(ALL(fLinearArray), ALL(value->fLinearArray), deacc);
    QwForEachInBoth(ALL(fHaloMonitor), ALL(value->fHaloMonitor), deacc);
  }
}

//*****************************************************************//
Bool_t QwBeamLine::Compare( VQwSubsystem* value ) {
  //  std::cout<<" Here in QwBeamLine::Compare \n";

  Bool_t res = kTRUE;
  if( typeid(*value) != typeid(*this) ) {
    res = kFALSE;
    //      std::cout<<" types are not ok \n";
    //      std::cout<<" this is bypassed just for now but should be fixed eventually \n";
  } else {
    auto* input = dynamic_cast<QwBeamLine*>(value);
    if( input->fStripline.size() != fStripline.size() ) {
      //	  std::cout<<" not the same number of striplines \n";
      res = kFALSE;
    } else if( input->fBCM.size() != fBCM.size() ) {
      res = kFALSE;
      //	  std::cout<<" not the same number of bcms \n";
    } else if( input->fHaloMonitor.size() != fHaloMonitor.size() ) {
      res = kFALSE;
      //	  std::cout<<" not the same number of halomonitors \n";
    } else if( input->fClock.size() != fClock.size() ) {
      res = kFALSE;
      //	  std::cout<<" not the same number of halomonitors \n";
    } else if( input->fBCMCombo.size() != fBCMCombo.size() ) {
      res = kFALSE;
    } else if( input->fBPMCombo.size() != fBPMCombo.size() ) {
      res = kFALSE;
    } else if( input->fLinearArray.size() != fLinearArray.size() ) {
      res = kFALSE;
    } else if( input->fECalculator.size() != fECalculator.size() ) {
      res = kFALSE;
    } else if( input->fCavity.size() != fCavity.size() ) {
      res = kFALSE;
    } else if( input->fQPD.size() != fQPD.size() ) {
      res = kFALSE;
    }

  }
  return res;
}


//*****************************************************************//
void QwBeamLine::ConstructHistograms( TDirectory* folder, TString& prefix ) {

  //  std::cout<<" here is QwBeamLine::ConstructHistogram with prefix ="<<prefix<<"\n";
  QwConstructHistograms construct(folder, prefix);

  for_each(ALL(fClock),       construct);
  for_each(ALL(fStripline),   construct);
  for_each(ALL(fQPD),         construct);
  for_each(ALL(fLinearArray), construct);
  for_each(ALL(fCavity),      construct);
  for_each(ALL(fBCM),         construct);
  for_each(ALL(fHaloMonitor), construct);
  for_each(ALL(fBCMCombo),    construct);
  for_each(ALL(fBPMCombo),    construct);
  for_each(ALL(fECalculator), construct);
}

//*****************************************************************//
void QwBeamLine::FillHistograms()
{
  QwFillHistograms fill;

  for_each(ALL(fClock),       fill);
  for_each(ALL(fStripline),   fill);
  for_each(ALL(fQPD),         fill);
  for_each(ALL(fLinearArray), fill);
  for_each(ALL(fCavity),      fill);
  for_each(ALL(fBCM),         fill);
  for_each(ALL(fHaloMonitor), fill);
  for_each(ALL(fBCMCombo),    fill);
  for_each(ALL(fBPMCombo),    fill);
  for_each(ALL(fECalculator), fill);
}


//*****************************************************************//
void QwBeamLine::ConstructBranchAndVector( TTree* tree, TString& prefix,
                                           std::vector<Double_t>& values )
{
  QwConstructBranchAndVector construct( tree, prefix, values );

  for_each(ALL(fClock),       construct);
  for_each(ALL(fStripline),   construct);
  for_each(ALL(fQPD),         construct);
  for_each(ALL(fLinearArray), construct);
  for_each(ALL(fCavity),      construct);
  for_each(ALL(fBCM),         construct);
  for_each(ALL(fHaloMonitor), construct);
  for_each(ALL(fBCMCombo),    construct);
  for_each(ALL(fBPMCombo),    construct);
  for_each(ALL(fECalculator), construct);
}

//*****************************************************************//
void QwBeamLine::ConstructBranch( TTree* tree, TString& prefix )
{
  QwConstructBranch construct(tree, prefix);

  for_each(ALL(fClock),       construct);
  for_each(ALL(fStripline),   construct);
  for_each(ALL(fQPD),         construct);
  for_each(ALL(fLinearArray), construct);
  for_each(ALL(fCavity),      construct);
  for_each(ALL(fBCM),         construct);
  for_each(ALL(fHaloMonitor), construct);
  for_each(ALL(fBCMCombo),    construct);
  for_each(ALL(fBPMCombo),    construct);
  for_each(ALL(fECalculator), construct);
}

//*****************************************************************//
void QwBeamLine::ConstructBranch( TTree* tree, TString& prefix, QwParameterFile& trim_file ) {
  TString tmp, varname, varvalue;
  tmp = "QwBCM";
  QwParameterFile* nextmodule;
  trim_file.RewindToFileStart();


  tmp = "QwBPMStripline";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fStripline.size(); i++ )
      fStripline[i].get()->ConstructBranch(tree, prefix, *nextmodule);

  }

  tmp = "QwQPD";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fQPD.size(); i++ )
      fQPD[i].ConstructBranch(tree, prefix, *nextmodule);
  }

  tmp = "QwLinearDiodeArray";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fLinearArray.size(); i++ )
      fLinearArray[i].ConstructBranch(tree, prefix, *nextmodule);
  }

  tmp = "QwBPMCavity";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fCavity.size(); i++ )
      fCavity[i].ConstructBranch(tree, prefix, *nextmodule);

  }

  tmp = "QwBCM";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fBCM.size(); i++ )
      fBCM[i].get()->ConstructBranch(tree, prefix, *nextmodule);
  }

  tmp = "QwClock";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fClock.size(); i++ )
      fClock[i].get()->ConstructBranch(tree, prefix, *nextmodule);
  }

  tmp = "QwHaloMonitor";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fHaloMonitor.size(); i++ )
      fHaloMonitor[i].ConstructBranch(tree, prefix, *nextmodule);
  }


  tmp = "QwCombinedBCM";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fBCMCombo.size(); i++ )
      fBCMCombo[i].get()->ConstructBranch(tree, prefix, *nextmodule);
  }


  tmp = "QwCombinedBPM";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fBPMCombo.size(); i++ )
      fBPMCombo[i].get()->ConstructBranch(tree, prefix, *nextmodule);
  }

  tmp = "QwEnergyCalculator";
  trim_file.RewindToFileStart();
  if( trim_file.FileHasModuleHeader(tmp) ) {
    nextmodule = trim_file.ReadUntilNextModule();//This section contains submodules and or channels to be included in the tree
    for( size_t i = 0; i < fECalculator.size(); i++ )
      fECalculator[i].ConstructBranch(tree, prefix, *nextmodule);
  }

}

//*****************************************************************//
void QwBeamLine::FillTreeVector( std::vector<Double_t>& values ) const
{
  QwFillTreeVector fill(values);

  for_each(ALL(fClock),       fill);
  for_each(ALL(fStripline),   fill);
  for_each(ALL(fQPD),         fill);
  for_each(ALL(fLinearArray), fill);
  for_each(ALL(fCavity),      fill);
  for_each(ALL(fBCM),         fill);
  for_each(ALL(fHaloMonitor), fill);
  for_each(ALL(fBCMCombo),    fill);
  for_each(ALL(fBPMCombo),    fill);
  for_each(ALL(fECalculator), fill);
}


//*****************************************************************//
void QwBeamLine::PrintInfo() const {
  std::cout << "Name of the subsystem =" << fSystemName << "\n";
  std::cout << "there are " << fClock.size() << " clock \n";
  std::cout << "there are " << fStripline.size() << " striplines \n";
  std::cout << "there are " << fQPD.size() << " QPDs \n";
  std::cout << "there are " << fLinearArray.size() << " LinearArrays \n";
  std::cout << "there are " << fCavity.size() << " cavities \n";
  std::cout << "there are " << fBCM.size() << " bcm \n";
  std::cout << "there are " << fHaloMonitor.size() << " halomonitors \n";
  std::cout << "there are " << fBCMCombo.size() << " combined bcms \n";
  std::cout << "there are " << fBPMCombo.size() << " combined bpms \n";
  std::cout << "there are " << fECalculator.size() << " energy calculators \n";
  std::cout << " Printing Running AVG and other channel info for BCMs" << std::endl;
  for( size_t i = 0; i < fBCM.size(); i++ )
    fBCM[i].get()->PrintInfo();
  for( size_t i = 0; i < fHaloMonitor.size(); i++ )
    fHaloMonitor[i].PrintInfo();
}


//*****************************************************************//
void QwBeamLine::PrintDetectorID() const {
  for( size_t i = 0; i < fBeamDetectorID.size(); i++ ) {
    std::cout << "=============================" << std::endl;
    std::cout << " Detector ID=" << i << std::endl;
    fBeamDetectorID[i].Print();
  }
}


//*****************************************************************//
void QwBeamLine::CopyTemplatedDataElements( const VQwSubsystem* source ) {
  const auto* input = dynamic_cast<const QwBeamLine*>(source);

  fClock.reserve(input->fClock.size());
  for( size_t i = 0; i < input->fClock.size(); i++ ) {
    fClock.push_back(VQwClock_ptr(VQwClock::Create(*(input->fClock[i].get()))));
  }

  fStripline.reserve(input->fStripline.size());
  for( size_t i = 0; i < input->fStripline.size(); i++ ) {
    fStripline.push_back(VQwBPM_ptr(
      VQwBPM::CreateStripline(*(input->fStripline[i].get()))));
  }

  fBCM.reserve(input->fBCM.size());
  for( size_t i = 0; i < input->fBCM.size(); i++ ) {
    fBCM.push_back(VQwBCM_ptr(
      VQwBCM::Create(*(input->fBCM[i].get()))));
  }

  fBCMCombo.reserve(input->fBCMCombo.size());
  for( size_t i = 0; i < input->fBCMCombo.size(); i++ ) {
    fBCMCombo.push_back(VQwBCM_ptr(
      VQwBCM::CreateCombo(*(
        input->fBCMCombo[i].get()))));
  }

  fBPMCombo.reserve(input->fBPMCombo.size());
  for( size_t i = 0; i < input->fBPMCombo.size(); i++ ) {
    fBPMCombo.push_back(VQwBPM_ptr(
      VQwBPM::CreateCombo(*(input->fBPMCombo[i].get()))));
  }
}

//*****************************************************************//
#ifdef __USE_DATABASE__
void QwBeamLine::FillDB(QwParityDB *db, TString datatype)
{

  Bool_t local_print_flag = false;

  if(local_print_flag) {
    QwMessage << " --------------------------------------------------------------- " << QwLog::endl;
    QwMessage << "                         QwBeamLine::FillDB                      " << QwLog::endl;
    QwMessage << " --------------------------------------------------------------- " << QwLog::endl;
  }

  std::vector<QwDBInterface> interface;
  std::vector<QwParitySSQLS::beam> entrylist;

  UInt_t analysis_id = db->GetAnalysisID();

  TString measurement_type_bcm;
  TString measurement_type_bpm;
  TString measurement_type_halo;

  measurement_type_bcm =
    QwDBInterface::DetermineMeasurementTypeID(datatype,"q");
  measurement_type_bpm =
    QwDBInterface::DetermineMeasurementTypeID(datatype,"p",kTRUE);
  measurement_type_halo =
    QwDBInterface::DetermineMeasurementTypeID(datatype);

  UInt_t i,j;
  i = j = 0;
  // try to access BCM mean and its error
  // there are 2 different types BCM data we have at the moment
  // Yield and Asymmetry
  if(local_print_flag)  QwMessage <<  QwColor(Qw::kGreen) << "Beam Current Monitors" <<QwLog::endl;

  for(i=0; i< fBCM.size(); i++) {
    interface.clear();
    interface = fBCM[i].get()->GetDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id );
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  ///   try to access BPM mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Beam Position Monitors" <<QwLog::endl;
  for(i=0; i< fStripline.size(); i++) {
    //    fStripline[i].MakeBPMList();
    interface.clear();
    interface = fStripline[i].get()->GetDBEntry();
    for (j=0; j<interface.size()-5; j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bpm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
    // effective charge (last 4 elements)  need to be saved as measurement_type_bcm
    for (j=interface.size()-5; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }


  ///   try to access CombinedBPM means and errors
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Combined Beam Position Monitors" <<QwLog::endl;
  for(i=0; i< fBPMCombo.size(); i++) {
    //    fBPMCombo[i].MakeBPMComboList();
    interface.clear();
    interface = fBPMCombo[i].get()->GetDBEntry();
    for (j=0; j<interface.size()-5; j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bpm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
    // effective charge (last element) need to be saved as measurement_type_bcm
    for (j=interface.size()-5; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  ///   try to access CombinedBCM means and errors
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Combined Beam Current Monitors" <<QwLog::endl;

  for(i=0; i< fBCMCombo.size(); i++) {
    interface.clear();
    interface = fBCMCombo[i].get()->GetDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id );
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  ///   try to access Energy Calculator mean and its error
  if(local_print_flag)  QwMessage <<  QwColor(Qw::kGreen) << "Energy Calculator" <<QwLog::endl;

  for(i=0; i< fECalculator.size(); i++) {
    interface.clear();
    interface = fECalculator[i].GetDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id );
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bpm );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }


  ///   try to access QPD mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Quadrant PhotoDiodes" <<QwLog::endl;
  for(i=0; i< fQPD.size(); i++) {
    //    fQPD[i].MakeQPDList();
    interface.clear();
    interface = fQPD[i].GetDBEntry();
    for (j=0; j<interface.size()-5; j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bpm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
    // effective charge need (last element) to be saved as measurement_type_bcm
    for (j=interface.size()-5; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  ///   try to access LinearArray mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Linear PhotoDiode Array" <<QwLog::endl;
  for(i=0; i< fLinearArray.size(); i++) {
    //    fLinearArray[i].MakeLinearArrayList();
    interface.clear();
    interface = fLinearArray[i].GetDBEntry();
    for (j=0; j<interface.size()-5; j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bpm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
    for (j=interface.size()-5; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  ///   try to access cavity bpm mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Cavity Monitors" <<QwLog::endl;
  for(i=0; i< fCavity.size(); i++) {
    //    fCavity[i].MakeBPMCavityList();
    interface.clear();
    interface = fCavity[i].GetDBEntry();
    for (j=0; j<interface.size()-5; j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bpm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
    for (j=interface.size()-5; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_bcm );
      interface.at(j).PrintStatus( local_print_flag);
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  // try to access halo mean and its error
  if(local_print_flag)  QwMessage <<  QwColor(Qw::kGreen) << "Halo Monitors" <<QwLog::endl;

  for(i=0; i< fHaloMonitor.size(); i++) {
    interface.clear();
    interface = fHaloMonitor[i].GetDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id );
      interface.at(j).SetMonitorID( db );
      interface.at(j).SetMeasurementTypeID( measurement_type_halo );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
  }

  if(local_print_flag){
    QwMessage << QwColor(Qw::kGreen)   << "Entrylist Size : "
        << QwColor(Qw::kBoldRed) << entrylist.size()
              << QwColor(Qw::kNormal)  << QwLog::endl;
  }

  db->Connect();
  // Check the entrylist size, if it isn't zero, start to query..
  if( entrylist.size() ) {
    mysqlpp::Query query= db->Query();
    query.insert(entrylist.begin(), entrylist.end());
    query.execute();
  }
  else {
    QwMessage << "QwBeamLine::FillDB :: This is the case when the entrlylist contains nothing in "<< datatype.Data() << QwLog::endl;
  }
  db->Disconnect();

  return;
}


void QwBeamLine::FillErrDB(QwParityDB *db, TString datatype)
{

  Bool_t local_print_flag = false;

  if(local_print_flag) {
    QwMessage << " --------------------------------------------------------------- " << QwLog::endl;
    QwMessage << "                      QwBeamLine::FillErrDB                      " << QwLog::endl;
    QwMessage << " --------------------------------------------------------------- " << QwLog::endl;
  }

  std::vector<QwErrDBInterface> interface;
  std::vector<QwParitySSQLS::beam_errors> entrylist;

  UInt_t analysis_id = db->GetAnalysisID();

  UInt_t i,j;
  i = j = 0;

  if(local_print_flag)  QwMessage <<  QwColor(Qw::kGreen) << "Beam Current Monitors" <<QwLog::endl;
  for(i=0; i< fBCM.size(); i++) {
    interface.clear();
    interface = fBCM[i].get()->GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id );
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  ///   try to access BPM mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Beam Position Monitors" <<QwLog::endl;
  for(i=0; i< fStripline.size(); i++) {
    interface.clear();
    interface = fStripline[i].get()->GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Combined Beam Position Monitors" <<QwLog::endl;
  for(i=0; i< fBPMCombo.size(); i++) {
    interface.clear();
    interface = fBPMCombo[i].get()->GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  ///   try to access CombinedBCM means and errors
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Combined Beam Current Monitors" <<QwLog::endl;
  for(i=0; i< fBCMCombo.size(); i++) {
    interface.clear();
    interface = fBCMCombo[i].get()->GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  ///   try to access Energy Calculator mean and its error
  if(local_print_flag)  QwMessage <<  QwColor(Qw::kGreen) << "Energy Calculator" <<QwLog::endl;
  for(i=0; i< fECalculator.size(); i++) {
    interface.clear();
    interface = fECalculator[i].GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  ///   try to access QPD mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Quadrant PhotoDiodes" <<QwLog::endl;
  for(i=0; i< fQPD.size(); i++) {
    interface.clear();
    interface = fQPD[i].GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  ///   try to access LinearArray mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Linear PhotoDiode Array" <<QwLog::endl;
  for(i=0; i< fLinearArray.size(); i++) {
    interface.clear();
    interface = fLinearArray[i].GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  ///   try to access cavity bpm mean and its error
  if(local_print_flag) QwMessage <<  QwColor(Qw::kGreen) << "Cavity Monitors" <<QwLog::endl;
  for(i=0; i< fCavity.size(); i++) {
    interface.clear();
    interface = fCavity[i].GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }

  // try to access halo mean and its error
  if(local_print_flag)  QwMessage <<  QwColor(Qw::kGreen) << "Halo Monitors" <<QwLog::endl;
  for(i=0; i< fHaloMonitor.size(); i++) {
    interface.clear();
    interface = fHaloMonitor[i].GetErrDBEntry();
    for (j=0; j<interface.size(); j++){
      interface.at(j).SetAnalysisID( analysis_id ) ;
      interface.at(j).SetMonitorID( db );
      interface.at(j).PrintStatus( local_print_flag );
      interface.at(j).AddThisEntryToList( entrylist );
    }
    if(local_print_flag) printf("\n");
  }


  if(local_print_flag){
    QwMessage << QwColor(Qw::kGreen)   << "Entrylist Size : "
          << QwColor(Qw::kBoldRed) << entrylist.size()
              << QwColor(Qw::kNormal)  << QwLog::endl;
  }

  db->Connect();
  // Check the entrylist size, if it isn't zero, start to query..
  if( entrylist.size() ) {
    mysqlpp::Query query= db->Query();
    query.insert(entrylist.begin(), entrylist.end());
    query.execute();
  }
  else {
    QwMessage << "QwBeamLine::FillErrDB :: This is the case when the entrlylist contains nothing in "<< datatype.Data() << QwLog::endl;
  }
  db->Disconnect();
  return;
}
#endif // __USE_DATABASE__

void QwBeamLine::WritePromptSummary( QwPromptSummary* ps, const TString& type ) {
  Bool_t local_print_flag = false;
  Bool_t local_add_element = type.Contains("yield");


  if( local_print_flag ) {
    printf("---------------------------------------------------------------\n");
    printf("QwBeamLine::WritePromptSummary()  Type : %12s\n", type.Data());
    printf("---------------------------------------------------------------\n");
  }

  const VQwHardwareChannel* tmp_channel = nullptr;
  TString element_name = "";
  Double_t element_value = 0.0;
  Double_t element_value_err = 0.0;
  Double_t element_value_width = 0.0;

  PromptSummaryElement* local_ps_element = nullptr;
  Bool_t local_add_these_elements = false;

  for( size_t i = 0; i < fBCM.size(); i++ ) {
    tmp_channel = GetChannel(kQwBCM, i, "");
    element_name = tmp_channel->GetElementName();
    element_value = 0.0;
    element_value_err = 0.0;
    element_value_width = 0.0;

    local_add_these_elements =
      element_name.EqualTo("bcm_an_us") || element_name.EqualTo("bcm_an_ds") || element_name.EqualTo("bcm_an_ds3") ||
      (element_name.Contains("cav4")); // Need to change this to add other BCMs in summary

    if( local_add_these_elements && local_add_element ) {
      ps->AddElement(new PromptSummaryElement(element_name));
    }


    local_ps_element = ps->GetElementByName(element_name);


    if( local_ps_element ) {
      element_value = tmp_channel->GetValue();
      element_value_err = tmp_channel->GetValueError();
      element_value_width = tmp_channel->GetValueWidth();

      local_ps_element->Set(type, element_value, element_value_err, element_value_width);
    }

    if( local_print_flag && local_ps_element ) {
      printf("Type %12s, Element %32s, value %12.4e error %8.4e  width %12.4e\n",
             type.Data(), element_name.Data(), element_value, element_value_err, element_value_width);
    }
  }

//Add BPM Cavity
  for( size_t i = 0; i < fCavity.size(); i++ ) {
    tmp_channel = GetChannel(kQwBPMCavity, i, "ef");
    element_name = tmp_channel->GetElementName();
    element_value = 0.0;
    element_value_err = 0.0;
    element_value_width = 0.0;

    local_add_these_elements = (element_name.EqualTo("cav4bQ") || element_name.EqualTo("cav4cQ") ||
                                element_name.EqualTo("cav4dQ")); // Need to change this to add other cavities in summary

    if( local_add_these_elements && local_add_element ) {
      ps->AddElement(new PromptSummaryElement(element_name));
    }


    local_ps_element = ps->GetElementByName(element_name);


    if( local_ps_element ) {
      element_value = tmp_channel->GetValue();
      element_value_err = tmp_channel->GetValueError();
      element_value_width = tmp_channel->GetValueWidth();

      local_ps_element->Set(type, element_value, element_value_err, element_value_width);
    }

    if( local_print_flag && local_ps_element ) {
      printf("Type %12s, Element %32s, value %12.4e error %8.4e  width %12.4e\n",
             type.Data(), element_name.Data(), element_value, element_value_err, element_value_width);
    }
  }

//////



  char property[2][6] = {"x", "y"};
  local_ps_element = nullptr;
  local_add_these_elements = false;


  for( size_t i = 0; i < fStripline.size(); i++ ) {
    for( Int_t j = 0; j < 2; j++ ) {
      tmp_channel = GetChannel(kQwBPMStripline, i, property[j]);
      element_name = tmp_channel->GetElementName();
      element_value = 0.0;
      element_value_err = 0.0;
      element_value_width = 0.0;

      local_add_these_elements =
        element_name.Contains("bpm4") || element_name.Contains("bpm18") || element_name.Contains("bpm14") ||
        element_name.Contains("bpm12"); //Need to change this to add other stripline BPMs in summary

      if( local_add_these_elements && local_add_element ) {
        ps->AddElement(new PromptSummaryElement(element_name));
      }

      local_ps_element = ps->GetElementByName(element_name);


      if( local_ps_element ) {
        element_value = tmp_channel->GetValue();
        element_value_err = tmp_channel->GetValueError();
        element_value_width = tmp_channel->GetValueWidth();
        local_ps_element->Set(type, element_value, element_value_err, element_value_width);
      }

      if( local_print_flag && local_ps_element ) {
        printf("Type %12s, Element %32s, value %12.4e error %8.4e  width %12.4e\n",
               type.Data(), element_name.Data(), element_value, element_value_err, element_value_width);
      }

    }
  }
}

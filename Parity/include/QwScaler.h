
#ifndef QWSCALER_H
#define QWSCALER_H

// System headers
#include <vector>

// Boost headers
#include <boost/shared_ptr.hpp>

// Qweak headers
#include "VQwSubsystemParity.h"
#include "QwScaler_Channel.h"


class QwScaler: public VQwSubsystemParity, public MQwSubsystemCloneable<QwScaler>
{
  public:
    /// No default constructor
    QwScaler() = delete;

    /// Constructor with name
    explicit QwScaler(const TString& name);
    /// Copy constructor
    QwScaler(const QwScaler& source)
      : VQwSubsystem(source)
      , VQwSubsystemParity(source)
      , fGoodEventCount{0}
    {
      fScaler.resize(source.fScaler.size());
      for( size_t i = 0; i < fScaler.size(); i++ ) {
        auto* scaler24 = dynamic_cast<QwSIS3801D24_Channel*>(source.fScaler.at(i));
        if( scaler24 ) {
          fScaler.at(i) = new QwSIS3801D24_Channel(*scaler24);
        } else {
          auto* scaler32 = dynamic_cast<QwSIS3801D32_Channel*>(source.fScaler.at(i));
          if( scaler32 ) {
            fScaler.at(i) = new QwSIS3801D32_Channel(*scaler32);
          }
        }
      }
    }

  /// Destructor
    virtual ~QwScaler();

    // Handle command line options
    static void DefineOptions(QwOptions &options);
    void ProcessOptions(QwOptions &options);

    Int_t LoadChannelMap( const TString& mapfile);
    Int_t LoadInputParameters( const TString& pedestalfile);

    void  ClearEventData();

    Int_t ProcessConfigurationBuffer(ROCID_t roc_id, BankID_t bank_id, UInt_t* buffer, UInt_t num_words);
    Int_t ProcessEvBuffer(ROCID_t roc_id, BankID_t bank_id, UInt_t *buffer, UInt_t num_words);
    void  ProcessEvent();

    using VQwSubsystem::ConstructHistograms;
    void  ConstructHistograms(TDirectory *folder, TString &prefix);
    void  FillHistograms();

    using VQwSubsystem::ConstructBranchAndVector;
    void  ConstructBranchAndVector(TTree *tree, TString &prefix, std::vector<Double_t> &values);
    void  ConstructBranch(TTree *tree, TString& prefix) { };
    void  ConstructBranch(TTree *tree, TString& prefix, QwParameterFile& trim_file) { };
    void  FillTreeVector(std::vector<Double_t> &values) const;

    Bool_t Compare(VQwSubsystem *source);

    VQwSubsystem& operator=(VQwSubsystem *value);
    VQwSubsystem& operator+=(VQwSubsystem *value);
    VQwSubsystem& operator-=(VQwSubsystem *value);
    void Sum(VQwSubsystem *value1, VQwSubsystem *value2);
    void Difference(VQwSubsystem *value1, VQwSubsystem *value2);
    void Ratio(VQwSubsystem *value1, VQwSubsystem  *value2);
    void Scale(Double_t factor);

    void AccumulateRunningSum(VQwSubsystem* value, Int_t count=0, Int_t ErrorMask=0xFFFFFFF);
    //remove one entry from the running sums for devices
    void DeaccumulateRunningSum(VQwSubsystem* value, Int_t ErrorMask=0xFFFFFFF);
    void CalculateRunningAverage();

    Int_t LoadEventCuts( const TString& filename);
    Bool_t SingleEventCuts();
    Bool_t ApplySingleEventCuts();

    Bool_t CheckForBurpFail(const VQwSubsystem *subsys){
        //QwError << "************* test inside Scaler *****************" << QwLog::endl;
        return kFALSE;
    };

    void IncrementErrorCounters();

    void PrintErrorCounters() const;
    UInt_t GetEventcutErrorFlag();
    //update the error flag in the subsystem level from the top level routines related to stability checks. This will uniquely update the errorflag at each channel based on the error flag in the corresponding channel in the ev_error subsystem
    void UpdateErrorFlag(const VQwSubsystem *ev_error){
    };

    void PrintValue() const;
    void PrintInfo() const;

    Double_t* GetRawChannelArray();

    Double_t GetDataForChannelInModule(Int_t modnum, Int_t channum) {
      Int_t index = fModuleChannel_Map[std::pair<Int_t,Int_t>(modnum,channum)];
      return fScaler.at(index)->GetValue();
    }

    Int_t GetChannelIndex(TString channelName, UInt_t module_number);

  private:

    // Number of good events
    Int_t fGoodEventCount;

    // Mapping from subbank to scaler channels
    typedef std::map< Int_t, std::vector< std::vector<Int_t> > > Subbank_to_Scaler_Map_t;
    Subbank_to_Scaler_Map_t fSubbank_Map;

    // Mapping from module and channel number to scaler channel
    typedef std::map< std::pair<Int_t,Int_t>, Int_t > Module_Channel_to_Scaler_Map_t;
    Module_Channel_to_Scaler_Map_t fModuleChannel_Map;

    // Mapping from name to scaler channel
    typedef std::map< TString, Int_t> Name_to_Scaler_Map_t;
    Name_to_Scaler_Map_t fName_Map;

    // Vector of scaler channels
    std::vector< VQwScaler_Channel* > fScaler; // Raw channels
    std::vector< UInt_t > fBufferOffset; // Offset in scaler buffer
    std::vector< std::pair< VQwScaler_Channel*, double > > fNorm;
};

#endif

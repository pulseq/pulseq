/** @file ExternalSequence.h */

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <map>
//#include <limits>	    // for std::numeric_limits<...>::quiet_NaN()

#ifndef _EXTERNAL_SEQUENCE_H_
#define _EXTERNAL_SEQUENCE_H_

//#define MASTER_SLAVE_FORMAT
#define EXTERNAL_GRADS

#define PI     3.1415926535897932384626433832795
#define TWO_PI 6.283185307179586476925286766558

#ifndef MIN
#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#endif

#ifndef MAX
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
#endif

// Define the path separator depending on the compile target
// it is different on host and scanner MPCU
#if defined(VXWORKS) 
#define PATH_SEPARATOR "/"
#elif defined (BUILD_PLATFORM_LINUX)
#define PATH_SEPARATOR "/"
#else
#define PATH_SEPARATOR "\\"
#endif

// Advance declaration(s)
class SeqBlock;
class ExternalSequence;

/**
 * @brief Output message types
 */
enum MessageType {
	ERROR_MSG=-2,
	WARNING_MSG,
	NORMAL_MSG,
	DEBUG_HIGH_LEVEL,
	DEBUG_MEDIUM_LEVEL,
	DEBUG_LOW_LEVEL
};

// Define the current level of messages to display
const MessageType MSG_LEVEL = NORMAL_MSG;
//const MessageType MSG_LEVEL = DEBUG_MEDIUM_LEVEL;
//const MessageType MSG_LEVEL = DEBUG_LOW_LEVEL;


/**
 * @brief Internal storage order
 */
enum Event {
//	DELAY,
	RF,
	GX,
	GY,
	GZ,
	ADC,
	EXT,
	LAST_UNUSED // this entry should be last in the list
};
const int NUM_EVENTS=LAST_UNUSED;
const int NUM_GRADS=ADC-GX;

/**
 * @brief RF event data
 *
 * Stores the amplitude, frequency offset and an index to basic shapes for the
 * magnitude and phase
 */
struct RFEvent
{
	float amplitude;     /**< @brief Amplitude of RF event (Hz) */
	int magShape;        /**< @brief ID of shape for magnitude */
	int phaseShape;      /**< @brief ID of shape for phase */
	int timeShape;       /**< @brief ID of shape for time sampling points */
	float center;        /**< @brief Effective RF center of the pulse shape measured from the start of the shape (us) */
    float freqPPM;       /**< @brief B0-dependent frequency offset of transmitter (ppm) */
    float phasePPM;      /**< @brief B0-dependent phase offset of transmitter (rad/MHz) */
    float freqOffset;    /**< @brief Constant frequency offset of transmitter (Hz) */
	float phaseOffset;   /**< @brief Phase offset of transmitter (rad) */
	int delay;           /**< @brief Delay prior to the pulse (us) */
	char use;            /**< @brief Single character indicating the intended use of the pulse, e.g. e,r,etc... */
};


/**
 * @brief Gradient event data
 *
 * Stores gradient amplitude and additional information
 * depending on the type:
 *  - **trapezoid:** ramp up, flat time, ramp down
 *  - **arbitrary:** index to basic shape (see CompressedShape)
 */
struct GradEvent
{
	float amplitude;      /**< @brief Amplitude of gradient (Hz/m) */
	int delay;
	// Trapezoid:
	long rampUpTime;      /**< @brief Ramp up time of trapezoid (us) */
	long flatTime;        /**< @brief Flat-top time of trapezoid (us) */
	long rampDownTime;    /**< @brief Ramp down time of trapezoid (us) */
	// Arbitrary:
	int waveShape;        /**< @brief whave shape ID for arbitrary gradient */
	int timeShape;        /**< @brief time shapeID for arbitrary gradient; 0 means regular sampling with the default step */
	float first;          /**< @brief amplitude at the start of the shape for arbitrary gradient */
	float last;           /**< @brief amplitude at the end of the shape for arbitrary gradient */
};

#define FLOAT_UNDEFINED 1.18e-38

/**
 * @brief ADC readout event
 *
 * Stores number of ADC samples, dwell time, delay and frequency/phase offset
 * information (used to shift the FOV in-plane)
 *
 */
struct ADCEvent
{
	int numSamples;           /**< @brief Number of samples */
	int dwellTime;            /**< @brief Dwell time of ADC readout (ns) */
	int delay;                /**< @brief Delay before first sample (us) */
    float freqPPM;            /**< @brief B0-dependent frequency offset of receiver (ppm) */
    float phasePPM;           /**< @brief B0-dependent phase offset of receiver (rad/MHz) */
    float freqOffset;         /**< @brief Constant frequency offset of receiver (Hz) */
	float phaseOffset;        /**< @brief Phase offset of receiver (rad) */
	int phaseModulationShape; /**< @brief Phase modulation shape of receiver (rad) */
};

/**
 * @brief Extension list entry
 *
 * Stores the information aout the extension events in a form of a single-linked list
 */
struct ExtensionListEntry
{
	int type;      /**< @brief extension type. see ExtType enum for details */
	int ref;       /**< @brief reference to the actual extension event */
	int next;      /**< @brief link to the next extension entry in the list. 0 indicates the end of the list */
};

/**
 * @brief Known extension types. Extenssions are recognized by their text ID and are mapped to these constants in the code
 */
enum ExtType {
	EXT_LIST=0,
	EXT_TRIGGER,
	EXT_ROTATION,
	EXT_LABELSET,
	EXT_LABELINC,
	EXT_DELAY,
	EXT_UNKNOWN /* marks the end of the enum, should always be the last */
};

/**
 * @brief Trigger event (extension)
 *
 * Stores trigger, type, duration
 */
struct TriggerEvent
{
	long duration;          /**< @brief Duration of trigger event (us) */
	long delay;             /**< @brief Delay prior to the trigger event (us) */
	int triggerType;        /**< @brief Type of trigger (system dependent). 0: undefined / unused */
	int triggerChannel;     /**< @brief channel of trigger (system dependent). 0: undefined / unused */
};


/**
 * @brief Trigger event (extension)
 *
 * Stores trigger, type, duration
 */
#define SOFT_DELAY_HINT_LENGTH 32
struct SoftDelayEvent
{
	int numID;                        /**< @brief Numeric index of the soft delay to help the intepreter (together with the hint string) to identify the delay and allocate it to the UI element*/
	int offset;                       /**< @brief Offset (positive or negative) added to the delay after the division by the factor (us) */
	int factor;                       /**< @brief Factor by which the value on the user interface needs to be divided for calculating the final delay applied to the sequence */
	char hint[SOFT_DELAY_HINT_LENGTH]; /**< @brief Text hint corresponding to this soft delay, e.g. TE */
};

/**
 * @brief Rotation event (extension)
 *
 * Stores rotation matrix
 */
struct RotationEvent
{
	bool defined;           /**< @brief Indicates whether a rotation object was defined in this block */
	double rotMatrix[9];    /**< @brief Gradient rotation matrix */
};

/**
 * @brief all supported labels
 */
enum Labels {
	SLC, 
	SEG, 
	ECO, 
	PHS, 
	SET,
	ACQ,
	LIN, 
	PAR, 
    AVG,
    REP, // LAST_ADC_RELEVANT_LABEL
    ONCE,
	LABEL_UNKNOWN // this entry should be the last in the list
};
const int NUM_LABELS=LABEL_UNKNOWN;
const int LAST_ADC_RELEVANT_LABEL = REP;
/**
 * @brief all supported flags
 */
enum Flags{
	NAV,
	REV,
	SMS, 
	REF,
	IMA,
	NOISE, // LAST_ADC_RELEVANT_FLAG
	PMC, 
	NOPOS,
	NOROT,
	NOSCL,
	FLAG_UNKNOWN
};
const int NUM_FLAGS=FLAG_UNKNOWN;
const int LAST_ADC_RELEVANT_FLAG = NOISE;

template <class T>
struct LabelStorage 
{
	std::vector<T>  val;							/**< @brief current label values */
	std::vector<bool> bValUpdated;					/**< @brief flag whether the numeric label value was updated in this block */
	std::vector<bool> bValUsed;						/**< @brief flag whether the numeric label value was ever updated in this sequence */
};

/**
 * @brief List of MDH IDs
 *
 * Stores IDs/Labels/Flags that relate MDH headers 
 */
struct LabelValueStorage 
{
	LabelStorage<int>  num;     /**< @brief current numeric label values and corresponding state flags */
    LabelStorage<bool> flag;    /**< @brief current flag value */
    std::vector<int>   getAdcCounters(bool bIgnoreREP, bool bAlsoIgnoreAVG); /**< @brief current ADC-relevant numeric label values as a simple vector */
};

struct MinMaxLabelStorage 
{
	std::vector<int>  numValMin;				/**< @brief minimum of the numeric label values */
	std::vector<int>  numValMax;				/**< @brief maximum of the numeric label values */
	std::vector<bool> bNumMinMaxValid;			/**< @brief flag whether the min/max values are valid, that is the numeric label value was ever updated in this sequence */
	// unsure whether we really need this
	std::vector<bool> flagValMin;				/**< @brief minimum of the flag label values */
	std::vector<bool> flagValMax;				/**< @brief maximum of the flag label values */
	std::vector<bool> bFlagMinMaxValid;			/**< @brief flag whether the min/max values are valid, that is the flag label value was ever updated in this sequence */
	// parallel imaging flags
	int minRefLin;
	int numRefLin;
	int numRefImaLin;
	int minRefPar;
	int numRefPar;
	int numRefImaPar;
};

class LabelStateAndBookkeeping
{
  public:
    LabelStateAndBookkeeping() { 
		initBookkeeping();
        initCurrState();
	}

    void initBookkeeping();
	void initBookkeepingADC();
    void initCurrState();

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  GetMinMaxDataLabelStorageADC						*
    // *                                                                    *
    // * Description :  Get minimum and maximum of the  values in the       *
    // *                LabelStorage that were used for an ADC event        *
    // * Return      :  MinMaxLabelStorage                                  *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    const MinMaxLabelStorage& GetMinMaxDataLabelStorageADC();

	// * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  updateLabelValues							        *
    // *                                                                    *
    // * Description :  Update the current global label values              *
    // * Return      :  bool                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    void updateLabelValues(SeqBlock* pBlock);
    
	bool checkLabelValuesADC();         // called for each ADC event during sequence execution
    void updateBookkeepingRecordsADC(); // called for each ADC event during sequence preparation
    
	bool isFirstScanInSlice();
    bool isLastScanInSlice();
    bool isLastScanInMeas();

	bool anyLabelsInUse();
    bool nonAdcLabelsInUse();
    bool adcLabelsInUse();

    int currVal(Labels lbl);
    int  maxVal(Labels lbl);
    int  minVal(Labels lbl);
    bool currFlag(Flags lbl);
    bool isFlagInUse(Flags lbl);

	void dump(const char* szMsg = NULL, bool bMinMax = true, bool bCurr = true);

  private:
    void dump_internal(const std::vector<int>& numVal, const std::vector<bool>& flagVal, const char* szMsg);
    void finalizeBookkeepingRecordsADC();

	bool m_bAdcLabelsInUse;
    bool m_bNonAdcLabelsInUse;

    LabelValueStorage m_currLabelValueStorage;
    MinMaxLabelStorage m_MinMaxLabelBookkeepingADC; // this only store min/max values in the context of the ADC, the conters
                                                    // may run out of counds temporarily but can be reset again before the ADC
    std::map<int, std::vector<int> > m_mapFirstInSlc;
    std::map<int, std::vector<int> > m_mapLastInSlc;
    std::vector<int>                m_lastInMeas;
    std::set<std::vector<int> > m_setFirstInSlc;
    std::set<std::vector<int> > m_setLastInSlc;
};

//inline std::vector<int> LabelValueStorage::getAdcCounters() { return std::vector<int>(num.val.begin(), num.val.begin()+LAST_ADC_RELEVANT_LABEL+1); }
inline int  LabelStateAndBookkeeping::currVal(Labels id) { return m_currLabelValueStorage.num.val[id]; }
inline bool LabelStateAndBookkeeping::currFlag(Flags id) { return m_currLabelValueStorage.flag.val[id]; }
inline bool LabelStateAndBookkeeping::isFlagInUse(Flags id) { return m_currLabelValueStorage.flag.bValUsed[id]; }
inline int LabelStateAndBookkeeping::maxVal(Labels lbl) { return m_MinMaxLabelBookkeepingADC.numValMax[lbl]; }
inline int LabelStateAndBookkeeping::minVal(Labels lbl) { return m_MinMaxLabelBookkeepingADC.numValMin[lbl]; }
inline const MinMaxLabelStorage& LabelStateAndBookkeeping::GetMinMaxDataLabelStorageADC() { return m_MinMaxLabelBookkeepingADC; }
inline bool LabelStateAndBookkeeping::anyLabelsInUse()    { return m_bAdcLabelsInUse || m_bNonAdcLabelsInUse; }
inline bool LabelStateAndBookkeeping::nonAdcLabelsInUse() { return m_bNonAdcLabelsInUse; }
inline bool LabelStateAndBookkeeping::adcLabelsInUse()    { return m_bAdcLabelsInUse; }
inline bool LabelStateAndBookkeeping::isLastScanInMeas(){ return m_lastInMeas == m_currLabelValueStorage.getAdcCounters(false,false); }


/**
 * @brief Label event (extension)
 *
 * Stores set/inc MDH headers and the corresponding target NameID (key) 
 */
struct LabelEvent
{
	std::pair<int,int >  numVal;		/**< @brief set/inc value of the target label */		
	std::pair<int,bool > flagVal;		/**< @brief set/unset bool value of the target flag */
};

struct LabelMap
{
	typedef std::pair<Labels,Flags> tP;
	typedef std::map<std::string, tP> tM; 
	tM mapStrToLabel;
	std::map<int,std::string> mapLabelIdToStr;
	std::map<int,std::string> mapFlagIdToStr;
};

/**
 * @brief List of event IDs
 *
 * Stores IDs that reference events stored in the event libraries.
 */
struct EventIDs
{
	int id[NUM_EVENTS];
};


/**
 * @brief Sequence block
 *
 * A class representing a custom sequence block. As a minimum the class stores the
 * indices of the events occurring in this block. However, the object can also
 * store a complete description of the block after it is decompressed.
 */
class SeqBlock
{
	friend class ExternalSequence;
public:
	/**
	 * @brief Constructor
	 */
	SeqBlock() { 
		gradWaveforms.resize(NUM_GRADS);
		gradExtTrapForms.resize(NUM_GRADS);
	}

	/**
	 * @brief Return `true` if block has RF event
	 */
	bool    isRF();

	/**
	 * @brief Return `true` if block has trapezoid event on given channel
	 */
	bool    isTrapGradient(int channel);

	/**
	 * @brief Return `true` if block has arbitrary gradient event on given channel
	 */
	bool    isArbitraryGradient(int channel);

	/**
	 * @brief Return `true` if block has arbitrary gradient event on given channel and it is oversampled by a factor of 2
	 */
	bool    isArbGradWithOversampling(int channel);


	/**
	 * @brief Return `true` if block has extended trapezoid event on given channel
	 */
	bool    isExtTrapGradient(int channel);

	/**
	 * @brief Return `true` if block has ADC readout event
	 */
	bool    isADC();

	/**
	 * @brief Return `true` if block has delay
	 */
	//bool    isDelay();

	/**
	 * @brief Return `true` if block has a trigger command
	 */
	bool    isTrigger();
	
	/**
	 * @brief Return `true` if block has a gradient rotation command
	 */
	bool    isRotation();

	/**
	 * @brief Return `true` if block has a soft delay extension 
	 */
	bool    isSoftDelay();
	
	/**
	 * @brief Return index of this block
	 */
	int     GetIndex();

	/**
	 * @brief Return ID of the corresponding given event type
	 * @param type Type of event
	 */
	int     GetEventIndex(Event type);

	/**
	 * @brief Return delay of block
	 */
	//long    GetDelay();

	/**
	 * @brief Return duration of block in units of us. If soft delay is defined in the block, it overrides the stored value
	 */
	double GetDuration();

	/**
	 * @brief Return duration of block in the units of duration raster. If soft delay is defined in the block, it overrides the stored value
	 */
	long    GetDuration_ru();

	/**
	 * @brief Return the stored duration of block in units of us, irrespective whether soft delay is defined
	 */
	double GetStoredDuration();

	/**
	 * @brief Return the stored duration of block in the units of duration raster, irrespective whether soft delay is defined
	 */
	long    GetStoredDuration_ru();

	/**
	 * @brief Return the number of samples of the arbitrary gradient on the given gradient channel.
	 * Only relevant for arbitrary gradients
	 */
	int     GetArbGradNumSamples(int channel);

	/**
	 * @brief Directly get a pointer to the samples of the arbitrary gradient
	 */
	float* GetArbGradShapePtr(int channel);

	/**
	 * @brief Return the timening and the shape of the ExtTrp grdient on the given gradient channel.
	 * Only relevant for ExtTrap gradients
	 */
	const std::vector<long>& GetExtTrapGradTimes(int channel);

	/**
	 * @brief Return the timening and the shape of the ExtTrp grdient on the given gradient channel.
	 * Only relevant for ExtTrap gradients
	 */
	const std::vector<float>& GetExtTrapGradShape(int channel);

	/**
	 * @brief Return the gradient event of the given channel
	 */
	GradEvent& GetGradEvent(int channel);

	/**
	 * @brief Return the RF event
	 */
	RFEvent&	GetRFEvent();

	/**
	 * @brief Return the number of samples of the RF shape
	 */
	int    GetRFLength();

	/**
	 * @brief Directly get a pointer to the samples of the RF amplitude shape
	 */
	float* GetRFAmplitudePtr();

	/**
	 * @brief Directly get a pointer to the samples of the RF phase shape
	 */
	float* GetRFPhasePtr();

	/**
	 * @brief Get dwell time for the RF amplitude and phase shapes (in us)
	 */
	float GetRFDwellTime();

	/**
	 * @brief Return the ADC event
	 */
	ADCEvent& GetADCEvent();

	/**
	 * @brief Return the trigger command event
	 */
	TriggerEvent& GetTriggerEvent();

	/**
	 * @brief Return the soft delay extension event
	 */
	SoftDelayEvent&  GetSoftDelayEvent();

	/**
	 * @brief set the actual soft delay for the current block calculated based on the extension event
	 */
	void setActualSoftDelay_ru(long softDelay_ru);
	
	/**
	 * @brief get the actual soft delay of the current block calculated based on the extension event
	 */
	long getActualSoftDelay_ru();
	
	/**
	 * @brief Return the trigger command event
	 */
	RotationEvent& GetRotationEvent();
	
	/**
	 * @brief Return a brief string description of the block
	 */
	std::string GetTypeString();

	/**
	 * @brief Free the memory associated with decompressed shapes of this block
	 */
	void   free();

	/**
	 * @brief Return the label command event array
	 */
	std::vector<LabelEvent>&  GetLabelSetEvents();
	/**
	 * @brief Return the label command event array
	 */
	std::vector<LabelEvent>&  GetLabelIncEvents();
	/**
	 * @brief Return `true` if block has a labelset/labelinc command
	 */
	bool     isLabel();

	static double getBlockDurationRaster();

protected:
	int index;          /**< @brief Index of this block */

	// Event array contains integer indices to events stored in the parent ExternalSequence object
	int events[NUM_EVENTS];	/**< @brief list of event indices (RF, GX, GY, GZ, ADC) */

	//long delay;         /**< @brief delay of this block (in us) */
	long duration_ru;      /**< @brief duration of this block in raster units */

	RFEvent rf;                 /**< @brief RF event */
	GradEvent grad[NUM_GRADS];  /**< @brief gradient events */
	ADCEvent adc;               /**< @brief ADC event  */
	TriggerEvent trigger;       /**< @brief trigger event (just one per block) */
	RotationEvent rotation;     /**< @brief optional rotation event */
	SoftDelayEvent softDelay;   /**< @brief optional soft delay event */
	long actualSoftDelay_ru;    /**< @brief actual soft delay to be applied in the block > */
	std::vector<LabelEvent> labelinc; /**< @brief labelinc event, can be more than one */ // MZ: TODO: check if we should switch to storing only IDs in the library
	std::vector<LabelEvent> labelset; /**< @brief labelset event, can be more than one */ // MZ: TODO: check if we should switch to storing only IDs in the library
	// Below is only valid once decompressed:

	// RF
	std::vector<float> rfAmplitude;    /**< @brief RF amplitude shape (uncompressed) */
	std::vector<float> rfPhase;        /**< @brief RF phase shape (uncompressed) */
	float              rfDwellTime_us; /**< @brief dwell time of the RF shapes (in us) */

	// Gradient waveforms
	std::vector< std::vector<float> > gradWaveforms;    /**< @brief Arbitrary gradient shapes for each channel (uncompressed) */

	// ExtTrap waveforms
	std::vector< std::pair< std::vector< long >, std::vector< float > > > gradExtTrapForms;    /**< @brief ExtTrap gradient shapes for each channel (uncompressed) */

	// static for the duraton raster
	static double s_blockDurationRaster;
};

// * ------------------------------------------------------------------ *
// * Inline functions                                                   *
// * ------------------------------------------------------------------ *
inline int       SeqBlock::GetIndex() { return index; }

inline bool      SeqBlock::isRF() { return (events[RF]>0); }
inline bool      SeqBlock::isTrapGradient(int channel) { return ((events[channel+GX]>0) && (grad[channel].waveShape==0)); }
inline bool      SeqBlock::isExtTrapGradient(int channel) { return ((events[channel+GX]>0) && (grad[channel].waveShape!=0) && (grad[channel].timeShape>0)); }
inline bool      SeqBlock::isArbitraryGradient(int channel) { return ((events[channel+GX]>0) && (grad[channel].waveShape!=0) && (grad[channel].timeShape<=0)); }
inline bool      SeqBlock::isArbGradWithOversampling(int channel) { return ((events[channel+GX]>0) && (grad[channel].waveShape!=0) && (grad[channel].timeShape==-1)); }
inline bool      SeqBlock::isADC() { return (events[ADC]>0); }
inline bool      SeqBlock::isRotation() { return (events[EXT]>0 && rotation.defined); }
inline bool      SeqBlock::isTrigger() { return (events[EXT]>0) && trigger.triggerType!=0; }
inline bool      SeqBlock::isSoftDelay() { return (events[EXT]>0) && softDelay.numID>=0; }
inline double    SeqBlock::GetDuration() { return ( (actualSoftDelay_ru<0) ? duration_ru : actualSoftDelay_ru) * SeqBlock::s_blockDurationRaster; }
inline long      SeqBlock::GetDuration_ru() { return (actualSoftDelay_ru<0) ? duration_ru : actualSoftDelay_ru; }
inline void      SeqBlock::setActualSoftDelay_ru(long softDelay_ru) {actualSoftDelay_ru=softDelay_ru;}
inline long      SeqBlock::getActualSoftDelay_ru() {return actualSoftDelay_ru;}
inline double    SeqBlock::GetStoredDuration() { return duration_ru * SeqBlock::s_blockDurationRaster; }
inline long      SeqBlock::GetStoredDuration_ru() { return duration_ru; }

inline int       SeqBlock::GetEventIndex(Event type) { return events[type]; }

inline GradEvent& SeqBlock::GetGradEvent(int channel) { return grad[channel]; }
inline RFEvent&   SeqBlock::GetRFEvent() { return rf; }
inline ADCEvent&  SeqBlock::GetADCEvent() { return adc; }
inline TriggerEvent&  SeqBlock::GetTriggerEvent() { return trigger; }
inline SoftDelayEvent&  SeqBlock::GetSoftDelayEvent() { return softDelay; }
inline RotationEvent&  SeqBlock::GetRotationEvent() { return rotation; }

inline std::vector<LabelEvent>&  SeqBlock::GetLabelSetEvents() { return labelset; }
inline std::vector<LabelEvent>&  SeqBlock::GetLabelIncEvents() { return labelinc; }
inline bool		 SeqBlock::isLabel() {return(events[EXT]>0) && !(labelset.empty()&&labelinc.empty()); }

inline std::string SeqBlock::GetTypeString() {
	std::string type;
	if (isRF()) type = type + "RF";
	if (isTrapGradient(0)) type += " TrapX";
	if (isTrapGradient(1)) type += " TrapY";
	if (isTrapGradient(2)) type += " TrapZ";
	if (isExtTrapGradient(0)) type += " ExtTrapX";
	if (isExtTrapGradient(1)) type += " ExtTrapY";
	if (isExtTrapGradient(2)) type += " ExtTrapZ";
	if (isArbitraryGradient(0)) type += " ArbX";
	if (isArbitraryGradient(1)) type += " ArbY";
	if (isArbitraryGradient(2)) type += " ArbZ";
	if (isADC()) type += " ADC";
	if (isRotation()) type = type + " Rot";
	if (isTrigger()) type = type + " Trig";
	if (isSoftDelay()) type = type + " SoftDelay";
	
	return type;
}

inline float*    SeqBlock::GetArbGradShapePtr(int channel) { return (gradWaveforms[channel].size()>0) ? &gradWaveforms[channel][0] : NULL; }
inline int       SeqBlock::GetArbGradNumSamples(int channel) {	return gradWaveforms[channel].size(); }

inline const std::vector<long>&  SeqBlock::GetExtTrapGradTimes(int channel) { return gradExtTrapForms[channel].first; }
inline const std::vector<float>& SeqBlock::GetExtTrapGradShape(int channel) { return gradExtTrapForms[channel].second; }

inline float*    SeqBlock::GetRFAmplitudePtr() { return &rfAmplitude[0]; }
inline float*    SeqBlock::GetRFPhasePtr() { return &rfPhase[0]; }
inline int       SeqBlock::GetRFLength() { return rfAmplitude.size(); }
inline float     SeqBlock::GetRFDwellTime() { return rfDwellTime_us; }
inline double SeqBlock::getBlockDurationRaster() {return SeqBlock::s_blockDurationRaster; }

inline void      SeqBlock::free() {
	// Force the memory to be freed
	std::vector<float>().swap(rfAmplitude);
	std::vector<float>().swap(rfPhase);
	std::vector<std::vector<float> >().swap(gradWaveforms);
	std::vector<std::pair<std::vector<long>,std::vector<float> > >().swap(gradExtTrapForms);
 }


/**
 * @brief Compressed shape data
 *
 * Stores an arbitrary shape compressed with run-length compression of the derivative.
 * The decompressed samples should be in the range [0,1].
 *
 */
struct CompressedShape
{
	int numUncompressedSamples;    /**< @brief Number of samples *after* decompression */
	bool isCompressed;             /**< @brief Flag whether the samples are compressed or not */
	std::vector<float> samples;    /**< @brief Compressed samples */
};


/**
 * @brief Data representing the entire MR sequence
 *
 * This class defines an abstract sequence consisting of arbitrary *blocks*
 * The basic structure of the sequence consists of three hierarchical levels:
 *     1. **Block** - list of integer indices pointing to simultaneous events for RF, gradients, ADC.
 *     2. **Event** - Basic event with one of the following types:
 *                   + *delay* - Simple delay
 *                   + *RF* - contains amplitude, frequency offset, and two indices pointing to basic shape for amp & phase
 *                   + *Trapezoid gradient* - contains amplitude, ramp up, flat top, ramp down times
 *                   + *Arbitrary gradient* - contains amplitude and index to basic shape
 *                   + *ADC readout* - contains number of samples, dwell time, delay, frequency & phase offsets
 *     3. **Shape** - List of compressed samples defining a basic shape (arbitrary RF or gradient)
 *
 * @author Kelvin Layton <kelvin.layton@uniklinik-freiburg.de>
 * @date May 2014
 */

class ExternalSequence
{
  public:

	/**
	 * @brief Constructor
	 */
	ExternalSequence();

	/**
	 * @brief Destructor
	 */
	~ExternalSequence();

	/**
	 * @brief Reset the sequence in memory 
	 *
	 * Resets the state of the class to an empty sequence
	 */
	void reset();

	/**
	 * @brief Load the sequence from file
	 *
	 * Reads the sequence files into the class members, the sequence is stored in
	 * compressed format. The given path refer to either:
	 *  1. A single file with all sequence definitions
	 *  2. A directory containing a single file (e.g. external.seq)
	 *  3. A directory containing three files (blocks.seq, events.seq, shapes.seq)
	 *
	 * @param  path location of file or directory
	 */
	bool load(std::string path);

	/**
	 * @brief Load the sequence from a single stream
	 *
	 * Reads the sequence from a text buffer into the class members, the sequence is stored in
	 * compressed format in a single file with all sequence definitions.
	 *
	 * @param  buffer
	 */

	bool load_from_buffer(char* buffer);

	/**
	 * @brief Load the sequence from a single stream
	 *
	 * Reads the sequence files into the class members, the sequence is stored in
	 * compressed format. The given path refer to either:
	 *  1. A single file with all sequence definitions
	 *  2. One of the three files (blocks.seq, events.seq, shapes.seq)
	 *
	 * @param  path location of file or directory
	 */

	enum load_mode {lm_singlefile=0, lm_shapes, lm_events, lm_blocks};
	bool load(std::istream &data_stream, load_mode loadMose = lm_singlefile);

	/**
	 * @brief Report the version of the loaded sequence
	 *
	 * Returns the version of the loaded PulSeq file combined to a single integer
	 * The format is MMmmmmrrrr where M is the major version number, m is the minor
	 * version number and r is the revision
	 *
	 */
	int GetVersion() { return version_combined; }


	/**
	 * @brief Display an output message
	 *
	 * Display a message only if the MSG_LEVEL is sufficiently high.
	 * This function calls the low-level output function, which can be overridden
	 * using SetPutMsgFunction().
	 *
	 * @param  level  type of message
	 * @param  ss     string stream containing the message
	 */
	static void print_msg(MessageType level, std::ostream& ss);

	/**
	 * @brief A pointer-type to the low-level print function
	 */
	typedef void (*PrintFunPtr)(const std::string &str);

	/**
	 * @brief Set the output print function
	 *
	 * Set the low-level output function to be used for printing or logging messages.
	 * If no function is specified messages will be printed to cout
	 *
	 * @param  fun  Pointer to function
	 */
	static void SetPrintFunction(PrintFunPtr fun);

	/**
	 * @brief Lookup the custom definition
	 *
	 * Search the list of user-specified definitions through the [DEFINITIONS] section.
	 * If the definition key is not found, an empty vector is returned.
	 *
	 * @param key  the definition name
	 * @return a list of values (or empty vector)
	 */
	std::vector<double> GetDefinition(std::string key);

	/**
	 * @brief Lookup the custom definition
	 *
	 * Search the list of user-specified definitions through the [DEFINITIONS] section.
	 * If the definition key is not found, an empty string is returned.
	 *
	 * @param key  the definition name
	 * @return definition value as a string (or empty string)
	 */
	std::string GetDefinitionStr(std::string key);

	/**
     * @brief Return a vector of all definitions (keys)
     *
     * Return the list of user-specified definitions through the [DEFINITIONS] section.
     * If no definitions present, an empty vector is returned.
     *
     * @return a list of strings (or empty vector)
     */
    std::vector<std::string> GetAllDefinitions();
    
	/**
	 * @brief Return number of sequence blocks
	 */
	int  GetNumberOfBlocks(void);

	/**
	 * @brief Construct a sequence block from the library events
	 *
	 * Events are loaded from the library. However, arbitrary waveforms are
	 * not decoded until decodeBlock() is called.
	 *
	 * @see decodeBlock()
	 */
	SeqBlock*  GetBlock(int blockIndex);

	/**
	 * @brief Decode a block by looking up indexed events
	 *
	 * This involves assigning the block's event objects from the libraries
	 * as well as decompressing arbitrary RF and gradient shapes.
	 *
	 * @return true if successful
	 */
	bool decodeBlock(SeqBlock *block);

	/**
	 * @brief Decode only the ext trap part of the block
	 *
	 * This involves decompressing extended gradient shapes.
	 *
	 * @return true if successful
	 */
	bool decodeExtTrapGradInBlock(SeqBlock *block);

	/**
	 * @brief Return `true` if block has a gradient which starts at a non-zero value on given channel (only possible for arbitrary or ExtTrap gradients)
	 */
	bool isGradientInBlockStartAtNonZero(SeqBlock *block, int channel);

	/**
	 * @brief Return `true` if block has no gradients which start at a non-zero value on any channel (only possible for arbitrary or ExtTrap gradients)
	 */
	bool isAllGradientsInBlockStartAtZero(SeqBlock *block);

	bool isSigned();
	std::string getSignature();
	std::string getSignatureType();

	bool isSignatureCheckSucceeded();

  private:

	static const int MAX_LINE_SIZE;	/**< @brief Maximum length of line */
	static const char COMMENT_CHAR;	/**< @brief Character defining the start of a comment line */

	// *** Private helper functions ***

	/**
	 * @brief Read a line from the input stream *independent of line ending*.
	 *
	 * Unlike istream::getline(), this function safely handles three
	 * different line endings:
	 *  - Unix/OSX (\\n)
	 *  - Windows (\\r\\n)
	 *  - Old Mac (\\r)
	 *
	 * @param stream the input file stream to read from
	 * @param buffer the output line buffer (null terminated)
	 * @param MAX_SIZE maximum size of the buffer
	 */
	// MZ: VC2008 / VD13 do not seem to properly handle the errors streams
	//static std::istream& 
	enum gl_ret {gl_false=0, gl_true, gl_truncated};
    static int getline(std::istream& stream, char* buffer, const int MAX_SIZE, bool bRaw=false);

	/**
	 * @brief Search the file stream for section headers e.g. [RF], [GRAD] etc
	 *
	 * Searches forward in the stream for sections enclosed in square brackets
	 * and writes to index
	 */
	void buildFileIndex(std::istream &stream);

	/**
	 * @brief Skip the comments and empty lines in the given input stream.
	 *
	 * @param stream the input file stream to process
	 * @param buffer return output buffer of next non-comment line
	 */
	void skipComments(std::istream &stream, char* buffer);

	/**
	 * @brief Strip any empty characters at th beginning or at the end of the string.
	 *
	 * @param buffer input line
	 * @return modified string
	 */
	char* stripWhiteSpace(char *buffer);

	/**
	 * @brief Decompress a run-length compressed shape
	 *
	 * @param encoded Compressed shape structure
	 * @param shape array of floating-point values (must be preallocated!)
	 */
	bool decompressShape(CompressedShape& encoded, float *shape);


	/**
	 * @brief Check the IDs contains references to valid events in the library
	 *
	 * @param  events The event IDs to check
	 * @return true if event references are ok
	 * @see checkRF(), checkGradient()
	 */
	bool checkBlockReferences(EventIDs& events);

	/**
	 * @brief Check the shapes defining the arbitrary gradient events (if present)
	 *
	 * Check the *decompressed* amplitude of the given block is
	 * between [0 1].
	 * @param  block The sequence block to check
	 * @see checkRF()
	 */
	void checkGradient(SeqBlock& block);

	/**
	 * @brief Check the shapes defining the RF event (if present)
	 *
	 * Check the *decompressed* RF amplitude of the given block is
	 * between [0 1] and the phase between [0 2pi].
	 * @param  block The sequence block to check
	 * @see checkGradient()
	 */
	void checkRF(SeqBlock& block);

	/**
	 * @brief Check the IDs contains references to valid labels in the library
	 *		  Transform order of Labels as in MHDheader
	 * @param  labels The labels IDs to check
	 * @return 0 if labels/values references are ok
	 * @return 1 if labels references are not recognized
	 * @return -1 if labels/values references are invalid
	 */
	int decodeLabel(ExtType, int&, char*, LabelEvent&);

	// *** Static helper function ***

	/**
	 * @brief Print string to standard output (with newline character)
	 */
	static void defaultPrint(const std::string &str);

	// *** Static members ***

	static PrintFunPtr print_fun;              /**< @brief Pointer to output print function */

	// *** Members ***

	int version_major;
	int version_minor;
	int version_revision;
	int version_combined;

	std::map<std::string,int> m_fileIndex;     /**< @brief File location of sections, [RF], [ADC] etc */
	std::set<int> m_fileSections;              /**< @brief File location of sections and EOF additionally */

	// Low level sequence blocks
	std::vector<EventIDs> m_blocks;            /**< @brief List of sequence blocks */

	std::vector<long> m_blockDurations_ru;     /**< @brief List of block durations expressed in duration raster units */

	// extension list storage
	//std::vector<ExtensionListEntry> m_extensions; /**< @brief the storage area of the extension list referenced by the enent table */

	// Global user-specified definitions
	std::map<std::string, std::string >m_definitions_str;  /**< @brief Custom definitions provided through [DEFINITIONS] section (stored as strings) */
	std::map<std::string, std::vector<double> >m_definitions;  /**< @brief Custom definitions provided through [DEFINITIONS] section (converted to doubles) */

	// pulseq file signature (read from file)
	std::map<std::string, std::string >m_signatureMap;  /**< @brief A hash of the internal sequence data structures stored in the [SIGNATURE] section along with auxilary pieces of information */
	bool m_bSignatureDefined;
	std::string m_strSignature;
	std::string m_strSignatureType;

	// pulseq file MD5 signature calculated during file loading
	std::string m_strCalculatedMD5Signature;
	bool m_bSignatureCheckSucceeded;

	// List of events (referenced by blocks)
	std::map<int,RFEvent>      m_rfLibrary;       /**< @brief Library of RF events */
	std::map<int,GradEvent>    m_gradLibrary;     /**< @brief Library of gradient events */
	std::map<int,ADCEvent>     m_adcLibrary;      /**< @brief Library of ADC readouts */
	std::map<int,long>         m_tmpDelayLibrary;    /**< @brief Library of delays, only used for loading older files and is cleaned immediately before load() is finished*/
	//std::map<int,ControlEvent> m_controlLibrary;  /**< @brief Library of control commands */
	std::map<int,ExtensionListEntry> m_extensionLibrary;  /**< @brief Library of extension list entries */
	std::map<int,std::pair<std::string,int> > m_extensionNameIDs; /**< @brief Map of extension IDs from the file to textIDs and internal known numeric IDs*/
	std::map<int, TriggerEvent>   m_triggerLibrary;   /**< @brief Library of trigger events */
	std::map<int, RotationEvent>  m_rotationLibrary;  /**< @brief Library of rotation events */
	std::map<int, LabelEvent>     m_labelsetLibrary;  /**< @brief Library of labelset events */
	std::map<int, LabelEvent>     m_labelincLibrary;  /**< @brief Library of labelinc events */
	std::map<int, SoftDelayEvent> m_softDelayLibrary; /**< @brief Library of soft delay events */
    LabelMap                      m_labelMap;        /**< @brief labelMap is useful for loading labels or damping/visualising values */

	// List of basic shapes (referenced by events)
	std::map<int,CompressedShape> m_shapeLibrary;    /**< @brief Library of compressed shapes */
	// raster times
	double m_dAdcRasterTime_us; // Siemens default: 1e-07s 
	double m_dGradientRasterTime_us; // Siemens default: 1e-05s 
	double m_dRadiofrequencyRasterTime_us; // Siemens default: 1e-06s 
	double m_dBlockDurationRaster_us; // Siemens default: 1e-05s 
};

// * ------------------------------------------------------------------ *
// * Inline functions                                                   *
// * ------------------------------------------------------------------ *

inline int ExternalSequence::GetNumberOfBlocks(void){return m_blocks.size();}

inline std::vector<double>	ExternalSequence::GetDefinition(std::string key){
	if (m_definitions.count(key)>0)
		return m_definitions[key];
	else
		return std::vector<double>();
}

inline std::string	ExternalSequence::GetDefinitionStr(std::string key){
	if (m_definitions_str.count(key)>0)
		return m_definitions_str[key];
	else
		return std::string();
}

inline void ExternalSequence::defaultPrint(const std::string &str) { std::cout << str << std::endl; }
inline void ExternalSequence::SetPrintFunction(PrintFunPtr fun) { print_fun=fun; }

inline bool ExternalSequence::isSigned() { return m_bSignatureDefined; }
inline std::string ExternalSequence::getSignature() { return m_strSignature; }
inline std::string ExternalSequence::getSignatureType() { return m_strSignatureType; }

inline bool ExternalSequence::isSignatureCheckSucceeded() { return m_bSignatureDefined && m_bSignatureCheckSucceeded; }

#endif	//_EXTERNAL_SEQUENCE_H_

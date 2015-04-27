/** @file ExternalSequence.h */

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

#ifndef CLASS_EXTERNAL_SEQUENCE
#define CLASS_EXTERNAL_SEQUENCE

//#define MASTER_SLAVE_FORMAT

#define PI     3.1415926535897932384626433832795
#define TWO_PI 6.283185307179586476925286766558


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
#define MSG_LEVEL NORMAL_MSG


/**
 * @brief Internal storage order
 */
enum Event {
	RF_EVENT,
	XGRAD_EVENT,
	YGRAD_EVENT,
	ZGRAD_EVENT,
	ADC_EVENT
};


/**
 * @brief RF event data
 *
 * Stores the amplitude, frequency offset and an index to basic shapes for the
 * magnitude and phase
 */
struct RFEvent
{
	float amplitude;	/**< @brief Amplitude of RF event (Hz) */
	int magShape;		/**< @brief ID of shape for magnitude */
	int phaseShape;		/**< @brief ID of shape for phase */
	float freqOffset;	/**< @brief Frequency offset of transmitter (Hz) */
	float phaseOffset;	/**< @brief Phase offset of transmitter (Hz) */
};


/**
 * @brief Gradient event data
 *
 * Stores gradient amplitude and additional information
 * depending on the type:
 *  - **trapezoid:** ramp up, flat time, ramp down
 *  - **arbirary:** index to basic shape (see CompressedShape)
 */
struct GradEvent
{
	float amplitude;	/**< @brief Amplitude of gradient (Hz/m) */
	// Trapezoid:
	long rampUpTime;	/**< @brief Ramp up time of trapezoid (us) */
	long flatTime;		/**< @brief Flat-top time of trapezoid (us) */
	long rampDownTime;	/**< @brief Ramp down time of trapezoid (us) */
	// Arbitrary:
	int shape;			/**< @brief ID of shape for arbitrary gradient */
};


/**
 * @brief ADC readout event
 *
 * Stores number of ADC samples, dwell time, delay and frequency/phase offset
 * information (used to shift the FOV in-plane)
 *
 */
struct ADCEvent
{
	int numSamples;		/**< @brief Number of samples */
	int dwellTime;		/**< @brief Dwell time of ADC readout (ns) */
	int delay;			/**< @brief Delay before first sample (us) */
	float freqOffset;	/**< @brief Frequency offset of receiver (Hz) */
	float phaseOffset;	/**< @brief Phase offset of receiver (Hz) */
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
	SeqBlock() {}

	/**
	 * @brief Return `true` if block has RF event
	 */
	bool    isRF();

	/**
	 * @brief Return `true` if block has trapezoid event on given channel
	 */
	bool    isTrapGradient(int channel);

	/**
	 * @brief Return `true` if block has arbitrary event on given channel
	 */
	bool    isArbitraryGradient(int channel);

	/**
	 * @brief Return `true` if block has ADC readout event
	 */
	bool    isADC();

	/**
	 * @brief Return `true` if block has delay
	 */
	bool    isDelay();

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
	long    GetDelay();

	/**
	 * @brief Return duration of block
	 */
	long    GetDuration();

	/**
	 * @brief Return the number of samples of the given gradient channel.
	 * Only relevant for arbitrary gradients
	 */
	int     GetGradientLength(int channel);
	
	/**
	 * @brief Return the gradient event of the given channel
	 */
	GradEvent GetGradEvent(int channel);

	/**
	 * @brief Directly get a pointer to the samples of the arbitrary X gradient
	 */
	float* GetXGradientPtr();

	/**
	 * @brief Directly get a pointer to the samples of the arbitrary Y gradient
	 */
	float* GetYGradientPtr();

	/**
	 * @brief Directly get a pointer to the samples of the arbitrary Z gradient
	 */
	float* GetZGradientPtr();

	/**
	 * @brief Return the RF event
	 */
	RFEvent	GetRFEvent();

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
	 * @brief Return the ADC event
	 */
	ADCEvent GetADCEvent();

	/**
	 * @brief Return a brief string description of the block
	 */
	std::string GetTypeString();

	/**
	 * @brief Free the memory associated with decompressed shapes of this block
	 */
	void   free();
protected:
	int index;			/**< @brief Index of this block */
	
	// Event array contains integer indices to events stored in the parent ExternalSequence object
	int events[5];		/**< @brief list of event indices (RF, GX, GY, GZ, ADC) */

	long delay;			/**< @brief delay of this block (in us) */
	long duration;		/**< @brief duration of this block (in us) used for error checking */

	RFEvent rf;			/**< @brief RF event */
	GradEvent grad[3];	/**< @brief gradient events */
	ADCEvent adc;		/**< @brief ADC event  */

	// Below is only valid once decompressed
		
	// RF
	std::vector<float> rfAmplitude;	/**< @brief RF amplitude shape (uncompressed) */
	std::vector<float> rfPhase;		/**< @brief RF phase shape (uncompressed) */

	// Grad
	std::vector<float> grad_x;		/**< @brief Arbitrary gradient shape for X channel (uncompressed) */
	std::vector<float> grad_y;		/**< @brief Arbitrary gradient shape for Y channel (uncompressed) */
	std::vector<float> grad_z;		/**< @brief Arbitrary gradient shape for Z channel (uncompressed) */

};

// * ------------------------------------------------------------------ *
// * Inline functions                                                   *
// * ------------------------------------------------------------------ *
inline int       SeqBlock::GetIndex() { return index; }

inline bool      SeqBlock::isRF() { return (events[RF_EVENT]>0); }
inline bool      SeqBlock::isTrapGradient(int channel) { return ((events[channel+1]>0) & (grad[channel].shape==0)); }
inline bool      SeqBlock::isArbitraryGradient(int channel) { return ((events[channel+1]>0) & (grad[channel].shape>0)); }
inline bool      SeqBlock::isADC() { return (events[ADC_EVENT]>0); }
inline bool      SeqBlock::isDelay() { return (delay>0); }

inline long      SeqBlock::GetDelay() { return delay; }
inline long      SeqBlock::GetDuration() { return duration; }	/* maybe not needed? */

inline GradEvent SeqBlock::GetGradEvent(int channel) { return grad[channel]; }
inline RFEvent   SeqBlock::GetRFEvent() { return rf; }
inline ADCEvent  SeqBlock::GetADCEvent() { return adc; }

inline int       SeqBlock::GetEventIndex(Event type) { return events[type]; }


inline int       SeqBlock::GetGradientLength(int channel) {
	if (channel==0)
		return grad_x.size();
	else if (channel==1)
		return grad_y.size();
	else if (channel==2)
		return grad_z.size();
	else
		return 0;
} 

inline std::string       SeqBlock::GetTypeString() {
	std::string type;
	if (isRF()) type = type + "RF";
	if (isTrapGradient(0)) type += " TrapX";
	if (isTrapGradient(1)) type += " TrapY";
	if (isTrapGradient(2)) type += " TrapZ";
	if (isArbitraryGradient(0)) type += " ArbX";
	if (isArbitraryGradient(1)) type += " ArbY";
	if (isArbitraryGradient(2)) type += " ArbZ";
	if (isADC()) type += " ADC";
	if (GetDelay()>0) type += " Delay";
	return type;
} 

inline float*    SeqBlock::GetXGradientPtr() { return (grad_x.size()>0) ? &grad_x[0] : NULL; }
inline float*    SeqBlock::GetYGradientPtr() { return (grad_y.size()>0) ? &grad_y[0] : NULL; }
inline float*    SeqBlock::GetZGradientPtr() { return (grad_z.size()>0) ? &grad_z[0] : NULL; }

inline float*    SeqBlock::GetRFAmplitudePtr() { return &rfAmplitude[0]; }
inline float*    SeqBlock::GetRFPhasePtr() { return &rfPhase[0]; }
inline int       SeqBlock::GetRFLength() { return rfAmplitude.size(); }

inline void      SeqBlock::free() {
	// Force the memory to be freed
	std::vector<float>().swap(rfAmplitude);
	std::vector<float>().swap(rfPhase);
	std::vector<float>().swap(grad_x);
	std::vector<float>().swap(grad_y);
	std::vector<float>().swap(grad_z);
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
	int numUncompressedSamples;		/**< @brief Number of samples *after* decompression */
	std::vector<float> samples;		/**< @brief Compressed samples */
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
	 * @brief Contructor
	 */
	ExternalSequence();
	
	/**
	 * @brief Destructor
	 */
	~ExternalSequence();
	
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
	 * @brief Display an output message
	 *
	 * Display a message only if the MSG_LEVEL is sufficiently high.
	 * This function calls the low-level output function, which can be overriden
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
	 * @brief Return the Scan ID set in the [DEFINITIONS] block
	 */
	int   GetScanID(void);

	/**
	 * @brief Return the FOV set in the [DEFINITIONS] block
	 */
	float GetFOV(int dim);

	/**
	 * @brief Return the FOV set in the [DEFINITIONS] block
	 */
	float GetSliceThickness(void);

	/**
	 * @brief Return the resoution set in the [DEFINITIONS] block
	 */
	int   GetBaseResolution(void);

	/**
	 * @brief Return the TE set in the [DEFINITIONS] block
	 */
	float GetTE(void);

	/**
	 * @brief Return the TR set in the [DEFINITIONS] block
	 */
	float GetTR(void);

	/**
	 * @brief Return an element of the rotation matrix set in the [DEFINITIONS] block
	 */
	double GetRotMatrix(int i, int j);

	/**
	 * @brief Return the gradient delay of a given channel as set in the [DEFINITIONS] block
	 */
	long  GetDelay(int channel);

	/**
	 * @brief Return the maximum gradient delay
	 */
	long  GetMaxDelay(void);

#ifdef MASTER_SLAVE_FORMAT
	/**
	 * @brief Return the Scan ID set in the [DEFINITIONS] block
	 */
	void  SetSlave(bool isSlave);
#endif
	
	/**
	 * @brief Return number of sequence blocks
	 */
	int  GetNumberOfBlocks(void);

	/**
	 * @brief Return a pointer to a sequence block
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

  private:
	
	static const int MAX_LINE_SIZE;	/**< @brief Maximum length of line */
	static const char COMMENT_CHAR;	/**< @brief Character defining the start of a comment line */

	// *** Private helper functions ***

	/**
	 * @brief Search the file stream for section headers e.g. [RF], [GRAD] etc
	 *
	 * Searches forward in the stream for sections enclosed in square brackets
	 * and writes to index
	 */
	void buildFileIndex(std::ifstream &stream);

	/**
	 * @brief Skip the comments and empty lines in the given input stream.
	 *
	 * @param stream the input file stream to process
	 * @param buffer return output buffer of next non-comment line
	 */
	void skipComments(std::ifstream &stream, char* buffer);
	
	//std::istream& safeGetLine(std::istream& is, std::string& t);

	/**
	 * @brief Decompress a run-length compressed shape
	 *
	 * @param encoded Compressed shape structure
	 * @param shape array of floating-point values (must be preallocated!)
	 */
	bool decompressShape(CompressedShape& encoded, float *shape);

	/**
	 * @brief Check the block contains references to defined events
	 *
	 * @param  block The sequence block to check
	 * @return true if block references are ok
	 * @see checkRF(), checkGradient()
	 */
	bool checkBlockReferences(SeqBlock& block);

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
	
	// *** Static helper function ***

	/**
	 * @brief Print string to standard output (with newline character)
	 */
	static void defaultPrint(const std::string &str);

	// *** Static members ***

	static PrintFunPtr print_fun;			/**< @brief Pointer to output print function */

	// *** Members ***

	std::map<std::string,int> m_fileIndex;	/**< @brief File location of sections, [RF], [ADC] etc */

	// Members to store parameters provided through the [DEFINITIONS] section
	int				m_BaseResolution[3];	/**< @brief Resolution (if provided through [DEFINITIONS] section) */
	float			m_FOV[3];				/**< @brief Field-of-view (if provided through [DEFINITIONS] section) */
	float			m_sliceThickness;		/**< @brief Slice thickness (if provided through [DEFINITIONS] section) */
	float			m_TE;					/**< @brief Echo time (if provided through [DEFINITIONS] section) */
	float			m_TR;					/**< @brief Repetition time (if provided through [DEFINITIONS] section) */
	int				m_scanID;				/**< @brief Scan ID (if provided through [DEFINITIONS] section) */
	double			m_RotMatrix[3][3];		/**< @brief Rotation matrix (if provided through [DEFINITIONS] section) */
	long			m_delays[3];			/**< @brief Gradient delays (if provided through [DEFINITIONS] section) */
	long			m_maxDelay;				/**< @brief Maximum gradient delay*/

#ifdef MASTER_SLAVE_FORMAT
	bool			m_isSlave;				/**< @brief True if running on the slave host/mpcu */
#endif

	// Low level sequence blocks
	std::vector<SeqBlock> m_blocks;	/**< @brief List of sequence blocks */

	// List of events (referenced by blocks)
	std::map<int,RFEvent>   m_rfLibrary;	/**< @brief Library of RF events */
	std::map<int,GradEvent> m_gradLibrary;	/**< @brief Library of gradient events */
	std::map<int,ADCEvent>  m_adcLibrary;	/**< @brief Library of ADC readouts */
	std::map<int,long>      m_delayLibrary;	/**< @brief Library of delays */
	
	// List of basic shapes (referenced by events)
	std::map<int,CompressedShape> m_shapeLibrary;	/**< @brief Library of compressed shapes */
};

// * ------------------------------------------------------------------ *
// * Inline functions                                                   *
// * ------------------------------------------------------------------ *
inline int		ExternalSequence::GetScanID(void){ return m_scanID; }
inline float	ExternalSequence::GetFOV(int dim){ return m_FOV[dim]; }
inline int		ExternalSequence::GetBaseResolution(void){ return m_BaseResolution[0]; }
inline float	ExternalSequence::GetSliceThickness(void){ return m_sliceThickness; }
inline float	ExternalSequence::GetTE(void){ return m_TE; }
inline float	ExternalSequence::GetTR(void){ return m_TR; }
inline int		ExternalSequence::GetNumberOfBlocks(void){return m_blocks.size();}
inline double	ExternalSequence::GetRotMatrix(int i, int j){return m_RotMatrix[i][j];}
inline long		ExternalSequence::GetDelay(int channel){return m_delays[channel];}
inline long		ExternalSequence::GetMaxDelay(void){return m_maxDelay;}
#ifdef MASTER_SLAVE_FORMAT
inline void		ExternalSequence::SetSlave(bool isSlave){ m_isSlave = isSlave; }
#endif

inline SeqBlock*	ExternalSequence::GetBlock(int index) { return &m_blocks[index]; }

inline void ExternalSequence::defaultPrint(const std::string &str) { std::cout << str << std::endl; }
inline void ExternalSequence::SetPrintFunction(PrintFunPtr fun) { print_fun=fun; }



#endif	//CLASS_EXTERNAL_SEQUENCE

#include "ExternalSequence.h"
extern "C" {
	#include "md5.h"
}

#include <stdio.h>		// sscanf
#include <cstring>		// strlen etc
#include <iomanip>		// std::setw etc

#include <algorithm>	// for std::max_element
#include <functional>	// for std::bind...

#include <math.h>		// fabs etc

#include <assert.h>		// assert (TODO: remove in the released version)

ExternalSequence::PrintFunPtr ExternalSequence::print_fun = &ExternalSequence::defaultPrint;
const int ExternalSequence::MAX_LINE_SIZE = 256;
const char ExternalSequence::COMMENT_CHAR = '#';
std::string& str_trim(std::string& str);
std::string str_tolower(std::string str);
double SeqBlock::s_blockDurationRaster = 10.0;

/***********************************************************/
ExternalSequence::ExternalSequence()
{
	m_blocks.clear();	// probably not needed.
	version_major=0;
	version_minor=0;
	version_revision=0;
	version_combined=0;
	m_bSignatureDefined=false;
	m_bSignatureCheckSucceeded=false;
}

/***********************************************************/
ExternalSequence::~ExternalSequence(){}

/***********************************************************/
void ExternalSequence::print_msg(MessageType level, std::ostream& ss) {
	if (MSG_LEVEL>=level) {
#if defined(VXWORKS) || defined (BUILD_PLATFORM_LINUX)
		// we skip messages on the scanner platforms due to performance limitations
		// we could trivially use UTRACE on newer scanners, but it is not compatible with older platforms
#else		
		std::ostringstream oss;
		oss.width(2*(level-1)); oss << "";
		oss << static_cast<std::ostringstream&>(ss).str();
		print_fun(oss.str().c_str());
#endif
	}
}

/***********************************************************/
void ExternalSequence::reset()
{
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "Resetting the sequence state");

	// reset the internal structures in case this is a repeated call to load()
	m_adcLibrary.clear();
	m_blockDurations_ru.clear();
	m_blocks.clear();
	m_bSignatureDefined=false;
	// not on vb17 // m_strSignature.clear();
	m_strSignature="";
	m_bSignatureCheckSucceeded=false;
	// not on vb17 // m_strCalculatedMD5Signature.clear();
	m_strCalculatedMD5Signature="";
	m_definitions.clear();
	m_definitions_str.clear();
	m_extensionLibrary.clear();
	m_extensionNameIDs.clear();
	m_fileIndex.clear();
	m_fileSections.clear();
	m_gradLibrary.clear();
	m_labelincLibrary.clear();
	m_labelsetLibrary.clear();
	m_rfLibrary.clear();
	m_rotationLibrary.clear();
	m_shapeLibrary.clear();
	m_signatureMap.clear();
	m_strSignature="";
	m_strSignatureType="";
	m_triggerLibrary.clear();
	version_combined=0;
}

/***********************************************************/
bool ExternalSequence::load(std::string path)
{
	//reset(); // moved to the stream loader

	// now time to read...
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "Reading external sequence files");

	// **********************************************************************************************************************
	// ************************ READ SHAPES ***********************************

	// Try single file mode (everything in a single .seq file)
	std::ifstream data_file;
	bool isSingleFileMode = true;
	std::string filepath = path;
	if (filepath.substr(filepath.size()-4) != std::string(".seq")) {
		filepath = path + PATH_SEPARATOR + "external.seq";
	}
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "Opening " << filepath);			
	// Open in binary mode to ensure all end-of-line characters are processed
	data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
	data_file.seekg(0, std::ios::beg);

	if (!data_file.is_open() || !data_file.good())
	{
		// Try separate file mode (blocks.seq, events.seq, shapes.seq)
		// not really logical, but the current code expects [VERSION] to be defined in every file, 
		// [DEFINITIONS] to be in the blocks.seq and only blocks.seq is signed... Not logical at all... TODO: FixMe!

		reset();

		filepath = path + PATH_SEPARATOR + "shapes.seq";
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);
		
		if (!data_file.is_open() || !data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}

		if (!load(data_file, lm_shapes)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to load shapes from file " << filepath);
			return false;
		}

		filepath = path + PATH_SEPARATOR + "events.seq";
		data_file.close();
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.is_open() || !data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}

		if (!load(data_file, lm_events)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to load events from file " << filepath);
			return false;
		}

		filepath = path + PATH_SEPARATOR + "blocks.seq";
		data_file.close();
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.is_open() || !data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}
		
		if (!load(data_file, lm_blocks)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to load blocks from file " << filepath);
			return false;
		}

		return true;
	}
	else
	{
		return load(data_file);
	}
}

bool ExternalSequence::load_from_buffer(char * buffer) 
{
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "Loading sequence from a text buffer");

	// Try single file mode (everything in a single .seq file)
	std::istringstream string_stream(buffer);
	string_stream.seekg(0, std::ios::beg);
	return load(string_stream);
}

bool ExternalSequence::load(std::istream& data_stream, load_mode loadMode /*=lm_singlefile*/)
{
	if (loadMode == lm_singlefile)
	{
		reset();
	}

	if (!data_stream.good())
	{
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Function load() failed to read from the stream provided");
		return false;
	}

	char buffer[MAX_LINE_SIZE];
	char tmpStr[MAX_LINE_SIZE];

	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Building index" );

	// Save locations of section tags
	buildFileIndex(data_stream);

	// Read version section
	if (m_fileIndex.find("[VERSION]") != m_fileIndex.end()) {
		print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "decoding VERSION section");
		// Version is a recommended but not a compulsory section
		// very basic reading code, repeated keywords will overwrite previous values, no serious error checking
		data_stream.seekg(m_fileIndex["[VERSION]"], std::ios::beg);
		skipComments(data_stream,buffer);			// load up some data and ignore comments & empty lines
		while (data_stream.good() && buffer[0]!='[')
		{
			//print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "buffer: \n" << buffer << std::endl );
			if (0==strncmp(buffer,"major",5)) {
				    if (1!=sscanf(buffer+5, "%d", &version_major)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode version_major");
					return false;
				}
			    print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "major=" << version_major);		
			} else if (0==strncmp(buffer,"minor",5)) {
				if (1!=sscanf(buffer+5, "%d", &version_minor)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode version_minor");
					return false;
				}
				print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "minor=" << version_minor);
			}
			else if (0==strncmp(buffer,"revision",8)) {
				if (1!=sscanf(buffer+8, "%d", &version_revision)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode version_revision \n" << buffer << std::endl );
					return false;
				}
				print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "revision=" << version_revision);
			}
			else
			{
				print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unknown field in the [VERSION] block");
				return false;
			}
			//getline(data_stream, buffer, MAX_LINE_SIZE);
			skipComments(data_stream,buffer);			// load up some data and ignore comments & empty lines
		}
		version_combined=version_major*1000000L+version_minor*1000L+version_revision;
	}
	else
	{
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: supported Pulseq files MUST contain the [VERSION] section");
		return false;
	}

	if (version_combined<1002000L)
	{
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: unsupported Pulseq file version " << version_combined << ". The oldest supported version is 1.2.0.");
		return false;
	}
	
	// **********************************************************************************************************************
	// ************************ READ SHAPES ***********************************

	if (loadMode == lm_singlefile || loadMode == lm_shapes) 
	{
		// Read shapes section
		// ------------------------
		m_shapeLibrary.clear();
		if (m_fileIndex.find("[SHAPES]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[SHAPES]"], std::ios::beg);
			skipComments(data_stream,buffer);			// Ignore comments & empty lines

			int shapeId, numSamples;
			float sample;

			while (data_stream.good() && buffer[0]=='s')
			{
				if (2!=sscanf(buffer, "%s%d", tmpStr, &shapeId)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode 'shapeId'\n" << buffer << std::endl );
					return false;
				}
				getline(data_stream, buffer, MAX_LINE_SIZE);
				if (2!=sscanf(buffer, "%s%d", tmpStr, &numSamples)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode 'numSamples'\n" << buffer << std::endl );
					return false;
				}

				//print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Reading shape " << shapeId );

				CompressedShape shape;
				shape.samples.clear();
				while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
					if (buffer[0]=='s' || strlen(buffer)==0) {
						break;
					}
					if (1!=sscanf(buffer, "%f", &sample)) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode 'sample'\n" << buffer << std::endl );
						return false;
					}
					shape.samples.push_back(sample);
				}
				// number of samples equal to the data length is used as a non-compressed flag
				// but only for v1.4.0 or above
				if (version_combined >= 1004000 && numSamples==shape.samples.size())
					shape.isCompressed=false;
				else 
					shape.isCompressed=true;
				shape.numUncompressedSamples=numSamples;

				print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape index " << shapeId << " has " << shape.samples.size()
					<< " compressed and " << shape.numUncompressedSamples << " uncompressed samples" );

				m_shapeLibrary[shapeId] = shape;

				skipComments(data_stream,buffer);			// Ignore comments & empty lines
			}
			data_stream.clear();	// In case EOF reached

			print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- SHAPES READ numShapes: " << m_shapeLibrary.size() );
		}
		else
		{
			print_msg(NORMAL_MSG, std::ostringstream().flush() << "-- No SHAPES section found, which is permisible but unusual" );
		}
	}

	// **********************************************************************************************************************
	// ************************ READ EVENTS ********************
	
	if (loadMode == lm_singlefile || loadMode == lm_events) 
	{
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "reading RF section");
		// Read RF section
		// ------------------------
		if (m_fileIndex.find("[RF]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[RF]"], std::ios::beg);

			int rfId;
			m_rfLibrary.clear();
			while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				RFEvent event;
				if (version_combined<1004000L)
				{
					// pre v1.4.0
					if (7!=sscanf(buffer, "%d%f%d%d%d%f%f", &rfId, &(event.amplitude),
								&(event.magShape),&(event.phaseShape), &(event.delay),
								&(event.freqOffset), &(event.phaseOffset)
								)) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode RF event\n" << buffer << std::endl );
						return false;
					}
                    event.freqPPM=0.0;
                    event.phasePPM = 0.0;
                    event.timeShape = 0;
					event.use='u'; // undefined use
					event.center=-1.0; // mark as invalid
				}
				else if (version_combined < 1005000) 
				{
					// 1.4.0
					if (8!=sscanf(buffer, "%d%f%d%d%d%d%f%f", &rfId, &(event.amplitude),
								&(event.magShape), &(event.phaseShape), &(event.timeShape),
								&(event.delay), &(event.freqOffset), &(event.phaseOffset)
								)) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode RF event\n" << buffer << std::endl );
						return false;
					}
                    event.freqPPM = 0.0;
                    event.phasePPM = 0.0;
                    event.use       = 'u';  // undefined use
					event.center=-1.0; // mark as invalid
				}
				else 
				{
					// 1.5.0
					if (12!=sscanf(buffer, "%d%f%d%d%d%f%d%f%f%f%f %c", &rfId, &(event.amplitude),
								&(event.magShape),&(event.phaseShape), &(event.timeShape), &(event.center),
								&(event.delay), &(event.freqPPM), &(event.phasePPM), &(event.freqOffset), &(event.phaseOffset), &(event.use)
								)) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode RF event\n" << buffer << std::endl );
						return false;
					}
				}
				m_rfLibrary[rfId] = event;
                ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "m_rfLibrary["<<rfId<<"].use="<<event.use);
			}
		}
		
		// Read *arbitrary* gradient section
		// -------------------------------
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "reading arbitrary gradient section");
		m_gradLibrary.clear();
		if (m_fileIndex.find("[GRADIENTS]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[GRADIENTS]"], std::ios::beg);

			while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "got line: " << buffer);
				int gradId;
				GradEvent event;
				if ( version_combined>=1005000L )
				{
					// v1.5.0
					if (7!=sscanf(buffer, "%d%f%f%f%d%d%d", &gradId, &(event.amplitude), &(event.first), &(event.last), &(event.waveShape), &(event.timeShape), &(event.delay))) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode v1.4.x gradient event\n" << buffer << std::endl );
						return false;
					}
				}
				else if ( version_combined>=1004000L )
				{
					// v1.4.0
					if (5!=sscanf(buffer, "%d%f%d%d%d", &gradId, &(event.amplitude), &(event.waveShape), &(event.timeShape), &(event.delay))) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode v1.4.x gradient event\n" << buffer << std::endl );
						return false;
					}
					event.first=FLOAT_UNDEFINED; // std::numeric_limits<float>::quiet_NaN(); <- did not work with older MSVC
					event.last=FLOAT_UNDEFINED; // std::numeric_limits<float>::quiet_NaN(); <- did not work with older MSVC
				}
				else
				{
					// pre v1.4.0
					event.timeShape=0;
					if (4!=sscanf(buffer, "%d%f%d%d", &gradId, &(event.amplitude), &(event.waveShape), &(event.delay))) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode v1.2.x gradient event\n" << buffer << std::endl );
						return false;
					}
					event.first=FLOAT_UNDEFINED; // std::numeric_limits<float>::quiet_NaN(); <- did not work with older MSVC
					event.last=FLOAT_UNDEFINED; // std::numeric_limits<float>::quiet_NaN(); <- did not work with older MSVC
				}
				print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "assigning the event to the library under the ID " << gradId);
				m_gradLibrary[gradId] = event;
				print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "done");
			}
		}

		// Read *trapezoid* gradient section
		// -------------------------------
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "reading trapezoids section");
		if (m_fileIndex.find("[TRAP]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[TRAP]"], std::ios::beg);

			while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				int gradId;
				GradEvent event;
				if (6!=sscanf(buffer, "%d%f%ld%ld%ld%d", &gradId, &(event.amplitude),
					&(event.rampUpTime),&(event.flatTime),&(event.rampDownTime),&(event.delay))) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode trapezoid gradient entry" << buffer << std::endl );
					return false;
				}					
				event.waveShape=0;
				event.timeShape=0;
				m_gradLibrary[gradId] = event;
			}
		}

		// Sort gradients based on index
		// -----------------------------
		//std::sort(m_gradLibrary.begin(),m_gradLibrary.end(),compareGradEvents);

		// Read ADC section
		// -------------------------------
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "reading ADC section");
		if (m_fileIndex.find("[ADC]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[ADC]"], std::ios::beg);

			int adcId;
			m_adcLibrary.clear();
			while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				ADCEvent event;
				if ( version_combined>=1005000L )
				{
					// v1.5.0
					if (9!=sscanf(buffer, "%d%d%d%d%f%f%f%f%d", &adcId, &(event.numSamples),
								&(event.dwellTime),&(event.delay),&(event.freqPPM),&(event.phasePPM),&(event.freqOffset),&(event.phaseOffset),&(event.phaseModulationShape))) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode ADC event\n" << buffer << std::endl );
						return false;
					}
				}
				else
				{
					// v1.4.0 and older
					if (6!=sscanf(buffer, "%d%d%d%d%f%f", &adcId, &(event.numSamples),
								&(event.dwellTime),&(event.delay),&(event.freqOffset),&(event.phaseOffset))) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode ADC event\n" << buffer << std::endl );
						return false;
					}
                    event.freqPPM=0.0; // no ppmOffset in older formats
                    event.phasePPM= 0.0; 
					event.phaseModulationShape=0; // no phase modulation shape provided 
				}
				
				m_adcLibrary[adcId] = event;
			}
		}

		// Read delays section (comatibility with Pulseq version prior to 1.4.0)
		// ---------------------------------------------------------------------
		//std::map<int,long> tmpDelayLibrary;
		m_tmpDelayLibrary.clear();
		if (m_fileIndex.find("[DELAYS]") != m_fileIndex.end()) {
			print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "reading DELAYS section (compatibility)");
			data_stream.seekg(m_fileIndex["[DELAYS]"], std::ios::beg);

			int delayId;
			long delay;
			while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				if (2!=sscanf(buffer, "%d%ld", &delayId, &delay)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode delay event\n" << buffer << std::endl );
					return false;
				}
				m_tmpDelayLibrary[delayId] = delay;
			}
		}

		// Read extensions section
		// -------------------------------
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "reading and processing extensions");
		m_extensionLibrary.clear();
		m_extensionNameIDs.clear();
		m_triggerLibrary.clear(); // clear also all known extension libraries
		m_labelsetLibrary.clear();
		m_labelincLibrary.clear();
		std::map<std::string,int>::iterator itFI = m_fileIndex.find("[EXTENSIONS]");
		if ( itFI != m_fileIndex.end()) {
			data_stream.seekg(itFI->second, std::ios::beg);
			std::set<int>::iterator itSFI = m_fileSections.find(itFI->second);
			if ( itSFI==m_fileSections.end() ||
				 (++itSFI)==m_fileSections.end() )
			{
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed find the end of the section while reading EXTENSIONS");
				return false;
			}
			int sectionEnd = *itSFI;
			// we first read in the extension list
			int nID;
			int nExtensionID=EXT_LIST; // EXT_LIST means we are reading the extension list
			while ( data_stream.tellg()<sectionEnd &&
					getline(data_stream, buffer, MAX_LINE_SIZE)) 
			{
				if (buffer[0]=='#' || buffer[0]=='[' || strlen(buffer)==0) {
					continue;
				}
				print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "input line: " << buffer);
				if (0==strncmp(buffer,"extension",9)) {
					// read new extension ID from the header
					char szStrID[MAX_LINE_SIZE];
					int nInternalID=0;
					int nKnownID=EXT_UNKNOWN;
					if (2!=sscanf(buffer, "extension %s %d", szStrID, &nInternalID)) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode extension header entry\n" << buffer << std::endl );
						return false;
					}
					// here is the list if extensions we currently recognize
					if (0==strcmp("TRIGGERS",szStrID))
						nKnownID=EXT_TRIGGER;
					else if (0==strcmp("ROTATIONS",szStrID))
						nKnownID=EXT_ROTATION;
					else if (0==strcmp("LABELSET",szStrID))
						nKnownID=EXT_LABELSET;
					else if (0==strcmp("LABELINC",szStrID))
						nKnownID=EXT_LABELINC;
					else if (0==strcmp("DELAYS",szStrID))
						nKnownID=EXT_DELAY;
					else if (0==strcmp("RF_SHIMS",szStrID))
						nKnownID=EXT_RF_SHIM;
					if (nKnownID!=EXT_UNKNOWN)
						m_extensionNameIDs[nInternalID]=std::make_pair(std::string(szStrID),nKnownID);
					else {
						print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unknown extension ignored\n" << buffer << std::endl );
					}
					nExtensionID=nKnownID;
				}
				else
				{
					ExtensionListEntry extEntry;
					TriggerEvent trigger;
					RotationEvent rotation;
					SoftDelayEvent delay;
                    RfShimmingEvent rfShim;
					int  nVal;					   // read label set/inc values from label set/inc extension
					int  nRet;                     // conversion result / return value
					char szLabelID[MAX_LINE_SIZE]; // read labels strings from label set/inc extension
					LabelEvent	label;			   // write label event
					switch (nExtensionID) {
						case EXT_LIST: 
							if (4!=sscanf(buffer, "%d%d%d%ld", &nID, &(extEntry.type), &(extEntry.ref), &(extEntry.next))) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode extension list entry\n" << buffer << std::endl );
								return false;
							}
							print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "decoding extension list entry " << buffer);
							m_extensionLibrary[nID] = extEntry;
							print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "nID:" << nID << " type:" << extEntry.type << " ref" << extEntry.ref << " next:" << extEntry.next);
							break;
						case EXT_TRIGGER: 
							if (5!=sscanf(buffer, "%d%d%d%ld%ld", &nID, &(trigger.triggerType), &(trigger.triggerChannel), &(trigger.delay), &(trigger.duration))) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode trigger event\n" << buffer << std::endl );
								return false;
							}
							m_triggerLibrary[nID] = trigger;
							break;
						case EXT_ROTATION: 
							if (5!=sscanf(buffer, "%d%lf%lf%lf%lf", &nID, &rotation.rotQuaternion[0], &rotation.rotQuaternion[1], &rotation.rotQuaternion[2], &rotation.rotQuaternion[3])) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode rotation event\n" << buffer << std::endl );
                                return false;
							}
                            {
								double dNorm=sqrt(rotation.rotQuaternion[0]*rotation.rotQuaternion[0]+rotation.rotQuaternion[1]*rotation.rotQuaternion[1]+rotation.rotQuaternion[2]*rotation.rotQuaternion[2]+rotation.rotQuaternion[3]*rotation.rotQuaternion[3]);
                                if (fabs(dNorm-1.0)>1e-3) {
									print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: rotation extension loaded a non-normalized quaternion " << buffer << std::endl );
									return false;
								}
								for (int i = 0; i < 4; ++i)
                                    rotation.rotQuaternion[i] /= dNorm; 
                            }
							rotation.defined=true;
							m_rotationLibrary[nID] = rotation; 
							break;
						case EXT_LABELSET: 
							if (3!=sscanf(buffer, "%d%d%s", &nID, &nVal, szLabelID)) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to load labelset event\n" << buffer << std::endl );
								return false;
							}
							nRet = decodeLabel(EXT_LABELSET,nVal,szLabelID,label);
							if (nRet<0) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode labelset event\n" << buffer << std::endl );
								return false;
							}else if(nRet>0) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** decoding labelset event returned 0\n" << buffer << std::endl );
							} 
							m_labelsetLibrary[nID] = label;
							break;
						case EXT_LABELINC: 
							if (3!=sscanf(buffer, "%d%d%s", &nID, &nVal, szLabelID)) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode labelinc event\n" << buffer << std::endl );
								return false;
							}
							nRet = decodeLabel(EXT_LABELINC,nVal,szLabelID,label);
							if (nRet<0) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode labelinc event\n" << buffer << std::endl );
								return false;
							}else if(nRet>0) {
								print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: decoding labelinc event returnd 0\n" << buffer << std::endl );
							}

							m_labelincLibrary[nID] = label;
							break;
						case EXT_DELAY: 
							{
								int n=0;
								if (4!=sscanf(buffer, "%d%d%d%d%n", &nID, &(delay.numID), &(delay.offset), &(delay.factor), &n)) {
									print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode soft delay event\n" << buffer << std::endl );
									return false;
								}
								strncpy(delay.hint,stripWhiteSpace(buffer+n),SOFT_DELAY_HINT_LENGTH);
								delay.hint[SOFT_DELAY_HINT_LENGTH-1]=0;
								print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "decoded soft delay " << delay.numID << " with the hint:" << delay.hint);
								m_softDelayLibrary[nID] = delay;
							}
							break;
                        case EXT_RF_SHIM:
                            {
								int nPos=0; 
								if (2!= sscanf(buffer, "%d%d%n", &rfShim.id, &(rfShim.nchan), &(nPos)))
                                {
                                    print_msg(
                                        ERROR_MSG,
                                        std::ostringstream().flush() << "*** ERROR: failed to decode RF shim event\n"
                                                                     << buffer << std::endl);
                                    return false;
                                }
								print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "decoded initial part of the RF shim event, id:" << rfShim.id << " nChan:" << rfShim.nchan << " curren position:" << nPos);
								rfShim.amplitudes.reserve(rfShim.nchan);
								rfShim.phases.reserve(rfShim.nchan);
								int n=nPos;
								for (int i=0;i<rfShim.nchan; ++i) 
								{
									float fa,fp;
									if (2!= sscanf(buffer+n, "%f%f%n", &(fa), &(fp), &(nPos)))
									{
										print_msg(
											ERROR_MSG,
											std::ostringstream().flush() << "*** ERROR: failed to decode RF shim event for channel " << i << " in\n"
																		 << buffer << std::endl);
										return false;
									}
									print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "decoded a further part of the RF shim event, fa:" << fa << " fp:" << fp << " curren position:" << n+nPos);
									rfShim.amplitudes.push_back(fa);
									rfShim.phases.push_back(fp);								
									n+=nPos;
								}
								print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "finished decoding RF shim event");
                                m_rfShimLibrary[rfShim.id] = rfShim;
                            }
                            break;
						case EXT_UNKNOWN:
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode unknown extension event\n" << buffer << std::endl );
							break; // just ignore unknown extensions
					}
				}
			}
		}		
		
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- EVENTS READ: "
			<<" RF: " << m_rfLibrary.size()
			<<" GRAD: " << m_gradLibrary.size()
			<<" ADC: " << m_adcLibrary.size()
			<<" EXTENSIONS: " << m_extensionLibrary.size() + m_triggerLibrary.size() + m_rotationLibrary.size() + m_labelsetLibrary.size() + m_labelincLibrary.size());
	}

	// **********************************************************************************************************************
	// ************************ READ BLOCKS ********************
	if (loadMode == lm_singlefile || loadMode == lm_blocks) 
	{
		unsigned int numBlocks = 0;
		
		// Read definition section
		// ------------------------
		if (m_fileIndex.find("[DEFINITIONS]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[DEFINITIONS]"], std::ios::beg);

			// Read each definition line
			m_definitions.clear();
			m_definitions_str.clear();
			int retgl=0;
			while (retgl=getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				// this was not compatible with Numaris4 VB line (MSVC6)
				/*std::istringstream ss(buffer);
				while(retgl==gl_truncated) { // if the line was truncated read in the rest of it into the stream/buffer
					retgl=getline(data_stream, buffer, MAX_LINE_SIZE);
					ss.str(ss.str()+buffer); 
				}*/
				// stupid compatible code
				std::string stmp(buffer);
				while(retgl==gl_truncated) { // if the line was truncated read in the rest of it into the temporary string
					retgl=getline(data_stream, buffer, MAX_LINE_SIZE);
					stmp+=buffer;
				}
				char* tmp_buff1= new char[stmp.length()+1]; // old compilers like MSVC6 require such stupid conversions
				strcpy(tmp_buff1,stmp.c_str());
				print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "--- reading definitions, tmp_buff1=`"<<tmp_buff1<<"'"<<std::endl);
				std::istringstream ss(tmp_buff1);
				// delete [] tmp_buff1; // moved few lines below
				// end of stupid compatible code (except for the delete line below)
				std::string key;
				ss >> key;
				std::string str_value;
				if (std::getline(ss,str_value)) {
					str_value = str_trim(str_value);
					m_definitions_str[key] = str_value; 
					// old compilers like MSVC6 require such stupid conversions
					char* tmp_buff2= new char[str_value.length()+1];
					strcpy(tmp_buff2,str_value.c_str());
					print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "--- reading definitions(2), tmp_buff2=`"<<tmp_buff2<<"'"<<std::endl);
					std::istringstream ssv(tmp_buff2);
					double value;
					std::vector<double> values;
					while (ssv >> value) {
						//print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "v["<<values.size()<<"]="<<value);
						values.push_back(value);
					}
					delete [] tmp_buff2;
					m_definitions[key] = values;
				}
				// this 'delete' should be here because some compilers pass the buffer by reference
				delete [] tmp_buff1;
			}

			std::ostringstream out;
			out << "-- " << "DEFINITIONS READ: " << m_definitions.size() << " : ";
			for (std::map<std::string,std::vector<double> >::iterator it=m_definitions.begin(); it!=m_definitions.end(); ++it)
			{
				out<< it->first << " ";
				for (int i=0; i<it->second.size(); i++)
					out << it->second[i] << " ";
			}

			print_msg(DEBUG_HIGH_LEVEL, out);

		} // if definitions exist

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "--- Starting to interpret definitions");

		if (version_combined<1004000L)
		{
			// initialize default raster times
			m_dAdcRasterTime_us=1e-1; // Siemens default: 1e-07 
			m_dGradientRasterTime_us=10.0; // Siemens default: 1e-05 
			m_dRadiofrequencyRasterTime_us=1.0; // Siemens default: 1e-06 
			m_dBlockDurationRaster_us = m_dGradientRasterTime_us;
		}
		else
		{
			// for v1.4.x and later we REQUIRE definitions to be present
			std::vector<double> def = GetDefinition("AdcRasterTime");
			if (def.empty()){
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: definition AdcRasterTime is not present in the file");
				return false;
			}
			m_dAdcRasterTime_us=1e6*def[0];
			def = GetDefinition("GradientRasterTime");
			if (def.empty()){
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: definition GradientRasterTime is not present in the file");
				return false;
			}
			m_dGradientRasterTime_us=1e6*def[0];
			def = GetDefinition("RadiofrequencyRasterTime");
			if (def.empty()){
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: definition RadiofrequencyRasterTime is not present in the file");
				return false;
			}
			m_dRadiofrequencyRasterTime_us=1e6*def[0];
			def = GetDefinition("BlockDurationRaster");
			if (def.empty()){
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: definition BlockDurationRaster is not present in the file");
				return false;
			}
			m_dBlockDurationRaster_us=1e6*def[0];
		}
		SeqBlock::s_blockDurationRaster=m_dBlockDurationRaster_us;

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "--- Finished reading definitions, reading blocks ...");

		// Read blocks section
		// ------------------------
		if (m_fileIndex.find("[BLOCKS]") == m_fileIndex.end()) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: [BLOCKS] section");
			return false;
		}
		data_stream.seekg(m_fileIndex["[BLOCKS]"], std::ios::beg);

		int blockIdx;
		EventIDs events;

		// Read blocks
		m_blocks.clear();
		while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}

			memset(events.id, 0, NUM_EVENTS*sizeof(int));
			long dur_ru =0;

			int ret=sscanf(buffer, "%d%d%d%d%d%d%d%d", &blockIdx,
					&dur_ru,                                        // block duration
					&events.id[RF],                                 // RF
					&events.id[GX],&events.id[GY],&events.id[GZ],   // Gradients
					&events.id[ADC],                                // ADCs
					&events.id[EXT]                                 // Extensions
					);
			if (7>ret
					) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode event table entry:\n" << buffer << std::endl );
						print_msg(ERROR_MSG, std::ostringstream().flush() << "***        number of fields read: " << ret << std::endl );
				return false;
			}

			if (!checkBlockReferences(events)) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Block " << blockIdx
					<< " contains references to undefined events" );
				print_msg(ERROR_MSG, std::ostringstream().flush() << "***        RF:" << events.id[RF] << " GX:" << events.id[GX] << " GY:" << events.id[GY] << " GZ:" << events.id[GZ] << " ADC:" << events.id[ADC] << " EXT:" << events.id[EXT]);
				return false;
			}
			// Add event IDs to list of blocks
			m_blocks.push_back(events);
			m_blockDurations_ru.push_back(dur_ru); // ATTENTION, for versions prior to 1.4.0 this will contain delayIDs, we fix it below
		}

		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- BLOCKS READ: " << m_blocks.size());
		// Num_Blocks definition (if defined) is used to check the correct number of blocks are read
		if (numBlocks>0 && m_blocks.size()!=numBlocks) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Expected " << numBlocks
				<< " blocks but read " << m_blocks.size() << " blocks");
			return false;
		}
		
		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "--- Reading signature...");
		// Read signature section
		// ------------------------
		if (m_fileIndex.find("[SIGNATURE]") != m_fileIndex.end()) {
			data_stream.seekg(m_fileIndex["[SIGNATURE]"], std::ios::beg);

			// Read each signature line
			m_signatureMap.clear();
			while (getline(data_stream, buffer, MAX_LINE_SIZE)) {
				if (buffer[0]=='[' || strlen(buffer)==0) {
					break;
				}
				if (buffer[0]=='#')
					continue;
				std::istringstream ss(buffer);
				std::string key;
				ss >> key;
				std::string str_value;
				if (std::getline(ss,str_value)) {
					str_value = str_trim(str_value);
					m_signatureMap[key] = str_value; // old compilers like MSVC6 require such stupid conversions
				}
			}

			std::ostringstream out;
			out << "-- " << "SIGNATURE SECTION READ with " << m_signatureMap.size() << " entries:" << std::endl;
			for (std::map<std::string, std::string >::iterator it=m_signatureMap.begin(); it!=m_signatureMap.end(); ++it)
				out << it->first << " : " << it->second << std::endl;		
			print_msg(DEBUG_HIGH_LEVEL, out);

			// convert the relevant field(s) into the internal structure(s)
			if (m_signatureMap.count("Hash")>0) {
				m_bSignatureDefined = true;
				m_strSignature = str_tolower(m_signatureMap["Hash"]);
				if (m_signatureMap.count("Type")>0) 
					m_strSignatureType = str_tolower(m_signatureMap["Type"]);
			}
		} // if signature exists

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "--- Finished reading signature");

		m_bSignatureCheckSucceeded=false;
		if (!m_bSignatureDefined) {
			print_msg(NORMAL_MSG, std::ostringstream().flush() << "signature block not found");
		}
		else
		{
			if (m_strSignatureType!="md5") {
				print_msg(NORMAL_MSG, std::ostringstream().flush() << "signature validation is not supported for " << m_strSignatureType);
			}
			else {
				m_bSignatureCheckSucceeded=(m_strSignature==m_strCalculatedMD5Signature);
				if (!m_bSignatureCheckSucceeded){
					print_msg(WARNING_MSG, std::ostringstream().flush() << "WARNING: signature validation failed! stored signature: " << m_strSignature << "; calculated hash: " << m_strCalculatedMD5Signature);
				}
			}
		}

		if (version_combined<1004000L) 
		{
			print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- converting blocks from version " << version_combined);
			// we need to calculate dutation of every block and save it in m_blockDurations_ru

			for (int b=0; b<m_blocks.size(); ++b) 
			{
				SeqBlock* block=GetBlock(b);
				// Calculate duration of block
				long duration = 0;
				// special processing of the delay objects (which are now eliminated)
				if (m_blockDurations_ru[b]) // non-zero means old delay library reference
				{
					if (m_tmpDelayLibrary.end()==m_tmpDelayLibrary.find(m_blockDurations_ru[b])) {
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: invalid delay library reference " << m_blockDurations_ru[b] << " in block " << b << " detected while convering the Pulseq file from older version");
						return false;
					}
					duration=m_tmpDelayLibrary[m_blockDurations_ru[b]]; // we know delay is still 0, see above
				}
				// fairly standard code, copied from the old version of GetBlock()
				if (block->isRF()) {
					RFEvent &rf = block->GetRFEvent();
					duration = MAX(duration, rf.delay+(long)m_shapeLibrary[rf.magShape].numUncompressedSamples); // in versions prior to v 1.4.0 RF raster was 1us and there was no RF time shape
				}
				for (int iC=0; iC<NUM_GRADS; iC++)
				{
					GradEvent &grad = block->GetGradEvent(iC);
					if (block->isArbitraryGradient(iC))
						duration = MAX(duration, (long)(m_dGradientRasterTime_us*m_shapeLibrary[grad.waveShape].numUncompressedSamples) + grad.delay); // in versions prior to v 1.4.0 there was no time shape
					else if (block->isTrapGradient(iC))
						duration = MAX(duration, grad.rampUpTime + grad.flatTime + grad.rampDownTime + grad.delay); 
					else if (block->isExtTrapGradient(iC)) {
						// in versions prior to 1.4.0 there were no extended trapezoids (no time shape IDs) so this should never happen
						print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: unexpected error while converting arbitrary gradients");
						return false;
					}
				}
				if (block->isADC()) {
					ADCEvent &adc = block->GetADCEvent();
					duration = MAX(duration, adc.delay + (adc.numSamples*adc.dwellTime)/1000);
				}
				if (block->isTrigger()) {
					TriggerEvent &trigger = block->GetTriggerEvent();
					duration = MAX(duration, trigger.delay+trigger.duration );
				}
				// clean up memory
				delete block;
				// convert duration to raster units and store it
				m_blockDurations_ru[b]=ceil(duration/SeqBlock::s_blockDurationRaster - 1e-12);
				// sanity check
				if (fabs(m_blockDurations_ru[b]*SeqBlock::s_blockDurationRaster-duration)>1e-9) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** WARNING: rounding up block duration for block" << b);
				}
			}
			m_tmpDelayLibrary.clear();
		}
	}

	//std::vector<double> def = GetDefinition("Scan_ID");
	//int scanID = def.empty() ? 0: (int)def[0];
	//print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );
	//print_msg(NORMAL_MSG, std::ostringstream().flush() << "===== EXTERNAL SEQUENCE #" << std::setw(5) << scanID << " ===========" );
	//print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );

	return true;
};


/***********************************************************/
void ExternalSequence::skipComments(std::istream &fileStream, char *buffer)
{
	while (getline(fileStream, buffer, MAX_LINE_SIZE)) {
		if (buffer[0]!=COMMENT_CHAR && strlen(buffer)>0) {
			break;
		}
	}
};


/***********************************************************/
char* ExternalSequence::stripWhiteSpace(char *buffer)
{
	while (*buffer && (*buffer==' ' || *buffer=='\t')) 
		buffer+=1;
	int n=strlen(buffer);
	while (n>0 && (buffer[n]==' ' || buffer[n]=='\t' || buffer[n]=='\r' || buffer[n]=='\n')) {
		buffer[n]=0;
		n-=1;
	}
	return buffer;
};


/***********************************************************/
void ExternalSequence::buildFileIndex(std::istream &fileStream)
{
	char buffer[MAX_LINE_SIZE];
    char* pBr;
    struct MD5Context mdc;
    unsigned char dg[16];
    MD5Init(&mdc);
    bool bSignatureSectionFound = false;
    std::string strippedLine;
    		
	while (getline(fileStream, buffer, MAX_LINE_SIZE, true)) {		
		//ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "buildFileIndex(): read line: [" << line << "]");
		if (buffer[0]=='[' ) {
            //print_msg(NORMAL_MSG, std::ostringstream().flush() << "opening bracket found in " << buffer);
            pBr=strchr(buffer+1,']');
			if (pBr) {
                std::string line  = std::string(buffer, pBr - buffer + 1);
                //print_msg(NORMAL_MSG, std::ostringstream().flush() << "closing bracket found, line " << line);
                m_fileIndex[line] = fileStream.tellg();
                m_fileSections.insert(fileStream.tellg());
                if (line == "[SIGNATURE]") {
                    bSignatureSectionFound = true;
                    // not on vb17 // strippedLine.clear();
					strippedLine="";
					//print_msg(NORMAL_MSG, std::ostringstream().flush() << "signature block found");
				}
			}
		}
        if (!bSignatureSectionFound) {
            if (!strippedLine.empty())
                MD5Update(&mdc, (unsigned char*)strippedLine.c_str(), strippedLine.length());
            if (0==strcmp(buffer, "\n") || 0==strcmp(buffer, "\r\n"))
                strippedLine = buffer;
			else {
                MD5Update(&mdc, (unsigned char*)buffer, strlen(buffer));
                // not on vb17 // strippedLine.clear();
				strippedLine="";
			}
		}			
	}
	m_fileSections.insert(fileStream.tellg()); // add the end-of-file (+1?)
	fileStream.clear();		// reset EOF flag
	fileStream.seekg(0, std::ios::beg);
    if (!strippedLine.empty())
        MD5Update(&mdc, (unsigned char*)strippedLine.c_str(), strippedLine.length());
    // finalize the MD5 hash calculation
	MD5Final(dg, &mdc);
    char hash[33];
    sprintf(hash,"%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",dg[0],dg[1],dg[2],dg[3],dg[4],dg[5],dg[6],dg[7],dg[8],dg[9],dg[10],dg[11],dg[12],dg[13],dg[14],dg[15]);
    print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "calculated md5 hash of the current pulseq file is: " << hash);
	m_strCalculatedMD5Signature=hash;    
};


std::vector<std::string> ExternalSequence::GetAllDefinitions()
{
    assert(m_definitions_str.size() == m_definitions.size());
    int n = m_definitions_str.size();
    std::vector<std::string> keys(n);
    std::map<std::string, std::string>::iterator it = m_definitions_str.begin();
    for (int i = 0; i < n; ++i)
    {
        keys[i] = it->first;
        ++it;
	}
    return keys;
}

/***********************************************************/
SeqBlock*	ExternalSequence::GetBlock(int index) {
	SeqBlock *block = new SeqBlock();

	// Copy event IDs
	EventIDs events = m_blocks[index];
	std::copy(events.id,events.id+NUM_EVENTS,&block->events[0]);

	//ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "GetBlock(" << index << ") : [ " << events.id[0] << " " << events.id[1] << " " << events.id[2] << " " << events.id[3] << " " << events.id[4] << " " << events.id[5] << " " << events.id[6] << "  ]");

	// Set some defaults
	block->index = index;
	block->duration_ru = m_blockDurations_ru[index];
	block->adc.numSamples = 0;
	block->trigger.triggerType=0;
	block->trigger.triggerChannel=0;
	block->rotation.defined=false;
	block->labelset.clear();
	block->labelinc.clear();
	block->softDelay.numID=-1;
	block->actualSoftDelay_ru=-1;
	block->rfShim.nchan=-1;
	block->rfShim.id=0;
	// Set event structures (if applicable) so e.g. gradient type can be determined
	if (events.id[RF]>0)     block->rf      = m_rfLibrary[events.id[RF]];
	if (events.id[ADC]>0)    block->adc     = m_adcLibrary[events.id[ADC]];
	for (unsigned int i=0; i<NUM_GRADS; i++)
		if (events.id[GX+i]>0) block->grad[i] = m_gradLibrary[events.id[GX+i]];
	// unpack (known) extension objects
	if (events.id[EXT]>0) {
		// oh yeah, the current data stuctures seem to be really ugly and slow...
		int nNextExtID=events.id[EXT];
		while (nNextExtID) {
			std::map<int,ExtensionListEntry>::iterator itEL = m_extensionLibrary.find(nNextExtID);
			if (itEL == m_extensionLibrary.end()) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "ERROR: could not find extension list entry " << nNextExtID);
				//return NULL;
				break;
			}
			// attempt to recognize the extension reference
			std::map<int,std::pair<std::string,int> >::iterator itEN=m_extensionNameIDs.find(itEL->second.type);
			if (itEN!=m_extensionNameIDs.end()) {
				// we have a known extension
				switch (itEN->second.second) {
					case EXT_TRIGGER:
						if (block->trigger.triggerType!=0) {
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: only one trigger per block is supported; error block: " << index );
						}
						else {
							// ok, lets find the trigger in the library
							block->trigger=m_triggerLibrary[itEL->second.ref]; // do we have to check whether it can be found?
						}
						break;
					case EXT_ROTATION:
						if (block->rotation.defined) {
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: only one rotation per block is supported; error block: " << index );
						}
						else {
							// ok, lets find the rotation in the library
							block->rotation=m_rotationLibrary[itEL->second.ref]; // do we have to check whether it can be found?
						}
						break;
					case EXT_LABELSET:
						//do we have to check anything ? //MZ: TODO: check that we find the evet in the library TODO: check for conflicts between set and inc
							// ok, lets find the labelset in the library
							block->labelset.push_back(m_labelsetLibrary[itEL->second.ref]); // do we have to check whether it can be found?
						break;
					case EXT_LABELINC:
						//do we have to check anything ? //MZ: TODO: check that we find the evet in the library TODO: check for conflicts between set and inc
							// ok, lets find the labelinc in the library
							block->labelinc.push_back(m_labelincLibrary[itEL->second.ref]); // do we have to check whether it can be found?
						break;
					case EXT_DELAY:
						if (block->softDelay.numID>=0) {
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: only one soft delay per block is supported; error block: " << index );
						}
						else {
							// ok, lets find the soft delay in the library
							block->softDelay=m_softDelayLibrary[itEL->second.ref]; // do we have to check whether it can be found?
						}
						break;
					case EXT_RF_SHIM:
						if (block->rfShim.nchan>0) {
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: only one soft delay per block is supported; error block: " << index );
						}
						else {
							// ok, lets find the RF shim event in the library
							block->rfShim=m_rfShimLibrary[itEL->second.ref]; // do we have to check whether it can be found?
						}
						break;
					default:
						print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unimplemented extension type " << itEN->second.first << " in block " << index );
				}
			}
			else
			{
				print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unrecognized extension type " << itEL->second.type << " in block " << index );
			}
			// update the next pointer
			nNextExtID=itEL->second.next;
		}
	}
	// Calculate duration of block
	/*long duration = 0;
	if (block->isRF()) {
		RFEvent &rf = block->GetRFEvent();
		duration = MAX(duration, rf.delay+(long)m_shapeLibrary[rf.magShape].numUncompressedSamples);
	}

	for (int iC=0; iC<NUM_GRADS; iC++)
	{
		GradEvent &grad = block->GetGradEvent(iC);
		if (block->isArbitraryGradient(iC))
			duration = MAX(duration, (long)(m_dGradientRasterTime_us*m_shapeLibrary[grad.waveShape].numUncompressedSamples) + grad.delay);
		else if (block->isTrapGradient(iC))
			duration = MAX(duration, grad.rampUpTime + grad.flatTime + grad.rampDownTime + grad.delay); 
		else if (block->isExtTrapGradient(iC)) {
			// we have to get the last time sample -- need to unpack the timeShape...
			std::vector<float> timepoints;
			// Decompress the shape for this channel
			CompressedShape& shape = m_shapeLibrary[grad.timeShape];
			timepoints.resize(shape.numUncompressedSamples);
			//if (!decompressShape(shape,&timepoints[0]))
			//	return false;
			//}
			decompressShape(shape,&timepoints[0]);
			duration = MAX(duration, timepoints[shape.numUncompressedSamples-1]*m_dGradientRasterTime_us + grad.delay);
		}
	}
	if (block->isADC()) {
		ADCEvent &adc = block->GetADCEvent();
		duration = MAX(duration, adc.delay + (adc.numSamples*adc.dwellTime)/1000);
	}
	if (block->isTrigger()) {
		TriggerEvent &trigger = block->GetTriggerEvent();
		duration = MAX(duration, trigger.delay+trigger.duration );
	}

	// TODO: conversion from previous versons
	assert(version_combined>=1004000L);
	//// handling of delays has changed in revision 1.2.0
	//if (version_combined<1002000L)
	//	block->duration = duration + block->delay;
	//else
	//	block->duration = MAX(duration, block->delay);
	*/

	//ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "block duration: " << block->duration);
    
	return block;
}

/***********************************************************/
bool ExternalSequence::decodeBlock(SeqBlock *block)
{
	int *events = &block->events[0];
	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Decoding block " << block->index << " events: "
		<< events[0]+1 << " " << events[1]+1 << " " << events[2]+1 << " "
		<< events[3]+1 << " " << events[4]+1 );
	
	std::vector<float> waveform;	
	// Decode RF
	if (block->isRF())
	{
		// Decompress the shape for this channel
		CompressedShape& shape = m_shapeLibrary[block->rf.magShape];
		waveform.resize(shape.numUncompressedSamples);
		if (!decompressShape(shape,&waveform[0]))
			return false;

		//MZ: original Kelvin's code follows
		CompressedShape& shapePhase = m_shapeLibrary[block->rf.phaseShape];
		std::vector<float> waveform_p;
		waveform_p.resize(shapePhase.numUncompressedSamples);
		if (!decompressShape(shapePhase,&waveform_p[0]))
			return false;
		// Scale phase by 2pi
		std::transform(waveform_p.begin(), waveform_p.end(), waveform_p.begin(), std::bind1st(std::multiplies<float>(),TWO_PI));

		float fDwellTime_us=0;
		// feature of v1.4.x
		// handle time shape 
		if (block->rf.timeShape) 
		{
			// new file format (v1.4.x)
			CompressedShape& shapeTime = m_shapeLibrary[block->rf.timeShape];
			// detect regular sampling 
			if (shapeTime.samples.size()!=shapeTime.numUncompressedSamples &&
				(shapeTime.samples.size()==3 || shapeTime.samples.size()==4)) 
			{
				// size=3: 1 1 N-2; size=4: t0 d d N-3 
				// in both cases the actuall dwell time is stored at shapeTime.samples[1]
				fDwellTime_us=m_dRadiofrequencyRasterTime_us*shapeTime.samples[1]; // no need to unpack 
			}
			else
			{
				std::vector<float> waveform_t;
				waveform_t.resize(shapeTime.numUncompressedSamples);
				if (!decompressShape(shapeTime,&waveform_t[0]))
					return false;
				// we resample the input on the fly 
				// for now we just use the RF raster time
				fDwellTime_us=m_dRadiofrequencyRasterTime_us;
				int nSamples=int(0.5+waveform_t.back());
				// for now we use nearest neighbour/right repetition interpolation
				// FIXME/TODO: convert to complex, use linear interpolation and convert back to magnitude&phase
				std::vector<float> wv_a(nSamples);
				std::vector<float> wv_p(nSamples);
				int tc=0;
				for(int c=0;c<nSamples;++c)
				{
					if(waveform_t[tc]<(c+1)) 
					{
						if (tc<waveform_t.size())
							++tc;
					}
					wv_a[c]=waveform[tc];
					wv_p[c]=waveform_p[tc];
				}
				// replace the waveforms
				waveform.swap(wv_a);
				waveform_p.swap(wv_p);
			}
		}
		else
		{
			if (version_combined >= 1004000)
				// use defined default RF raster time
                fDwellTime_us = m_dRadiofrequencyRasterTime_us;
            else
				// compatibility mode with the older pulseq versions
				fDwellTime_us=1.0; // old Pulseq's predefined RF raster time
		}
		//
		block->rfAmplitude = std::vector<float>(waveform);
		block->rfPhase = std::vector<float>(waveform_p);
		block->rfDwellTime_us = fDwellTime_us;
	}

	// Decode gradients
	for (int iC=GX; iC<ADC; iC++)
	{
		if (block->isArbitraryGradient(iC-GX))	// is arbitrary gradient?
		{
			// Decompress the arbitrary shape for this channel
			CompressedShape& shape = m_shapeLibrary[block->grad[iC-GX].waveShape];

			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Loaded shape with "
				<< shape.samples.size() << " compressed samples" );

			waveform.resize(shape.numUncompressedSamples);
			if (!decompressShape(shape,&waveform[0]))
				return false;

			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape uncompressed to "
				<< shape.numUncompressedSamples << " samples" );

			if (fabs(m_dGradientRasterTime_us-10)>1e-3)
			{
				print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape is on a raster that is different from the system raster, exitting... (will try resampling in the future versions...)" );
				// TODO: !!!
				// PROBLEM: we need 'first' and 'last' to be able to interpolate correctly...
				return false;
			}

			if (block->isArbGradWithOversampling(iC-GX))	// oversampling?
			{
				block->gradWaveforms[iC-GX] = std::vector<float>((waveform.size()+1)/2);
				std::vector<float>::iterator it_os=waveform.begin();
				for (std::vector<float>::iterator it=block->gradWaveforms[iC-GX].begin(); it !=block->gradWaveforms[iC-GX].end(); ++it){
					*it=*it_os;
					// std::advance(it_os,2); // this doen't work because of the odd number of elements 
					++it_os;
					if (it_os!=waveform.end())
						++it_os;
				}
			}
			else
				block->gradWaveforms[iC-GX] = std::vector<float>(waveform);
		}
	}

	if (!decodeExtTrapGradInBlock(block))
		return false;

	// Decode ADC
	/*if (block->isADC())
	{
		block->adc = m_adcLibrary[events[ADC]];
	}

	// Decode Delays
	if (block->isDelay())
	{
		block->delay = m_delayLibrary[events[DELAY]];
	}*/

	checkGradient(*block);
	checkRF(*block);

	return true;
}

bool ExternalSequence::decodeExtTrapGradInBlock(SeqBlock *block)
{
	int *events = &block->events[0];
	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Decoding ext gradient in block " << block->index << " events: "
		<< events[0]+1 << " " << events[1]+1 << " " << events[2]+1 << " "
		<< events[3]+1 << " " << events[4]+1 );

	std::vector<float> waveform;
	block->gradExtTrapForms.clear();
	block->gradExtTrapForms.resize(NUM_GRADS);
	// Decode gradients
	for (int iC=GX; iC<ADC; iC++)
	{
		if (block->isExtTrapGradient(iC-GX))	// is arbitrary gradient?
		{
			// Decompress the ExtTrap shapes for this channel
			// time shape first
			CompressedShape& tshape = m_shapeLibrary[block->grad[iC-GX].timeShape];
			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Loaded time shape " << block->grad[iC-GX].timeShape << " with " << tshape.samples.size() << " compressed samples" );
			//for (int a=0; a<tshape.samples.size(); ++a) {
			//	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << tshape.samples[a] );
			//}
			waveform.resize(tshape.numUncompressedSamples);
			if (!decompressShape(tshape,&waveform[0])) return false;
			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Time shape uncompressed to " << tshape.numUncompressedSamples << " samples" );
			//block->gradExtTrapForms[iC-GX].first = std::vector<float>(waveform); 
			block->gradExtTrapForms[iC-GX].first.resize(waveform.size());
			for (int i=0;i<waveform.size();++i)
				block->gradExtTrapForms[iC-GX].first[i]=long(0.5+m_dGradientRasterTime_us*waveform[i]); // convert to long usec from grad rasters 
			// now wave amplitude shape
			CompressedShape& wshape = m_shapeLibrary[block->grad[iC-GX].waveShape];
			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Loaded wave shape " << block->grad[iC-GX].waveShape << " with " << wshape.samples.size() << " compressed samples" );
			waveform.resize(wshape.numUncompressedSamples);
			if (!decompressShape(wshape,&waveform[0])) return false;
			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Wave shape uncompressed to " << wshape.numUncompressedSamples << " samples" );
			block->gradExtTrapForms[iC-GX].second = std::vector<float>(waveform);
			if (block->gradExtTrapForms[iC-GX].first.size() != block->gradExtTrapForms[iC-GX].second.size()) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "ERROR: uncompressed extended trapezoid time and wave shape lengths do not match" );
				return false;
			}
		}
	}
	return true;
}

/***********************************************************/
bool ExternalSequence::decompressShape(CompressedShape& encoded, float *shape)
{
	if (!encoded.isCompressed) {
		memcpy(shape,&encoded.samples.front(),sizeof(float)*encoded.numUncompressedSamples);
		return true;
	}
	// need to uncompress
	float *packed = &encoded.samples[0];
	int numPacked = encoded.samples.size();
	int numSamples = encoded.numUncompressedSamples;

	int countPack=1;
	int countUnpack=1;
	while (countPack<numPacked)
	{
		if (packed[countPack-1]!=packed[countPack])
		{
			shape[countUnpack-1]=((float)packed[countPack-1]);
			countPack++; countUnpack++;
		}
		else
		{
			int rep = ((int)packed[countPack+1])+2;
			if (fabs(packed[countPack+1]+2-rep)>1e-6) // MZ: detect format error present in some Pulseq Matlab toolbox versions
			{
				print_msg(ERROR_MSG, std::ostringstream().flush() << "ERROR: compressed shape format error detected \n"
																	 "  packed[countPack-1]=" << packed[countPack-1] << "  packed[countPack]=" << packed[countPack] << std::endl <<
																	 "  packed[countPack+1]=" << packed[countPack+1] << "  rep=" << rep << "  countPack=" << countPack );
				return false;
			}
			for (int i=countUnpack-1; i<=countUnpack+rep-2; i++)
				shape[i]= ((float)packed[countPack-1]);
			countPack += 3;
			countUnpack += rep;
		}

	}
	if (countPack==numPacked) {
		shape[countUnpack-1]=((float)packed[countPack-1]);
	}

	// Cumulative sum
	for (int i=1; i<numSamples; i++) {
		shape[i]=shape[i]+shape[i-1];
	}

	return true;
};

/***********************************************************/
bool ExternalSequence::checkBlockReferences(EventIDs& events)
{
	bool error;
	error = (events.id[RF]>0    && m_rfLibrary.count(events.id[RF])==0);
	error|= (events.id[GX]>0    && m_gradLibrary.count(events.id[GX])==0);
	error|= (events.id[GY]>0    && m_gradLibrary.count(events.id[GY])==0);
	error|= (events.id[GZ]>0    && m_gradLibrary.count(events.id[GZ])==0);
	error|= (events.id[ADC]>0   && m_adcLibrary.count(events.id[ADC])==0);
	//error|= (events.id[DELAY]>0 && m_delayLibrary.count(events.id[DELAY])==0);
	//error|= (events.id[CTRL]>0  && m_controlLibrary.count(events.id[CTRL])==0); // TODO: currently all error checking is done in getBlock(); it needs to be done here 
	
	return (!error);
}

/***********************************************************/
void ExternalSequence::checkGradient(SeqBlock& block)
{
	for (unsigned int i=0; i<NUM_GRADS; i++)
	{
		std::vector<float> &waveform = block.gradWaveforms[i];
		for (unsigned j=0; j<waveform.size(); j++)
		{
			if (waveform[j]>1.0)  waveform[j]= 1.0;
			if (waveform[j]<-1.0) waveform[j]=-1.0;
		}
		// Ensure last point is zero // MZ: no, its wrong! trapezoid gradients have a non-zero at the end!
		// if (waveform.size()>0) waveform[waveform.size()-1]=0.0;
	}
}


/***********************************************************/
void ExternalSequence::checkRF(SeqBlock& block)
{
	for (unsigned int i=0; i<block.rfAmplitude.size(); i++)
	{
		if (block.rfAmplitude[i]>1.0) block.rfAmplitude[i]=1.0;
		if (block.rfAmplitude[i]<0.0) block.rfAmplitude[i]=0.0;
		if (block.rfPhase[i]>TWO_PI-1.e-4) block.rfPhase[i]=(float)(TWO_PI-1.e-4);
		if (block.rfPhase[i]<0) block.rfPhase[i]=0.0;
	}
}


/***********************************************************/
int ExternalSequence::getline(std::istream& is, char *buffer, int MAX_SIZE, bool bRaw)
{
	//std::cout << "ExternalSequence::getline()" << std::endl;

	int i=0;	// Current position in buffer
	for(;;) {
		char c = (char)is.get();

		if (!bRaw) {
            switch (c) {
                case '\n':
                    buffer[i] = '\0';
                    // std::cout << "ExternalSequence::getline() EOL i:" << i << std::endl;
                    // std::cout << "buffer: " << buffer << std::endl;
                    return gl_true;
                case '\r':
                    if (is.peek() == '\n')
                    {
                        is.get(); // Discard character
                    }
                    buffer[i] = '\0';
                    // std::cout << "ExternalSequence::getline() CR i:" << i << std::endl;
                    // std::cout << "buffer: " << buffer << std::endl;
                    return gl_true;
            }
		}
        switch (c) {
			case EOF:
				if(i==0)
					is.seekg(-1);   // Create error on stream
				buffer[i]='\0';     // In case last line has no line ending
				//std::cout << "ExternalSequence::getline() EOF i:" << i << std::endl;
				return (i!=0)?gl_true:gl_false;
			default:
				buffer[i++] = c;
				if (i>=MAX_SIZE-1 || // -1 because we have to reserve one place for the \0 symbol
					c=='\n' )        // we are in the *raw* mode
				{
					buffer[i]='\0';
                    if (c == '\n')
                        return gl_true;
                    else
						return gl_truncated;
				}
				
		}
	}
}

#define LABELMAP_COUNTER(LBL) \
	m_labelMap.mapLabelIdToStr[LBL]=#LBL;\
	m_labelMap.mapStrToLabel[#LBL]=std::make_pair(LBL,FLAG_UNKNOWN);

#define LABELMAP_FLAG(LBL) \
	m_labelMap.mapFlagIdToStr[LBL]=#LBL;\
	m_labelMap.mapStrToLabel[#LBL]=std::make_pair(LABEL_UNKNOWN,LBL);

int ExternalSequence::decodeLabel(ExtType exttype, int& nVal, char* szLabelID, LabelEvent& label)
{
	//initialize label map if needed
	if (m_labelMap.mapStrToLabel.empty())
	{
		LABELMAP_COUNTER(SLC);
		LABELMAP_COUNTER(SEG);
		LABELMAP_COUNTER(REP);
		LABELMAP_COUNTER(AVG);
		LABELMAP_COUNTER(ECO);
		LABELMAP_COUNTER(PHS);
		LABELMAP_COUNTER(SET);
		LABELMAP_COUNTER(ACQ);
		LABELMAP_COUNTER(LIN);
		LABELMAP_COUNTER(PAR);
		LABELMAP_COUNTER(ONCE);		
		LABELMAP_FLAG(NAV);
		LABELMAP_FLAG(REV);
		LABELMAP_FLAG(SMS);
		LABELMAP_FLAG(REF);
		LABELMAP_FLAG(IMA);
		LABELMAP_FLAG(NOISE);
		LABELMAP_FLAG(PMC);
		LABELMAP_FLAG(NOPOS);
		LABELMAP_FLAG(NOROT);
		LABELMAP_FLAG(NOSCL);
		// check if all labels/flags have been added to the map
		//print_msg(WARNING_MSG, std::ostringstream().flush() << "*** m_labelMap.mapLabelIdToStr.size()= " << m_labelMap.mapLabelIdToStr.size());
		//print_msg(WARNING_MSG, std::ostringstream().flush() << "*** m_labelMap.mapFlagIdToStr.size()= " << m_labelMap.mapFlagIdToStr.size());
		//print_msg(WARNING_MSG, std::ostringstream().flush() << "*** m_labelMap.mapStrToLabel.size()= " << m_labelMap.mapStrToLabel.size());
		assert(m_labelMap.mapLabelIdToStr.size()==NUM_LABELS);
		assert(m_labelMap.mapFlagIdToStr.size()==NUM_FLAGS);
		assert(m_labelMap.mapStrToLabel.size()==NUM_LABELS+NUM_FLAGS);
	}

	// now search
	LabelMap::tM::const_iterator it = m_labelMap.mapStrToLabel.find(szLabelID);
	if (it==m_labelMap.mapStrToLabel.end())
	{
		//label.defined=false; // when no label is founded, reset label.defined to false
		print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unknown label specification\n");
		return 1;
	}

	int nKnownLBL=it->second.first;
	int nKnownFG=it->second.second;
	
	//assemble LabelEvent
	if (nKnownFG !=FLAG_UNKNOWN){
		label.flagVal=std::make_pair(nKnownFG,bool(nVal));
		label.numVal=std::make_pair(nKnownLBL,0);
	}else if (nKnownLBL !=LABEL_UNKNOWN){
		label.flagVal=std::make_pair(nKnownFG,false);
		label.numVal=std::make_pair(nKnownLBL,nVal);
	}

	//here we check if the labels/flags are valid //No boundary check, boundary check moves to prep()
	if (exttype==EXT_LABELSET){				//here we check if the values are valid
		if (nKnownLBL!=LABEL_UNKNOWN){
			/*if (nVal<0){
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Extension specification LABELSET for int-type MDH Headers is incorrect\n");
				return -1;
			}else*/
				return 0;
		}else if (nKnownFG!=FLAG_UNKNOWN){
			if ((nVal!=0)&&(nVal!=1)){
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Extension specification LABELSET for bool-type MDH Headers is incorrect\n");
				return -1;
			}else
				return 0;
		}else{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: EXT_LABELSET only support LABEL&&FLAG\n");
			return -1;					 
		}
	}else if (exttype==EXT_LABELINC){
		if (nKnownLBL!=LABEL_UNKNOWN){			
			return 0;				 
		}else if (nKnownFG!=FLAG_UNKNOWN){				// EXT_LABELINC should NOT be used for bool type MDH Headers
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Extension specification LABELINC is NOT for bool-type MDH Headers\n");
			return -1;
		}else {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: EXT_LABELINC only support LABEL&&FLAG\n");
			return -1;	
		}		
	}else{											// No ExtType recognized
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Extension specification is NOT recognized\n");
		return -1;					 
	}
}

std::string ExternalSequence::getCounterIdAsString(int nID) {
	std::map<int,std::string>::iterator it = m_labelMap.mapLabelIdToStr.find(nID);
	if (it==m_labelMap.mapLabelIdToStr.end())
		return "";
	return it->second;
}

std::string ExternalSequence::getFlagIdAsString(int nID) {
	std::map<int,std::string>::iterator it = m_labelMap.mapFlagIdToStr.find(nID);
	if (it==m_labelMap.mapFlagIdToStr.end())
		return "";
	return it->second;
}
bool ExternalSequence::isGradientInBlockStartAtNonZero(SeqBlock *block, int channel) {
	if (!block->isExtTrapGradient(channel) && !block->isArbitraryGradient(channel)) {
		//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "isGradientInBlockStartAtNonZero() returns FALSE because there is no arb grad in channel " << channel);
		return false;
	}
	// only Arbitrary or ExtTrap gradients reach here
	if (block->grad[channel].delay>0) {
		//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "isGradientInBlockStartAtNonZero() returns FALSE because there is delay in channel " << channel);
		return false;
	}
	// new since v1.5.0 : first/last
    if (version_combined >= 1005000 && block->grad[channel].first!=FLOAT_UNDEFINED )
    {
		return fabs(block->grad[channel].first)>0;
	}
	// older formats
	if (!block->gradWaveforms[channel].empty()) { 
		//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "isGradientInBlockStartAtNonZero() uses decompressed shape and returns " << (fabs(block->gradWaveforms[channel].front())>0));
		return fabs(block->gradWaveforms[channel].front())>0;
	}
	// we could decode the block's shapes at this point, but we can also just look up the first sample of the compressed shape
	//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "isGradientInBlockStartAtNonZero() uses compressed shape and returns " << (fabs(m_shapeLibrary[block->grad[channel].waveShape].samples.front())>0));
	return fabs(m_shapeLibrary[block->grad[channel].waveShape].samples.front())>0;
}

bool ExternalSequence::isAllGradientsInBlockStartAtZero(SeqBlock *block) {
	for (int i=0; i<NUM_GRADS; ++i)
		if (isGradientInBlockStartAtNonZero(block,i))
			return false;
	//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "isAllGradientsStartAtZero() returns TRUE");
	return true;
}

const std::string& trim_chars = "\t\n\v\f\r ";
std::string& str_ltrim(std::string& str)
{
    str.erase(0, str.find_first_not_of(trim_chars));
    return str;
}
 
std::string& str_rtrim(std::string& str)
{
    str.erase(str.find_last_not_of(trim_chars) + 1);
    return str;
}
 
std::string& str_trim(std::string& str)
{
    return str_ltrim(str_rtrim(str));
}

std::string str_tolower(std::string str)
{
	for (std::string::iterator it=str.begin(); it!=str.end(); ++it)
		if ((*it)<='Z' && (*it)>='A')
			*it += 'a'-'A';
    return str;
}


void LabelStateAndBookkeeping::initBookkeeping()
{
    m_bAdcLabelsInUse=false;
    m_bNonAdcLabelsInUse=false;
    initBookkeepingADC();
    m_currLabelValueStorage.flag.bValUsed.assign(NUM_FLAGS, false);
    m_currLabelValueStorage.num.bValUsed.assign(NUM_LABELS, false);
}

void LabelStateAndBookkeeping::initBookkeepingADC()
{
    // TRACE_PUT0(TC_INFO,TF_SEQ,">>> PULSEQ-DEBUG: initializing the MinMaxLabelStorage...");
    m_MinMaxLabelBookkeepingADC.flagValMin.assign(NUM_FLAGS, false);
    m_MinMaxLabelBookkeepingADC.flagValMax.assign(NUM_FLAGS, false);
    m_MinMaxLabelBookkeepingADC.bFlagMinMaxValid.assign(NUM_FLAGS, false);
    m_MinMaxLabelBookkeepingADC.numValMin.assign(NUM_LABELS, 0);
    m_MinMaxLabelBookkeepingADC.numValMax.assign(NUM_LABELS, 0);
    m_MinMaxLabelBookkeepingADC.bNumMinMaxValid.assign(NUM_LABELS, false);
    // parallel imaging flags
    m_MinMaxLabelBookkeepingADC.minRefLin    = -1;
    m_MinMaxLabelBookkeepingADC.numRefLin    = 0;
    m_MinMaxLabelBookkeepingADC.numRefImaLin = 0;
    m_MinMaxLabelBookkeepingADC.minRefPar    = -1;
    m_MinMaxLabelBookkeepingADC.numRefPar    = 0;
    m_MinMaxLabelBookkeepingADC.numRefImaPar = 0;
	// first/last
    m_mapFirstInSlc.clear();
    m_mapLastInSlc.clear();
    m_lastInMeas.clear();
    m_setFirstInSlc.clear();
    m_setLastInSlc.clear();
}

void LabelStateAndBookkeeping::initCurrState()
{
    m_currLabelValueStorage.flag.val.assign(NUM_FLAGS, false);
    m_currLabelValueStorage.flag.bValUpdated.assign(NUM_FLAGS, false);
    //m_currLabelValueStorage.flag.bValUsed.assign(NUM_FLAGS, false);
    m_currLabelValueStorage.num.val.assign(NUM_LABELS, 0);
    m_currLabelValueStorage.num.bValUpdated.assign(NUM_LABELS, false);
    //m_currLabelValueStorage.num.bValUsed.assign(NUM_LABELS, false);
}

// * ------------------------------------------------------------------ *
// *                                                                    *
// * Name        :  UpdateDataLabelStorage
// *
// * Description :  Update the current/maximum global DataLabelStorage
// *
// * Return      :  bool
// *
// * ------------------------------------------------------------------ *
void LabelStateAndBookkeeping::updateLabelValues(SeqBlock* pBlock)
{
    // update m_currLabelValueStorage according to labelinc/labelset from the block
    // for each block, we first check labelset then labelinc, which means labelinc will affect labelset, not the other
    // way around.
    static const char*      ptModule     = {"ARBITRARY_SBB::updateLabelValueStorage"};
    std::vector<LabelEvent> Tmp_labelset = pBlock->GetLabelSetEvents();
    int                     id;
    m_currLabelValueStorage.flag.bValUpdated.assign(NUM_FLAGS, false);
    m_currLabelValueStorage.num.bValUpdated.assign(NUM_LABELS, false);
    for (id = 0; id < Tmp_labelset.size(); ++id)
    {
        // if (Tmp_labelset[id].defined){
        if (Tmp_labelset[id].numVal.first != LABEL_UNKNOWN)
        {
            m_currLabelValueStorage.num.val[Tmp_labelset[id].numVal.first] = Tmp_labelset[id].numVal.second;
            if (Tmp_labelset[id].numVal.second)
                m_currLabelValueStorage.num.bValUsed[Tmp_labelset[id].numVal.first]= true; // only mark as used if non-zero
            m_currLabelValueStorage.num.bValUpdated[Tmp_labelset[id].numVal.first] = true;
            if (Tmp_labelset[id].numVal.first <= LAST_ADC_RELEVANT_LABEL)
                m_bAdcLabelsInUse = true;
            else
                m_bNonAdcLabelsInUse = true;
        }
        if (Tmp_labelset[id].flagVal.first != FLAG_UNKNOWN)
        {
            m_currLabelValueStorage.flag.val[Tmp_labelset[id].flagVal.first] = Tmp_labelset[id].flagVal.second;
            if (Tmp_labelset[id].flagVal.second)
                m_currLabelValueStorage.flag.bValUsed[Tmp_labelset[id].flagVal.first]
                    = true; // only mark as used if non-zero
            m_currLabelValueStorage.flag.bValUpdated[Tmp_labelset[id].flagVal.first] = true;
            if (Tmp_labelset[id].numVal.first <= LAST_ADC_RELEVANT_FLAG)
                m_bAdcLabelsInUse = true;
            else
                m_bNonAdcLabelsInUse = true;
        }
    }
    std::vector<LabelEvent> Tmp_labelinc = pBlock->GetLabelIncEvents();
    for (id = 0; id < Tmp_labelinc.size(); ++id)
    {
        if (Tmp_labelinc[id].numVal.first != LABEL_UNKNOWN)
        {
            m_currLabelValueStorage.num.val[Tmp_labelinc[id].numVal.first] += Tmp_labelinc[id].numVal.second;
            m_currLabelValueStorage.num.bValUsed[Tmp_labelinc[id].numVal.first]
                = true; // always mark as used because it is always non-zero
            m_currLabelValueStorage.num.bValUpdated[Tmp_labelinc[id].numVal.first] = true;
        }
        if (Tmp_labelinc[id].numVal.first <= LAST_ADC_RELEVANT_LABEL)
            m_bAdcLabelsInUse = true;
        else
            m_bNonAdcLabelsInUse = true;
    }
}

bool LabelStateAndBookkeeping::checkLabelValuesADC()
{
    if (m_currLabelValueStorage.flag.val[NOISE] || m_currLabelValueStorage.flag.val[NAV]) // noise scans and navigator scans are not included in first/last/min/max and therefore cannot be checked
		return true;

    // label boundary check
    int id;
    for (id = 0; id < m_currLabelValueStorage.num.val.size(); ++id)
    {
        if (m_MinMaxLabelBookkeepingADC.bNumMinMaxValid[id])
        {
            if ((m_currLabelValueStorage.num.val[id] > m_MinMaxLabelBookkeepingADC.numValMax[id])
                || (m_currLabelValueStorage.num.val[id] < m_MinMaxLabelBookkeepingADC.numValMin[id]))
            {
                ExternalSequence::print_msg(
                    NORMAL_MSG,
                    std::ostringstream().flush()
                        << "ERROR: Current data label (counter) " << id << " with value "
                        << m_currLabelValueStorage.num.val[id] << " is either below the min boundary "
                        << m_MinMaxLabelBookkeepingADC.numValMin[id] << " or above the max boundary "
                        << m_MinMaxLabelBookkeepingADC.numValMax[id]);
                dump_internal(m_currLabelValueStorage.num.val, m_currLabelValueStorage.flag.val, "current: ");
                dump_internal(m_MinMaxLabelBookkeepingADC.numValMin, m_MinMaxLabelBookkeepingADC.flagValMin, "min: ");
                dump_internal(m_MinMaxLabelBookkeepingADC.numValMax, m_MinMaxLabelBookkeepingADC.flagValMax, "max: ");
                return false;
            }
        }
    }
    for (id = 0; id < m_currLabelValueStorage.flag.val.size(); ++id)
    {
        if (m_MinMaxLabelBookkeepingADC.bFlagMinMaxValid[id])
        {
            if (m_currLabelValueStorage.flag.val[id] < m_MinMaxLabelBookkeepingADC.flagValMin[id]
                || m_currLabelValueStorage.flag.val[id] > m_MinMaxLabelBookkeepingADC.flagValMax[id])
            {
                ExternalSequence::print_msg(
                    NORMAL_MSG,
                    std::ostringstream().flush()
                            << "ERROR: Current data label (flag) " << id << " with value " << m_currLabelValueStorage.flag.val[id]
                            << " is either below the min boundary "
							<< m_MinMaxLabelBookkeepingADC.flagValMin[id] << " or above the max boundary " << m_MinMaxLabelBookkeepingADC.flagValMax[id]);
                dump_internal(m_currLabelValueStorage.num.val, m_currLabelValueStorage.flag.val, "current: ");
                dump_internal(m_MinMaxLabelBookkeepingADC.numValMin, m_MinMaxLabelBookkeepingADC.flagValMin, "min: ");
                dump_internal(m_MinMaxLabelBookkeepingADC.numValMax, m_MinMaxLabelBookkeepingADC.flagValMax, "max: ");
                return false;
            }
        }
    }
    return true;
}

// some older compilers (or environments) do not know vec.cbegin()...
// ve12u has GCC version 5.3, but still no cbegin() and cend(), so checking for versions doenät bring us there, therefor an uglz combo at the end...
#if defined(VXWORKS) || (defined(_MSC_VER) && _MSC_VER<=1600) || (defined(__GNUC__) && (__GNUC__<4 || (__GNUC__==4 && __GNUC_MINOR__<=6))) || !defined(COMPAT_nX)
#define cbegin begin
#define cend end 
#endif //VXWORKS //_MSC_VER //__GNUC__

std::string vec2str(const std::vector<int>& vec) 
{
	std::ostringstream os;
	for (std::vector<int>::const_iterator it = vec.cbegin(); it != vec.cend(); ++it) {
		os << *it;
		if (it + 1 != vec.cend()) os << ' ';
	}
	return os.str();
}

void LabelStateAndBookkeeping::updateBookkeepingRecordsADC()
{
    // label boundary evaluation
    // data flags for LastLine/LastSlice/LastPar, etc are analyzed further below
    if (!m_currLabelValueStorage.flag.val[NOISE] && !m_currLabelValueStorage.flag.val[NAV]) // noise scans and navigator scans are not included in first/last/min/max QC: FIRSTSCANINSLICE==false for navigator scans. 2025.01.31
    {
		int id;
		// flags
		for (id = 0; id < m_currLabelValueStorage.flag.val.size(); ++id)
		{
			/*if (m_currLabelValueStorage.bFlagUsed[id])*/ {
				if (m_MinMaxLabelBookkeepingADC.bFlagMinMaxValid[id]) // QC: if the min & max for the current flag is already set, then check and update the min and max. 2025.01.31
				{
					if (m_MinMaxLabelBookkeepingADC.flagValMin[id] > m_currLabelValueStorage.flag.val[id])
						m_MinMaxLabelBookkeepingADC.flagValMin[id] = m_currLabelValueStorage.flag.val[id];
					if (m_MinMaxLabelBookkeepingADC.flagValMax[id] < m_currLabelValueStorage.flag.val[id])
						m_MinMaxLabelBookkeepingADC.flagValMax[id] = m_currLabelValueStorage.flag.val[id];
				}
				else // QC: otherwise, set the min & max to the current value and set the flag up. 2025.01.31
				{
					ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "initializing bookkeeping for flag " << id << " with " << m_currLabelValueStorage.flag.val[id]);
					m_MinMaxLabelBookkeepingADC.flagValMin[id]       = m_currLabelValueStorage.flag.val[id];
					m_MinMaxLabelBookkeepingADC.flagValMax[id]       = m_currLabelValueStorage.flag.val[id];
					m_MinMaxLabelBookkeepingADC.bFlagMinMaxValid[id] = true;
				}
			}
		}
		// numeric values
		for (id = 0; id < m_currLabelValueStorage.num.val.size(); ++id)
		{
			/*if (m_currLabelValueStorage.bNumValUsed[id])*/ {
				if (m_MinMaxLabelBookkeepingADC.bNumMinMaxValid[id])
				{
					if (m_MinMaxLabelBookkeepingADC.numValMin[id] > m_currLabelValueStorage.num.val[id])
						m_MinMaxLabelBookkeepingADC.numValMin[id] = m_currLabelValueStorage.num.val[id];
					if (m_MinMaxLabelBookkeepingADC.numValMax[id] < m_currLabelValueStorage.num.val[id])
						m_MinMaxLabelBookkeepingADC.numValMax[id] = m_currLabelValueStorage.num.val[id];
				}
				else
				{
					ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "initializing bookkeeping for counter " << id << " with " << m_currLabelValueStorage.num.val[id]);
					m_MinMaxLabelBookkeepingADC.numValMin[id]       = m_currLabelValueStorage.num.val[id];
					m_MinMaxLabelBookkeepingADC.numValMax[id]       = m_currLabelValueStorage.num.val[id];
					m_MinMaxLabelBookkeepingADC.bNumMinMaxValid[id] = true;
				}
			}
		}
		// counter & flags tracking for parallel imaging
		// QC: always correct by using bookkeeping? e.g. in unregular or random samping case. 2024.12.20
		if (m_currLabelValueStorage.flag.bValUsed[REF])
		{
			if (m_currLabelValueStorage.flag.val[REF])
			{
				if (m_MinMaxLabelBookkeepingADC.minRefLin < 0
					|| m_MinMaxLabelBookkeepingADC.minRefLin > m_currLabelValueStorage.num.val[LIN])
				{
					m_MinMaxLabelBookkeepingADC.minRefLin = m_currLabelValueStorage.num.val[LIN];
					// reset counters because we only want to count once (e.g. for minPar and minLin)
					m_MinMaxLabelBookkeepingADC.numRefLin    = 0;
					m_MinMaxLabelBookkeepingADC.numRefImaLin = 0;
					m_MinMaxLabelBookkeepingADC.numRefPar    = 0;
					m_MinMaxLabelBookkeepingADC.numRefImaPar = 0;
				}
				if (m_MinMaxLabelBookkeepingADC.minRefPar < 0
					|| m_MinMaxLabelBookkeepingADC.minRefPar > m_currLabelValueStorage.num.val[PAR])
				{
					m_MinMaxLabelBookkeepingADC.minRefPar = m_currLabelValueStorage.num.val[PAR];
					// reset counters because we only want to count once (e.g. for minPar and minLin)
					m_MinMaxLabelBookkeepingADC.numRefLin    = 0;
					m_MinMaxLabelBookkeepingADC.numRefImaLin = 0;
					m_MinMaxLabelBookkeepingADC.numRefPar    = 0;
					m_MinMaxLabelBookkeepingADC.numRefImaPar = 0;
				}
				// this is a tricky part: which labels do we require to be 0 for the bookkeeping?
				// SEG and AVG can be non-zero
				if (m_currLabelValueStorage.num.val[SLC] == 0 && m_currLabelValueStorage.num.val[REP] == 0
					&& m_currLabelValueStorage.num.val[AVG] == 0 && m_currLabelValueStorage.num.val[ECO] == 0
					&& m_currLabelValueStorage.num.val[PHS] == 0 && m_currLabelValueStorage.num.val[SET] == 0)
				{
					if (m_MinMaxLabelBookkeepingADC.minRefLin == m_currLabelValueStorage.num.val[LIN])
					{
						++m_MinMaxLabelBookkeepingADC.numRefPar;
						if (m_currLabelValueStorage.flag.val[IMA])
							++m_MinMaxLabelBookkeepingADC.numRefImaPar;
					}
					if (m_MinMaxLabelBookkeepingADC.minRefPar == m_currLabelValueStorage.num.val[PAR])
					{
						++m_MinMaxLabelBookkeepingADC.numRefLin;
						if (m_currLabelValueStorage.flag.val[IMA])
							++m_MinMaxLabelBookkeepingADC.numRefImaLin;
					}
				}
			}
		}
		// track first/last scan in slice, etc
        int nCurSlc = m_currLabelValueStorage.num.val[SLC]; 
		// last tracking is easy: we just use the current as the last candidate and the finalize call will make it right
		m_lastInMeas = m_currLabelValueStorage.getAdcCounters(false,false);
		//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "updating m_lastInMeas to " << vec2str(m_lastInMeas));
        m_mapLastInSlc[nCurSlc] = m_currLabelValueStorage.getAdcCounters(true, false); // fixing the lastInSlice behaviour for REP!=0
		// first tracking: only update if this slice has no record
        if (m_mapFirstInSlc.find(nCurSlc) == m_mapFirstInSlc.end()) //QC: if firstScanInSlice is empty, then set it to the current ADC counter. 2025.01.31
        {
            m_mapFirstInSlc[nCurSlc] = m_mapLastInSlc[nCurSlc]; // m_mapLastInSlc[nCurSlc] actually contains the current counters, see above
        }
    }
}

std::vector<int> LabelValueStorage::getAdcCounters(bool bIgnoreREP, bool bAlsoIgnoreAVG)
{
	//ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "LabelValueStorage::getAdcCounters() called for " << vec2str(num.val));
    if (bAlsoIgnoreAVG && !bIgnoreREP)
    {
        ExternalSequence::print_msg(WARNING_MSG, std::ostringstream().flush() << "WARNING: LabelValueStorage::getAdcCounters() was called with bAlsoIgnoreAVG, fixing bIgnoreREP");
        bIgnoreREP = bAlsoIgnoreAVG;
	}
    std::vector<int> r(num.val.begin(), num.val.begin() + LAST_ADC_RELEVANT_LABEL + 1 - int(bIgnoreREP) - int(bAlsoIgnoreAVG));
    //ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "returning " << vec2str(r));
	return r;
}


void LabelStateAndBookkeeping::finalizeBookkeepingRecordsADC()
{
	ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "Entering LabelStateAndBookkeeping::finalizeBookkeepingRecordsADC()...");
	std::map<int, std::vector<int> >::iterator it;
    for (it = m_mapFirstInSlc.begin(); it != m_mapFirstInSlc.end(); ++it)
    {
		ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "Adding " << vec2str(it->second) << " as FirstInSlc");
        m_setFirstInSlc.insert(it->second);
    }
    for (it = m_mapLastInSlc.begin(); it != m_mapLastInSlc.end(); ++it)
    {
		ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "Adding " << vec2str(it->second) << " as LastInSlc");
		m_setLastInSlc.insert(it->second);
    }
}

bool LabelStateAndBookkeeping::isFirstScanInSlice() 
{
	if (m_setFirstInSlc.empty() && !m_mapFirstInSlc.empty())
		finalizeBookkeepingRecordsADC();
	return m_setFirstInSlc.end() != m_setFirstInSlc.find(m_currLabelValueStorage.getAdcCounters(true,false));
}

bool LabelStateAndBookkeeping::isLastScanInSlice() 
{
	if (m_setLastInSlc.empty() && !m_mapLastInSlc.empty())
		finalizeBookkeepingRecordsADC();
	return m_setLastInSlc.end() != m_setLastInSlc.find(m_currLabelValueStorage.getAdcCounters(true,false));
}


void LabelStateAndBookkeeping::dump(const char* szMsg /*=NULL*/, bool bMinMax /*=true*/, bool bCurr /*=true*/)
{
    if (bMinMax)
    {
        ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "Min Values:");
        dump_internal(m_MinMaxLabelBookkeepingADC.numValMin, m_MinMaxLabelBookkeepingADC.flagValMin, szMsg);
        ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "Max Values:");
        dump_internal(m_MinMaxLabelBookkeepingADC.numValMax, m_MinMaxLabelBookkeepingADC.flagValMax, szMsg);
    }
    if (bCurr)
    {
        ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << "Current Values:");
        dump_internal(m_currLabelValueStorage.num.val, m_currLabelValueStorage.flag.val, szMsg);
    }
}

void LabelStateAndBookkeeping::dump_internal(const std::vector<int>& numVal, const std::vector<bool>& flagVal, const char* szMsg)
{
    const char* szSpe = ""; 
	char* szSpe0=NULL;
    if (szMsg == NULL)
        szMsg ="";
    else
    {
        int nLen = strlen(szMsg);
        szSpe0  = new (char[nLen + 1]);
        memset(szSpe0, ' ', nLen);
        szSpe0[nLen] = '\0';
        szSpe = szSpe0;
    }
    ExternalSequence::print_msg(
        NORMAL_MSG,
        std::ostringstream().flush() << szMsg << " SLC " << numVal[SLC] << " SEG " << numVal[SEG] << " REP "
                                     << numVal[REP]);
    ExternalSequence::print_msg(
        NORMAL_MSG,
        std::ostringstream().flush() << szSpe << " AVG " << numVal[AVG] << " ECO " << numVal[ECO] << " PHS "
                                     << numVal[PHS]);
    ExternalSequence::print_msg(
        NORMAL_MSG,
        std::ostringstream().flush() << szSpe << " SET " << numVal[SET] << " LIN " << numVal[LIN] << " PAR "
                                     << numVal[PAR]);
    ExternalSequence::print_msg(
        NORMAL_MSG, std::ostringstream().flush() << szSpe << " NAV " << flagVal[NAV] << " REV " << flagVal[REV]);
    ExternalSequence::print_msg(
        NORMAL_MSG, std::ostringstream().flush() << szSpe << " SMS " << flagVal[SMS] << " PMC " << flagVal[PMC]);
    ExternalSequence::print_msg(
        NORMAL_MSG,
        std::ostringstream().flush() << szSpe << " REF " << flagVal[REF] << " IMA " << flagVal[IMA] << " NOISE "
                                     << flagVal[NOISE]);
    ExternalSequence::print_msg(
        NORMAL_MSG,
        std::ostringstream().flush() << szSpe << " NOPOS " << flagVal[NOPOS] << " NOROT " << flagVal[NOROT] << " NOSCL "
                                     << flagVal[NOSCL]);
    ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << szSpe << " ONCE " << numVal[ONCE]);
    /*int c=0, r=0;
    int i;
    for (l=0;i<NUM_LABELS;++i)
    {
        if (c==0) {
            ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << (r ? szSpe : szMsg) );
        }
        ExternalSequence::print_msg(NORMAL_MSG, std::ostringstream().flush() << " " << ... <-- we have a problem here
    because this is a static function if (++c>3) {
            ++r;
            c=0;
            <-- we have another problem here because print_msg always adds a newline...
        }
    }
    if (c)
        <-- print new line */
    delete szSpe0;
}

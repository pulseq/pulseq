#include "ExternalSequence.h"

#include <stdio.h>		// sscanf
#include <cstring>		// strlen etc
#include <iomanip>		// std::setw etc

#include <algorithm>	// for std::sort
#include <functional>	// for std::sort


#ifndef MAX
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
#endif

ExternalSequence::PrintFunPtr ExternalSequence::print_fun = &ExternalSequence::defaultPrint;
const int ExternalSequence::MAX_LINE_SIZE = 256;
const char ExternalSequence::COMMENT_CHAR = '#';
	
/***********************************************************/
ExternalSequence::ExternalSequence()
{
	m_scanID			= 1;
	m_FOV[0]			= 256.0;
	m_FOV[1]			= 256.0;
	m_BaseResolution[0]	= 64;
	m_BaseResolution[1]	= 64;
	m_TE				= 15.0;
	m_TR				= 50.0;
	m_sliceThickness	= 3.0;
#ifdef MASTER_SLAVE_FORMAT
	m_isSlave			= false;
#endif

	m_RotMatrix[0][0] = 1.0; m_RotMatrix[1][0] = 0.0; m_RotMatrix[2][0] = 0.0;
	m_RotMatrix[0][1] = 0.0; m_RotMatrix[1][1] = 1.0; m_RotMatrix[2][1] = 0.0;
	m_RotMatrix[0][2] = 0.0; m_RotMatrix[1][2] = 0.0; m_RotMatrix[2][2] = 1.0;
	
	m_delays[0] = 0; m_delays[1] = 0; m_delays[2] = 0;

	m_blocks.clear();
	
}


/***********************************************************/
ExternalSequence::~ExternalSequence(){}

/***********************************************************/
void ExternalSequence::print_msg(MessageType level, std::ostream& ss) {

	if (MSG_LEVEL>=level) {
		std::ostringstream oss;
		oss.width(2*(level-1)); oss << "";
		oss << static_cast<std::ostringstream&>(ss).str();
		print_fun(oss.str().c_str());
	}
}

/***********************************************************/
bool ExternalSequence::load(std::string path)
{
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "Reading external sequence files");
#ifdef MASTER_SLAVE_FORMAT
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "Expecting 6-channel file format");
#endif

	char buffer[MAX_LINE_SIZE];
	char tmpStr[MAX_LINE_SIZE];
	
	// **********************************************************************************************************************
	// ************************ READ SHAPES ***********************************

	// Try single file mode (everything in .seq file)
	std::ifstream data_file;
	bool isSingleFileMode = true;
	std::string filepath = path;
	if (filepath.substr(filepath.size()-4) != std::string(".seq")) {
		filepath = path + "/external.seq";
	}
	data_file.open(filepath.c_str(), std::ifstream::in);
	data_file.seekg(0, std::ios::beg);

	if (!data_file.good())
	{
		// Try separate file mode (blocks.seq, events.seq, shapes.seq)
		filepath = path + "/shapes.seq";
		data_file.open(filepath.c_str(), std::ifstream::in);
		data_file.seekg(0, std::ios::beg);
		isSingleFileMode = false;
	}
	if (!data_file.good())
	{
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
		return false;
	}
	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Buildingindex" );
	
	// Save locations of section tags
	buildFileIndex(data_file);
	
	// Read shapes section
	// ------------------------
	if (m_fileIndex.find("[SHAPES]") == m_fileIndex.end()) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: [SHAPES] section");
		return false;
	}
	data_file.seekg(m_fileIndex["[SHAPES]"], std::ios::beg);
	skipComments(data_file,buffer);			// Ignore comments & empty lines

	int shapeId, numSamples;
	float sample;
	
	m_shapeLibrary.clear();
	while (data_file.good() && buffer[0]=='s')
	{
		sscanf(buffer, "%s%d", tmpStr, &shapeId);
		data_file.getline(buffer, MAX_LINE_SIZE);
		sscanf(buffer, "%s%d", tmpStr, &numSamples);

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Reading shape " << shapeId );

		CompressedShape shape;
		shape.samples.clear();
		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='s' || strlen(buffer)==0) {
				break;
			}
			sscanf(buffer, "%f", &sample);
			shape.samples.push_back(sample);
		}
		shape.numUncompressedSamples=numSamples;

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape index " << shapeId << " has " << shape.samples.size()
			<< " compressed and " << shape.numUncompressedSamples << " uncompressed samples" );

		m_shapeLibrary[shapeId] = shape;
		
		skipComments(data_file,buffer);			// Ignore comments & empty lines
	}
	data_file.clear();	// In case EOF reached

	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- SHAPES READ numShapes: " << m_shapeLibrary.size() );
	

	// **********************************************************************************************************************
	// ************************ READ EVENTS ********************
	if (!isSingleFileMode) {
		filepath = path + "/events.seq";
		data_file.close();
		data_file.open(filepath.c_str(), std::ifstream::in);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}
		buildFileIndex(data_file);
	}

	// Read RF section
	// ------------------------
	if (m_fileIndex.find("[RF]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[RF]"], std::ios::beg);

		int rfId;
		m_rfLibrary.clear();
		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			RFEvent event;
			sscanf(buffer, "%d%f%d%d%f%f", &rfId, &(event.amplitude),
				&(event.magShape),&(event.phaseShape),
				&(event.freqOffset), &(event.phaseOffset)
				);

			m_rfLibrary[rfId] = event;
		}
	}

	// Read *arbitrary* gradient section
	// -------------------------------
	if (m_fileIndex.find("[GRADIENTS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[GRADIENTS]"], std::ios::beg);

		m_gradLibrary.clear();
		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			int gradId;
			GradEvent event;
			sscanf(buffer, "%d%f%d", &gradId, &(event.amplitude), &(event.shape) );

			m_gradLibrary[gradId] = event;
		}
	}
	
	// Read *trapezoid* gradient section
	// -------------------------------
	if (m_fileIndex.find("[TRAP]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[TRAP]"], std::ios::beg);

		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			int gradId;
			GradEvent event;
			sscanf(buffer, "%d%f%ld%ld%ld", &gradId, &(event.amplitude),
				&(event.rampUpTime),&(event.flatTime),&(event.rampDownTime)
				);
			event.shape=0;
			m_gradLibrary[gradId] = event;
		}
	}

	// Sort gradients based on index
	// -----------------------------
	//std::sort(m_gradLibrary.begin(),m_gradLibrary.end(),compareGradEvents);

	
	// Read ADC section
	// -------------------------------
	if (m_fileIndex.find("[ADC]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[ADC]"], std::ios::beg);

		int adcId;
		m_adcLibrary.clear();
		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			ADCEvent event;
			sscanf(buffer, "%d%d%d%d%f%f", &adcId, &(event.numSamples),
				&(event.dwellTime),&(event.delay),&(event.freqOffset),&(event.phaseOffset)
				);
			m_adcLibrary[adcId] = event;
		}
	}
	
	// Read delays section
	// -------------------------------
	if (m_fileIndex.find("[DELAYS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[DELAYS]"], std::ios::beg);

		int delayId;
		long delay;
		m_delayLibrary.clear();
		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			sscanf(buffer, "%d%ld", &delayId, &delay);
			m_delayLibrary[delayId] = delay;
		}
	}

	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- EVENTS READ: "
		<<" RF: " << m_rfLibrary.size() 
		<<" GRAD: " << m_gradLibrary.size()
		<<" ADC: " << m_adcLibrary.size()
		<<" DELAY: " << m_delayLibrary.size());


	// **********************************************************************************************************************
	// ************************ READ BLOCKS ********************
	if (!isSingleFileMode) {
		filepath = path + "/blocks.seq";
		data_file.close();
		data_file.open(filepath.c_str(), std::ifstream::in);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}
		buildFileIndex(data_file);
	}

	// File storage order
	// -----------------
	enum Fields {
		DELAY, RF,
		GX, GY, GZ,
#ifdef MASTER_SLAVE_FORMAT
		GA, GB, GC,
#endif
		ADC
	};
	const int NUM_FIELDS = ADC+1;

	// Read definition section
	// ------------------------
	unsigned int numBlocks = 0;
	if (m_fileIndex.find("[DEFINITIONS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[DEFINITIONS]"], std::ios::beg);
		
	#ifdef MASTER_SLAVE_FORMAT
		std::vector<long> delays(6,0);
	#else
		std::vector<long> delays(3,0);
	#endif

		// Define the definition structure: a mapping of string keys to
		// a vector of pointers (destination variables) and the corresponding types for scanf.
		// The vector allows multiple values per line
		std::map<std::string,std::pair<std::string,std::vector<void*> > > definitions;
		definitions["Num_Blocks"]        = std::make_pair(std::string("%d"),  std::vector<void*>(1,&numBlocks));
		definitions["Scan_ID"]           = std::make_pair(std::string("%d"),  std::vector<void*>(1,&m_scanID));
		definitions["SliceThickness_mm"] = std::make_pair(std::string("%f"),  std::vector<void*>(1,&m_sliceThickness));
		definitions["TE_ms"]             = std::make_pair(std::string("%f"),  std::vector<void*>(1,&m_TE));
		definitions["TR_ms"]             = std::make_pair(std::string("%f"),  std::vector<void*>(1,&m_TR));
		definitions["FOV_mm"]            = std::make_pair(std::string("%f %f"),std::vector<void*>(1,&m_FOV[0]));
		definitions["FOV_mm"].second.push_back(&m_FOV[1]);
		void *ptr[9] = {&m_RotMatrix[0][0],&m_RotMatrix[0][1],&m_RotMatrix[0][2],
						&m_RotMatrix[1][0],&m_RotMatrix[1][1],&m_RotMatrix[1][2],
						&m_RotMatrix[2][0],&m_RotMatrix[2][1],&m_RotMatrix[2][2] };
		definitions["Rot_Matrix"]        = std::make_pair(std::string("%lf %lf %lf %lf %lf %lf %lf %lf %lf"),std::vector<void*>());
		unsigned int i;
		for (i=0; i<9; i++) {
			definitions["Rot_Matrix"].second.push_back(ptr[i]);
		}
		definitions["Delays_us"]         = std::make_pair(std::string(""),std::vector<void*>());
		for (i=0; i<delays.size(); i++) {
			definitions["Delays_us"].first += std::string(" %ld");
			definitions["Delays_us"].second.push_back(&delays[i]);
		}

		// Read each line and search for the key.
		// If found read into corresponding pointer
		while (data_file.getline(buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			std::istringstream ss_buf(buffer);
			std::string key, val, type;
			ss_buf >> key;
			if (definitions.find(key) == definitions.end()) {
				print_msg(WARNING_MSG, std::ostringstream().flush() << "WARNING: Definition " << key << " ignored");
			} else {
				// Found key, so parse value according to type string
				std::istringstream ss_type((char*)definitions[key].first.c_str());
				std::vector<void*> ptrs = definitions[key].second;
				for (unsigned int i=0; i<ptrs.size(); i++) {
					ss_buf >> val;		// get the value as a string
					ss_type >> type;	// get the scanf type
					sscanf(val.c_str(), type.c_str(), ptrs[i]);
				}
			}
		}
	
	
		print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- DEFINITIONS READ: "
			<<"ID: "<< m_scanID << " "
			<<"Num blocks: "<< numBlocks << " "
			<<"FOV: "<< m_FOV[0] << " " << m_FOV[1] << " "
			<<"ST: " << m_sliceThickness << " "
			<<"TE: "<< m_TE << " TR: " << m_TR << " "
			<<"Rot:[" << m_RotMatrix[0][0] << " " << m_RotMatrix[0][1] << " " << m_RotMatrix[0][2] << "; "
			<< m_RotMatrix[1][0] << " " << m_RotMatrix[1][1] << " " << m_RotMatrix[1][2] << "; "
			<< m_RotMatrix[2][0] << " " << m_RotMatrix[2][1] << " " << m_RotMatrix[2][2] << "] " 
			<< "Delays: " << delays[0] << " " << delays[1] << " " << delays[2] << " "
#ifdef MASTER_SLAVE_FORMAT
			<< delays[3] << " " << delays[4] << " " << delays[5] 
#endif
			);

		m_maxDelay = *((std::max_element)(delays.begin(), delays.end()));
		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Max delay: " << m_maxDelay);

		m_delays[0]=delays[0];
		m_delays[1]=delays[1];
		m_delays[2]=delays[2];
		
#ifdef MASTER_SLAVE_FORMAT
		if (m_isSlave) {
			m_delays[0]=delays[3];
			m_delays[1]=delays[4];
			m_delays[2]=delays[5];
		}
#endif
	} // if definitions exist
	
	// Read blocks section
	// ------------------------
	if (m_fileIndex.find("[BLOCKS]") == m_fileIndex.end()) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: [BLOCKS] section");
		return false;
	}
	data_file.seekg(m_fileIndex["[BLOCKS]"], std::ios::beg);

	// System timing parameters
	long COIL_LEAD_TIME=100;
	long lCoilLeadTime=0;
	long lFreqResetTime=10;
	
	int blockIdx;
	int eventId[NUM_FIELDS];
	
	// Read blocks
	m_blocks.clear();
	while (data_file.getline(buffer, MAX_LINE_SIZE)) {
		if (buffer[0]=='[' || strlen(buffer)==0) {
			break;
		}

#ifdef MASTER_SLAVE_FORMAT
		sscanf(buffer, "%d%d%d%d%d%d%d%d%d%d", &blockIdx, 
			&eventId[DELAY],							// Delay
			&eventId[RF],								// RF
			&eventId[GX],&eventId[GY],&eventId[GY],		// Linear
			&eventId[GA],&eventId[GB],&eventId[GC],		// Patloc
			&eventId[ADC]								// ADCs
			);
#else
		sscanf(buffer, "%d%d%d%d%d%d%d", &blockIdx, 
			&eventId[DELAY],							// Delay
			&eventId[RF],								// RF
			&eventId[GX],&eventId[GY],&eventId[GZ],		// Gradients
			&eventId[ADC]								// ADCs
			);
#endif

		SeqBlock block;
		block.index = blockIdx-1;

		// Calculate duration of block
		long duration = 0;
		lCoilLeadTime=0;
		if (eventId[RF]>0) {
			RFEvent &rf = m_rfLibrary[eventId[RF]];
			lCoilLeadTime=COIL_LEAD_TIME;
			//duration = MAX(duration, lFreqResetTime+lCoilLeadTime+(long)m_shapeLibrary[rf.magShape].numUncompressedSamples);
			duration = MAX(duration, lCoilLeadTime+(long)m_shapeLibrary[rf.magShape].numUncompressedSamples);
		}
		for (int iC=GX; iC<ADC; iC++)
		{
			if (eventId[iC]>0) {	// has gradient
				GradEvent &grad = m_gradLibrary[eventId[iC]];
				if (grad.shape>0)
					duration = MAX(duration, lCoilLeadTime + (long)(10*m_shapeLibrary[grad.shape].numUncompressedSamples));
				else
					duration = MAX(duration, lCoilLeadTime + grad.rampUpTime + grad.flatTime + grad.rampDownTime);
			}
		}
		if (eventId[ADC]>0) {
			ADCEvent &adc = m_adcLibrary[eventId[ADC]];
			duration = MAX(duration, lFreqResetTime + adc.delay + adc.numSamples*(adc.dwellTime/1000));
		}

		block.events[RF_EVENT]    = eventId[RF];		// RF
		block.events[ADC_EVENT]   = eventId[ADC];		// ADC
		block.events[XGRAD_EVENT] = eventId[GX];		// Linear gradients: X, Y, Z
		block.events[YGRAD_EVENT] = eventId[GY];
		block.events[ZGRAD_EVENT] = eventId[GZ];
		
#ifdef MASTER_SLAVE_FORMAT
		// Assign events based on patloc mode
		if (m_isSlave) {
			block.events[RF_EVENT]    = 0;				// No RF on PatLoc slave
			block.events[ADC_EVENT]   = 0;				// No ADC on PatLoc slave
			block.events[XGRAD_EVENT] = eventId[GA];	// PatLoc gradients: A, B, C
			block.events[YGRAD_EVENT] = eventId[GB];
			block.events[ZGRAD_EVENT] = eventId[GC];
		}
#endif
		// Set some defaults
		block.delay = 0;
		block.adc.numSamples = 0;
		
		// Set event structures (if applicable) so e.g. gradient type can be determined
		if (block.events[XGRAD_EVENT]>0)  block.grad[0] = m_gradLibrary[block.events[XGRAD_EVENT]];
		if (block.events[YGRAD_EVENT]>0)  block.grad[1] = m_gradLibrary[block.events[YGRAD_EVENT]];
		if (block.events[ZGRAD_EVENT]>0)  block.grad[2] = m_gradLibrary[block.events[ZGRAD_EVENT]];
		if (block.events[RF_EVENT]>0)     block.rf      = m_rfLibrary[block.events[RF_EVENT]];
		if (block.events[ADC_EVENT]>0)    block.adc     = m_adcLibrary[block.events[ADC_EVENT]];
		if (eventId[DELAY]>0)             block.delay   = m_delayLibrary[eventId[DELAY]];
		
#ifdef MASTER_SLAVE_FORMAT
		// Adjust timing for coil lead time when nonlinear gradients are present during RF pulse
		if (m_isSlave && eventId[RF]>0 && (eventId[GX]>0 || eventId[GY]>0 || eventId[GZ]>0))
			block.delay += lCoilLeadTime;
#endif

		block.duration = duration + block.delay;
	
		//print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Block duration: " << block->duration );
		
		if (!checkBlockReferences(block)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Block " << blockIdx
				<< " contains references to undefined events" );
			return false;
		}
		// Add block to list
		m_blocks.push_back(block);
	}
	data_file.close();

	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- BLOCKS READ: " << m_blocks.size());

	// Num_Blocks definition (if defined) is used to check the correct number of blocks are read
	if (numBlocks>0 && m_blocks.size()!=numBlocks) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Expected " << numBlocks
			<< " blocks but read " << m_blocks.size() << " blocks");
		return false;
	}

	print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "===== EXTERNAL SEQUENCE #" << std::setw(5) << m_scanID << " ===========" );
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );
	
	return true;
};


/***********************************************************/
void ExternalSequence::skipComments(std::ifstream &fileStream, char *buffer)
{
	while (fileStream.getline(buffer, MAX_LINE_SIZE)) {
		if (buffer[0]!=COMMENT_CHAR && strlen(buffer)>0) {
			break;
		}
	}
};


/***********************************************************/
void ExternalSequence::buildFileIndex(std::ifstream &fileStream)
{
	char buffer[MAX_LINE_SIZE];
	int pos=0;
	while (fileStream.getline(buffer, MAX_LINE_SIZE)) {
		// tellp() has bugs in ancient compilers, so we accumulate gcount()
		pos += fileStream.gcount();
		std::string line = std::string(buffer);
		
		if (line[0]=='[' && line[line.length()-1]==']') {
			m_fileIndex[line] = pos;
		}
	}
	fileStream.clear();		// reset EOF flag
	fileStream.seekg(0, std::ios::beg);
};

/***********************************************************/
bool ExternalSequence::decodeBlock(SeqBlock *block)
{

	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Decoding block " << block->index << " events: "
		<< block->events[0]+1 << " " << block->events[1]+1 << " " << block->events[2]+1 << " "
		<< block->events[3]+1 << " " << block->events[4]+1 );

	std::vector<float> waveform;

	// Decode RF
	if (block->events[RF_EVENT]>0)
	{		
		// Decompress the shape for this channel
		CompressedShape& shape = m_shapeLibrary[block->rf.magShape];
		waveform.resize(shape.numUncompressedSamples);
		decompressShape(shape,&waveform[0]);
		block->rfAmplitude = std::vector<float>(waveform);

		CompressedShape& shapePhase = m_shapeLibrary[block->rf.phaseShape];
		waveform.resize(shapePhase.numUncompressedSamples);
		decompressShape(shapePhase,&waveform[0]);
		
		// Scale phase by 2pi
		std::transform(waveform.begin(), waveform.end(), waveform.begin(), std::bind1st(std::multiplies<float>(),TWO_PI));
		block->rfPhase = std::vector<float>(waveform);
	}

	// Decode gradients
	for (int iC=XGRAD_EVENT; iC<=ZGRAD_EVENT; iC++)
	{
		if (block->events[iC]>0)	// has a gradient?
		{
			if (block->grad[iC-1].shape>0)	// is arbitrary gradient?
			{
				// Decompress the arbitrary shape for this channel
				CompressedShape& shape = m_shapeLibrary[block->grad[iC-1].shape];

				print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Loaded shape with "
					<< shape.samples.size() << " compressed samples" );

				waveform.resize(shape.numUncompressedSamples);
				decompressShape(shape,&waveform[0]);
		
				print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape uncompressed to "
					<< shape.numUncompressedSamples << " samples" );

				if (iC==XGRAD_EVENT)
					block->grad_x = std::vector<float>(waveform);
				else if (iC==YGRAD_EVENT)
					block->grad_y = std::vector<float>(waveform);
				else if (iC==ZGRAD_EVENT)
					block->grad_z = std::vector<float>(waveform);
			} 
		}
	}

	checkGradient(*block);
	checkRF(*block);

	return true;
}

/***********************************************************/
bool ExternalSequence::decompressShape(CompressedShape& encoded, float *shape)
{
	float *packed = &encoded.samples[0];
	int numPacked = encoded.samples.size();
	int numSamples = encoded.numUncompressedSamples;

	/*for (int j=0; j<numPacked; j++)
		std::cout << "packed: " << packed[j] << std::endl;*/

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
			for (int i=countUnpack-1; i<=countUnpack+rep-2; i++)
				shape[i]= ((float)packed[countPack-1]);
			countPack += 3;
			countUnpack += rep;
		}

	}
	if (countPack==numPacked) {
		shape[countUnpack-1]=((float)packed[countPack-1]);
	}

	/*for ( j=numSamples-10; j<numSamples; j++) {
		std::cout << "unpacked " << shape[j] << std::endl;	// Prev problem with the last decoded sampled
	}*/

	// Cumulative sum
	for (int i=1; i<numSamples; i++) {
		shape[i]=shape[i]+shape[i-1];
	}

	return true;
};

/***********************************************************/
bool ExternalSequence::checkBlockReferences(SeqBlock& block)
{
	bool error;
	error = (block.events[RF_EVENT]>0    && m_rfLibrary.count(block.events[RF_EVENT])==0);
	error|=	(block.events[XGRAD_EVENT]>0 && m_gradLibrary.count(block.events[XGRAD_EVENT])==0);
	error|=	(block.events[YGRAD_EVENT]>0 && m_gradLibrary.count(block.events[YGRAD_EVENT])==0);
	error|=	(block.events[ZGRAD_EVENT]>0 && m_gradLibrary.count(block.events[ZGRAD_EVENT])==0);
	error|=	(block.events[ADC_EVENT]>0   && m_adcLibrary.count(block.events[ADC_EVENT])==0);
	
	return (!error);
}

/***********************************************************/
void ExternalSequence::checkGradient(SeqBlock& block)
{
	unsigned int i;
	for (i=0; i<block.grad_x.size(); i++)
	{
		if (block.grad_x[i]>1.0)  block.grad_x[i]= 1.0;
		if (block.grad_x[i]<-1.0) block.grad_x[i]=-1.0;
	}
	for (i=0; i<block.grad_y.size(); i++)
	{
		if (block.grad_y[i]>1.0)  block.grad_y[i]= 1.0;
		if (block.grad_y[i]<-1.0) block.grad_y[i]=-1.0;
	}
	for (i=0; i<block.grad_z.size(); i++)
	{
		if (block.grad_z[i]>1.0)  block.grad_z[i]= 1.0;
		if (block.grad_z[i]<-1.0) block.grad_z[i]=-1.0;
	}

	// Ensure last point is zero
	if (block.grad_x.size()>0) block.grad_x[block.grad_x.size()-1]=0.0;
	if (block.grad_y.size()>0) block.grad_y[block.grad_y.size()-1]=0.0;
	if (block.grad_z.size()>0) block.grad_z[block.grad_z.size()-1]=0.0;
}

/***********************************************************/
void ExternalSequence::checkRF(SeqBlock& block)
{
	unsigned int i;
	for (i=0; i<block.rfAmplitude.size(); i++)
	{
		if (block.rfAmplitude[i]>1.0) block.rfAmplitude[i]=1.0;
		if (block.rfAmplitude[i]<0.0) block.rfAmplitude[i]=0.0;
		if (block.rfPhase[i]>TWO_PI) block.rfPhase[i]=(float)TWO_PI;
		if (block.rfPhase[i]<0)      block.rfPhase[i]=0.0;
	}
}



// * ------------------------------------------------------------------ *
// * Name        :  safeGetLine
// *
// * Description :  Read a line from the input stream handling three different line endings: 
// *                UNIX/OSX (\n), Windows (\r\n), Old Mac (\r). 
// *
// *                This allows robust reading of files on the host and MPCU
// *
// * Return      :  The stream
// *
// * ------------------------------------------------------------------ *
/*std::istream& ExternalSequence::safeGetLine(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

	return is;
    //std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            //if(t.empty())
            //    is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}*/

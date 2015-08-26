#include "ExternalSequence.h"

#include <stdio.h>		// sscanf
#include <cstring>		// strlen etc
#include <iomanip>		// std::setw etc

#include <algorithm>	// for std::max_element

#ifndef MAX
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
#endif

ExternalSequence::PrintFunPtr ExternalSequence::print_fun = &ExternalSequence::defaultPrint;
const int ExternalSequence::MAX_LINE_SIZE = 256;
const char ExternalSequence::COMMENT_CHAR = '#';

/***********************************************************/
ExternalSequence::ExternalSequence()
{
	m_blocks.clear();	// probably not needed.
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
	// Open in binary mode to ensure all end-of-line characters are processed
	data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
	data_file.seekg(0, std::ios::beg);

	if (!data_file.good())
	{
		// Try separate file mode (blocks.seq, events.seq, shapes.seq)
		filepath = path + "/shapes.seq";
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);
		isSingleFileMode = false;
	}
	if (!data_file.good())
	{
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
		return false;
	}
	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Building index" );

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
		getline(data_file, buffer, MAX_LINE_SIZE);
		sscanf(buffer, "%s%d", tmpStr, &numSamples);

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Reading shape " << shapeId );

		CompressedShape shape;
		shape.samples.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
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
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
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
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
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
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
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

		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
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
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
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
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
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
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}
		buildFileIndex(data_file);
	}


	// Read definition section
	// ------------------------
	unsigned int numBlocks = 0;
	if (m_fileIndex.find("[DEFINITIONS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[DEFINITIONS]"], std::ios::beg);

		// Read each definition line
		m_definitions.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			std::istringstream ss(buffer);
			std::string key;
			ss >> key;
			double value;
			std::vector<double> values;
			while (ss >> value) {
				values.push_back(value);
			}
			m_definitions[key] = values;
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

	// Read blocks section
	// ------------------------
	if (m_fileIndex.find("[BLOCKS]") == m_fileIndex.end()) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: [BLOCKS] section");
		return false;
	}
	data_file.seekg(m_fileIndex["[BLOCKS]"], std::ios::beg);

	int blockIdx;
	EventIDs events;

	// Read blocks
	m_blocks.clear();
	while (getline(data_file, buffer, MAX_LINE_SIZE)) {
		if (buffer[0]=='[' || strlen(buffer)==0) {
			break;
		}

		sscanf(buffer, "%d%d%d%d%d%d%d", &blockIdx,
			&events.id[DELAY],                              // Delay
			&events.id[RF],                                 // RF
			&events.id[GX],&events.id[GY],&events.id[GZ],   // Gradients
			&events.id[ADC]                                 // ADCs
			);

		if (!checkBlockReferences(events)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Block " << blockIdx
				<< " contains references to undefined events" );
			return false;
		}
		// Add event IDs to list of blocks
		m_blocks.push_back(events);
	}
	data_file.close();

	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- BLOCKS READ: " << m_blocks.size());

	// Num_Blocks definition (if defined) is used to check the correct number of blocks are read
	if (numBlocks>0 && m_blocks.size()!=numBlocks) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Expected " << numBlocks
		    << " blocks but read " << m_blocks.size() << " blocks");
		return false;
	}
	std::vector<double> def = GetDefinition("Scan_ID");
	int scanID = def.empty() ? 0: (int)def[0];
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "===== EXTERNAL SEQUENCE #" << std::setw(5) << scanID << " ===========" );
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );

	return true;
};


/***********************************************************/
void ExternalSequence::skipComments(std::ifstream &fileStream, char *buffer)
{
	while (getline(fileStream, buffer, MAX_LINE_SIZE)) {
		if (buffer[0]!=COMMENT_CHAR && strlen(buffer)>0) {
			break;
		}
	}
};


/***********************************************************/
void ExternalSequence::buildFileIndex(std::ifstream &fileStream)
{
	char buffer[MAX_LINE_SIZE];
	
	while (getline(fileStream, buffer, MAX_LINE_SIZE)) {
		std::string line = std::string(buffer);
		if (line[0]=='[' && line[line.length()-1]==']') {
			m_fileIndex[line] = fileStream.tellg();
			
		}
	}
	fileStream.clear();		// reset EOF flag
	fileStream.seekg(0, std::ios::beg);
};

/***********************************************************/
SeqBlock*	ExternalSequence::GetBlock(int index) {
	SeqBlock *block = new SeqBlock();

	// Copy event IDs
	EventIDs events = m_blocks[index];
	std::copy(events.id,events.id+NUM_EVENTS,&block->events[0]);

	// Set some defaults
	block->index = index;
	block->delay = 0;
	block->adc.numSamples = 0;

	// Set event structures (if applicable) so e.g. gradient type can be determined
	if (events.id[RF]>0)     block->rf     = m_rfLibrary[events.id[RF]];
	if (events.id[ADC]>0)    block->adc    = m_adcLibrary[events.id[ADC]];
	if (events.id[DELAY]>0)  block->delay  = m_delayLibrary[events.id[DELAY]];
	for (unsigned int i=0; i<NUM_GRADS; i++)
		if (events.id[GX+i]>0) block->grad[i] = m_gradLibrary[events.id[GX+i]];

	// Calculate duration of block
	long duration = 0;
	if (block->isRF()) {
		RFEvent &rf = block->GetRFEvent();
		duration = MAX(duration, (long)m_shapeLibrary[rf.magShape].numUncompressedSamples);
	}
	for (int iC=0; iC<NUM_GRADS; iC++)
	{
		GradEvent &grad = block->GetGradEvent(iC);
		if (block->isArbitraryGradient(iC))
			duration = MAX(duration, (long)(10*m_shapeLibrary[grad.shape].numUncompressedSamples));
		else if (block->isTrapGradient(iC))
			duration = MAX(duration, grad.rampUpTime + grad.flatTime + grad.rampDownTime);
	}
	if (block->isADC()) {
		ADCEvent &adc = block->GetADCEvent();
		duration = MAX(duration, adc.delay + (adc.numSamples*adc.dwellTime)/1000);
	}

	block->duration = duration + block->delay;

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
	for (int iC=GX; iC<ADC; iC++)
	{

		if (block->isArbitraryGradient(iC-GX))	// is arbitrary gradient?
		{
			// Decompress the arbitrary shape for this channel
			CompressedShape& shape = m_shapeLibrary[block->grad[iC-GX].shape];

			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Loaded shape with "
				<< shape.samples.size() << " compressed samples" );

			waveform.resize(shape.numUncompressedSamples);
			decompressShape(shape,&waveform[0]);

			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape uncompressed to "
				<< shape.numUncompressedSamples << " samples" );

			block->gradWaveforms[iC-GX] = std::vector<float>(waveform);
		}
	}

	// Decode ADC
	if (block->isADC())
	{
		block->adc = m_adcLibrary[events[ADC]];
	}

	// Decode Delays
	if (block->isDelay())
	{
		block->delay = m_delayLibrary[events[DELAY]];
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
	error|= (events.id[DELAY]>0 && m_delayLibrary.count(events.id[DELAY])==0);

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
		// Ensure last point is zero
		if (waveform.size()>0) waveform[waveform.size()-1]=0.0;
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
std::istream& ExternalSequence::getline(std::istream& is, char *buffer, int MAX_SIZE)
{
	int i=0;	// Current position in buffer
	for(;;) {
		char c = (char)is.get();

		switch (c) {
			case '\n':
				buffer[i]='\0';
				return is;
			case '\r':
				if(is.peek() == '\n') {
					is.get();       // Discard character
				}
				buffer[i]='\0';
				return is;
			case EOF:
				if(i==0)
					is.seekg(-1);   // Create error on stream
				buffer[i]='\0';     // In case last line has no line ending
				return is;
			default:
				buffer[i++] = c;
				if (i>=MAX_SIZE)
					return is;
		}
	}
}
